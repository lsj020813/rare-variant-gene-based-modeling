
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import argparse
import csv
import gzip
import json
import math
import os
import sys
import time

import numpy as np
import pandas as pd

from hcf import align as halign
from hcf import fastpath
from hcf import phase as hphase
from hcf import state as hstate
from hcf import tiles as htiles
from hcf import vcfio

ORIG = _config_path("${PROJECT_ROOT}/work/ref/orig_index/chr%s.vcf.gz")
RUN_ID = "HCF-20260922-v1/run_20260922T052403Z"
PROTOCOL_VERSION = "HCF-20260922-v1"

MARKERS = 16
STEP = 8
MAX_GAP = 20000
MAX_VARIANTS_PER_TILE = 4096
MIN_DISTINCT_A = 50
MAX_STATES = 32
MAX_PAIRS_PER_TILE = 64
MIN_DH, MIN_CIS, MIN_TRANS = 100, 25, 25
SPECTRUM_MIN_N = 5000
SPECTRUM_MAX_N = 20000

def _ridge_residual_ratio(Y, X, lam=1e-3):
    Y = np.asarray(Y, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    if Y.ndim == 1:
        Y = Y[:, None]
    ok = np.all(np.isfinite(X), axis=1) & np.all(np.isfinite(Y), axis=1)
    if ok.sum() < 50 or Y.shape[1] == 0:
        return float("nan"), 0
    Xo, Yo = X[ok], Y[ok]
    Xc = Xo - Xo.mean(axis=0, keepdims=True)
    sd = Xc.std(axis=0, ddof=1)
    keep = sd > 0
    if keep.sum() == 0:
        return float("nan"), int(ok.sum())
    Xc = Xc[:, keep] / sd[keep]
    Yc = Yo - Yo.mean(axis=0, keepdims=True)
    G = Xc.T.dot(Xc)
    G[np.diag_indices_from(G)] += lam * Xc.shape[0]
    B = np.linalg.solve(G, Xc.T.dot(Yc))
    R = Yc - Xc.dot(B)
    denom = (Yc ** 2).sum(axis=0)
    num = (R ** 2).sum(axis=0)
    good = denom > 0
    if not good.any():
        return float("nan"), int(ok.sum())
    return float(np.mean(num[good] / denom[good])), int(ok.sum())

def _univariate_design(Gsub):
    G = np.asarray(Gsub, dtype=np.float64)
    return np.concatenate([G, G ** 2, (G == 1.0).astype(np.float64)], axis=1)

def _finite_rows(mat):
    mat = np.asarray(mat, dtype=np.float64)
    if mat.ndim == 1:
        mat = mat[:, None]
    if mat.size == 0:
        return mat
    return mat[np.all(np.isfinite(mat), axis=1)]

def _spectrum(mat, max_n, seed=20260922):
    mat = _finite_rows(mat)
    n = mat.shape[0]
    if n > max_n:
        rs = np.random.RandomState(seed)
        sel = rs.choice(n, size=max_n, replace=False)
        mat = mat[sel]
    er = hstate.effective_rank(mat, standardize=True)
    er_raw = hstate.effective_rank(mat, standardize=False)
    sv = er_raw["singular_values"]
    return {
        "n_rows_used": int(mat.shape[0]),
        "M": int(mat.shape[1]),
        "M_eff_standardized": er["effective_rank"],
        "M_eff_raw": er_raw["effective_rank"],
        "numerical_rank": er["numerical_rank"],
        "status": er["status"],
        "n_constant_columns": er["n_constant_columns"],
        "sv_top5": ";".join("%.6g" % s for s in sv[:5]),
        "sv_sum": float(np.sum(sv)) if sv else float("nan"),
        "sv_max": float(sv[0]) if sv else float("nan"),
    }

def select_tile_pairs(subwins, pos, var_keys, max_pairs=MAX_PAIRS_PER_TILE):
    cand = {}
    for wi, (a, b) in enumerate(subwins):
        for x in range(a, b):
            for y in range(x + 1, b):
                if (x, y) not in cand:
                    cand[(x, y)] = wi
    if not cand:
        return [], 0, float("nan")
    items = sorted(cand.keys())
    dist = np.array([pos[y] - pos[x] for (x, y) in items], dtype=np.int64)
    n_possible = len(items)
    if n_possible <= max_pairs:
        chosen = items
    else:
        qs = np.quantile(dist, [0.25, 0.5, 0.75])
        bucket = np.digitize(dist, qs, right=True)
        per = max_pairs // 4
        chosen = []
        for q in range(4):
            idx = [i for i in range(n_possible) if bucket[i] == q]
            idx.sort(key=lambda i: hstate.key_hash("%s|%s" % (var_keys[items[i][0]], var_keys[items[i][1]])))
            chosen.extend(items[i] for i in idx[:per])
        if len(chosen) < max_pairs:
            rest = [it for it in items if it not in set(chosen)]
            rest.sort(key=lambda it: hstate.key_hash("%s|%s" % (var_keys[it[0]], var_keys[it[1]])))
            chosen.extend(rest[: max_pairs - len(chosen)])
    return [(x, y, cand[(x, y)]) for (x, y) in chosen], n_possible, float(len(chosen)) / n_possible

def process_tile(chrom, tile_start, vt, splits, batch_proxy, out, verify=True,
                 n_people_expected=None):
    t_all0 = time.time()
    status = {"run_id": RUN_ID, "protocol_version": PROTOCOL_VERSION, "build": "GRCh37",
              "chrom": str(chrom), "tile_id": "chr%s:%d" % (chrom, tile_start),
              "tile_start": tile_start, "n_primary_variants": int(len(vt)),
              "status": "OK", "limitation": ""}
    if len(vt) > MAX_VARIANTS_PER_TILE:
        status["status"] = "DENSE_WINDOW_RESOURCE_LIMIT"
        status["limitation"] = "n_primary=%d > cap %d; tile listed, not pruned" % (
            len(vt), MAX_VARIANTS_PER_TILE)
        out["tile_status"].append(status)
        return status
    if len(vt) < MARKERS:
        status["status"] = "TOO_FEW_MARKERS"
        status["limitation"] = "n_primary=%d < markers_per_subwindow=%d" % (len(vt), MARKERS)
        out["tile_status"].append(status)
        return status

    keep = set(zip(vt["pos"].tolist(), vt["ref"].tolist(), vt["alt"].tolist()))
    t0 = time.time()
    rd = vcfio.read_region(ORIG % chrom, "%s:%d-%d" % (chrom, tile_start + 1, tile_start + htiles.TILE_BP),
                           n_samples=n_people_expected, keep=keep)
    t_read = time.time() - t0
    if rd["hap"].shape[0] < MARKERS:
        status["status"] = "TOO_FEW_MARKERS_AFTER_READ"
        status["limitation"] = "read %d of %d primary variants" % (rd["hap"].shape[0], len(vt))
        out["tile_status"].append(status)
        return status

    pos = rd["pos"]
    var_keys = [halign.variant_key(chrom, int(pos[i]), rd["ref"][i], rd["alt"][i])
                for i in range(len(pos))]
    if not np.all(np.diff(pos) >= 0):
        raise halign.AlignmentError(
            "tile %s:%d positions not sorted (silent reordering is forbidden)" % (chrom, tile_start))
    n_dup_pos = int(np.sum(np.diff(pos) == 0))
    vkey = list(zip(vt["pos"].tolist(), vt["ref"].tolist(), vt["alt"].tolist()))
    afmap = dict(zip(vkey, vt["af"].tolist()))
    r2map = dict(zip(vkey, vt["r2"].tolist()))
    tymap = dict(zip(vkey, vt["typed"].tolist()))
    trip = [(int(pos[i]), rd["ref"][i], rd["alt"][i]) for i in range(len(pos))]
    af = np.array([afmap[t] for t in trip], dtype=np.float64)
    r2v = np.array([r2map[t] for t in trip], dtype=np.float64)
    tyv = np.array([tymap[t] for t in trip], dtype=np.float64)

    hap = vcfio.to_person_major(rd["hap"])
    n_people = hap.shape[0]
    idxA, idxB, idxC = splits["A"], splits["B"], splits["C"]

    subwins = htiles.subwindows(pos, markers=MARKERS, step=STEP, max_gap_bp=MAX_GAP)
    status["n_variants_read"] = int(len(pos))
    status["n_duplicate_positions"] = n_dup_pos
    status["n_subwindows"] = len(subwins)
    status["read_wall_s"] = round(t_read, 3)
    if not subwins:
        status["status"] = "NO_SUBWINDOW_AFTER_GAP_RULE"
        status["limitation"] = "20 kb gap rule left no run of %d markers" % MARKERS
        out["tile_status"].append(status)
        return status

    Gtile = (hap[idxA, 0, :].astype(np.float64) + hap[idxA, 1, :].astype(np.float64))
    Gtile[(hap[idxA, 0, :] < 0) | (hap[idxA, 1, :] < 0)] = np.nan
    nspec = int(min(SPECTRUM_MAX_N, max(SPECTRUM_MIN_N, 2 * len(pos))))
    sp = _spectrum(Gtile, nspec)
    out["repr_info"].append(dict(
        run_id=RUN_ID, protocol_version=PROTOCOL_VERSION, evidence_class="[실측-신규]",
        scope_id="A", build="GRCh37", chrom=str(chrom), tile_id=status["tile_id"],
        subwindow_id="", level="tile", representation_id="DS",
        n_people=int(len(idxA)), n_haplotypes=int(2 * len(idxA)), split_id="A",
        residual_ratio_vs_univariate="", residual_n="",
        spectrum_n_lt_M=int(sp["n_rows_used"] < sp["M"]), **sp))

    t_state = t_ctr = t_pair = t_spec = 0.0
    tile_pairs, n_possible_pairs, incl_p = select_tile_pairs(subwins, pos, var_keys)
    pair_rows_by_sw = {}
    for (x, y, wi) in tile_pairs:
        pair_rows_by_sw.setdefault(wi, []).append((x, y))

    for wi, (a, b) in enumerate(subwins):
        sw_id = "chr%s:%d:w%03d" % (chrom, tile_start, wi)
        sub = hap[:, :, a:b]
        t0 = time.time()
        keys = fastpath.encode_keys(sub)
        if verify and wi == 0:
            fastpath.assert_equivalent_encoding(sub, keys)
        keysA = keys[idxA]
        dic = hstate.StateDictionary.fit(np.asarray(keysA, dtype=object),
                                         min_distinct_people=MIN_DISTINCT_A,
                                         max_states=MAX_STATES)
        ZA, cols = fastpath.transform_counts(dic, keysA)
        if verify and wi == 0:
            fastpath.assert_equivalent_transform(dic, keysA, ZA)
        t_state += time.time() - t0

        allref = hstate.all_reference_key(MARKERS)
        flatA = keysA.ravel()
        uniqA, cntA = np.unique(flatA, return_counts=True)
        a_key_set = set(str(u) for u in uniqA)
        pA = cntA.astype(np.float64) / cntA.sum()
        entropy = float(-np.sum(pA * np.log(pA)))
        allref_copies = float(cntA[uniqA == allref].sum()) if allref in a_key_set else 0.0
        nonref_mask = uniqA != allref
        nonref_tot = float(cntA[nonref_mask].sum())
        rep_rate = float(np.mean(cntA[nonref_mask] >= 2)) if nonref_mask.any() else float("nan")

        row = dict(
            run_id=RUN_ID, protocol_version=PROTOCOL_VERSION, evidence_class="[실측-신규]",
            scope_id="A_dictionary", build="GRCh37", chrom=str(chrom),
            tile_id=status["tile_id"], subwindow_id=sw_id,
            variant_filter_id="noncoding_primary_r2ge0.8_or_TYPED",
            representation_id="S_EXACT", input_route="Route_P_phased_GT",
            phase_provenance="VCF_ASSUMED_CONTIG_PHASE",
            marker_start_pos=int(pos[a]), marker_end_pos=int(pos[b - 1]),
            physical_span_bp=int(pos[b - 1] - pos[a]), n_markers=MARKERS,
            maf_min=float(np.min(np.minimum(af[a:b], 1 - af[a:b]))),
            maf_median=float(np.median(np.minimum(af[a:b], 1 - af[a:b]))),
            maf_max=float(np.max(np.minimum(af[a:b], 1 - af[a:b]))),
            r2_median=float(np.median(r2v[a:b])),
            typed_fraction=float(np.mean(tyv[a:b])),
            n_people_A=int(len(idxA)), n_haplotypes_A=int(2 * len(idxA)),
            n_distinct_states_A=int(dic.diagnostics["n_distinct_observed_states"]),
            n_states_with_missing_markers=int(dic.diagnostics["n_states_with_missing_markers"]),
            n_supported_states_A=int(dic.diagnostics["n_supported_states"]),
            n_retained_states=int(dic.diagnostics["n_retained_states"]),
            n_dropped_unsupported=int(dic.diagnostics["n_states_dropped_unsupported"]),
            n_dropped_over_cap=int(dic.diagnostics["n_states_dropped_over_cap"]),
            reference_state_is_all_reference=int(dic.reference_state == allref),
            all_reference_retained=int(dic.diagnostics["all_reference_key_retained"]),
            all_reference_copy_share_A=allref_copies / float(cntA.sum()),
            top_nonreference_share_A=(float(cntA[nonref_mask].max()) / nonref_tot)
            if (nonref_mask.any() and nonref_tot > 0) else float("nan"),
            nonreference_state_repeat_rate_A=rep_rate,
            entropy_A_nats=entropy, exp_entropy_A=float(math.exp(entropy)),
            entropy_is_not_information="entropy/exp(entropy) is not independent sample size nor phenotype information",
            min_distinct_people_threshold=MIN_DISTINCT_A, max_supported_states=MAX_STATES,
            support_threshold_is_power_proof=False,
            batch_proxy_assoc="[미정]",
            status="OK", limitation="",
        )
        for nm, idx in (("A", idxA), ("B", idxB), ("C", idxC)):
            if len(idx) == 0:
                continue
            kk = keys[idx]
            Z, _ = fastpath.transform_counts(dic, kk)
            tot = float(Z.sum())
            row["n_people_%s" % nm] = int(len(idx))
            row["other_copy_share_%s" % nm] = float(Z[:, cols.index(dic.other_label)].sum() / tot)
            flat = kk.ravel()
            uq, ct = np.unique(flat, return_counts=True)
            unseen = float(sum(int(ct[i]) for i in range(uq.size) if str(uq[i]) not in a_key_set))
            row["unseen_in_A_copy_share_%s" % nm] = unseen / float(flat.size)
            row["coverage_by_dictionary_%s" % nm] = 1.0 - row["other_copy_share_%s" % nm]
            row["all_reference_copy_share_%s" % nm] = float(
                ct[uq == allref].sum()) / float(flat.size) if (uq == allref).any() else 0.0
            row["n_missing_bearing_copies_%s" % nm] = int(
                sum(int(ct[i]) for i in range(uq.size) if hstate.has_missing(str(uq[i]))))

        if batch_proxy is not None:
            t0 = time.time()
            lab = np.argmax(ZA, axis=1)
            ct = pd.crosstab(lab, batch_proxy[idxA])
            if ct.shape[0] > 1 and ct.shape[1] > 1:
                from scipy.stats import chi2_contingency
                try:
                    chi2, pv, dof, _ = chi2_contingency(ct.values)
                    n = ct.values.sum()
                    cramv = math.sqrt(chi2 / (n * (min(ct.shape) - 1)))
                    row["batch_proxy_assoc"] = "cramers_v=%.4g;p=%.3g;dof=%d" % (cramv, pv, dof)
                except Exception as exc:
                    row["batch_proxy_assoc"] = "PROXY_TEST_FAILED:%s" % type(exc).__name__
            row["batch_proxy_note"] = ("sample-ID-prefix/order proxy only; true cohort(CT/NC) "
                                       "and PC association = [미정] (phenotype file not read)")
            t_ctr += time.time() - t0
        out["state_inventory"].append(row)

        t0 = time.time()
        gsub = sub[idxA, 0, :].astype(np.int16) + sub[idxA, 1, :].astype(np.int16)
        miss = np.any(sub[idxA] < 0, axis=(1, 2))
        gb = np.ascontiguousarray((gsub + 1).astype(np.uint8)).view("S%d" % MARKERS).ravel()
        k1 = keysA[:, 0]
        k2 = keysA[:, 1]
        lo = np.where(k1 <= k2, k1, k2)
        hi = np.where(k1 <= k2, k2, k1)
        dip = np.char.add(np.char.add(lo, "|"), hi)
        df = pd.DataFrame({"g": gb[~miss], "d": dip[~miss]})
        if len(df):
            nun = df.groupby("g")["d"].nunique()
            sz = df.groupby("g").size()
            multi = nun > 1
            out["diplotype_contrast"].append(dict(
                run_id=RUN_ID, protocol_version=PROTOCOL_VERSION, evidence_class="[실측-신규]",
                scope_id="A", build="GRCh37", chrom=str(chrom), tile_id=status["tile_id"],
                subwindow_id=sw_id, representation_id="S_EXACT", split_id="A",
                n_markers=MARKERS, n_people=int(len(df)),
                n_people_excluded_missing=int(miss.sum()),
                n_unphased_genotype_groups=int(len(nun)),
                n_groups_multi_diplotype=int(multi.sum()),
                n_people_in_multi_diplotype_groups=int(sz[multi].sum()),
                share_people_in_multi_diplotype_groups=float(sz[multi].sum()) / float(len(df)),
                n_singleton_genotype_groups=int((sz == 1).sum()),
                share_people_in_singleton_groups=float(sz[sz == 1].sum()) / float(len(df)),
                max_group_size=int(sz.max()), median_group_size=float(sz.median()),
                max_diplotypes_in_a_group=int(nun.max()),
                trivial_one_marker=0,
                identifiability_note="rare direct matched contrast != proven non-identifiability (§10.2)",
                status="OK", limitation=""))
        t_ctr += time.time() - t0

        t0 = time.time()
        for (x, y) in pair_rows_by_sw.get(wi, []):
            diag = hphase.pair_diagnostics(hap[idxA], x, y,
                                           min_double_het=MIN_DH, min_cis=MIN_CIS,
                                           min_trans=MIN_TRANS)
            q = hphase.phase_contrast_q(hap[idxA], x, y)
            fin = np.isfinite(q)
            mafx = min(af[x], 1 - af[x])
            mafy = min(af[y], 1 - af[y])
            out["phase_pairs"].append(dict(
                run_id=RUN_ID, protocol_version=PROTOCOL_VERSION, evidence_class="[실측-신규]",
                scope_id="A", build="GRCh37", chrom=str(chrom), tile_id=status["tile_id"],
                subwindow_id=sw_id, representation_id="Q_PHASE_CONTRAST", split_id="A",
                pair_key="%s__%s" % (var_keys[x], var_keys[y]),
                marker_j=int(x), marker_k=int(y),
                pos_j=int(pos[x]), pos_k=int(pos[y]), distance_bp=int(pos[y] - pos[x]),
                maf_j=float(mafx), maf_k=float(mafy),
                maf_min_of_pair=float(min(mafx, mafy)),
                n_people=int(len(idxA)), n_missing=int(diag["n_missing"]),
                n_double_het=int(diag["n_double_het"]), n_cis=int(diag["n_cis"]),
                n_trans=int(diag["n_trans"]),
                q_variance=float(np.var(q[fin], ddof=1)) if fin.sum() > 1 else float("nan"),
                q_mean=float(np.mean(q[fin])) if fin.any() else float("nan"),
                haplotype_r2=diag["haplotype_r2"],
                supported=int(diag["supported"]),
                flags=";".join(diag["flags"]),
                unsupported_reasons=";".join(diag["unsupported_reasons"]),
                n_possible_pairs_in_tile=int(n_possible_pairs),
                inclusion_probability=float(incl_p),
                selection_rule="within-subwindow pairs; distance-quartile balanced; sha256 hash order; phenotype-blind",
                support_threshold_is_power_proof=False,
                status="OK", limitation=""))
        t_pair += time.time() - t0

        t0 = time.time()
        Gsw = gsub.astype(np.float64)
        Gsw[np.any(sub[idxA] < 0, axis=1)] = np.nan
        U = _univariate_design(Gsw)
        Zdes, dcols = dic.design_matrix(np.asarray(keysA, dtype=object))
        r_state, n_state = _ridge_residual_ratio(Zdes, U)
        sp_ds = _spectrum(Gsw, SPECTRUM_MAX_N)
        sp_st = _spectrum(Zdes, SPECTRUM_MAX_N)
        base = dict(run_id=RUN_ID, protocol_version=PROTOCOL_VERSION, evidence_class="[실측-신규]",
                    scope_id="A", build="GRCh37", chrom=str(chrom), tile_id=status["tile_id"],
                    subwindow_id=sw_id, level="subwindow", split_id="A",
                    n_people=int(len(idxA)), n_haplotypes=int(2 * len(idxA)))
        out["repr_info"].append(dict(base, representation_id="DS",
                                     residual_ratio_vs_univariate="", residual_n="",
                                     spectrum_n_lt_M=int(sp_ds["n_rows_used"] < sp_ds["M"]), **sp_ds))
        out["repr_info"].append(dict(base, representation_id="S_EXACT",
                                     residual_ratio_vs_univariate=r_state, residual_n=n_state,
                                     spectrum_n_lt_M=int(sp_st["n_rows_used"] < sp_st["M"]), **sp_st))
        pl = pair_rows_by_sw.get(wi, [])
        if pl:
            Qm = np.column_stack([hphase.phase_contrast_q(hap[idxA], x, y, check_identity=False)
                                  for (x, y) in pl])
            r_q, n_q = _ridge_residual_ratio(Qm, U)
            sp_q = _spectrum(Qm, SPECTRUM_MAX_N)
            out["repr_info"].append(dict(base, representation_id="Q_PHASE_CONTRAST",
                                         residual_ratio_vs_univariate=r_q, residual_n=n_q,
                                         spectrum_n_lt_M=int(sp_q["n_rows_used"] < sp_q["M"]),
                                         **sp_q))
        t_spec += time.time() - t0

    status.update({"state_wall_s": round(t_state, 2), "contrast_wall_s": round(t_ctr, 2),
                   "pair_wall_s": round(t_pair, 2), "spectrum_wall_s": round(t_spec, 2),
                   "total_wall_s": round(time.time() - t_all0, 2),
                   "n_pairs_selected": len(tile_pairs),
                   "n_possible_pairs": int(n_possible_pairs)})
    out["tile_status"].append(status)
    return status

def load_variants(chrom, vardir):
    p = os.path.join(vardir, "chr%s.primary.tsv.gz" % chrom)
    return pd.read_csv(p, sep="\t")

def main(argv=None):
    ap = argparse.ArgumentParser("hcf.build_states")
    ap.add_argument("--tilelist", required=True)
    ap.add_argument("--vardir", required=True)
    ap.add_argument("--split", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tmpdir", required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-verify", action="store_true")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.tmpdir, exist_ok=True)
    sp = pd.read_csv(args.split, sep="\t")
    splits = dict((k, sp.index[sp["split"] == k].values) for k in ("A", "B", "C"))
    n_people = len(sp)
    batch_proxy = sp["batch_proxy"].values if "batch_proxy" in sp.columns else None

    tl = pd.read_csv(args.tilelist, sep="\t")
    tl = tl.iloc[args.shard::args.nshards].reset_index(drop=True)
    if args.limit:
        tl = tl.iloc[: args.limit]

    out = {"state_inventory": [], "diplotype_contrast": [], "phase_pairs": [],
           "repr_info": [], "tile_status": [], "tmpdir": args.tmpdir}
    cur_chrom, vt_all = None, None
    for i in range(len(tl)):
        ch = str(tl["chrom"].iloc[i])
        ts = int(tl["tile_start"].iloc[i])
        if ch != cur_chrom:
            vt_all = load_variants(ch, args.vardir)
            cur_chrom = ch
        vt = vt_all[(vt_all["pos"] >= ts) & (vt_all["pos"] < ts + htiles.TILE_BP)]
        try:
            st = process_tile(ch, ts, vt.reset_index(drop=True), splits, batch_proxy, out,
                              verify=not args.no_verify, n_people_expected=n_people)
        except Exception as exc:
            out["tile_status"].append({"run_id": RUN_ID, "chrom": ch,
                                       "tile_id": "chr%s:%d" % (ch, ts), "tile_start": ts,
                                       "status": "TILE_FAILED",
                                       "limitation": "%s: %s" % (type(exc).__name__, exc)})
            st = {"status": "TILE_FAILED"}
        sys.stderr.write("[shard %d] %d/%d chr%s:%d %s\n" % (args.shard, i + 1, len(tl), ch, ts,
                                                             st.get("status")))
        sys.stderr.flush()

    for name in ("state_inventory", "diplotype_contrast", "phase_pairs", "repr_info", "tile_status"):
        rows = out[name]
        path = os.path.join(args.outdir, "%s.shard%02d.csv" % (name, args.shard))
        if rows:
            pd.DataFrame(rows).to_csv(path, index=False)
        else:
            open(path, "w").write("")
    return 0

if __name__ == "__main__":
    sys.exit(main())
