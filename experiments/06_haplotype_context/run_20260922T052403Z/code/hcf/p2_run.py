
import argparse
import errno
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

from hcf import align as halign
from hcf import fastpath
from hcf import p2
from hcf import phase as hphase
from hcf import state as hstate
from hcf import tiles as htiles
from hcf import vcfio

HAPLA = None

def log(msg):
    sys.stderr.write("[%s] %s\n" % (time.strftime("%H:%M:%S"), msg))
    sys.stderr.flush()

def sh(cmd, timeout=None):
    r = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError("cmd failed (%d): %s\n%s" % (r.returncode, cmd,
                                                        r.stdout.decode()[-2000:]))
    return r.stdout.decode()

class Ctx(object):
    def __init__(self, args):
        self.args = args
        run = args.run
        self.run = run
        sp = pd.read_csv(os.path.join(run, "private/sample_split_manifest.private.tsv"), sep="\t")
        ph = pd.read_csv(args.pheno, sep="\t")
        cols = ["sample_id", "y"] + p2.COV_COLS
        ph = ph[cols]
        sp_ab = sp[sp["split"].isin(["A", "B"])]
        mg = sp_ab.merge(ph, on="sample_id", how="inner")
        mg = mg.dropna(subset=["y"] + p2.COV_COLS)
        mg = mg.sort_values("sample_index").reset_index(drop=True)
        self.n = len(mg)
        self.sample_ids = mg["sample_id"].tolist()
        self.split = mg["split"].values
        self.y = mg["y"].values.astype(np.float32)
        self.cov = mg[p2.COV_COLS].values.astype(np.float32)
        self.idxA = np.where(self.split == "A")[0]
        self.idxB = np.where(self.split == "B")[0]
        self.nA, self.nB = len(self.idxA), len(self.idxB)
        self.s2A = float(np.var(self.y[self.idxA], ddof=1))
        rng = np.random.RandomState(p2.SEED)
        perm = rng.permutation(self.nA)
        self.folds = [self.idxA[perm[k::p2.INNER_FOLDS]] for k in range(p2.INNER_FOLDS)]
        self.samples_file = os.path.join(args.tmp, "p2_samples.txt")
        with open(self.samples_file, "w") as fh:
            fh.write("\n".join(self.sample_ids) + "\n")
        pp = pd.read_csv(os.path.join(run, "states/phase_pair_support.csv"),
                         usecols=["chrom", "tile_id", "pair_key", "pos_j", "pos_k", "supported"])
        pp = pp[pp["supported"] == True]
        self.sup_pairs = {}
        for t, g in pp.groupby("tile_id"):
            self.sup_pairs[t] = list(zip(g["pos_j"].astype(int), g["pos_k"].astype(int),
                                         g["pair_key"]))
        self.vardir = os.path.join(run, "variants")
        self._vt_chrom, self._vt = None, None

    def variants(self, chrom):
        if self._vt_chrom != chrom:
            self._vt = pd.read_csv(os.path.join(self.vardir, "chr%s.primary.tsv.gz" % chrom),
                                   sep="\t")
            self._vt_chrom = chrom
        return self._vt

def hapla_blocks(ctx, chrom, tile_start, pos_keep, train_rows, tag, tmpd):
    try:
        return _hapla_blocks(ctx, chrom, tile_start, pos_keep, train_rows, tag, tmpd)
    except Exception as exc:
        msg = "%s: %s" % (type(exc).__name__, str(exc)[:200])
        clog = os.path.join(tmpd, "hap_%s.clog" % tag)
        if os.path.exists(clog):
            try:
                msg += " | clog: " + open(clog).read()[-300:].replace("\n", " ")
            except OSError:
                pass
        return None, None, "H_CLUSTER_FAILED(%s)" % msg

def _hapla_blocks(ctx, chrom, tile_start, pos_keep, train_rows, tag, tmpd):
    if HAPLA is None:
        return None, None, "H_CLUSTER_DISABLED"
    pref = os.path.join(tmpd, "hap_%s" % tag)
    allbcf = os.path.join(tmpd, "tile_all.bcf")
    if not os.path.exists(allbcf):
        posf = os.path.join(tmpd, "tile_pos.tsv")
        with open(posf, "w") as fh:
            for pth in pos_keep:
                fh.write("%s\t%d\n" % (chrom, pth))
        sh("%s view -r %s:%d-%d -T %s -S %s -Ou %s | %s annotate -x FORMAT/GP,FORMAT/DS,INFO "
           "-Ob -o %s" % (p2.BCFTOOLS, chrom, tile_start + 1, tile_start + htiles.TILE_BP,
                          posf, ctx.samples_file, p2.ORIG % chrom, p2.BCFTOOLS, allbcf))
        sh("%s index -f %s" % (p2.BCFTOOLS, allbcf))
    trf = os.path.join(tmpd, "train_%s.txt" % tag)
    with open(trf, "w") as fh:
        fh.write("\n".join(ctx.sample_ids[i] for i in sorted(train_rows)) + "\n")
    trbcf = os.path.join(tmpd, "train_%s.bcf" % tag)
    sh("%s view -S %s -Ob -o %s %s" % (p2.BCFTOOLS, trf, trbcf, allbcf))
    sh("%s index -f %s" % (p2.BCFTOOLS, trbcf))
    sh("nice -n 19 ionice -c3 %s cluster -g %s -f 16 -s 8 -p 0.1 --max-clusters 32 "
       "--min-mac 100 --min-freq 0.000001 --medians -t 1 -o %s > %s.clog 2>&1"
       % (HAPLA, trbcf, pref, pref))
    sh("nice -n 19 ionice -c3 %s predict -g %s -r %s -t 1 -o %s.prj > %s.plog 2>&1"
       % (HAPLA, allbcf, pref, pref, pref))
    bca = pref + ".prj.bca"
    if not os.path.exists(bca):
        return None, None, "H_CLUSTER_NO_PROJECTION"
    raw = np.fromfile(bca, dtype=np.uint8)
    nhap = 2 * ctx.n
    hdr = len(raw) % nhap
    if hdr not in (0, 1, 2, 3, 4, 8):
        return None, None, "H_CLUSTER_BCA_UNEXPECTED_SIZE"
    arr = raw[hdr:].reshape(-1, nhap)
    Z, names = p2.cluster_block(arr, ctx.n)
    for f in (trbcf, trbcf + ".csi", trf):
        try:
            os.remove(f)
        except OSError:
            pass
    return Z, names, "OK"

def build_tile_inputs(ctx, chrom, tile_start, tile_id):
    info = {}
    vt_all = ctx.variants(chrom)
    vt = vt_all[(vt_all["pos"] >= tile_start) & (vt_all["pos"] < tile_start + htiles.TILE_BP)]
    vt = vt.reset_index(drop=True)
    info["n_primary_tile"] = int(len(vt))
    if len(vt) < p2.MARKERS:
        return None, dict(info, status="TOO_FEW_MARKERS")
    keep = set(zip(vt["pos"].tolist(), vt["ref"].tolist(), vt["alt"].tolist()))
    region = "%s:%d-%d" % (chrom, tile_start + 1, tile_start + htiles.TILE_BP)
    t0 = time.time()
    rd = p2.read_gt(chrom, region, ctx.samples_file, ctx.n, keep)
    rdd = p2.read_ds(chrom, region, ctx.samples_file, ctx.n, keep=keep)
    info["read_wall_s"] = round(time.time() - t0, 2)
    if rd["hap"].shape[0] < p2.MARKERS:
        return None, dict(info, status="TOO_FEW_MARKERS_AFTER_READ")
    if list(rd["pos"]) != list(rdd["pos"]) or rd["ref"] != rdd["ref"] or rd["alt"] != rdd["alt"]:
        raise halign.AlignmentError("GT/DS variant order mismatch in tile %s" % tile_id)
    pos = rd["pos"]
    if not np.all(np.diff(pos) >= 0):
        raise halign.AlignmentError("tile %s positions not sorted" % tile_id)
    hap = vcfio.to_person_major(rd["hap"])
    DS = np.ascontiguousarray(rdd["ds"].T)
    var_keys = [halign.variant_key(chrom, int(pos[i]), rd["ref"][i], rd["alt"][i])
                for i in range(len(pos))]

    posmap = {}
    for i, pth in enumerate(pos):
        posmap.setdefault(int(pth), i)
    pairs = []
    for (pj, pk, pkey) in ctx.sup_pairs.get(tile_id, []):
        if pj in posmap and pk in posmap:
            pairs.append((posmap[pj], posmap[pk], pkey))
    info["n_supported_pairs"] = len(pairs)

    m = len(pos)
    info["prune_applied"] = 0
    info["prune_r2"] = ""
    if m > p2.MAX_TILE_MARKERS:
        force = sorted(set([x for pr in pairs for x in pr[:2]]))
        R = p2.corr_matrix(DS, ctx.idxA)
        sel = None
        for r2 in p2.TILE_PRUNE_R2_LADDER:
            cand = p2.greedy_ld_prune(R, r2, force_keep=force)
            if len(cand) <= p2.MAX_TILE_MARKERS:
                sel = (cand, r2)
                break
        if sel is None:
            cand = p2.greedy_ld_prune(R, p2.TILE_PRUNE_R2_LADDER[-1], force_keep=force,
                                      max_keep=p2.MAX_TILE_MARKERS)
            sel = (cand, p2.TILE_PRUNE_R2_LADDER[-1])
        idxk = np.asarray(sel[0], dtype=np.int64)
        remap = dict((int(o), i) for i, o in enumerate(idxk))
        pairs = [(remap[a], remap[b], k) for (a, b, k) in pairs if a in remap and b in remap]
        pos = pos[idxk]
        hap = np.ascontiguousarray(hap[:, :, idxk])
        DS = np.ascontiguousarray(DS[:, idxk])
        var_keys = [var_keys[i] for i in idxk]
        info["prune_applied"] = 1
        info["prune_r2"] = sel[1]
        info["n_markers_before_prune"] = m
        del R
    info["n_markers"] = int(len(pos))

    lo = max(0, tile_start - p2.NUISANCE_FLANK_BP)
    hi = tile_start + htiles.TILE_BP + p2.NUISANCE_FLANK_BP
    vf = vt_all[((vt_all["pos"] >= lo) & (vt_all["pos"] < tile_start)) |
                ((vt_all["pos"] >= tile_start + htiles.TILE_BP) & (vt_all["pos"] < hi))]
    vf = vf[vf["af"].between(0.01, 0.99)].reset_index(drop=True)
    info["n_flank_candidates_all"] = int(len(vf))
    if len(vf) > p2.MAX_FLANK_CANDIDATES:
        take = np.linspace(0, len(vf) - 1, p2.MAX_FLANK_CANDIDATES).round().astype(int)
        vf = vf.iloc[np.unique(take)].reset_index(drop=True)
    NUIS = np.zeros((ctx.n, 0), dtype=np.float32)
    nuis_keys = []
    if len(vf):
        regf = os.path.join(ctx.args.tmp, "flank_%s.tsv" % tile_id.replace(":", "_"))
        with open(regf, "w") as fh:
            for pth in vf["pos"].tolist():
                fh.write("%s\t%d\n" % (chrom, pth))
        fkeep = set(zip(vf["pos"].tolist(), vf["ref"].tolist(), vf["alt"].tolist()))
        rf = p2.read_ds(chrom, None, ctx.samples_file, ctx.n, keep=fkeep, regions_file=regf)
        os.remove(regf)
        F = np.ascontiguousarray(rf["ds"].T)
        if F.shape[1]:
            Rf = p2.corr_matrix(F, ctx.idxA)
            kf = p2.greedy_ld_prune(Rf, p2.NUIS_LD_R2, max_keep=p2.MAX_NUISANCE_COLS)
            NUIS = np.ascontiguousarray(F[:, kf])
            nuis_keys = [halign.variant_key(chrom, int(rf["pos"][i]), rf["ref"][i], rf["alt"][i])
                         for i in kf]
            del Rf, F
    info["n_nuisance_cols"] = int(NUIS.shape[1])
    info["nuisance_coding_variants"] = 0
    info["nuisance_scope"] = "tile+-500kb primary noncoding, MAF>=0.01, LD-pruned r2<%.2f in A" % p2.NUIS_LD_R2

    G = (hap[:, 0, :].astype(np.float32) + hap[:, 1, :].astype(np.float32))
    miss = (hap[:, 0, :] < 0) | (hap[:, 1, :] < 0)
    info["n_missing_calls"] = int(miss.sum())
    if info["n_missing_calls"]:
        G[miss] = np.nan
    PAIR = np.zeros((ctx.n, len(pairs)), dtype=np.float32)
    Q = np.zeros((ctx.n, len(pairs)), dtype=np.float32)
    for c, (a, b, _k) in enumerate(pairs):
        PAIR[:, c] = G[:, a] * G[:, b]
        Q[:, c] = (hap[:, 0, a].astype(np.float32) - hap[:, 1, a].astype(np.float32)) * \
                  (hap[:, 0, b].astype(np.float32) - hap[:, 1, b].astype(np.float32))
    if len(pairs):
        rng = np.random.RandomState(p2.SEED)
        sel = rng.choice(ctx.n, size=min(256, ctx.n), replace=False)
        hsel = np.ascontiguousarray(hap[sel])
        for c, (a, b, _k) in enumerate(pairs[:4]):
            ref = np.asarray(hphase.phase_contrast_q(hsel, a, b), dtype=np.float32)
            if not np.allclose(ref, Q[sel, c], equal_nan=True):
                raise AssertionError("vectorised Q differs from hcf.phase.phase_contrast_q")

    subwins = htiles.subwindows(pos, markers=p2.MARKERS, step=p2.STEP, max_gap_bp=p2.MAX_GAP)
    info["n_subwindows"] = len(subwins)
    if not subwins:
        return None, dict(info, status="NO_SUBWINDOW_AFTER_GAP_RULE")
    keys_by_sw = p2.encode_subwindow_keys(hap, subwins, verify=True)

    return dict(pos=pos, hap=hap, DS=DS, G=G, PAIR=PAIR, Q=Q, NUIS=NUIS,
                pairs=pairs, var_keys=var_keys, nuis_keys=nuis_keys, subwins=subwins,
                keys_by_sw=keys_by_sw), dict(info, status="OK")

def run_tile(ctx, chrom, tile_start, tmpd):
    tile_id = "chr%s:%d" % (chrom, tile_start)
    t_start = time.time()
    blocks, info = build_tile_inputs(ctx, chrom, tile_start, tile_id)
    if blocks is None:
        return [], None, info
    n_pairs = len(blocks["pairs"])
    arms = list(p2.ARMS)
    arm_status = dict((a, "OK") for a in arms)
    if n_pairs == 0:
        arm_status["B_PAIR"] = "NOT_IDENTIFIABLE_NO_SUPPORTED_PAIR"
        arm_status["H_PHASE"] = "NOT_IDENTIFIABLE_NO_SUPPORTED_PAIR"

    nfix = blocks["NUIS"].shape[1] + 4 * blocks["DS"].shape[1] + 2 * n_pairs
    cands = [(a, l) for a in p2.ALPHAS for l in p2.LAMBDAS]
    ncand = len(cands)
    oof = dict((a, np.full((ncand, ctx.nA), np.nan, dtype=np.float64)) for a in arms)
    fold_mse = dict((a, np.full((ncand, p2.INNER_FOLDS), np.nan)) for a in arms)
    conv = dict((a, np.ones((ncand, p2.INNER_FOLDS), dtype=bool)) for a in arms)
    posA = dict((int(v), i) for i, v in enumerate(ctx.idxA))
    notes = {"hapla": "OK"}
    sel_final = {}
    timing = {"gram": 0.0, "solve": 0.0, "state": 0.0, "hapla": 0.0}

    def pass_setup(fk, train):
        t0 = time.time()
        ZEX, zex_names, dicts = p2.state_block(blocks["keys_by_sw"], train,
                                               verify_first=(fk == 0))
        timing["state"] += time.time() - t0
        t0 = time.time()
        ZCL, zcl_names, cst = hapla_blocks(ctx, chrom, tile_start, blocks["pos"], train,
                                           "f%d" % fk, tmpd)
        timing["hapla"] += time.time() - t0
        if ZCL is None:
            ZCL = np.zeros((ctx.n, 0), dtype=np.float32)
            zcl_names = []
            notes["hapla"] = cst
            arm_status["H_CLUSTER"] = cst
        m_tile = blocks["DS"].shape[1]
        ptot = nfix + ZEX.shape[1] + ZCL.shape[1]
        if fk == 0:
            info["n_design_cols_union"] = int(ptot)
            info["n_state_cols"] = int(ZEX.shape[1])
            info["n_cluster_cols"] = int(ZCL.shape[1])
            if ptot > p2.MAX_DESIGN_COLS:
                info["limitation"] = "union design %d > cap %d" % (ptot, p2.MAX_DESIGN_COLS)
        widths = [("NUIS", blocks["NUIS"].shape[1]), ("DS", m_tile), ("IND", 3 * m_tile),
                  ("PAIR", blocks["PAIR"].shape[1]), ("ZEX", ZEX.shape[1]),
                  ("ZCL", ZCL.shape[1]), ("Q", blocks["Q"].shape[1])]
        off, colmap = 0, {}
        for name, w in widths:
            colmap[name] = np.arange(off, off + w)
            off += w
        ncov = ctx.cov.shape[1]
        U = np.empty((ctx.n, 1 + ncov + ptot + 1), dtype=np.float32)
        U[:, 0] = 1.0
        U[:, 1:1 + ncov] = ctx.cov
        base = 1 + ncov
        for name, blk in (("NUIS", blocks["NUIS"]), ("DS", blocks["DS"]),
                          ("PAIR", blocks["PAIR"]), ("ZEX", ZEX), ("ZCL", ZCL),
                          ("Q", blocks["Q"])):
            if blk.shape[1]:
                U[:, base + colmap[name][0]: base + colmap[name][-1] + 1] = blk
        io = base + colmap["IND"][0]
        np.multiply(blocks["DS"], blocks["DS"], out=U[:, io:io + m_tile])
        U[:, io + m_tile:io + 2 * m_tile] = blocks["G"]
        np.multiply(blocks["G"], blocks["G"], out=U[:, io + 2 * m_tile:io + 3 * m_tile])
        U[:, -1] = ctx.y
        t0 = time.time()
        ff = p2.FoldFit(U, ncov, train)
        timing["gram"] += time.time() - t0
        return ff, colmap, U, (zex_names, zcl_names, dicts)

    for fk in range(p2.INNER_FOLDS):
        val = ctx.folds[fk]
        train = np.setdiff1d(ctx.idxA, val, assume_unique=False)
        ff, colmap, U, _names = pass_setup(fk, train)
        t0 = time.time()
        rrows = np.array([posA[int(v)] for v in val])
        for arm in arms:
            if arm_status[arm] != "OK":
                continue
            cols = (np.concatenate([colmap[b] for b in p2.ARM_BLOCKS[arm]])
                    if p2.ARM_BLOCKS[arm] else np.zeros(0, dtype=np.int64))
            w0, prev_alpha = None, None
            for ci, (a, l) in enumerate(cands):
                if prev_alpha != a:
                    w0, prev_alpha = None, a
                fit = ff.solve(cols, a, l, w0=w0)
                w0 = fit["beta_std"] if a > 0 else None
                yh = ff.predict(fit, val)
                oof[arm][ci, rrows] = yh
                fold_mse[arm][ci, fk] = float(np.mean((ctx.y[val] - yh) ** 2))
                conv[arm][ci, fk] = bool(fit["converged"])
        timing["solve"] += time.time() - t0
        del U, ff
        log("  %s fold%d done (%.0fs cum)" % (tile_id, fk, time.time() - t_start))

    sel, cv_grid = {}, []
    for arm in arms:
        if arm_status[arm] != "OK":
            continue
        fm = fold_mse[arm]
        cvm = np.nanmean(fm, axis=1)
        cvse = np.nanstd(fm, axis=1, ddof=1) / np.sqrt(p2.INNER_FOLDS)
        elig = np.where(conv[arm].all(axis=1))[0]
        n_nonconv = int(ncand - elig.size)
        all_nonconv = bool(elig.size == 0)
        if all_nonconv:
            elig = np.arange(ncand)
        cvm_e = np.where(np.isin(np.arange(ncand), elig), cvm, np.inf)
        jmin = int(np.nanargmin(cvm_e))
        pick, _ = p2.one_se_select(cvm_e, cvse[jmin], cands)
        sel[arm] = dict(n_nonconverged_candidates=n_nonconv, all_candidates_nonconverged=all_nonconv,
                        alpha=cands[pick][0], lam=cands[pick][1], cv_mse=float(cvm[pick]),
                        cv_mse_min=float(cvm[jmin]), cv_se_at_min=float(cvse[jmin]),
                        alpha_min=cands[jmin][0], lam_min=cands[jmin][1],
                        oof_sse=float(np.nansum((ctx.y[ctx.idxA] - oof[arm][pick]) ** 2)),
                        oof_sse_cvmin=float(np.nansum((ctx.y[ctx.idxA] - oof[arm][jmin]) ** 2)))
        cv_grid.extend(dict(tile_id=tile_id, chrom=str(chrom), representation_id=arm,
                            alpha=cands[ci][0], lam=cands[ci][1], cv_mse=float(cvm[ci]),
                            cv_se=float(cvse[ci]), converged=bool(conv[arm][ci].all()),
                            is_cv_min=(ci == jmin), is_one_se=(ci == pick))
                       for ci in range(ncand))

    ff, colmap, U, (zex_names, zcl_names, dicts) = pass_setup(p2.INNER_FOLDS, ctx.idxA)
    yB = ctx.y[ctx.idxB].astype(np.float64)
    ssB_tot = float(np.sum((yB - yB.mean()) ** 2))
    predB, fitinfo = {}, {}
    t0 = time.time()
    for arm in arms:
        if arm_status[arm] != "OK":
            continue
        cols = (np.concatenate([colmap[b] for b in p2.ARM_BLOCKS[arm]])
                if p2.ARM_BLOCKS[arm] else np.zeros(0, dtype=np.int64))
        s = sel[arm]
        fit = ff.solve(cols, s["alpha"], s["lam"], want_df=True)
        predB[arm] = ff.predict(fit, ctx.idxB)
        fitinfo[arm] = dict(n_features=int(fit["keep"].size), df=float(fit["df"]),
                            nnz=int(fit["nnz"]))
        if (s["alpha_min"], s["lam_min"]) == (s["alpha"], s["lam"]):
            predB[arm + "__cvmin"] = predB[arm]
            fitinfo[arm]["df_cvmin"] = float(fit["df"])
            fitinfo[arm]["nnz_cvmin"] = int(fit["nnz"])
        else:
            f2 = ff.solve(cols, s["alpha_min"], s["lam_min"], want_df=True)
            predB[arm + "__cvmin"] = ff.predict(f2, ctx.idxB)
            fitinfo[arm]["df_cvmin"] = float(f2["df"])
            fitinfo[arm]["nnz_cvmin"] = int(f2["nnz"])
    timing["solve"] += time.time() - t0
    del U, ff

    rows = []
    cov_sse = float(np.sum((yB - predB["B_COV"]) ** 2)) if "B_COV" in predB else float("nan")
    for arm in arms:
        base = dict(run_id=p2.RUN_ID, protocol_version=p2.PROTOCOL_VERSION,
                    evidence_class="[실측-신규]", scope_id="A_dev_B_select", build="GRCh37",
                    chrom=str(chrom), tile_id=tile_id, subwindow_id="",
                    variant_filter_id="primary_noncoding_r2ge0.8" +
                    ("_LDpruned" if info.get("prune_applied") else ""),
                    representation_id=arm, input_route="P_phased_GT",
                    phase_provenance="VCF_ASSUMED_CONTIG_PHASE",
                    split_id="A_train/B_select", n_people_A=ctx.nA, n_people_B=ctx.nB,
                    n_haplotypes_A=2 * ctx.nA, family_aware=False,
                    n_markers=info.get("n_markers"), n_supported_pairs=n_pairs,
                    source_path=p2.ORIG % chrom, software_version="numpy%s/sklearn-gram" % np.__version__,
                    status=arm_status[arm], limitation=info.get("limitation", ""))
        if arm_status[arm] != "OK":
            rows.append(dict(base, alpha_selected="", lambda_selected="", n_features="",
                             df_effective="", cv_mse_selected="", cv_mse_min="",
                             cv_se_at_min="", oof_r2_A="", sse_B="", mse_B="", r2_B="",
                             gain_vs_BCOV_normalized=""))
            continue
        s = sel[arm]
        sseB = float(np.sum((yB - predB[arm]) ** 2))
        sseB_cv = float(np.sum((yB - predB[arm + "__cvmin"]) ** 2))
        ssA_tot = float(np.sum((ctx.y[ctx.idxA] - ctx.y[ctx.idxA].mean()) ** 2))
        degenerate = bool(np.allclose(predB[arm], predB["B_COV"], atol=1e-12))
        rows.append(dict(
            base, alpha_selected=s["alpha"], lambda_selected=s["lam"],
            n_features=fitinfo[arm]["n_features"], df_effective=round(fitinfo[arm]["df"], 3),
            n_nonzero=fitinfo[arm]["nnz"],
            cv_mse_selected=s["cv_mse"], cv_mse_min=s["cv_mse_min"],
            cv_se_at_min=s["cv_se_at_min"], alpha_at_min=s["alpha_min"],
            lambda_at_min=s["lam_min"],
            n_nonconverged_candidates=s["n_nonconverged_candidates"],
            all_candidates_nonconverged=s["all_candidates_nonconverged"],
            oof_r2_A=1.0 - s["oof_sse"] / ssA_tot, oof_sse_A=s["oof_sse"],
            sse_B=sseB, mse_B=sseB / ctx.nB, r2_B=1.0 - sseB / ssB_tot,
            gain_vs_BCOV_normalized=(cov_sse - sseB) / (ctx.nB * ctx.s2A),
            sse_B_cvmin=sseB_cv, r2_B_cvmin=1.0 - sseB_cv / ssB_tot,
            df_effective_cvmin=round(fitinfo[arm]["df_cvmin"], 3),
            gain_vs_BCOV_normalized_cvmin=(cov_sse - sseB_cv) / (ctx.nB * ctx.s2A),
            shrunk_to_covariate_only=degenerate))
    payload = dict(tile_id=tile_id, chrom=str(chrom), tile_start=int(tile_start),
                   info=info, sel=sel, fitinfo=fitinfo, arm_status=arm_status,
                   hapla_note=notes["hapla"], timing=dict((k, round(v, 1)) for k, v in timing.items()),
                   wall_s=round(time.time() - t_start, 1),
                   n_zex_cols=len(zex_names), n_zcl_cols=len(zcl_names),
                   pair_keys=[k for (_a, _b, k) in blocks["pairs"]],
                   var_keys_sha256=p2.sha256_text("|".join(blocks["var_keys"])),
                   nuis_keys_sha256=p2.sha256_text("|".join(blocks["nuis_keys"])),
                   state_dict_sha256=p2.sha256_text(
                       "|".join(",".join(d.states) for d in dicts)),
                   predB=dict((a, predB[a].astype(np.float32)) for a in predB),
                   cv_grid=cv_grid, s2A=ctx.s2A)
    return rows, payload, info

def claim(queue_dir, tile):
    lock = os.path.join(queue_dir, tile.replace(":", "_") + ".claim")
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except OSError as e:
        if e.errno == errno.EEXIST:
            return False
        raise
    os.write(fd, ("%d %s\n" % (os.getpid(), time.strftime("%FT%TZ", time.gmtime()))).encode())
    os.close(fd)
    return True

def main(argv=None):
    global HAPLA
    ap = argparse.ArgumentParser("hcf.p2_run")
    ap.add_argument("--run", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--tilelist", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--queue", required=True)
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--hapla", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--deadline-s", type=float, default=0.0)
    ap.add_argument("--order", default="file", choices=["file", "desc"])
    args = ap.parse_args(argv)
    HAPLA = args.hapla or None
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.queue, exist_ok=True)
    args.tmp = os.path.join(args.tmp, "w%02d" % args.worker)
    os.makedirs(args.tmp, exist_ok=True)

    ctx = Ctx(args)
    log("ctx n=%d A=%d B=%d s2A=%.5f folds=%s" %
        (ctx.n, ctx.nA, ctx.nB, ctx.s2A, [len(f) for f in ctx.folds]))
    tl = pd.read_csv(args.tilelist, sep="\t")
    if args.order == "desc":
        tl = tl.sort_values("n_primary", ascending=False).reset_index(drop=True)
    t_begin = time.time()
    ndone = 0
    for i in range(len(tl)):
        ch = str(tl["chrom"].iloc[i])
        ts = int(tl["tile_start"].iloc[i])
        tile_id = "chr%s:%d" % (ch, ts)
        if args.deadline_s and (time.time() - t_begin) > args.deadline_s:
            log("deadline reached before claiming %s, stopping" % tile_id)
            break
        if not claim(args.queue, tile_id):
            continue
        tmpd = os.path.join(args.tmp, tile_id.replace(":", "_"))
        os.makedirs(tmpd, exist_ok=True)
        try:
            rows, payload, info = run_tile(ctx, ch, ts, tmpd)
            if rows:
                pd.DataFrame(rows).to_csv(
                    os.path.join(args.outdir, "bc_%s.csv" % tile_id.replace(":", "_")),
                    index=False)
                np.savez_compressed(
                    os.path.join(args.outdir, "pred_%s.npz" % tile_id.replace(":", "_")),
                    **payload.pop("predB"))
                with open(os.path.join(args.outdir, "lock_%s.json" % tile_id.replace(":", "_")),
                          "w") as fh:
                    json.dump(payload, fh, ensure_ascii=False, default=str)
                log("%s OK wall=%.0fs %s" % (tile_id, payload["wall_s"], payload["timing"]))
            else:
                with open(os.path.join(args.outdir, "skip_%s.json" % tile_id.replace(":", "_")),
                          "w") as fh:
                    json.dump(info, fh, ensure_ascii=False, default=str)
                log("%s SKIP %s" % (tile_id, info.get("status")))
        except Exception as exc:
            import traceback
            with open(os.path.join(args.outdir, "fail_%s.json" % tile_id.replace(":", "_")),
                      "w") as fh:
                json.dump({"tile_id": tile_id, "error": "%s: %s" % (type(exc).__name__, exc),
                           "traceback": traceback.format_exc()[-4000:]}, fh, ensure_ascii=False)
            log("%s FAILED %s: %s" % (tile_id, type(exc).__name__, exc))
        finally:
            shutil.rmtree(tmpd, ignore_errors=True)
        ndone += 1
        if args.limit and ndone >= args.limit:
            break
    log("worker %d finished, %d tiles" % (args.worker, ndone))
    return 0

if __name__ == "__main__":
    sys.exit(main())
