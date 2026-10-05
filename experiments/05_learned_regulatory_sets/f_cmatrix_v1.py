#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re
import math as _cfg_math

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)

def _config_number(name, cast=int, positive=False):
    value = cast(_cfg_os.environ[name])
    if not _cfg_math.isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value

import os, sys, json, glob, gzip, time
import numpy as np, pandas as pd
ROOT = _config_path("${PROJECT_ROOT}/work"); P = f"{ROOT}/fset/primary"; Q = f"{P}/qc"; os.makedirs(Q, exist_ok=True)
t0 = time.time()
S = pd.read_csv(f"{P}/primary_sample.tsv", sep="\t", dtype={"chr": str, "key": str, "domain": str})
S["ub_member"] = ((S.has_ccre > 0) | (S.n_tf > 0) | (S.n_re2g > 0)).astype(int)
N = len(S); assert N == _config_number("N_MATRIX_ROWS", int, True), N
feat = pd.concat([pd.read_csv(f, sep="\t", dtype={"chr": str, "key": str, "domain": str, "ccre_class": str, "rep_class": str},
                              keep_default_na=False, na_values=[""], quoting=3) for f in sorted(glob.glob(f"{P}/feat_chr*.tsv.gz"))])
NK = S.key.nunique(); assert len(feat) == NK and feat.key.is_unique, (len(feat), NK)
feat = feat.set_index("key").loc[S.key.values].reset_index()
assert len(feat) == N
S["key_multiplicity"] = S.groupby("key")["key"].transform("size").values
assert (feat.pos38.values == S.pos38.values).all()
g = pd.read_csv(f"{P}/gpn_scores.tsv.gz", sep="\t", dtype={"key": str, "chr": str}).drop_duplicates("key").set_index("key").loc[S.key.values]
feat["gpn_score"] = g.gpn_score.values; gpn_match = g.gpn_match.values
del g
vre = pd.concat([pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False, quoting=3) for f in sorted(glob.glob(f"{P}/vset_re2g_chr*.tsv.gz"))]).drop_duplicates("key").set_index("key").loc[S.key.values]
vst = pd.concat([pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False, quoting=3) for f in sorted(glob.glob(f"{P}/vset_chr*.tsv.gz"))])
vst["kd"] = vst.key + "|" + vst.domain; vst = vst.drop_duplicates("kd").set_index("kd").loc[(S.key + "|" + S.domain).values]
assert len(vre) == N and len(vst) == N
print("loaded", round(time.time() - t0), "s", flush=True)

CCRE_LV = ["PLS", "pELS", "dELS", "CA-CTCF", "CA-H3K4me3", "CA-TF", "CA", "TF", "none"]
REP_MAIN = ["SINE", "LINE", "LTR", "DNA", "Simple_repeat", "Low_complexity", "Satellite", "Retroposon"]
REP_LV = REP_MAIN + ["other", "none"]
def rep_collapse(v):
    v = str(v).rstrip("?")
    return v if v in REP_MAIN else ("none" if v == "none" else "other")
feat["rep_class_c"] = feat.rep_class.map(rep_collapse)
TFC = [c for c in feat.columns if c.startswith("tf_") and c != "tf_n"]; assert len(TFC) == 60
rand701 = (S.forced.values == 0)
masks = {"ub_inub": S.in_ub.values == 1, "ub_member": S.ub_member.values == 1, "ua": np.ones(N, bool)}

def robust_z(x, fit_mask):
    v = x[fit_mask]; v = v[~np.isnan(v)]
    med = float(np.median(v)); iqr = float(np.subtract(*np.percentile(v, [75, 25])))
    if iqr > 0: sc, rule = iqr / 1.349, "IQR/1.349"
    else:
        mad = float(np.median(np.abs(v - med))) * 1.4826
        if mad > 0: sc, rule = mad, "MAD*1.4826"
        else:
            sd = float(np.std(v)); sc, rule = (sd, "std") if sd > 0 else (1.0, "unit(zero-var)")
    return (x - med) / sc, {"median": med, "scale": sc, "rule": rule}

CONT = {"ccre_dist": "log10(d+1)", "rep_dist": "log10(d+1)", "cpg_dist": "log10(d+1)", "tf_n": "log10(n+1)", "map_k36": "identity", "gpn_score": "identity"}
def build(fit_mask):
    cols, mats, params = [], [], {}
    for lv in CCRE_LV: cols.append(f"ccre_{lv}"); mats.append((feat.ccre_class.values == lv).astype(np.float32))
    for lv in REP_LV: cols.append(f"rep_{lv}"); mats.append((feat.rep_class_c.values == lv).astype(np.float32))
    cols.append("cpg"); mats.append(feat.cpg.values.astype(np.float32))
    for c in TFC: cols.append(c); mats.append(feat[c].values.astype(np.float32))
    for c, tr in CONT.items():
        x = feat[c].values.astype(float)
        if tr.startswith("log10"): x = np.log10(x + 1)
        z, prm = robust_z(x, fit_mask); prm["transform"] = tr; params[c] = prm
        miss = np.isnan(z); z = np.where(miss, 0.0, z)
        cols.append(f"{c}_rz"); mats.append(z.astype(np.float32))
        if miss[fit_mask].any(): cols.append(f"{c}_missing"); mats.append(miss.astype(np.float32))
    X = np.column_stack(mats).astype(np.float32)
    return X, cols, params

def pca_summary(X):
    Xc = X.astype(np.float64); Xc = Xc - Xc.mean(0)
    C = (Xc.T @ Xc) / (len(Xc) - 1)
    ev = np.linalg.eigvalsh(C)[::-1]; ev = np.clip(ev, 0, None)
    tot = ev.sum(); p = ev / tot; pn = p[p > 0]
    er = float(np.exp(-(pn * np.log(pn)).sum()))
    tol = ev[0] * 1e-8; nz = int((ev <= tol).sum()); evp = ev[ev > tol]
    cum = np.cumsum(p)
    pr = float(tot ** 2 / (ev ** 2).sum())
    return {"eig": ev, "evr": p, "cum": cum, "eff_rank": er, "participation_ratio": pr, "cond": float(evp[0] / evp[-1]),
            "n_zero_eig": nz, "n90": int(np.searchsorted(cum, 0.90) + 1), "n95": int(np.searchsorted(cum, 0.95) + 1),
            "pc1_evr": float(p[0])}

def nzv(v):
    v = v[~pd.isna(v)]
    if len(v) == 0: return True, np.nan, 0
    vc = pd.Series(v).value_counts()
    fr = vc.iloc[0] / vc.iloc[1] if len(vc) > 1 else np.inf
    pu = 100 * len(vc) / len(v)
    return bool(fr > 19 and pu < 10), float(fr), int(len(vc))

MISS_MEANING = {
    "ccre_class": "none = biological absence (no cCRE overlap); never NA", "ccre_dist": "0 when overlapping; NA only if no cCRE on chromosome (technical)",
    "rep_class": "none = biological absence (no RepeatMasker element); never NA", "rep_dist": "0 when overlapping; NA only if no element on chromosome (technical)",
    "cpg": "0 = biological absence", "cpg_dist": "0 when overlapping; NA only if no island on chromosome (technical)",
    "map_k36": "NA = bigWig has no value at position (technical: unmappable/gap)", "tf_n": "0 = biological absence (no ReMap peak)",
    "gpn_score": "NA = no GPN-MSA record (alignment gap) or hg19->hg38 allele mismatch (technical)",
}
qc_rows = []
def add_row(name, axis, source, ftype, vals, scaling, meaning, raw_of=""):
    vals = np.asarray(vals, dtype=object) if ftype == "categorical" else np.asarray(vals, dtype=float)
    def stats(m):
        v = vals[m]; miss = pd.isna(v).mean() if ftype != "categorical" else 0.0
        z, fr, nu = nzv(v)
        zv = (nu <= 1)
        return miss, zv, z, fr, nu
    r = {"feature": name, "axis": axis, "source": source, "type": ftype, "scaling_rule": scaling, "missing_meaning": meaning, "raw_parent": raw_of}
    for tag, m in [("all1000", np.ones(N, bool)), ("rand701", rand701)]:
        miss, zv, z, fr, nu = stats(m)
        r.update({f"missing_rate_{tag}": round(float(miss), 6), f"zero_variance_{tag}": zv, f"near_zero_var_{tag}": z,
                  f"freq_ratio_{tag}": (None if not np.isfinite(fr) else round(fr, 2)), f"n_unique_{tag}": nu})
    for tag, m in masks.items():
        r[f"n_rows_{tag}"] = int(m.sum())
        r[f"n_nonmissing_{tag}"] = int((~pd.isna(vals[m])).sum()) if ftype != "categorical" else int(m.sum())
        r[f"n_rows_{tag}_rand701"] = int((m & rand701).sum())
    qc_rows.append(r)

SRC = {"ccre": "ENCODE cCRE (ref/b6_cards/ccre.s.bed)", "rep": "UCSC RepeatMasker rmsk", "cpg": "UCSC cpgIslandExt", "map": "Umap k36 MultiTrackMappability",
       "tf": "ReMap2022 NR (top-60 TF by peak count, fset/out/remap_top_tfs.txt)", "gpn": "GPN-MSA scores.tsv.bgz (tabix -R, allele matched)"}
add_row("ccre_class", "C-set", SRC["ccre"], "categorical", feat.ccre_class.values, f"one-hot {len(CCRE_LV)} levels", MISS_MEANING["ccre_class"])
add_row("ccre_dist", "C-set", SRC["ccre"], "continuous", feat.ccre_dist.values, "log10(d+1) -> robust z", MISS_MEANING["ccre_dist"])
add_row("rep_class", "C-set", SRC["rep"], "categorical", feat.rep_class.values, f"collapse '?'-suffix, minor classes->other; one-hot {len(REP_LV)} levels", MISS_MEANING["rep_class"])
add_row("rep_dist", "C-set", SRC["rep"], "continuous", feat.rep_dist.values, "log10(d+1) -> robust z", MISS_MEANING["rep_dist"])
add_row("cpg", "C-set", SRC["cpg"], "binary", feat.cpg.values, "0/1 as is", MISS_MEANING["cpg"])
add_row("cpg_dist", "C-set", SRC["cpg"], "continuous", feat.cpg_dist.values, "log10(d+1) -> robust z", MISS_MEANING["cpg_dist"])
add_row("map_k36", "C-set", SRC["map"], "continuous", feat.map_k36.values, "robust z (+missing indicator)", MISS_MEANING["map_k36"])
add_row("tf_n", "C-set", SRC["tf"], "continuous(count)", feat.tf_n.values, "log10(n+1) -> robust z", MISS_MEANING["tf_n"])
for c in TFC: add_row(c, "C-set", SRC["tf"], "binary", feat[c].values, "0/1 as is", "0 = biological absence (no peak of this TF)")
add_row("gpn_score", "C-set", SRC["gpn"], "continuous", feat.gpn_score.values, "robust z (+missing indicator)", MISS_MEANING["gpn_score"])
add_row("maf", "diagnostic-only(not a feature)", "primary_sample.tsv", "continuous", S.maf.values, "NOT scaled, NOT in any matrix", "n/a")
def tofloat(a): return pd.to_numeric(pd.Series(a).replace("", np.nan), errors="coerce").values
bios = pd.read_csv(f"{P}/vset_re2g_biosamples.tsv", sep="\t").set_index("index").biosample.to_dict()
for i in range(10):
    src_i = "ENCODE rE2G %s (ref/re2g, last bed field = Score)" % bios.get(i, i)
    add_row(f"re2g{i}_score", "V-set", src_i, "continuous", tofloat(vre[f"re2g{i}_score"].values), "held out (F3 only, within rE2G-overlapping variants, A1-3iii)", "NA = no rE2G element of this biosample overlaps (biological absence for this biosample)")
    add_row(f"re2g{i}_target", "V-set", src_i, "categorical", np.where(vre[f"re2g{i}_target"].values == "", "none", vre[f"re2g{i}_target"].values), "held out", "none = no element overlap")
    add_row(f"re2g{i}_dist_tss", "V-set", src_i, "continuous", tofloat(vre[f"re2g{i}_dist_tss"].values), "held out", "NA = no element overlap")
re2g_any = np.zeros(N, bool)
for i in range(10): re2g_any |= ~np.isnan(tofloat(vre[f"re2g{i}_score"].values))
add_row("re2g_any_biosample", "V-set(derived)", "ENCODE rE2G 10 biosamples", "binary", re2g_any.astype(float), "held out; A1-3(ii): NOT a validation metric (U-b eligibility axis)", "0 = no element in any biosample")
add_row("n_re2g_sample_col", "diagnostic(eligibility axis)", "primary_sample.tsv n_re2g (f_annflag)", "continuous(count)", S.n_re2g.values.astype(float), "not a feature", "0 = no element")
for c in ["gn_lof_oe", "gn_loeuf", "gn_pli", "gn_lof_z", "gn_mis_z", "gn_syn_z"]:
    add_row(c, "V-set", "gnomAD v4.1 constraint (gene-level, joined on domain ENSG)", "continuous(domain-level)", tofloat(vst[c].values), "held out (F3 only)", "NA = domain gene absent from constraint table (technical)")
inband = (S.maf.values >= 0.001) & (S.maf.values < 0.01)
for c in ["st_cons", "st_epi_active", "st_epi_repr", "st_epi_trans", "st_tf", "st_linsight"]:
    add_row(c, "V-set", "band-limited STAAR annotation (ref/annot/extract, MAF 0.1-1% band only)", "continuous", tofloat(vst[c].values), "held out (F3 only, band-restricted)", "NA outside 0.1-1% MAF band = technical (band extraction); NA inside band = source NA")
st_rows_present = vst["st_cons"].values != ""
for c, lab in [("st_cage_prom", "CAGE promoter id (text)"), ("st_genehancer", "GeneHancer record (text; Name=score)")]:
    pres = (vst[c].values != "").astype(float); pres[~st_rows_present] = np.nan
    add_row(c + "_present", "V-set", "band-limited STAAR annotation, %s" % lab, "binary(sparse text->presence)", pres, "held out (F3 only, band-restricted)", "NA = row absent from band annotation (technical); 0 = no element (biological absence); 1 = element present")
gh = pd.Series(vst["st_genehancer"].values).str.extract(r"Name=([0-9.]+)")[0].astype(float).values
gh[np.isnan(gh) & st_rows_present] = 0.0
add_row("st_genehancer_score", "V-set(derived)", "GeneHancer Name= score parsed from st_genehancer", "continuous", gh, "held out", "NA = row absent from band annotation; 0 = no GeneHancer element")
for c in ["cadd_raw", "cadd_phred"]:
    add_row(c, "V-set", "CADD v1.6 (ref/features, MAF 0.1-1% band only)", "continuous", tofloat(vst[c].values), "held out (F3 only, band-restricted)", "NA outside 0.1-1% MAF band = technical (band extraction)")
qc = pd.DataFrame(qc_rows)
for c in [c for c in qc.feature if (c.startswith("st_") or c.startswith("cadd_")) and c in vst.columns and c not in ("st_cage_prom", "st_genehancer")]:
    v = tofloat(vst[c].values)
    qc.loc[qc.feature == c, "missing_rate_inband_maf0.1-1pct"] = round(float(np.isnan(v[inband]).mean()), 6)
    qc.loc[qc.feature == c, "n_inband_rows"] = int(inband.sum())
qc["gpn_match_counts"] = ""
qc.loc[qc.feature == "gpn_score", "gpn_match_counts"] = json.dumps(pd.Series(gpn_match).value_counts().to_dict())

built = {}
FN = {"ub_member": "Cset_ub", "ub_inub": "Cset_ub_inub", "ua": "Cset_ua"}
for name in ["ub_member", "ub_inub", "ua"]:
    m = masks[name]; X, cols, prm = build(m); built[name] = (X[m], cols, prm)
    print(name, X[m].shape, flush=True)
proc_rows = []
for name in ["ub_member", "ub_inub", "ua"]:
    Xm, cols, prm = built[name]; m = masks[name]; r7 = rand701[m]
    for j, cname in enumerate(cols):
        col = Xm[:, j]
        for tag, mm in [("all1000", np.ones(len(col), bool)), ("rand701", r7)]:
            v = col[mm]; nu = len(np.unique(v)); z, fr, _ = nzv(v.astype(float))
            proc_rows.append({"matrix": FN[name], "column": cname, "subset": tag, "n_rows": int(mm.sum()), "mean": float(v.mean()), "sd": float(v.std()),
                              "n_unique": nu, "zero_variance": nu <= 1, "near_zero_var": z, "freq_ratio": (None if not np.isfinite(fr) else round(fr, 2)),
                              "minority_frac_if_binary": (float(min(v.mean(), 1 - v.mean())) if nu <= 2 else None)})
proc = pd.DataFrame(proc_rows)
ub_proc = proc[(proc.matrix == "Cset_ub")]
for cname in built["ub_member"][1]:
    a = ub_proc[(ub_proc.column == cname) & (ub_proc.subset == "all1000")].iloc[0]
    b = ub_proc[(ub_proc.column == cname) & (ub_proc.subset == "rand701")].iloc[0]
    parent = cname.split("_rz")[0].split("_missing")[0] if ("_rz" in cname or "_missing" in cname) else ("ccre_class" if cname.startswith("ccre_") else "rep_class" if cname.startswith("rep_") else cname)
    qc_rows.append({"feature": cname, "axis": "C-set(processed column)", "source": "derived", "type": ("binary" if a.n_unique <= 2 else "continuous(scaled)"),
                    "scaling_rule": "see raw parent", "missing_meaning": "imputed 0 (=median) with indicator" if "_rz" in cname else "", "raw_parent": parent,
                    "missing_rate_all1000": 0.0, "zero_variance_all1000": a.zero_variance, "near_zero_var_all1000": a.near_zero_var, "freq_ratio_all1000": a.freq_ratio, "n_unique_all1000": a.n_unique,
                    "missing_rate_rand701": 0.0, "zero_variance_rand701": b.zero_variance, "near_zero_var_rand701": b.near_zero_var, "freq_ratio_rand701": b.freq_ratio, "n_unique_rand701": b.n_unique,
                    "n_rows_ub_member": int(masks["ub_member"].sum()), "n_rows_ub_inub": int(masks["ub_inub"].sum()), "n_rows_ua": N,
                    "minority_frac_ub_member": a.minority_frac_if_binary})
qc = pd.DataFrame(qc_rows)
for c in [c for c in qc.feature if (c.startswith("st_") or c.startswith("cadd_")) and c in vst.columns and c not in ("st_cage_prom", "st_genehancer")]:
    v = tofloat(vst[c].values); qc.loc[qc.feature == c, "missing_rate_inband_maf0.1-1pct"] = round(float(np.isnan(v[inband]).mean()), 6); qc.loc[qc.feature == c, "n_inband_rows"] = int(inband.sum())
qc["gpn_match_counts"] = ""; qc.loc[qc.feature == "gpn_score", "gpn_match_counts"] = json.dumps(pd.Series(gpn_match).value_counts().to_dict())
qc["row_unit"] = "variant x domain instance; features are per-variant"
qc["ub_member_definition"] = "has_ccre>0 | n_tf>0 | n_re2g>0 (A2); ub_inub = in_ub==1 (A2 primary analysis rows)"
qc.to_csv(f"{Q}/regulatory_feature_matrix_qc.csv", index=False)
proc.to_csv(f"{Q}/processed_column_qc_long.csv", index=False)

red_rows = []
def redundancy(name, Xm, cols):
    sd = Xm.std(0); keep = sd > 0
    R = np.corrcoef(Xm[:, keep].astype(np.float64), rowvar=False); cn = [c for c, k in zip(cols, keep) if k]
    R = np.nan_to_num(R); np.fill_diagonal(R, 0)
    pairs = [(cn[i], cn[j], float(R[i, j])) for i in range(len(cn)) for j in range(i + 1, len(cn)) if abs(R[i, j]) > 0.8]
    parent = {c: c for c in cn}
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a, b, _ in pairs: parent[find(a)] = find(b)
    groups = {}
    for cc in cn: groups.setdefault(find(cc), []).append(cc)
    groups = {k: v for k, v in groups.items() if len(v) > 1}
    rep_of = {}
    for gid, (k, mem) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1]))):
        idx = [cn.index(m) for m in mem]; sub = np.abs(R[np.ix_(idx, idx)]).sum(1)
        rep = mem[int(np.argmax(sub))]
        for m in mem: rep_of[m] = (gid + 1, rep, len(mem))
    for a, b, r in sorted(pairs, key=lambda t: -abs(t[2])):
        red_rows.append({"matrix": name, "section": "pair_abs_r_gt_0.8", "feature_a": a, "feature_b": b, "r": round(r, 4), "group_id": rep_of[a][0]})
    for m, (gid, rep, sz) in sorted(rep_of.items(), key=lambda kv: (kv[1][0], kv[0])):
        red_rows.append({"matrix": name, "section": "group_membership", "feature_a": m, "group_id": gid, "group_size": sz, "representative": rep, "is_representative": m == rep})
    keep_all = list(cols); rep_only = [c for c in cols if (c not in rep_of) or rep_of[c][1] == c]
    zero_var = [c for c, k in zip(cols, keep) if not k]
    for c in keep_all: red_rows.append({"matrix": name, "section": "feature_list_keep_all_groups", "feature_a": c, "zero_variance": c in zero_var})
    for c in rep_only: red_rows.append({"matrix": name, "section": "feature_list_representative_only", "feature_a": c})
    ps_all = pca_summary(Xm); jr = [cols.index(c) for c in rep_only]; ps_rep = pca_summary(Xm[:, jr])
    for lab, ps, nf in [("keep_all_groups", ps_all, len(cols)), ("representative_only", ps_rep, len(rep_only))]:
        red_rows.append({"matrix": name, "section": "pca_summary", "condition": lab, "n_features": nf, "n_rows": len(Xm), "effective_rank_exp_entropy": round(ps["eff_rank"], 4),
                         "participation_ratio": round(ps["participation_ratio"], 4), "condition_number_nonzero_eig": round(ps["cond"], 2), "n_numerically_zero_eig": ps["n_zero_eig"],
                         "n_comp_90pct": ps["n90"], "n_comp_95pct": ps["n95"], "pc1_evr": round(ps["pc1_evr"], 4), "max_abs_r": round(max(abs(r) for _, _, r in pairs), 4) if pairs else None,
                         "n_pairs_abs_r_gt_0.8": len(pairs), "n_groups": len(groups)})
    for i, (e, p, cm) in enumerate(zip(ps_all["eig"], ps_all["evr"], ps_all["cum"])):
        red_rows.append({"matrix": name, "section": "pca_full_spectrum_keep_all", "component": i + 1, "eigenvalue": float(e), "explained_var_ratio": float(p), "cumulative": float(cm)})
    return rep_only
rep_only_ub = redundancy("Cset_ub", *built["ub_member"][:2])
redundancy("Cset_ub_inub", *built["ub_inub"][:2])
redundancy("Cset_ua", *built["ua"][:2])
pd.DataFrame(red_rows).to_csv(f"{Q}/feature_redundancy.csv", index=False)

er_rows = []
for name, m in masks.items():
    for tag, mm in [("all1000", m), ("rand701", m & rand701)]:
        X, cols, prm = build(mm); Xs = X[mm]
        for cond, cc in [("keep_all_groups", cols), ("representative_only", [c for c in cols if c in rep_only_ub])]:
            ps = pca_summary(Xs[:, [cols.index(c) for c in cc]])
            er_rows.append({"universe": name, "domains": tag, "condition": cond, "row_type": "summary", "n_rows": int(mm.sum()), "n_features": len(cc),
                            "effective_rank_exp_entropy": round(ps["eff_rank"], 4), "participation_ratio": round(ps["participation_ratio"], 4),
                            "condition_number_nonzero_eig": round(ps["cond"], 2), "n_numerically_zero_eig": ps["n_zero_eig"], "n_comp_90pct": ps["n90"], "n_comp_95pct": ps["n95"], "pc1_evr": round(ps["pc1_evr"], 4)})
            if cond == "keep_all_groups":
                for i, (e, p, cm) in enumerate(zip(ps["eig"], ps["evr"], ps["cum"])):
                    er_rows.append({"universe": name, "domains": tag, "condition": cond, "row_type": "component", "component": i + 1, "eigenvalue": float(e), "explained_var_ratio": float(p), "cumulative": float(cm)})
        print("eff rank", name, tag, flush=True)
pd.DataFrame(er_rows).to_csv(f"{Q}/feature_effective_rank.csv", index=False)

meta = {"n_rows_total": N, "created": time.strftime("%Y-%m-%dT%H:%M:%S"), "matrices": {}}
for name in ["ub_member", "ub_inub", "ua"]:
    Xm, cols, prm = built[name]; m = masks[name]; fn = FN[name]
    np.savez_compressed(f"{P}/{fn}.npz", X=Xm, columns=np.array(cols), key=S.key.values[m].astype(str), domain=S.domain.values[m].astype(str), chr=S.chr.values[m].astype(str),
                        pos38=S.pos38.values[m], in_ub=S.in_ub.values[m], key_multiplicity=S.key_multiplicity.values[m], ub_member=S.ub_member.values[m], forced=S.forced.values[m], diag_maf=S.maf.values[m].astype(np.float32), diag_r2=S.r2.values[m].astype(np.float32))
    with gzip.open(f"{P}/{fn}.keys.txt.gz", "wt") as o:
        for k in S.key.values[m]: o.write(k + "\n")
    json.dump({"columns": cols, "scaling": prm, "row_mask": name, "n_rows": int(m.sum())}, open(f"{P}/{fn}.scaling.json", "w"), indent=1)
    meta["matrices"][fn] = {"rows": int(m.sum()), "cols": len(cols), "row_rule": {"ub_member": "has_ccre>0 | n_tf>0 | n_re2g>0 (A2 membership formula; PRIMARY per task rule; NOTE: in capped domains exceeds the S5 400/domain cap because in_ub=0 extras that are U-b members are included)", "ub_inub": "in_ub==1 (S5-capped seeded U-b sampling path; all rows are ub_member)", "ua": "all rows (U-a parallel; over-samples U-b: ub_member frac %.3f)" % S.ub_member.mean()}[name]}
meta["n_distinct_variants_total"] = int(NK); meta["n_variant_x_domain_duplicate_rows"] = int(N - NK)
for name in ["ub_member", "ub_inub", "ua"]:
    fn = FN[name]; meta["matrices"][fn]["distinct_variants"] = int(S.key[masks[name]].nunique()); meta["matrices"][fn]["max_rows_per_domain"] = int(S[masks[name]].groupby("domain").size().max())
meta["ub_member_rows"] = int(S.ub_member.sum()); meta["rand701_rows"] = int(rand701.sum()); meta["seconds"] = round(time.time() - t0, 1)
json.dump(meta, open(f"{Q}/cset_meta.json", "w"), indent=1)
open(f"{P}/cset.done", "w").write(json.dumps(meta) + "\n")
print(json.dumps(meta))
