
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import json, glob, os, numpy as np, pandas as pd
OUTD=os.environ.get("OUTD",_config_path("${PROJECT_ROOT}/work/gate1/out/tables"))
os.makedirs(OUTD, exist_ok=True)

files = sorted(glob.glob(os.environ.get("CTXGLOB",_config_path("${PROJECT_ROOT}/work/gate1/out/ctx_chr*.json"))))
meta = None; genes = []
for fn in files:
    d = json.load(open(fn))
    meta = {k: d[k] for k in ("seed","var_cap","n_sub_individuals","radii_hamming","radii_quantile","noncoding_def")}
    genes += [g for g in d["genes"]]
ok = [g for g in genes if g.get("status") == "ok"]

inv = []
for g in ok:
    for u, d in g["universes"].items():
        if not d.get("n_var"): continue
        inv.append(dict(gene=g["gene"], chr=g["chr"], universe=u, n_var=d["n_var"],
                        nB1=g["n_by_bin_full"]["B1"], nB2=g["n_by_bin_full"]["B2"],
                        nB3=g["n_by_bin_full"]["B3"], nB4=g["n_by_bin_full"]["B4"],
                        capped=g["capped"], window_bp=g["window_bp"], n_intervals=g["n_intervals"],
                        src_dist_only=g["src_counts"]["dist_only"], src_re2g_only=g["src_counts"]["re2g_only"],
                        src_both=g["src_counts"]["both"],
                        median_maf=d["median_maf"], iqr_maf=d["iqr_maf"], median_r2=d["median_r2"],
                        median_dosage_var=d["median_dosage_var"],
                        median_hard_carrier=d["median_hard_carrier"],
                        median_expected_allele=d["median_expected_allele"],
                        M_eff=d.get("M_eff"), M_eff_over_M=d.get("M_eff_over_M"),
                        condition_number=d.get("condition_number"), eff_rank_entropy=d.get("eff_rank_entropy"),
                        n_blocks_r2_0_8=d.get("n_blocks_r2_0.8"), highLD_pair_frac=d.get("highLD_pair_frac"),
                        n_ccre_units=d.get("n_ccre_units"), n_var_in_ccre=d.get("n_var_in_ccre")))
inv = pd.DataFrame(inv); inv.to_csv(OUTD+"/gene_variant_inventory.csv", index=False)

kd = []
for g in ok:
    for u, d in g["universes"].items():
        if not d.get("n_var"): continue
        kd.append(dict(gene=g["gene"], chr=g["chr"], universe=u, n_var=d["n_var"],
                       P_K0=d["P_K0"], P_K1=d["P_K1"], P_Kge2=d["P_Kge2"], P_Kge3=d["P_Kge3"],
                       P_Kge2_given_ge1=d["P_Kge2_given_ge1"], P_Kge3_given_ge1=d["P_Kge3_given_ge1"],
                       median_K=d["median_K"], median_A=d["median_A"]))
kd = pd.DataFrame(kd); kd.to_csv(OUTD+"/individual_gene_K_distribution.csv", index=False)

co = []
for g in ok:
    for u, d in g["universes"].items():
        for pair, s in (d.get("cocarrier") or {}).items():
            co.append(dict(gene=g["gene"], chr=g["chr"], universe=u, bin_pair=pair, **s))
co = pd.DataFrame(co); co.to_csv(OUTD+"/frequency_cocarrier_matrix.csv", index=False)

rc = []
for g in ok:
    for u, d in g["universes"].items():
        for defn in ("ctxA","ctxB","ctxC","ctxD"):
            s = d.get(defn)
            if not s: continue
            rc.append(dict(gene=g["gene"], chr=g["chr"], universe=u, definition=defn, **s))
rc = pd.DataFrame(rc); rc.to_csv(OUTD+"/context_recurrence.csv", index=False)

nb = []
for g in ok:
    for u, d in g["universes"].items():
        for key, defn in (("nbhd_A","A_exact"),("nbhd_C","C_ccre_burden"),("nbhd_D","D_freq_aware")):
            s = d.get(key)
            if not s: continue
            nb.append(dict(gene=g["gene"], chr=g["chr"], universe=u, definition=defn, **s))
nb = pd.DataFrame(nb); nb.to_csv(OUTD+"/context_neighborhood.csv", index=False)

rows = []
for u in ("U1","U2","U3"):
    k = kd[kd.universe == u]; i = inv[inv.universe == u]
    r = rc[rc.universe == u]
    row = dict(universe=u, n_genes=int(k.gene.nunique()),
               median_n_var=float(i.n_var.median()) if len(i) else np.nan,
               median_P_Kge2_given_ge1=float(k.P_Kge2_given_ge1.median()) if len(k) else np.nan,
               median_P_Kge3_given_ge1=float(k.P_Kge3_given_ge1.median()) if len(k) else np.nan,
               frac_genes_P_Kge2_given_ge1_ge0_2=float((k.P_Kge2_given_ge1 >= 0.2).mean()) if len(k) else np.nan,
               median_M_eff_over_M=float(i.M_eff_over_M.median()) if len(i) else np.nan)
    for defn in ("ctxA","ctxB","ctxC","ctxD"):
        rr = r[r.definition == defn]
        row[f"{defn}_median_singleton_frac"] = float(rr.singleton_frac.median()) if len(rr) else np.nan
        row[f"{defn}_median_eff_ctx"] = float(rr.eff_ctx_count.median()) if len(rr) else np.nan
        row[f"{defn}_median_frac_ind_ctx_ge50"] = float(rr.frac_ind_in_ctx_ge50.median()) if len(rr) else np.nan
    cu = co[co.universe == u]
    row["cocarrier_pairs_ge50"] = int(cu.n_ge50.sum()) if len(cu) else 0
    row["cocarrier_pairs_ge50_lowLD"] = int(cu.n_ge50_lowLD.sum()) if len(cu) else 0
    rows.append(row)
fx = pd.DataFrame(rows); fx.to_csv(OUTD+"/frequency_expansion_curve.csv", index=False)
print(fx.to_string(index=False))

print("meta", json.dumps(meta))
print("genes_ok", len(ok), "genes_total", len(genes))
from collections import Counter
print("status", json.dumps(dict(Counter(g.get("status") for g in genes))))
