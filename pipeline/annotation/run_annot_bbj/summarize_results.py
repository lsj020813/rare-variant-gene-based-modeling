#!/usr/bin/env python3
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import json, glob, os, sys
E6 = _config_path("${PROJECT_ROOT}/work/run_band15/model_v10_out/e6_runs")
rows = []; lad = []; meta = []
def rnd(x, k):
    return round(x, k) if isinstance(x, (int, float)) else x
for d in sorted(glob.glob(f"{E6}/*")):
    name = os.path.basename(d)
    p = f"{d}/RESULT.json"; partial = False
    if not os.path.exists(p):
        p = f"{d}/RESULT.partial.json"; partial = True
    if not os.path.exists(p):
        meta.append(dict(run=name, status="no result", rc=open(f"{d}/rc.txt").read().strip() if os.path.exists(f"{d}/rc.txt") else "running"))
        continue
    r = json.load(open(p))
    fc = r.get("feature_contract", {})
    meta.append(dict(run=name, status=("partial" if partial else "complete"), fold=r.get("fold"), train_bands=r.get("train_bands"),
                     phi_columns=fc.get("phi_columns"), n_selected=len(fc.get("selected", [])), phi_lambda=r.get("final_fit", {}).get("phi_lambda"),
                     phi_lambda_key=r.get("selection", {}).get("phi", {}).get("chosen"), C2_lambda_key=r.get("selection", {}).get("C2", {}).get("chosen"),
                     split_counts=r.get("split_counts"), verdict=(r.get("verdict") or {}).get("category"), verdict_reason=(r.get("verdict") or {}).get("reason"),
                     labelled_pairs=r.get("join_audit", {}).get("labelled_pairs"), annotation_rows=r.get("join_audit", {}).get("annotation_rows"),
                     traits=len(r.get("join_audit", {}).get("joins", [])), lc_points=len((r.get("learning_curves") or {}).get("rows", [])),
                     sim="yes" if r.get("simulation") else "no"))
    for axis, ev in r.get("evaluation", {}).items():
        if ev.get("status") == "not_requested": continue
        for band, pb in ev.get("bands", {}).items():
            if pb.get("status") != "ok":
                rows.append(dict(run=name, axis=axis, band=band, n=pb.get("n"), model="-", status=pb.get("status"))); continue
            for m, v in pb["models"].items():
                if m.startswith("C1:"): continue
                w = v["weighted_ce"]; d2 = v.get("delta_ce_c2") or {}; sp = v.get("cs_spearman") or {}
                rows.append(dict(run=name, axis=axis, band=band, n=pb["n"], model=m, ce=rnd(w["observed"], 5), ce_null_mean=rnd(w.get("null_mean", float("nan")), 5), ce_z=rnd(w.get("z", float("nan")), 2),
                                 dCE_vs_C2=rnd(d2.get("observed", float("nan")), 5), dCE_z=rnd(d2.get("z", float("nan")), 2), dCE_p=d2.get("p_greater"),
                                 spearman=rnd(sp.get("observed", float("nan")), 4) if isinstance(sp.get("observed"), (int, float)) else sp.get("observed"), spearman_z=rnd(sp.get("z", float("nan")), 2) if isinstance(sp.get("z"), (int, float)) else None,
                                 cs_defined=v.get("cs_defined")))
            for step, s in pb["ladder"].items():
                lad.append(dict(run=name, axis=axis, band=band, step=step, observed=rnd(s.get("observed", float("nan")), 5), null_mean=rnd(s.get("null_mean", float("nan")), 5), null_sd=rnd(s.get("null_sd", float("nan")), 5), z=rnd(s.get("z", float("nan")), 2), p_greater=s.get("p_greater"), status=s.get("status")))
def dump(name, rs):
    if not rs: return
    keys = []
    for r in rs:
        for k in r:
            if k not in keys: keys.append(k)
    with open(f"{E6}/{name}.tsv", "w") as o:
        o.write("\t".join(keys) + "\n")
        for r in rs: o.write("\t".join(json.dumps(r.get(k), ensure_ascii=False) if isinstance(r.get(k), (dict, list)) else str(r.get(k)) for k in keys) + "\n")
dump("summary_meta", meta); dump("summary_models", rows); dump("summary_ladder", lad)
print(json.dumps(meta, ensure_ascii=False, indent=0))
