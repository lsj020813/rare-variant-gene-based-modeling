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


import csv, glob, json, math, os, random, statistics as st
W=_config_path("${PROJECT_ROOT}/work"); D=f"{W}/run_l3b/out"
sv={}
for f in glob.glob(f"{D}/arm_none/part*.singleAssoc.txt"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        try:
            k=f'{r["CHR"]}:{r["POS"]}:{r["Allele1"]}:{r["Allele2"]}'
            af=float(r["AF_Allele2"]); maf=min(af,1-af); mac=float(r["AC_Allele2"]); mac=min(mac, 2*float(r["N"])-mac)
            sv[k]=(float(r["p.value"]), maf, mac)
        except (ValueError,KeyError): pass
print("singleAssoc:", len(sv))
genes={}; cur=None
for line in open(f"{D}/G_phi.txt"):
    t=line.split()
    if len(t)<3: continue
    g,kind=t[0],t[1]
    genes.setdefault(g,{})
    if kind=="var": genes[g]["var"]=t[2:]
    elif kind=="weight": genes[g]["phi"]=[float(x) for x in t[2:]]
print("genes:", len(genes)); sample_var=next(iter(genes.values()))["var"][:2]; print("var key sample:", sample_var)
def norm_key(v):
    parts=v.replace("-",":").replace("_",":").split(":")
    return ":".join(parts[:4])
skat={}
for r in csv.DictReader(open(f"{D}/none.merged", encoding="utf-8", errors="replace"), delimiter="\t"):
    if r.get("Group")=="all":
        try: skat[r["Region"]]=(float(r["Pvalue"]), float(r["Pvalue_SKAT"]))
        except (ValueError,KeyError): pass
def beta125(m): return 25.0*(1-m)**24
EPS=1e-15
def acatv(ps, ws):
    ws=[w for w in ws]; S=sum(ws)
    if S<=0: return float("nan")
    T=sum(w/S*math.tan((0.5-min(max(p,EPS),1-EPS))*math.pi) for p,w in zip(ps,ws))
    return 0.5-math.atan(T)/math.pi
def burden_p_placeholder(): return None
res={}; nmac_low=0; nmatch=0; ntot=0; missing_genes=0
rngs={s: random.Random(s) for s in (20260909,20260910,20260911)}
for g,d in genes.items():
    if "var" not in d or "phi" not in d or len(d["var"])!=len(d["phi"]): missing_genes+=1; continue
    rows=[]
    for v,ph in zip(d["var"], d["phi"]):
        k=norm_key(v); ntot+=1
        if k in sv:
            pj,maf,mac=sv[k]; nmatch+=1
            if mac<=10: nmac_low+=1; continue
            rows.append((pj, maf, ph))
    if len(rows)<2: continue
    ps=[r[0] for r in rows]; base=[beta125(r[1])**2*r[1]*(1-r[1]) for r in rows]; phi=[r[2] for r in rows]
    out={"n":len(rows),
         "flat": acatv(ps, base),
         "equal": acatv(ps, [1.0]*len(rows)),
         "phi": acatv(ps, [b*f for b,f in zip(base,phi)])}
    for s,rg in rngs.items():
        pp=phi[:]; rg.shuffle(pp); out[f"perm{s%100-8}"]=acatv(ps,[b*f for b,f in zip(base,pp)])
    res[g]=out
print(f"genes scored {len(res)} | variants matched {nmatch}/{ntot} | MAC<=10: {nmac_low} | skipped genes {missing_genes}")
arms=["flat","equal","phi","perm1","perm2","perm3"]
THR=[2.5e-6,1e-5,1e-4,1e-3]
print("  "+f"{'팔':8s}"+"".join(f"{f'<{t:.0e}':>9s}" for t in THR))
for a in arms: print("  "+f"{a:8s}"+"".join(f"{sum(1 for g in res if res[g][a]<t):9d}" for t in THR))
print("  "+f"{'SKATnone':8s}"+"".join(f"{sum(1 for g in res if g in skat and skat[g][1]<t):9d}" for t in THR))
print("  "+f"{'omni':8s}"+"".join(f"{sum(1 for g in res if g in skat and skat[g][0]<t):9d}" for t in THR))
def sign(x_get, y_get, label):
    d=[]; 
    for g in res:
        try: a=x_get(g); b=y_get(g)
        except KeyError: continue
        if a is None or b is None or math.isnan(a) or math.isnan(b): continue
        d.append(math.log10(b)-math.log10(a))
    win=sum(1 for v in d if v>1e-12); lose=sum(1 for v in d if v<-1e-12); tie=len(d)-win-lose; n=win+lose
    z=(win-n/2)/math.sqrt(n/4) if n else float("nan")
    print(f"  {label:26s} 중위차 {st.median(d):+.4f} | 승 {win} 패 {lose} 동률 {tie} | z {z:+.2f}")
    return {"median":st.median(d),"win":win,"lose":lose,"tie":tie,"z":z}
S={}
print("  === 전수 부호검정 (동률 제외; 양수/양 z = 앞쪽 우세)")
S["phi_vs_flat"]=sign(lambda g:res[g]["phi"], lambda g:res[g]["flat"], "ACATV_phi vs ACATV_flat")
for k in ("perm1","perm2","perm3"): S[f"phi_vs_{k}"]=sign(lambda g:res[g]["phi"], lambda g,k=k:res[g][k], f"ACATV_phi vs ACATV_{k}")
for k in ("perm1","perm2","perm3"): S[f"{k}_vs_flat"]=sign(lambda g,k=k:res[g][k], lambda g:res[g]["flat"], f"ACATV_{k} vs ACATV_flat")
S["flat_vs_skat"]=sign(lambda g:res[g]["flat"], lambda g:skat[g][1], "ACATV_flat vs SKAT_none")
S["flat_vs_omni"]=sign(lambda g:res[g]["flat"], lambda g:skat[g][0], "ACATV_flat vs omni_none")
S["equal_vs_flat"]=sign(lambda g:res[g]["equal"], lambda g:res[g]["flat"], "ACATV_equal vs ACATV_flat")
tracked=["ENSG00000129353.15","ENSG00000213892.12","ENSG00000186567.14","ENSG00000130202.10","ENSG00000130204.13","ENSG00000104856.15","ENSG00000069399.15","ENSG00000079805.19","ENSG00000142453.13","ENSG00000127616.22","ENSG00000129354.12"]
print("  === 추적 11개 (omni / SKAT / ACATV flat / phi / perm1-3)")
for g in tracked:
    if g in res: print(f"  {g:22s} omni={skat[g][0]:.2e} skat={skat[g][1]:.2e} | flat={res[g]['flat']:.2e} phi={res[g]['phi']:.2e} | " + " ".join(f"{res[g][k]:.2e}" for k in ("perm1","perm2","perm3")))
json.dump({"n_genes":len(res),"matched":nmatch,"total":ntot,"mac_le10":nmac_low,"sign_tests":S,
           "threshold_counts":{a:[sum(1 for g in res if res[g][a]<t) for t in THR] for a in arms},
           "per_gene":res}, open(f"{D}/acatv_axis.json","w"), indent=1)
print("ACATV_DONE")
