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


import gzip, zipfile, io, csv, json, statistics as st, sys
W=_config_path("${PROJECT_ROOT}/work")
ours={}
for line in open(f"{W}/run_l3b/out/G_phi.txt"):
    t=line.split()
    if len(t)>=3 and t[1]=="var":
        for v in t[2:]: ours.setdefault(v, t[0])
tracked={"ENSG00000129353.15","ENSG00000213892.12","ENSG00000186567.14","ENSG00000130202.10","ENSG00000130204.13",
         "ENSG00000104856.15","ENSG00000069399.15","ENSG00000079805.19","ENSG00000142453.13","ENSG00000127616.22","ENSG00000129354.12"}
print(f"our chr19 band variants: {len(ours)} | tracked-gene variants: {sum(1 for k,g in ours.items() if g in tracked)}", flush=True)
Z=f"{W}/ref/bbj/hum0197.v3.BBJ.TC.v1.zip"
zf=zipfile.ZipFile(Z); mem=[n for n in zf.namelist() if n.endswith(".gz")][0]
print("member:", mem, flush=True)
hit={}; n_lines=0; n_chr19=0
with zf.open(mem) as raw, gzip.open(raw, "rt") as fh:
    rd=csv.DictReader(fh, delimiter="\t")
    for r in rd:
        n_lines+=1
        if r["CHR"]!="19": continue
        n_chr19+=1
        k=f'{r["CHR"]}:{r["POS"]}:{r["Allele1"]}:{r["Allele2"]}'
        k2=f'{r["CHR"]}:{r["POS"]}:{r["Allele2"]}:{r["Allele1"]}'
        for kk,flip in ((k,False),(k2,True)):
            if kk in ours:
                try:
                    af=float(r["AF_Allele2"]); maf=min(af,1-af)
                    hit[kk]=dict(gene=ours[kk], beta=float(r["BETA"])*(-1 if flip else 1), se=float(r["SE"]),
                                 maf=maf, info=float(r["imputationInfo"]), n=float(r["N"]), flip=flip)
                except (ValueError,KeyError): pass
                break
        if n_lines%2000000==0: print(f"  scanned {n_lines/1e6:.0f}M lines, chr19 {n_chr19}, hits {len(hit)}", flush=True)
print(f"total lines {n_lines} | chr19 rows {n_chr19} | matched {len(hit)}", flush=True)
def summ(sub, label):
    if not sub: print(f"  {label}: 0"); return
    z=[abs(v["beta"]/v["se"]) for v in sub.values() if v["se"]>0]
    print(f"  {label}: n={len(sub)} | |z| 중위 {st.median(z):.2f} 최대 {max(z):.2f} | |z|>=2 {sum(1 for x in z if x>=2)} | |z|>=4 {sum(1 for x in z if x>=4)} "
          f"| SE 중위 {st.median([v['se'] for v in sub.values()]):.4f} | info 중위 {st.median([v['info'] for v in sub.values()]):.3f} | N 중위 {st.median([v['n'] for v in sub.values()]):.0f}")
print("=== 커버리지")
print(f"  우리 대역 변이 중 BBJ TC 보유: {len(hit)}/{len(ours)} = {100*len(hit)/len(ours):.1f}%")
summ(hit, "전체 매칭")
summ({k:v for k,v in hit.items() if v["gene"] in tracked}, "추적 11유전자")
summ({k:v for k,v in hit.items() if v["maf"]<0.01}, "BBJ MAF<1%")
summ({k:v for k,v in hit.items() if v["maf"]<0.01 and v["info"]>=0.7}, "BBJ MAF<1% & info>=0.7")
byg={}
for k,v in hit.items(): byg.setdefault(v["gene"],[]).append(v)
tg={g:len(v) for g,v in byg.items() if g in tracked}
print("  추적 유전자별 매칭 수:", tg)
print("  추적 유전자 중 매칭 3개 이상:", sum(1 for g,c in tg.items() if c>=3), "/ 11")
json.dump({k:v for k,v in hit.items()}, open(f"{W}/run_l3b/out/bbj_tc_chr19_beta.json","w"))
print("BBJCOV_DONE")
