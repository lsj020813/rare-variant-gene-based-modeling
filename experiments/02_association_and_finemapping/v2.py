
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, json, os
W=_config_path("${PROJECT_ROOT}/work/run_cond")
R=_config_path("${PROJECT_ROOT}/work/ref")
sym={}
with gzip.open(f"{R}/deductive/gencode.sorted.gtf.gz","rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f=line.split("\t",9)
        if len(f)<9 or f[2]!="gene": continue
        a=f[8]
        if 'gene_name "' in a:
            sym[a.split('gene_id "')[1].split('"')[0].split(".")[0]]=a.split('gene_name "')[1].split('"')[0]
THR=2.5e-6
rows=[]
for L in ("APOE19","LDLR19","SORT1_1","CETP16","CHR5","APOB2","APO11"):
    fp=f"{W}/{L}/tchl.cond"
    if L=="SORT1_1": fp=f"{W}/{L}/tchl.cond3"
    if not os.path.exists(fp): print(f"{L}: MISSING {fp}"); continue
    with open(fp) as fh:
        h=fh.readline().rstrip("\n").split("\t")
        has_cond = "Pvalue_cond" in h
        print(f"{L}: has Pvalue_cond = {has_cond}  (file {os.path.basename(fp)})")
        if not has_cond: continue
        gi=h.index("Region"); pi=h.index("Pvalue"); ci=h.index("Pvalue_cond")
        bi=h.index("Pvalue_Burden_cond")
        for line in fh:
            fl=line.rstrip("\n").split("\t")
            try:
                pu, pc, bc = float(fl[pi]), float(fl[ci]), float(fl[bi])
            except (ValueError,IndexError): continue
            rows.append((L, sym.get(fl[gi].split(".")[0], fl[gi]), pu, pc, bc))
print()
print(f"{'locus':<9}{'gene':<13}{'uncond':>10}{'cond':>10}{'burden_c':>10}  판정")
surv=0; lost=0
for L,g,pu,pc,bc in sorted(rows, key=lambda x:(x[0],x[2])):
    was, still = pu<THR, pc<THR
    if was:
        v = "유지" if still else "소실"
        surv += still; lost += (not still)
    else:
        v = "-"
    print(f"{L:<9}{g:<13}{pu:>10.1e}{pc:>10.1e}{bc:>10.1e}  {v}")
print(f"\n유의 43개 중: 생존 {surv} / 소실 {lost}")
json.dump([{"locus":L,"gene":g,"uncond":pu,"cond":pc,"burden_cond":bc} for L,g,pu,pc,bc in rows],
          open(f"{W}/conditional_verdict2.json","w"), indent=1)
print("VERDICT2_WRITTEN")
