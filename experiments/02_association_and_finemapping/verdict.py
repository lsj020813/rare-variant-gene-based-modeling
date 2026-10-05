
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, json, glob, os
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
def read(fp):
    out={}
    if not os.path.exists(fp): return out
    with open(fp) as fh:
        h=fh.readline().rstrip("\n").split("\t")
        if "Pvalue" not in h: return out
        gi=h.index("Region"); pi=h.index("Pvalue")
        bu=h.index("Pvalue_Burden") if "Pvalue_Burden" in h else None
        for line in fh:
            fl=line.rstrip("\n").split("\t")
            try: out[fl[gi].split(".")[0]]=(float(fl[pi]), float(fl[bu]) if bu is not None else None)
            except (ValueError,IndexError): pass
    return out
THR=2.5e-6
allrows=[]
for L in ("APOE19","LDLR19","SORT1_1","CETP16","CHR5","APOB2","APO11"):
    lead=json.load(open(f"{W}/{L}/lead.json"))
    un=read(f"{W}/{L}/tchl.uncond"); co=read(f"{W}/{L}/tchl.cond")
    print(f"\n===== {L}  lead {lead['lead']}  lead_p={lead['lead_p']:.1e} =====")
    print(f"  {'gene':<12}{'uncond p':>11}{'cond p':>11}{'ratio':>9}  판정")
    for g,(pu,_) in sorted(un.items(), key=lambda x:x[1][0]):
        pc = co.get(g,(None,None))[0]
        if pc is None: continue
        was = pu<THR; still = pc<THR
        verdict = ("유지" if still else "소실") if was else ("신규?" if still else "-")
        allrows.append((L,sym.get(g,g),pu,pc,verdict))
        print(f"  {sym.get(g,g):<12}{pu:>11.1e}{pc:>11.1e}{pc/pu:>9.0f}x  {verdict}")
surv=[r for r in allrows if r[4]=="유지"]
lost=[r for r in allrows if r[4]=="소실"]
print(f"\n=== 총괄: 유의였던 것 {len(surv)+len(lost)}개 중 조건부 생존 {len(surv)} / 소실 {len(lost)} ===")
for L,g,pu,pc,_ in surv: print(f"  생존: {L} {g}  {pu:.1e} -> {pc:.1e}")
json.dump([{"locus":L,"gene":g,"uncond":pu,"cond":pc,"verdict":v} for L,g,pu,pc,v in allrows],
          open(f"{W}/conditional_verdict.json","w"), indent=1)
print("WRITTEN")
