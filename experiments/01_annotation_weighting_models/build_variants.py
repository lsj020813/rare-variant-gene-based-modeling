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


import numpy as np, glob, csv, math, json, hashlib, sys
W=_config_path("${PROJECT_ROOT}/work"); O=_config_path("${PROJECT_ROOT}/work/run_trackB")
GP=f"{W}/run_l3b/out/G_phi.txt"
T="ENSG00000129353.15 ENSG00000213892.12 ENSG00000186567.14 ENSG00000130202.10 ENSG00000130204.13 ENSG00000104856.15 ENSG00000069399.15 ENSG00000079805.19 ENSG00000142453.13 ENSG00000127616.22 ENSG00000129354.12".split()
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda: f.read(1<<20), b""): h.update(b)
    return h.hexdigest()
gk={}
for line in open(GP):
    t=line.rstrip("\n").split()
    if len(t)>2 and t[1]=="var": gk[t[0]]=t[2:]
rng=np.random.default_rng(20260910)
others=sorted(g for g in gk if g not in T)
ctrl=[str(x) for x in rng.choice(others, 20, replace=False)]
sel=T+ctrl
sz={}; files=sorted(glob.glob(f"{W}/run_l3b/out/arm_none/part*.singleAssoc.txt"))
for f in files:
    for r in csv.DictReader(open(f), delimiter="\t"):
        try:
            v=float(r["var"])
            if v<=0: continue
            k=f'{r["CHR"]}:{r["POS"]}:{r["Allele1"]}:{r["Allele2"]}'
            sz[k]=(abs(float(r["Tstat"])/math.sqrt(v)), float(r["BETA"]), float(r["SE"]), float(r["p.value"]))
        except (ValueError,KeyError): pass
rows={}
for g in sel:
    for k in gk[g]:
        d=rows.setdefault(k, {"genes":[], "is_target":0})
        d["genes"].append(g)
        if g in T: d["is_target"]=1
n_nozz=0
with open(f"{O}/variants.tsv","w") as fo:
    fo.write("key\tchr\tpos\tref\talt\tgenes\tn_genes_sel\tis_target\tz_abs\tbeta\tse\tp\tis_indel\n")
    for k in sorted(rows, key=lambda k:int(k.split(":")[1])):
        c,p,ref,alt=k.split(":")
        if k not in sz: n_nozz+=1; continue
        z,b,s,pv=sz[k]
        fo.write(f"{k}\t{c}\t{p}\t{ref}\t{alt}\t{';'.join(rows[k]['genes'])}\t{len(rows[k]['genes'])}\t{rows[k]['is_target']}\t{z:.6g}\t{b:.6g}\t{s:.6g}\t{pv:.6g}\t{int(len(ref)!=1 or len(alt)!=1)}\n")
n_pairs=sum(len(gk[g]) for g in sel)
uniq=len(rows); n_indel=sum(1 for k in rows if len(k.split(":")[2])!=1 or len(k.split(":")[3])!=1)
man={"G_phi_sha256":sha(GP),"singleAssoc_files":len(files),"seed_ctrl":20260910,"target_genes":T,"ctrl_genes":ctrl,
     "n_variant_gene_pairs":n_pairs,"n_unique_variants":uniq,"n_unique_with_z":uniq-n_nozz,"n_no_z":n_nozz,"n_indel":n_indel,
     "per_gene_n":{g:len(gk[g]) for g in sel},"build":"GRCh37 (G_phi keys == fm_all.npz key37, 115237/115237)"}
json.dump(man, open(f"{O}/variants_manifest.json","w"), indent=1)
print(json.dumps({k:v for k,v in man.items() if k not in ("per_gene_n","ctrl_genes","target_genes")}))
print("per-gene n: min", min(man["per_gene_n"].values()), "median", int(np.median(list(man["per_gene_n"].values()))), "max", max(man["per_gene_n"].values()))
print("BUILD_DONE")
