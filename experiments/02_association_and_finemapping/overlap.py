import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import glob, itertools, json
G=_config_path('${PROJECT_ROOT}/work/ref/groupfiles_chunks')
GENES={"ENSG00000130203":"APOE","ENSG00000130204":"TOMM40","ENSG00000130208":"APOC1",
       "ENSG00000234906":"APOC2","ENSG00000267467":"APOC4","ENSG00000224916":"APOC4-APOC2",
       "ENSG00000104853":"CLPTM1","ENSG00000130202":"NECTIN2","ENSG00000142273":"CBLC",
       "ENSG00000187244":"BCAM","ENSG00000069399":"BCL3","ENSG00000073008":"PVR",
       "ENSG00000130164":"LDLR"}
V={}
for fp in glob.glob(f'{G}/chr19.part*.txt'):
    for line in open(fp):
        f=line.rstrip("\n").split()
        if len(f)<3 or f[1]!="var": continue
        base=f[0].split(".")[0]
        if base in GENES: V.setdefault(GENES[base], set()).update(f[2:])
print("gene sizes:", {k:len(v) for k,v in sorted(V.items(), key=lambda x:-len(x[1]))})
print()
print("pairwise Jaccard (|A∩B| / |A∪B|) and shared count:")
names=[g for g in ("APOE","TOMM40","APOC1","APOC2","APOC4","NECTIN2","CLPTM1","BCAM","BCL3","CBLC","PVR","LDLR") if g in V]
print("        " + "".join(f"{n[:7]:>9}" for n in names))
J={}
for a in names:
    row=f"{a[:7]:<8}"
    for b in names:
        if a==b: row+=f"{'-':>9}"; continue
        inter=len(V[a]&V[b]); uni=len(V[a]|V[b])
        j=inter/uni if uni else 0
        J[f"{a}|{b}"]=[inter,round(j,3)]
        row+=f"{j:>9.2f}"
    print(row)
print()
print("APOE vs TOMM40 detail:", "shared", len(V["APOE"]&V["TOMM40"]),
      "| APOE only", len(V["APOE"]-V["TOMM40"]), "| TOMM40 only", len(V["TOMM40"]-V["APOE"]))
print("LDLR overlaps:", {b: len(V["LDLR"]&V[b]) for b in names if b!="LDLR"})
json.dump({"sizes":{k:len(v) for k,v in V.items()}, "jaccard":J}, open("/tmp/overlap.json","w"))
