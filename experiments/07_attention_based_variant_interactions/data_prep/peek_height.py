#!/usr/bin/env python
import os, glob, csv, json
R=os.environ["PHENOTYPE_DATA_ROOT"]
for f in glob.glob(os.environ["HEIGHT_SCAN_GLOB"]):
    h=open(f,errors="ignore").readline().rstrip("\n").split("\t")
    hh=[x for x in h if any(k in x.upper() for k in ("HEIGHT","HTCM","WEIGHT","BMI"))]
    if hh: print("source", os.path.basename(f)[:50], hh[:6])
for rel,col in json.loads(os.environ["HEIGHT_SOURCES_JSON"]):
    f=os.path.join(R,rel)
    rd=csv.reader(open(f,errors="ignore"),delimiter="\t"); h=next(rd); i=h.index(col); n=0; ok=0
    for row in rd:
        n+=1
        try:
            v=float(row[i]); ok+= (100<v<220)
        except: pass
    print(os.path.basename(f), "rows", n, "valid_height", ok)
