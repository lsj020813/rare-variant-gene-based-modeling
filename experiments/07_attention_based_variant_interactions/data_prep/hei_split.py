#!/usr/bin/env python
import os as _os
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
import csv, json, collections, os
W=(_os.environ["PROJECT_ROOT"] + '/work/prs'); os.makedirs(W+"/private/hei",exist_ok=True)
H={r["sample_id"] for r in csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_hei/hei.tsv')),delimiter="\t")}
rows=[l.rstrip("\n").split("\t") for l in open(W+"/private/split.tsv")]
hdr,rows=rows[0],rows[1:]; keep=[r for r in rows if r[0] in H]
with open(W+"/private/hei/split_hei.tsv","w") as g:
    g.write("\t".join(hdr)+"\n"); [g.write("\t".join(r)+"\n") for r in keep]
c0=collections.Counter(r[1] for r in rows); c1=collections.Counter(r[1] for r in keep)
json.dump(dict(before=dict(c0),after=dict(c1)),open(W+"/out/hei/split_hei_summary.json","w")); print(dict(c0),dict(c1))
