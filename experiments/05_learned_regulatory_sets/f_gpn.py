#!/usr/bin/env python
import os as _cfg_os
import re as _cfg_re

def _config_path(value):
    names = _cfg_re.findall(r"\$(?:\{([A-Za-z_]\w*)\}|([A-Za-z_]\w*))", value)
    for pair in names:
        name = pair[0] or pair[1]
        if not _cfg_os.environ.get(name, "").strip():
            raise ValueError("Set required path environment variable: " + name)
    return _cfg_os.path.expandvars(value)


import gzip, os, sys, json, collections
ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/fset/out"
SRC=f"{ROOT}/ref/gpnmsa/scores.tsv.bgz"
need=collections.defaultdict(set)
with open(sys.argv[1]) as f:
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[0]=="key" or len(p)<3: continue
        need[p[1]].add(int(p[2]))
tot=sum(len(v) for v in need.values())
hdr=None; n_lines=0; n_hit=0
tmp=f"{OUT}/gpn_scores.tsv.gz.tmp"
with gzip.open(SRC,"rt") as f, gzip.open(tmp,"wt") as out:
    hdr=f.readline().rstrip("\n").split("\t")
    out.write("chr\tpos38\t"+"\t".join(hdr[2:])+"\n")
    ci,pi=0,1
    for line in f:
        n_lines+=1
        p=line.rstrip("\n").split("\t")
        ch=p[ci].replace("chr","")
        s=need.get(ch)
        if not s: continue
        try: po=int(p[pi])
        except ValueError: continue
        if po in s:
            n_hit+=1
            out.write(ch+"\t"+str(po)+"\t"+"\t".join(p[2:])+"\n")
os.replace(tmp,f"{OUT}/gpn_scores.tsv.gz")
json.dump({"header":hdr,"lines_scanned":n_lines,"positions_requested":tot,"positions_found":n_hit},
          open(f"{OUT}/gpn_scores.meta.json","w"))
print(json.dumps({"lines_scanned":n_lines,"requested":tot,"found":n_hit,"header":hdr[:8]}))
