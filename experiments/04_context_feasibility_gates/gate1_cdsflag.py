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


import gzip, sys, os, re, subprocess, bisect, json
CHR=sys.argv[1]; ROOT=_config_path("${PROJECT_ROOT}/work"); OUT=f"{ROOT}/gate1/out"
BCF="bcftools"
GENCODE=f"{ROOT}/ref/deductive/gencode.nochr.gtf.gz"
KEYED=f"{ROOT}/ref/lift38_keyed/chr{CHR}.keyed38.vcf.gz"
WJ=f"{OUT}/windows_chr{CHR}.json.gz"

need=set()
with gzip.open(WJ,"rt") as f:
    W=json.load(f)
for g,ent in W.items():
    for b,lst in ent["bins"].items():
        for rec in lst: need.add(rec[0])

iv=[]
with gzip.open(GENCODE,"rt") as f:
    for line in f:
        if line[0]=="#": continue
        p=line.split("\t",9)
        if p[0]!=CHR or p[2]!="CDS": continue
        iv.append((int(p[3]),int(p[4])))
iv.sort(); merged=[]
for s,e in iv:
    if merged and s<=merged[-1][1]: merged[-1][1]=max(merged[-1][1],e)
    else: merged.append([s,e])
starts=[x[0] for x in merged]
def coding(p):
    i=bisect.bisect_right(starts,p)-1
    return 1 if (i>=0 and merged[i][1]>=p) else 0

n=0; nc=0
with gzip.open(f"{OUT}/vpos_chr{CHR}.tsv.gz.tmp","wt") as out:
    out.write("hg19_id\thg38_pos\tis_coding\n")
    proc=subprocess.Popen([BCF,"query","-f","%CHROM\t%POS\t%ID\n",KEYED],stdout=subprocess.PIPE,text=True,bufsize=1<<20)
    for line in proc.stdout:
        ch,pos,vid=line.rstrip("\n").split("\t")
        if vid not in need: continue
        if ch not in (f"chr{CHR}",CHR): continue
        p=int(pos); cd=coding(p); n+=1; nc+= (cd==0)
        out.write(f"{vid}\t{p}\t{cd}\n")
    proc.wait()
os.replace(f"{OUT}/vpos_chr{CHR}.tsv.gz.tmp", f"{OUT}/vpos_chr{CHR}.tsv.gz")
print(json.dumps({"chr":CHR,"assigned_variants_needed":len(need),"written":n,"noncoding":nc,
                  "cds_intervals":len(merged)}))
