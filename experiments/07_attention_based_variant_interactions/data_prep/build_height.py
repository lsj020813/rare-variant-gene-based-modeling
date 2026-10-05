#!/usr/bin/env python
import os as _os
if not _os.environ["PHENOTYPE_DATA_ROOT"].strip():
    raise ValueError("Set nonblank PHENOTYPE_DATA_ROOT")
import math as _number_math
def _required_number(name, cast, positive=False):
    raw = _os.environ.get(name, "")
    if not raw.strip():
        raise ValueError(name + " must be set and nonblank")
    try:
        value = cast(raw)
    except (ValueError, OverflowError):
        raise ValueError(name + " has an invalid numeric value") from None
    if isinstance(value, float) and not _number_math.isfinite(value):
        raise ValueError(name + " must be finite")
    if positive and value <= 0:
        raise ValueError(name + " must be positive")
    return value
if not _os.environ.get("PROJECT_ROOT", "").strip():
    raise ValueError("PROJECT_ROOT must be set and nonblank")
N_SAMPLES = _required_number("N_SAMPLES", int, True)
import os, csv, json, subprocess, collections
BASE=_os.environ["PHENOTYPE_DATA_ROOT"]; OUT=(_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_hei'); os.makedirs(OUT,exist_ok=True)
BCF=_os.environ.get("BCFTOOLS_BIN", "bcftools"); VCF22=(_os.environ["PROJECT_ROOT"] + '/work/ref/band_vcf/chr22.band.vcf.gz')
raw=subprocess.run([BCF,"query","-l",VCF22],capture_output=True,text=True,check=True).stdout.split()
def norm(s):
    h=s.split("_"); return h[0] if len(h)==2 and h[0]==h[1] else s
dup_of={norm(s):s for s in raw}; assert len(dup_of)==N_SAMPLES
MISS={"","NA","NaN","nan",".","-9","99999","77777","88888","9999","999"}
def read_col(path,col):
    d={}
    with open(path,newline="",errors="replace") as fh:
        hdr=fh.readline().rstrip("\n").split("\t"); assert col in hdr, (path,col); i=hdr.index(col)
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)!=len(hdr): continue
            k,v=p[0],p[i].strip()
            if k in dup_of and v not in MISS: d.setdefault(k,v)
    return d
H=json.loads(_os.environ["HEIGHT_SOURCES_JSON"])
indicator_columns=json.loads(_os.environ["COHORT_INDICATOR_COLUMNS_JSON"])
assert len(indicator_columns) == 2
hei={}; log={}
for rel,col in H:
    d=read_col(f"{BASE}/{rel}",col); n0=len(hei); bad=0
    for k,v in d.items():
        try: x=float(v)
        except ValueError: bad+=1; continue
        if 100<x<220: hei.setdefault(k,x)
        else: bad+=1
    log[col]=dict(rows=len(d),added=len(hei)-n0,out_of_range_or_bad=bad)
n=0; byc=collections.Counter(); ys=[]
with open(OUT+"/hei.tsv.tmp","w") as g:
    rd=csv.DictReader(open((_os.environ["PROJECT_ROOT"] + '/work/ref/pheno_v3/tchl_v3.tsv')),delimiter="\t"); cols=rd.fieldnames
    g.write("\t".join(cols)+"\n")
    for r in rd:
        nk=norm(r["sample_id"])
        if nk not in hei: continue
        r["y"]=f"{hei[nk]:.2f}"; g.write("\t".join(r[c] for c in cols)+"\n"); n+=1; ys.append(hei[nk])
        byc["cohort_1" if r[indicator_columns[0]]=="1" else "cohort_2" if r[indicator_columns[1]]=="1" else "cohort_other"]+=1
os.replace(OUT+"/hei.tsv.tmp",OUT+"/hei.tsv")
import statistics
summ=dict(sources=log,pooled=len(hei),n_with_covariates=n,by_cohort=dict(byc),mean_cm=round(statistics.mean(ys),2),sd_cm=round(statistics.pstdev(ys),2))
json.dump(summ,open(OUT+"/hei_summary.json","w"),indent=1); print(json.dumps(summ))
