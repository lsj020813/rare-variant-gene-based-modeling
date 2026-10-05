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


import subprocess, glob, gzip, os, sys, json, collections

BCF="bcftools"
R=_config_path("${PROJECT_ROOT}/work/ref")
RAW=_config_path("${GENOTYPE_DIR}")
W=_config_path("${PROJECT_ROOT}/work/run_cond")
UD=(_config_path("PROOT_NO_SECCOMP=1 UDOCKER_DIR=${CONTAINER_STATE_DIR} "
    "${UDOCKER_BIN} run --volume=/data:/data ${SAIGE_IMAGE} "))

LOCI={
 "APOE19":  ("19", ["APOC4","APOC2","APOC4-APOC2","CLPTM1","TOMM40","APOE","APOC1","NECTIN2",
                    "CBLC","BCAM","BCL3","PVR","RELB","CEACAM16","CEACAM19"]),
 "LDLR19":  ("19", ["DNM2","AP1M2","ILF3","SLC44A2","CARM1","ATG4D","KRI1","CDKN2D","QTRT1",
                    "C19orf38","TMED1","SMARCA4","S1PR5","YIPF2","TIMM29","LDLR","ELOF1"]),
 "SORT1_1": ("1",  ["MYBPHL","SORT1"]),
 "CETP16":  ("16", ["SLC12A3","CETP","HERPUD1","NUP93","NLRC5"]),
 "CHR5":    ("5",  ["ANKDD1B","POLK","GCNT4","POC5"]),
 "APOB2":   ("2",  ["APOB"]),
 "APO11":   ("11", ["APOC3","APOA5","APOA1"]),
}
LOCUS=sys.argv[1]
CH, SYMS = LOCI[LOCUS]
OUT=f"{W}/{LOCUS}"; os.makedirs(OUT, exist_ok=True)

sym2id={}
with gzip.open(f"{R}/deductive/gencode.sorted.gtf.gz","rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f=line.split("\t",9)
        if len(f)<9 or f[2]!="gene" or f[0].replace("chr","")!=CH: continue
        a=f[8]
        if 'gene_name "' in a:
            nm=a.split('gene_name "')[1].split('"')[0]
            if nm in SYMS: sym2id[nm]=a.split('gene_id "')[1].split('"')[0]
missing=[s for s in SYMS if s not in sym2id]
assert not missing, f"GATE FAIL: symbols not in GTF chr{CH}: {missing}"
ids=set(sym2id.values())

lines=[]; pos=[]
for gf in sorted(glob.glob(f"{R}/groupfiles_chunks/chr{CH}.part*.txt")):
    for line in open(gf):
        f=line.split()
        if len(f)>2 and f[0] in ids:
            lines.append(line.rstrip("\n"))
            if f[1]=="var": pos += [int(k.split(":")[1]) for k in f[2:]]
assert lines and pos, "GATE FAIL: no group-file lines for locus genes"
got={l.split()[0] for l in lines}
assert got==ids, f"GATE FAIL: group-file genes {len(got)} != wanted {len(ids)}: missing {ids-got}"
lo, hi = min(pos)-200000, max(pos)+200000
open(f"{OUT}/group.txt","w").write("\n".join(lines)+"\n")
print(f"[{LOCUS}] genes {len(ids)}  window {CH}:{lo:,}-{hi:,}  var-positions {len(pos):,}")

rawf=f"{RAW}/chr{CH}.vcf.gz"
cand=f"{OUT}/common.vcf.gz"
subprocess.run(f"{BCF} view -r {CH}:{lo}-{hi} -i 'INFO/MAF>=0.01 && INFO/R2>=0.8' {rawf} -Oz -o {cand} && {BCF} index -c {cand}",
               shell=True, check=True)
ncand=int(subprocess.run(f"{BCF} index -n {cand}", shell=True, capture_output=True, text=True).stdout.strip() or 0)
assert ncand>0, "GATE FAIL: zero common candidates in window"
print(f"[{LOCUS}] common candidates: {ncand:,}")

sv=f"{OUT}/sv.tchl"
cmd=(UD + "step2_SPAtests.R "
     f"--vcfFile={cand} --vcfFileIndex={cand}.csi --vcfField=DS --chrom={CH} "
     "--AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE "
     f"--GMMATmodelFile={R}/saige_step1_v4/tchl_v4.rda "
     f"--varianceRatioFile={R}/saige_step1_v4/tchl_v4.varianceRatio.txt "
     f"--SAIGEOutputFile={sv}")
rc=subprocess.run(cmd+f" > {sv}.log 2>&1", shell=True).returncode
nrow=sum(1 for _ in open(sv))-1 if os.path.exists(sv) else 0
assert rc==0 and nrow>0, f"GATE FAIL: single-variant run rc={rc} rows={nrow}"
best=None
with open(sv) as fh:
    h=fh.readline().rstrip("\n").split("\t")
    ci,pi_,ri,ai,pv = h.index("CHR"),h.index("POS"),h.index("Allele1"),h.index("Allele2"),h.index("p.value")
    for line in fh:
        fl=line.rstrip("\n").split("\t")
        try: p=float(fl[pv])
        except ValueError: continue
        if best is None or p<best[0]: best=(p, fl[ci], fl[pi_], fl[ri], fl[ai])
assert best, "GATE FAIL: no parseable single-variant p"
p_lead, cch, cpos, cref, calt = best
cond=f"{cch}:{cpos}_{cref}/{calt}"
print(f"[{LOCUS}] lead SNP {cond}  p={p_lead:.2e}  (of {nrow:,} tested)")
json.dump({"locus":LOCUS,"lead":cond,"lead_p":p_lead,"n_candidates":nrow,
           "window":[lo,hi]}, open(f"{OUT}/lead.json","w"))

band=f"{R}/band_vcf/chr{CH}.band.vcf.gz"
slice_=f"{OUT}/band_slice.vcf.gz"; leadv=f"{OUT}/lead.vcf.gz"; mini=f"{OUT}/mini.vcf.gz"
subprocess.run(f"{BCF} view -r {CH}:{lo}-{hi} {band} -Oz -o {slice_} && {BCF} index -c {slice_}", shell=True, check=True)
subprocess.run(f"{BCF} view -r {CH}:{cpos}-{cpos} {rawf} -i 'REF==\"{cref}\" && ALT==\"{calt}\"' -Oz -o {leadv} && {BCF} index -c {leadv}", shell=True, check=True)
nlead=int(subprocess.run(f"{BCF} index -n {leadv}", shell=True, capture_output=True, text=True).stdout.strip() or 0)
assert nlead==1, f"GATE FAIL: lead extraction got {nlead} records"
subprocess.run(f"{BCF} concat -a {slice_} {leadv} 2>/dev/null | {BCF} sort -Oz -o {mini} && {BCF} index -c {mini}", shell=True, check=True)
nmini=int(subprocess.run(f"{BCF} index -n {mini}", shell=True, capture_output=True, text=True).stdout.strip() or 0)
nslice=int(subprocess.run(f"{BCF} index -n {slice_}", shell=True, capture_output=True, text=True).stdout.strip() or 0)
assert nmini==nslice+1, f"GATE FAIL: mini {nmini} != slice {nslice}+1"

for tag, extra in (("uncond",""), ("cond", f"--condition={cond} --weights_for_condition=1 ")):
    o=f"{OUT}/tchl.{tag}"
    cmd=(UD + "step2_SPAtests.R "
         f"--vcfFile={mini} --vcfFileIndex={mini}.csi --vcfField=DS --chrom={CH} "
         "--AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE "
         f"--GMMATmodelFile={R}/saige_step1_v4/tchl_v4.rda "
         f"--varianceRatioFile={R}/saige_step1_v4/tchl_v4.varianceRatio.txt "
         f"--groupFile={OUT}/group.txt --annotation_in_groupTest=all --maxMAF_in_groupTest=0.01 "
         + extra + f"--SAIGEOutputFile={o}")
    rc=subprocess.run(cmd+f" > {o}.log 2>&1", shell=True).returncode
    n=sum(1 for _ in open(o))-1 if os.path.exists(o) else 0
    assert rc==0 and n==len(ids), f"GATE FAIL: {tag} rc={rc} rows={n}/{len(ids)}"
    if tag=="cond":
        logtxt=open(f"{o}.log", errors="ignore").read()
        assert "condition" in logtxt.lower(), "GATE FAIL: condition marker not acknowledged in log"
print(f"[{LOCUS}] LOCUS_DONE")
