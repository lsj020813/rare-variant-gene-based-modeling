#!/usr/bin/env bash
: "${CONDA_INIT_SCRIPT:?Set CONDA_INIT_SCRIPT}"
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
N=$1; D="${GENOTYPE_PREFIX:?Set GENOTYPE_PREFIX}"
M=${D}.map; P=${D}.ped
R=${GENOTYPE_DIR}/chr$N.vcf.gz
O=${PROJECT_ROOT}/work/run_cmp_hexa; mkdir -p $O; T=$O/chr$N
awk -v c=$N -F'\t' '$1==c{print NR"\t"$2"\t"$4}' $M > $T.map.tsv
NM=$(wc -l < $T.map.tsv)
source ${CONDA_INIT_SCRIPT}; conda activate "${CONDA_ENV:?Set CONDA_ENV}" 2>/dev/null
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
bcftools view -i 'INFO/TYPED=1' -G $R 2>/dev/null | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/AF\n' > $T.vcf_typed.tsv
NT=$(wc -l < $T.vcf_typed.tsv)
bcftools view -G $R 2>/dev/null | bcftools query -f '%POS\n' | wc -l > $T.vcf_all.n
python3 - "$P" "$T" <<'EOF'
import sys, collections
P, T = sys.argv[1], sys.argv[2]
idx = [int(l.split("\t")[0]) for l in open(T + ".map.tsv")]          # marker index within map (1-based)
cols = [6 + (i - 1) for i in idx]                                     # 0-based token index of "A B" pair in ped (6 leading cols)
cnt = [collections.Counter() for _ in idx]; n = 0
with open(P) as f:
    for line in f:
        t = line.rstrip("\n").split("\t"); n += 1
        for j, c in enumerate(cols):
            g = t[c]
            if g != "0 0": cnt[j][g[0]] += 1; cnt[j][g[2]] += 1
with open(T + ".ped_af.tsv", "w") as o:
    for j, l in enumerate(open(T + ".map.tsv")):
        rsid, pos = l.split("\t")[1], l.split("\t")[2].strip()
        c = cnt[j]; tot = sum(c.values())
        if tot == 0 or len(c) > 2: o.write(f"{pos}\t{rsid}\tNA\tNA\tNA\n"); continue
        al = sorted(c, key=lambda k: -c[k]); a1 = al[0]; a2 = al[1] if len(al) > 1 else "."
        o.write(f"{pos}\t{rsid}\t{a1}\t{a2}\t{c[a2]/tot if a2!='.' else 0.0:.6f}\n")   # minor allele, its freq
print(f"PED_SAMPLES {n}")
EOF
python3 - "$T" "$N" "$NM" "$NT" <<'EOF'
import sys
T, N, NM, NT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
vcf = {}
for l in open(T + ".vcf_typed.tsv"):
    pos, ref, alt, af = l.rstrip("\n").split("\t")
    if "," in alt: continue
    try: vcf[pos] = (ref, alt, float(af))
    except ValueError: pass
flip = {"A":"T","T":"A","C":"G","G":"C"}
n_pos = n_allele = n_af01 = n_af001 = n_afexact = n_strand = 0; n_ped_ok = 0
for l in open(T + ".ped_af.tsv"):
    pos, rsid, a1, a2, mf = l.rstrip("\n").split("\t")
    if a1 == "NA": continue
    n_ped_ok += 1
    if pos not in vcf: continue
    n_pos += 1
    ref, alt, af = vcf[pos]; maf_v = min(af, 1 - af)
    s = {a1, a2} - {"."}; v = {ref, alt}
    if s <= v: n_allele += 1
    elif {flip.get(x, x) for x in s} <= v: n_allele += 1; n_strand += 1
    else: continue
    d = abs(float(mf) - maf_v)
    if d < 0.01: n_af01 += 1
    if d < 0.001: n_af001 += 1
    if d < 1e-6: n_afexact += 1
vall = int(open(T + ".vcf_all.n").read().strip())
print(f"chr{N} | 새map 마커 {NM} (빈도계산가능 {n_ped_ok}) | 우리VCF 전체 {vall} · TYPED {NT}")
print(f"chr{N} | 좌표 일치 {n_pos} | 대립유전자 일치 {n_allele} (그중 가닥반전 {n_strand})")
print(f"chr{N} | 빈도 일치: |ΔMAF|<0.01 {n_af01} · <0.001 {n_af001} · <1e-6 {n_afexact}")
print(f"chr{N} CMP_DONE")
EOF
