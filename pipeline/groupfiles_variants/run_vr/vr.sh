#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_vr/.lock
flock -n 9 || { echo "another instance running"; exit 0; }

BCF=${BCFTOOLS:-bcftools}
P2=${PLINK2:-plink2}
RAW=${GENOTYPE_DIR:?Set GENOTYPE_DIR}
PRUNED=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
OUT=${PROJECT_ROOT}/work/ref/vr_plink_v3
W=${PROJECT_ROOT}/work/run_vr
mkdir -p "$OUT" "$W/chunks"
export TMPDIR=${PROJECT_ROOT}/work/tmp

LO=${MAF_LO}
HI=${MAF_HI}
R2MIN=0.3
TARGET_PER_CHR=200

for f in "$BCF" "$P2" "$PRUNED.bed" "$PRUNED.bim" "$PRUNED.fam"; do
  [ -e "$f" ] || { echo "PREFLIGHT FAIL missing $f"; exit 3; }
done
echo "preflight OK  pruned markers=$(wc -l < $PRUNED.bim)  samples=$(wc -l < $PRUNED.fam)"

one() {
  local N=$1
  local V=$RAW/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}${N}${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
  local O=$W/chunks/chr${N}.rare
  [ -s "$O.done" ] && { echo "  chr$N cached"; return 0; }
  [ -e "$V" ] || { echo "  chr$N SOURCE MISSING"; return 1; }
  $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MAF\t%INFO/R2\n' "$V" 2>/dev/null | \
    awk -v lo=$LO -v hi=$HI -v r2=$R2MIN -v want=$TARGET_PER_CHR \
      'BEGIN{n=0} {m=$5+0; r=$6+0;
        if(m>lo && m<=hi && r>=r2){print $1"\t"$2"\t"$3"\t"$4; n++; if(n>=want) exit}}' \
    > "$O.tsv"
  local k=$(wc -l < "$O.tsv")
  if [ "$k" -eq 0 ]; then echo "  chr$N GATE FAIL: 0 rare markers"; return 1; fi
  touch "$O.done"; echo "  chr$N rare markers: $k"
}
export -f one; export BCF RAW W LO HI R2MIN TARGET_PER_CHR
seq 1 22 | xargs -P 8 -I{} bash -c 'one {}'

cat $W/chunks/chr*.rare.tsv > $W/rare_sites.tsv
NRARE=$(wc -l < $W/rare_sites.tsv)
echo "rare sites collected: $NRARE"
[ "$NRARE" -ge 1000 ] || { echo "GATE FAIL: only $NRARE rare sites (<1000)"; exit 6; }
echo "VR_STAGE1_COMPLETE"
