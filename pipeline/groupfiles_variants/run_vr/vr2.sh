#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_vr/.lock4
flock -n 9 || { echo "another instance running"; exit 0; }

BCF=${BCFTOOLS:-bcftools}
P2=${PLINK2:-plink2}
P1=${PROJECT_ROOT}/work/ref/plink19_env/bin/plink
RAW=${GENOTYPE_DIR:?Set GENOTYPE_DIR}
PRUNED=${PRUNED_PLINK_PREFIX:?Set PRUNED_PLINK_PREFIX}
OUT=${PROJECT_ROOT}/work/ref/vr_plink_v3
W=${PROJECT_ROOT}/work/run_vr
mkdir -p "$OUT" "$W/g"
export TMPDIR=${PROJECT_ROOT}/work/tmp

one() {
  local N=$1; local V=$RAW/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}${N}${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
  local S=$W/hc/chr${N}.tsv; local O=$W/g/chr${N}
  [ -s "$O.bed" ] && [ -s "$O.ok" ] && [ -s "$O.idfix" ] && { echo "  chr$N cached"; return 0; }
  awk '{print $1"\t"$2}' "$S" > "$O.reg"
  $BCF view -R "$O.reg" -Oz -o "$O.vcf.gz" "$V" 2>"$O.err"
  local nrec=$($BCF view -H "$O.vcf.gz" 2>/dev/null | wc -l)
  [ "$nrec" -eq 0 ] && { echo "  chr$N GATE FAIL: 0 records"; tail -2 "$O.err"; return 1; }
  $P2 --vcf "$O.vcf.gz" --make-bed --const-fid 0 \
      --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 60 missing \
      --out "$O" >"$O.p2.log" 2>&1
  [ -s "$O.bed" ] || { echo "  chr$N GATE FAIL: no bed"; tail -3 "$O.p2.log"; return 1; }
  touch "$O.ok" "$O.idfix"; echo "  chr$N records=$nrec markers=$(wc -l < $O.bim) samples=$(wc -l < $O.fam)"
}
export -f one; export BCF P2 RAW W
seq 1 22 | xargs -P 6 -I{} bash -c 'one {}'
nd=$(ls $W/g/chr*.ok 2>/dev/null | wc -l); echo "extracted: $nd/22"
[ "$nd" -eq 22 ] || { echo "GATE FAIL: $nd/22"; exit 6; }

: > $W/ml.txt; for N in $(seq 2 22); do echo "$W/g/chr${N}" >> $W/ml.txt; done
$P1 --bfile $W/g/chr1 --merge-list $W/ml.txt --keep-allele-order \
    --make-bed --out $W/rare_all >$W/m1.log 2>&1
[ -s $W/rare_all.bed ] || { echo "GATE FAIL: rare merge"; tail -6 $W/m1.log; exit 6; }
echo "rare merged: $(wc -l < $W/rare_all.bim) markers"

echo "$W/rare_all" > $W/ml2.txt
$P1 --bfile "$PRUNED" --bmerge $W/rare_all --keep-allele-order \
    --make-bed --out $OUT/vr_v3 >$W/m2.log 2>&1
[ -s $OUT/vr_v3.bed ] || { echo "  merge error:"; grep -m3 -i error $W/m2.log; }
[ -s $OUT/vr_v3.bed ] || { echo "GATE FAIL: final merge"; tail -8 $W/m2.log; exit 6; }
NM=$(wc -l < $OUT/vr_v3.bim); NS=$(wc -l < $OUT/vr_v3.fam)
echo "FINAL: $NM markers  $NS samples"
EXP=$(( $(wc -l < ${PRUNED}.bim) + $(wc -l < $W/rare_all.bim) ))
[ "$NM" -eq "$EXP" ] || { echo "GATE FAIL: markers $NM != pruned+rare $EXP"; exit 6; }
echo "  identity OK: $(wc -l < ${PRUNED}.bim) + $(wc -l < $W/rare_all.bim) = $NM"
DUP=$(awk '{print $2}' $OUT/vr_v3.bim | sort | uniq -d | wc -l)
[ "$DUP" -eq 0 ] || { echo "GATE FAIL: $DUP duplicate variant IDs"; exit 6; }

$P2 --bfile $OUT/vr_v3 --freq counts --out $W/ff >/dev/null 2>&1
awk 'NR>1{c=$6+0; if(c>10&&c<=20.5)a++; else if(c>20.5)b++; else z++}
     END{printf "MAC<=10: %d\nMAC 10-20.5: %d\nMAC>20.5: %d\n",z,a,b;
         if(a<200){print "GATE FAIL: cat1 <200"; exit 6}
         if(b<200){print "GATE FAIL: cat2 <200"; exit 6}}' $W/ff.acount || exit 6
[ "$NS" -eq ${N_SAMPLES} ] || { echo "GATE FAIL: samples $NS != ${N_SAMPLES}"; exit 6; }
echo "VR_PLINK_V3_COMPLETE"
