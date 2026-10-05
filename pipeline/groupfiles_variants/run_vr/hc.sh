#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_vr/.lock3
flock -n 9 || { echo "another instance running"; exit 0; }

BCF=${BCFTOOLS:-bcftools}
RAW=${GENOTYPE_DIR:?Set GENOTYPE_DIR}
W=${PROJECT_ROOT}/work/run_vr
mkdir -p "$W/hc"
export TMPDIR=${PROJECT_ROOT}/work/tmp

WANT=${WANT:-200}
R2MIN=0.3

one() {
  local N=$1
  local V=$RAW/${GENOTYPE_FILENAME_PREFIX:?Set GENOTYPE_FILENAME_PREFIX}${N}${GENOTYPE_FILENAME_SUFFIX:?Set GENOTYPE_FILENAME_SUFFIX}
  local O=$W/hc/chr${N}
  [ -s "$O.done" ] && { echo "  chr$N cached ($(wc -l < $O.tsv))"; return 0; }
  [ -e "$V" ] || { echo "  chr$N SOURCE MISSING"; return 1; }
  $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/R2\t[%GT]\n' "$V" 2>/dev/null | \
    awk -v want=$WANT -v r2=$R2MIN 'BEGIN{n=0}
      {
        if ($5+0 < r2) next
        s=$6; ac=gsub(/1/,"1",s)
        if (ac > 10 && ac <= 20) { print $1"\t"$2"\t"$3"\t"$4"\t"ac; n++; if (n>=want) exit }
      }' > "$O.tsv"
  local k=$(wc -l < "$O.tsv")
  [ "$k" -eq 0 ] && { echo "  chr$N GATE FAIL: 0 hard-call MAC 11-20 markers"; return 1; }
  touch "$O.done"; echo "  chr$N hard-call markers: $k"
}
export -f one; export BCF RAW W WANT R2MIN
seq 1 22 | xargs -P 8 -I{} bash -c 'one {}'

ndone=$(ls $W/hc/chr*.done 2>/dev/null | wc -l)
cat $W/hc/chr*.tsv > $W/hc_sites.tsv
NR_=$(wc -l < $W/hc_sites.tsv)
echo "chromosomes: $ndone/22   hard-call sites: $NR_"
[ "$ndone" -eq 22 ] || { echo "GATE FAIL: only $ndone/22 chromosomes"; exit 6; }
[ "$NR_" -ge 1000 ] || { echo "GATE FAIL: $NR_ sites < 1000"; exit 6; }
awk '{s+=$5; if($5<mn||mn==0)mn=$5; if($5>mx)mx=$5} END{printf "  MAC range %d-%d, mean %.1f\n", mn, mx, s/NR}' $W/hc_sites.tsv
echo "HC_STAGE1_COMPLETE"
