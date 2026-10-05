#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_proc/.proc.lock
flock -n 9 || { echo "another procurement running"; exit 0; }
R=${PROJECT_ROOT}/work/ref
mkdir -p $R

get(){
  local NAME=$1 URL=$2 DIR=$R/$3 MIN=${4:-100000}
  mkdir -p "$DIR"
  local F="$DIR/$(basename "${URL%%\?*}")"
  [ -e "$DIR/.done" ] && { echo "[$NAME] skip"; return 0; }
  echo "[$NAME] downloading $(basename "$F")"
  curl -sSL --retry 3 --retry-delay 5 -m 36000 -o "$F.part" "$URL" || { echo "[$NAME] FAIL download"; return 1; }
  local SZ=$(stat -c %s "$F.part" 2>/dev/null || echo 0)
  if [ "$SZ" -lt "$MIN" ]; then echo "[$NAME] FAIL size $SZ < $MIN"; rm -f "$F.part"; return 1; fi
  if head -c 200 "$F.part" | grep -qiE '<!doctype|<html'; then
    echo "[$NAME] FAIL html error page"; rm -f "$F.part"; return 1; fi
  mv "$F.part" "$F"
  echo "[$NAME] OK $(du -h "$F" | cut -f1)"
  touch "$DIR/.done"
}

U=https://hgdownload.soe.ucsc.edu
get phyloP447    $U/goldenPath/hg38/phyloP447way/hg38.phyloP447way.bw        phylop     1000000000 &
get phastCons470 $U/goldenPath/hg38/phastCons470way/hg38.phastCons470way.bw  phastcons  1000000000 &
get CpG_island   $U/goldenPath/hg38/database/cpgIslandExt.txt.gz             cpg        100000 &
get ReMap2022    https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/remap2022_nr_macs2_hg38_v1_0.bed.gz  remap  100000000 &
wait
get rmsk         $U/goldenPath/hg38/database/rmsk.txt.gz                     repeats    10000000 &
get mappability  $U/gbdb/hg38/hoffmanMappability/k36.Umap.MultiTrackMappability.bw  mappability 100000000 &
get JASPAR2024   https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt  jaspar 100000 &
get UTRannotator https://raw.githubusercontent.com/ImperialCardioGenetics/UTRannotator/master/uORF_5UTR_GRCh38_PUBLIC.txt  utrannotator 1000000 &
get TargetScan   https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Target_Locations.default_predictions.hg19.bed.zip  targetscan 1000000 &
get gnomad_constraint https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv  gnomad_constraint 10000000 &
wait

echo "=== PROCUREMENT SUMMARY ==="
for d in phylop phastcons cpg remap repeats mappability jaspar utrannotator targetscan gnomad_constraint; do
  if [ -e "$R/$d/.done" ]; then printf '  HAVE %-20s %s\n' "$d" "$(du -sh $R/$d 2>/dev/null | cut -f1)";
  else printf '  FAIL %-20s\n' "$d"; fi
done
echo PROCUREMENT_COMPLETE
