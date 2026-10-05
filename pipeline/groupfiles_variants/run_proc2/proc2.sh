#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_proc2/.proc2.lock
flock -n 9 || { echo "another run"; exit 0; }
R=${PROJECT_ROOT}/work/ref

get(){
  local NAME=$1 URL=$2 DIR=$R/$3 MIN=$4
  mkdir -p "$DIR"; local F="$DIR/$(basename "${URL%%\?*}")"
  [ -e "$DIR/.done" ] && { echo "[$NAME] skip"; return 0; }
  echo "[$NAME] start"
  curl -sSL --retry 5 --retry-delay 10 -m 86400 -o "$F.part" "$URL" || { echo "[$NAME] FAIL curl"; return 1; }
  local SZ=$(stat -c %s "$F.part" 2>/dev/null || echo 0)
  [ "$SZ" -lt "$MIN" ] && { echo "[$NAME] FAIL size $SZ<$MIN"; rm -f "$F.part"; return 1; }
  head -c 200 "$F.part" | grep -qiE '<!doctype|<html' && { echo "[$NAME] FAIL html"; rm -f "$F.part"; return 1; }
  mv "$F.part" "$F"; touch "$DIR/.done"; echo "[$NAME] OK $(du -h "$F" | cut -f1)"
}

get Pangolin https://zenodo.org/api/records/15649338/files/Pangolin_hg38_snvs_masked.zip/content pangolin 10000000000 &

DV=https://dataverse.harvard.edu/api/access/datafile
mkdir -p $R/favor
(for spec in "1:${DATASET_ITEM_06}" "2:${DATASET_ITEM_01}" "3:${DATASET_ITEM_02}" "4:${DATASET_ITEM_21}" "5:${DATASET_ITEM_11}" "6:${DATASET_ITEM_16}" "7:${DATASET_ITEM_05}" "8:${DATASET_ITEM_13}" "9:${DATASET_ITEM_67}" "10:${DATASET_ITEM_07}" "11:${DATASET_ITEM_17}" "12:${DATASET_ITEM_20}" "13:${DATASET_ITEM_03}" "14:${DATASET_ITEM_09}" "15:${DATASET_ITEM_15}" "16:${DATASET_ITEM_18}" "17:${DATASET_ITEM_10}" "18:${DATASET_ITEM_08}" "19:${DATASET_ITEM_14}" "20:${DATASET_ITEM_12}" "21:${DATASET_ITEM_19}" "22:${DATASET_ITEM_04}"; do
   N=${spec%%:*}; ID=${spec##*:}
   [ -e "$R/favor/chr${N}.done" ] && continue
   curl -sSL --retry 5 -m 86400 -o "$R/favor/chr${N}.tar.gz.part" "$DV/$ID" \
     && mv "$R/favor/chr${N}.tar.gz.part" "$R/favor/chr${N}.tar.gz" \
     && touch "$R/favor/chr${N}.done" && echo "[FAVOR chr$N] OK $(du -h $R/favor/chr${N}.tar.gz | cut -f1)" \
     || echo "[FAVOR chr$N] FAIL"
 done) &

get CATlas_meta https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184462/suppl/GSE184462_metadata.tsv.gz catlas 30000000 &
wait
echo "=== SUMMARY ==="
for d in gpnmsa pangolin favor catlas remap repeats mappability jaspar utrannotator targetscan gnomad_constraint cpg; do
  if compgen -G "$R/$d/.done" > /dev/null || compgen -G "$R/$d/*.done" > /dev/null; then
    printf '  HAVE %-20s %s\n' "$d" "$(du -sh $R/$d 2>/dev/null | cut -f1)"
  else printf '  PEND %-20s %s\n' "$d" "$(du -sh $R/$d 2>/dev/null | cut -f1)"; fi
done
echo PROC2_COMPLETE
