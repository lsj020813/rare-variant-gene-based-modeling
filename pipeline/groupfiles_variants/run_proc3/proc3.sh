#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_proc3/.proc3.lock
flock -n 9 || { echo "another run"; exit 0; }
R=${PROJECT_ROOT}/work/ref

fetch(){
  local URL=$1 OUT=$2 LBL=$3 tries=0
  local WANT=$(curl -sIL -m 30 "$URL" | grep -i '^content-length' | tail -1 | tr -d '\r' | awk '{print $2}')
  [ -z "$WANT" ] && { echo "[$LBL] FAIL no content-length"; return 1; }
  while :; do
    local HAVE=$(stat -c %s "$OUT" 2>/dev/null || echo 0)
    if [ "$HAVE" -ge "$WANT" ]; then echo "[$LBL] COMPLETE $((HAVE/1073741824))G"; return 0; fi
    tries=$((tries+1))
    [ "$tries" -gt 40 ] && { echo "[$LBL] GIVEUP after $tries tries ($HAVE/$WANT)"; return 1; }
    echo "[$LBL] resume try $tries at $((HAVE/1048576))M / $((WANT/1048576))M"
    curl -sSL -C - --retry 5 --retry-delay 15 -m 86400 -o "$OUT" "$URL" || true
    sleep 5
  done
}

mkdir -p $R/pangolin
( fetch https://zenodo.org/api/records/15649338/files/Pangolin_hg38_snvs_masked.zip/content \
        $R/pangolin/content.part Pangolin \
    && mv $R/pangolin/content.part $R/pangolin/Pangolin_hg38_snvs_masked.zip \
    && touch $R/pangolin/.done ) &

mkdir -p $R/favor
( declare -A ID=( [1]=${DATASET_ITEM_06} [2]=${DATASET_ITEM_01} [3]=${DATASET_ITEM_02} [4]=${DATASET_ITEM_21} [5]=${DATASET_ITEM_11} [6]=${DATASET_ITEM_16}
                  [7]=${DATASET_ITEM_05} [8]=${DATASET_ITEM_13} [9]=${DATASET_ITEM_67} [10]=${DATASET_ITEM_07} [11]=${DATASET_ITEM_17} [12]=${DATASET_ITEM_20}
                  [13]=${DATASET_ITEM_03} [14]=${DATASET_ITEM_09} [15]=${DATASET_ITEM_15} [16]=${DATASET_ITEM_18} [17]=${DATASET_ITEM_10} [18]=${DATASET_ITEM_08}
                  [19]=${DATASET_ITEM_14} [20]=${DATASET_ITEM_12} [21]=${DATASET_ITEM_19} [22]=${DATASET_ITEM_04} )
  for N in 21 22 19 20 18 16 17 15 14 13 9 10 12 11 8 7 6 5 4 3 1 2; do
    [ -e "$R/favor/chr${N}.done" ] && { echo "[FAVOR chr$N] skip"; continue; }
    fetch https://dataverse.harvard.edu/api/access/datafile/${ID[$N]} \
          $R/favor/chr${N}.tar.gz.part "FAVOR chr$N" \
      && mv $R/favor/chr${N}.tar.gz.part $R/favor/chr${N}.tar.gz \
      && touch $R/favor/chr${N}.done
  done ) &
wait
echo "=== SUMMARY ==="
printf '  gpnmsa   %s\n' "$(du -sh $R/gpnmsa 2>/dev/null | cut -f1)"
printf '  pangolin %s (done=%s)\n' "$(du -sh $R/pangolin 2>/dev/null | cut -f1)" "$([ -e $R/pangolin/.done ] && echo Y || echo N)"
printf '  favor    %s (%s/22 chroms)\n' "$(du -sh $R/favor 2>/dev/null | cut -f1)" "$(ls $R/favor/*.done 2>/dev/null | wc -l)"
echo PROC3_COMPLETE
