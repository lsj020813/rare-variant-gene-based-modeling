#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
exec 9>${PROJECT_ROOT}/work/run_proc4/.proc4.lock
flock -n 9 || { echo "another run"; exit 0; }
R=${PROJECT_ROOT}/work/ref

size_of(){ curl -sL -m 40 -r 0-0 -D - -o /dev/null "$1" 2>/dev/null \
           | grep -i '^content-range' | tail -1 | tr -d '\r' | sed 's#.*/##'; }

fetch_norsm(){
  local URL=$1 OUT=$2 LBL=$3
  local WANT=$(size_of "$URL"); case "$WANT" in ''|*[!0-9]*) echo "[$LBL] FAIL size unknown"; return 1;; esac
  for i in $(seq 1 8); do
    rm -f "$OUT.part"
    curl -sSL --retry 10 --retry-delay 15 -m 86400 -o "$OUT.part" "$URL" || true
    local S=$(stat -c %s "$OUT.part" 2>/dev/null || echo 0)
    if [ "$S" -eq "$WANT" ]; then mv "$OUT.part" "$OUT"; echo "[$LBL] COMPLETE_EXACT $S"; return 0; fi
    echo "[$LBL] try$i got $S want $WANT — restart from zero"
  done
  echo "[$LBL] GIVEUP"; return 1
}

fetch_rsm(){
  local URL=$1 OUT=$2 LBL=$3
  local WANT=$(size_of "$URL"); case "$WANT" in ''|*[!0-9]*) echo "[$LBL] FAIL size unknown"; return 1;; esac
  for i in $(seq 1 40); do
    local H=$(stat -c %s "$OUT.part" 2>/dev/null || echo 0)
    if [ "$H" -eq "$WANT" ]; then mv "$OUT.part" "$OUT"; echo "[$LBL] COMPLETE_EXACT $H"; return 0; fi
    if [ "$H" -gt "$WANT" ]; then echo "[$LBL] OVERSIZE $H>$WANT — discarding (append corruption)"; rm -f "$OUT.part"; continue; fi
    curl -sSL -C - --retry 5 --retry-delay 15 -m 86400 -o "$OUT.part" "$URL" || true
  done
  echo "[$LBL] GIVEUP"; return 1
}

mkdir -p $R/pangolin
( [ -e $R/pangolin/.done ] || { rm -f $R/pangolin/content.part $R/pangolin/P.part
    fetch_norsm https://zenodo.org/api/records/15649338/files/Pangolin_hg38_snvs_masked.zip/content \
                $R/pangolin/Pangolin_hg38_snvs_masked.zip Pangolin && touch $R/pangolin/.done; } ) &

( declare -A ID=( [1]=${DATASET_ITEM_06} [2]=${DATASET_ITEM_01} [3]=${DATASET_ITEM_02} [4]=${DATASET_ITEM_21} [5]=${DATASET_ITEM_11} [6]=${DATASET_ITEM_16}
                  [7]=${DATASET_ITEM_05} [8]=${DATASET_ITEM_13} [9]=${DATASET_ITEM_67} [10]=${DATASET_ITEM_07} [11]=${DATASET_ITEM_17} [12]=${DATASET_ITEM_20}
                  [13]=${DATASET_ITEM_03} [14]=${DATASET_ITEM_09} [15]=${DATASET_ITEM_15} [16]=${DATASET_ITEM_18} [17]=${DATASET_ITEM_10} [18]=${DATASET_ITEM_08}
                  [19]=${DATASET_ITEM_14} [20]=${DATASET_ITEM_12} [21]=${DATASET_ITEM_19} [22]=${DATASET_ITEM_04} )
  for N in 18 17 15 14 13 9 10 12 11 8 7 6 5 4 3 1 2; do
    [ -e "$R/favor/chr${N}.done" ] && continue
    fetch_rsm https://dataverse.harvard.edu/api/access/datafile/${ID[$N]} \
              $R/favor/chr${N}.tar.gz "FAVOR chr$N" && touch $R/favor/chr${N}.done
  done ) &
wait
echo "=== SUMMARY (marker-based only) ==="
for d in pangolin favor gpnmsa remap; do
  printf '  %-10s %-8s markers=%s\n' "$d" "$(du -sh $R/$d 2>/dev/null | cut -f1)" \
         "$(ls $R/$d/.done $R/$d/*.done 2>/dev/null | wc -l)"
done
echo PROC4_COMPLETE
