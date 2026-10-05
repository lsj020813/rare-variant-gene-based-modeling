#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
TAG=$1; GF=$2
W=${PROJECT_ROOT}/work; B=$W/run_l3b; T=tchl; N=19
export TMPDIR=$W/tmp TMP=$W/tmp TEMP=$W/tmp; mkdir -p "$TMPDIR" "$B/out"
S1=$W/ref/saige_step1_v4; V=$W/ref/band_vcf/chr$N.band.vcf.gz
O=$B/out/$TAG
for f in "$GF" "$V" "$V.csi" $S1/${T}_v4.rda $S1/${T}_v4.varianceRatio.txt; do
  [ -s "$f" ] || { echo "GATE FAIL: missing $f"; exit 2; }
done
[ -s "$O.done" ] && { echo "$TAG already complete"; exit 0; }
export UDOCKER_DIR=${CONTAINER_STATE_DIR}
UD=$(grep -m1 '^UD=' $W/run_bwg2/assoc.sh | cut -d= -f2- | tr -d '"')
[ -n "$UD" ] || { echo "GATE FAIL: could not read the udocker invocation"; exit 3; }
echo "[$TAG] start $(date '+%m-%d %H:%M:%S')"
eval nice -n 19 ionice -c3 env $UD step2_SPAtests.R \
  --vcfFile=\"$V\" --vcfFileIndex=\"$V.csi\" \
  --vcfField=\"DS\" --chrom=\"$N\" --AlleleOrder=ref-first \
  --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
  --GMMATmodelFile=\"$S1/${T}_v4.rda\" --varianceRatioFile=\"$S1/${T}_v4.varianceRatio.txt\" \
  --groupFile=\"$GF\" --annotation_in_groupTest=\"all\" --maxMAF_in_groupTest=0.01 \
  --SAIGEOutputFile=\"$O\" > "$O.log" 2>&1
rc=$?
exp=$(awk '$2=="var"{print $1}' "$GF" | sort -u | wc -l)
nrow=0; [ -s "$O" ] && nrow=$(( $(wc -l < "$O") - 1 ))
echo "[$TAG] end $(date '+%m-%d %H:%M:%S') rc=$rc rows=$nrow/$exp"
if [ "$rc" -ne 0 ] || [ "$nrow" -ne "$exp" ]; then
  echo "[$TAG] FAIL"; tail -3 "$O.log"; exit 1
fi
printf '{"status":"complete","tag":"%s","rows":%d,"groupfile_sha256":"%s"}\n' \
  "$TAG" "$nrow" "$(sha256sum "$GF" | cut -d" " -f1)" > "$O.done"
echo "[$TAG] OK"
