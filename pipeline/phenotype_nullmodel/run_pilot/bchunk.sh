#!/usr/bin/env bash
set -uo pipefail
: "${PROJECT_ROOT:?Set PROJECT_ROOT to a generic workspace root}"
RUN_DIR=${RUN_DIR:-"${PROJECT_ROOT}/work/run_pilot"}
mkdir -p "$RUN_DIR"
exec 9>"$RUN_DIR/.lock3"
flock -n 9 || { echo "another instance running"; exit 0; }

R=${REFERENCE_DIR:-"${PROJECT_ROOT}/work/ref"}
S1=$R/saige_step1_v4
GP=$R/groupfiles_pilot
BAND=$R/band_vcf
OUT=$R/saige_step2_pilot
CHK=$GP/chunks
mkdir -p "$OUT" "$CHK"
SAIGE_RUNNER=${SAIGE_RUNNER:-Rscript}
SAIGE_STEP2_SCRIPT=${SAIGE_STEP2_SCRIPT:?Set SAIGE_STEP2_SCRIPT}
NCHUNK=${NCHUNK:-6}
SRC=$GP/chr19.B_3kb_re2g.txt

[ -s "$SRC" ] || { echo "PREFLIGHT FAIL missing $SRC"; exit 3; }
NG=$(awk '$2=="var"{n++} END{print n}' "$SRC")
PER=$(( (NG + NCHUNK - 1) / NCHUNK ))
rm -f $CHK/B.part*.txt
awk -v per=$PER -v d="$CHK" '
  $2=="var" { g++; part = int((g-1)/per) + 1 }
  { printf "%s\n", $0 >> sprintf("%s/B.part%03d.txt", d, part) }
' "$SRC"

SG=0; SV=0
for f in $CHK/B.part*.txt; do
  SG=$(( SG + $(awk '$2=="var"{n++} END{print n+0}' $f) ))
  SV=$(( SV + $(awk '$2=="var"{for(i=3;i<=NF;i++) n++} END{print n+0}' $f) ))
done
TV=$(awk '$2=="var"{for(i=3;i<=NF;i++) n++} END{print n+0}' "$SRC")
echo "split: $(ls $CHK/B.part*.txt | wc -l) chunks  genes $SG/$NG  var-slots $SV/$TV"
[ "$SG" -eq "$NG" ] || { echo "GATE FAIL: genes $SG != $NG"; exit 6; }
[ "$SV" -eq "$TV" ] || { echo "GATE FAIL: var slots $SV != $TV"; exit 6; }
for f in $CHK/B.part*.txt; do
  a=$(awk '$2=="var"{n++} END{print n+0}' $f); b=$(awk '$2=="anno"{n++} END{print n+0}' $f)
  [ "$a" -eq "$b" ] || { echo "GATE FAIL: $f var=$a anno=$b"; exit 6; }
  lbl=$(awk '$2=="anno"{print $3; exit}' $f)
  [ "$lbl" = "all" ] || { echo "GATE FAIL: $f label=$lbl (expect all)"; exit 6; }
done
echo "split gates OK"

for p in $(pgrep -f 'groupFile.*B_3kb_re2g'); do kill -TERM $p 2>/dev/null; done
sleep 3

run_chunk() {
  local P=$1
  local G=$CHK/$P.txt
  local O=$OUT/tchl.chr19.B_3kb_re2g.$P
  [ -s "$O.done" ] && { echo "  skip $P"; return 0; }
  "$SAIGE_RUNNER" "$SAIGE_STEP2_SCRIPT" \
    --vcfFile="$BAND/chr19.band.vcf.gz" --vcfFileIndex="$BAND/chr19.band.vcf.gz.csi" \
    --vcfField="DS" --chrom="19" --AlleleOrder=ref-first \
    --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=TRUE \
    --GMMATmodelFile="$S1/tchl_v4.rda" --varianceRatioFile="$S1/tchl_v4.varianceRatio.txt" \
    --groupFile="$G" --annotation_in_groupTest="all" --maxMAF_in_groupTest=0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --SAIGEOutputFile="$O" > "$O.log" 2>&1
  local rc=$?
  local exp=$(awk '$2=="var"{print $1}' "$G" | sort -u | wc -l)
  local nrow=0; [ -s "$O" ] && nrow=$(( $(wc -l < "$O") - 1 ))
  if [ "$rc" -ne 0 ] || [ "$nrow" -le 0 ]; then
    echo "  FAIL $P rc=$rc rows=$nrow/$exp"; tail -3 "$O.log"; return 1
  fi
  touch "$O.done"
  [ "$nrow" -eq "$exp" ] && echo "  OK $P $nrow" || echo "  PARTIAL $P $nrow/$exp"
}
export -f run_chunk; export CHK BAND OUT S1 SAIGE_RUNNER SAIGE_STEP2_SCRIPT

ls $CHK/B.part*.txt | sed 's@.*/@@; s@\.txt$@@' | xargs -P "$NCHUNK" -I{} bash -c 'run_chunk {}'
n=$(ls $OUT/tchl.chr19.B_3kb_re2g.B.part*.done 2>/dev/null | wc -l)
echo "chunks done: $n/$(ls $CHK/B.part*.txt | wc -l)"
echo "B_CHUNKED_COMPLETE"
