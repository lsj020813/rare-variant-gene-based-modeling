#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -o pipefail
BASE=${PROJECT_ROOT}/work/fset/settest
PREP=$BASE/prep; SUB=$BASE/sub; LOG=$BASE/logs
mkdir -p "$SUB" "$LOG"
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
BT=bcftools
ORIG=${PROJECT_ROOT}/work/ref/orig_index

one() {
  local N=$1
  local O=$SUB/chr$N.sub.vcf.gz
  [ -s "$O.done" ] && { echo "skip chr$N"; return 0; }
  local t0=$(date +%s)
  nice -n 19 ionice -c3 $BT view -R "$PREP/chr$N.pos.bed" -Ou "$ORIG/chr$N.vcf.gz" 2>"$LOG/ext_chr$N.err" \
    | nice -n 19 ionice -c3 $BT annotate -x FORMAT/GP -Oz -o "$O.tmp" --threads 2 2>>"$LOG/ext_chr$N.err"
  local rc=$?
  if [ $rc -ne 0 ]; then echo "FAIL chr$N rc=$rc"; return 1; fi
  local eof=$(tail -c 28 "$O.tmp" | od -An -tx1 | tr -d ' \n')
  if [ "$eof" != "1f8b08040000000000ff0600424302001b0003000000000000000000" ]; then
    echo "FAIL chr$N no-EOF"; mv "$O.tmp" "$O.trunc"; return 1; fi
  mv "$O.tmp" "$O"
  nice -n 19 ionice -c3 $BT index -f --csi "$O" || { echo "FAIL chr$N index"; return 1; }
  local nrec=$($BT index -n "$O")
  local ntgt=$(wc -l < "$PREP/chr$N.pos.bed")
  local t1=$(date +%s)
  echo "chr$N nrec=$nrec npos=$ntgt secs=$((t1-t0))" > "$O.done"
  echo "OK chr$N nrec=$nrec npos=$ntgt secs=$((t1-t0))"
}
export -f one; export BASE PREP SUB LOG BT ORIG

P=${P:-4}
seq 22 -1 1 | xargs -P "$P" -n1 bash -c 'one "$@"' _
nd=$(ls $SUB/*.done 2>/dev/null | wc -l)
echo "extract_done: $nd / 22"
[ "$nd" -eq 22 ] && echo EXTRACT_COMPLETE > $BASE/extract.done
