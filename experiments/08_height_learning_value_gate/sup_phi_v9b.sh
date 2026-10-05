#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u
SOURCE_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WATCHDOG="$SOURCE_DIR/../07_attention_based_variant_interactions/infra/wd_prs_big.sh"
B=${PROJECT_ROOT}/work/phi_gate; P=${PROJECT_ROOT}/work/prs; PY=python3
WDLOG=$P/logs/wd_phigate9.log; SL=$B/logs/sup_v9.log; cd $B
say(){ echo "$(date -Is) $*" >> $SL; }
pid_of(){ pgrep -f "^$PY final_phi_gate_v9[.]py chr $1\$"; }
wd_up(){
  pgrep -f 'wd_prs_big.sh .*phigate9' >/dev/null && return 0
  tmux kill-session -t wd_phigate9 2>/dev/null
  tmux new-session -d -s wd_phigate9 "bash '$WATCHDOG' '[f]inal_phi_gate_v9[.]py' phigate9"; sleep 6
  pgrep -f 'wd_prs_big.sh .*phigate9' >/dev/null && [ -s $WDLOG ] && { say "watchdog (re)started"; return 0; }
  say "WATCHDOG START FAILED -> abort (no heavy job without watchdog)"; exit 9; }
quarantine(){ N=$1; [ -s out/chr$N/chr.done ] && return 0
  Q=_v9_aborted_$(date +%Y%m%dT%H%M%S)_chr$N; mkdir -p $Q
  [ -e out/chr$N ] && mv out/chr$N $Q/out_chr$N; [ -e private/chr$N ] && mv private/chr$N $Q/private_chr$N
  say "quarantined partial chr$N -> $Q"; }
run_one(){ N=$1; [ -s out/chr$N/chr.done ] && { say "chr$N already done, skip"; return 0; }
  quarantine $N; wd_up
  say "SEQ start chr $N"
  setsid $PY final_phi_gate_v9.py chr $N > logs/v9_seq_chr$N.log 2>&1; rc=$?
  echo "exit=$rc $(date -Is)" >> logs/v9_seq_chr$N.log; say "SEQ end chr $N rc=$rc"
  [ $rc -eq 0 ] || { say "chr $N failed in sequential mode -> stop (manual review)"; exit 1; }; }
fallback(){ why="$1"; say "FALLBACK to sequential: $why"
  tmux kill-session -t phigate9_chain 2>/dev/null; say "chain session stopped"
  if [ "$why" = "long_stop" ]; then
    p=$(pid_of 2); [ -n "$p" ] && { kill -CONT $p; kill -TERM $p; sleep 15; kill -KILL $p 2>/dev/null; say "chr2 pid $p terminated to relieve load; chr12 continues"; }
    while [ -n "$(pid_of 12)" ]; do sleep 60; done
    rc12=$(grep -o 'exit=[0-9]*' logs/v9_chr12.log | tail -1 | cut -d= -f2); say "concurrent chr12 finished rc=$rc12"
    [ "$rc12" = "0" ] || quarantine 12
  else
    for N in 12 2; do while [ -n "$(pid_of $N)" ]; do sleep 10; done; done
  fi
  run_one 12; run_one 2
  wd_up; say "SEQ start pool"
  $PY final_phi_gate_v9.py pool 2 12 > logs/v9_pool.log 2>&1; rc=$?; say "SEQ end pool rc=$rc"
  [ $rc -eq 0 ] && echo "chain ok (sequential fallback) $(date -Is)" > out/chain_v9.done
  exit $rc; }
say "supervisor start (concurrent mode)"
KILLS0=$(grep -c 'KILL TRIGGERED' $WDLOG 2>/dev/null); KILLS0=${KILLS0:-0}
while true; do
  [ -s out/chain_v9.done ] && { say "chain done -> supervisor exit"; exit 0; }
  grep -qE 'CHAIN STOPPED|END pool rc=[1-9]' logs/chain_v9.log 2>/dev/null && { say "chain reported genuine failure -> no retry, exit"; exit 1; }
  K=$(grep -c 'KILL TRIGGERED' $WDLOG 2>/dev/null); K=${K:-0}
  [ "$K" -gt "$KILLS0" ] && fallback killed
  if [ "$(tail -60 $WDLOG | grep -c 'stopped=1')" -ge 60 ] && [ -n "$(pid_of 2)" ] && [ -n "$(pid_of 12)" ]; then fallback long_stop; fi
  sleep 60
done
