#!/usr/bin/env bash
set -u
PATTERN="${1:?}"; WD="${2:?}"; LOAD="${3:?}"; RSS="${4:?}"; WAITS="${5:?}"; CHAIN="${6:?}"
LOG="$WD/wd_attach.log"
while ! grep -qE 'CHAIN_DONE|CHAIN ABORTED' "$CHAIN" 2>/dev/null; do
  bash "$WD/wd_attach.sh" "$PATTERN" "$WD" "$LOAD" "$RSS" "$WAITS"
  sleep 10
done
echo "$(date '+%Y-%m-%d %H:%M:%S') wd_keep: chain finished, watchdog keeper stopping" >> "$LOG"
