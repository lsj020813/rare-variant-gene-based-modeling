import io
p = "run_arm_seq.sh"
s = open(p).read()
old = (
    'L1=$(cut -d\' \' -f1 /proc/loadavg)\n'
    'if [ "$(echo "$L1 > $MAXSTART" | bc -l)" = "1" ]; then\n'
    '  echo "[$TAG] HOLD: 1-min load $L1 exceeds start ceiling $MAXSTART"; exit 3\n'
    'fi\n'
)
new = (
    'L1=$(cut -d\' \' -f1 /proc/loadavg)\n'
    'WAITED=0\n'
    'while [ "$(echo "$L1 > $MAXSTART" | bc -l)" = "1" ]; do\n'
    '  echo "[$TAG] WAIT $(date \'+%H:%M:%S\') 1-min load $L1 above start ceiling $MAXSTART '
    '(waited ${WAITED}s)"\n'
    '  sleep 60; WAITED=$((WAITED+60)); L1=$(cut -d\' \' -f1 /proc/loadavg)\n'
    '  if [ "$WAITED" -ge 7200 ]; then echo "[$TAG] HOLD: load stayed above $MAXSTART for 2h"; '
    'exit 3; fi\n'
    'done\n'
    '[ "$WAITED" -gt 0 ] && echo "[$TAG] WAIT_END $(date \'+%H:%M:%S\') load $L1 after ${WAITED}s"\n'
)
assert s.count(old) == 1, s.count(old)
open(p, "w").write(s.replace(old, new))
print("PATCH_OK")
