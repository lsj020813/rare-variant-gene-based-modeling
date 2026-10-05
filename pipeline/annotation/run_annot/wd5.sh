#!/usr/bin/env bash
set -u
PATTERN="${1:?pkill pattern required}"
WD="${2:-.}"
LOAD_MAX="${3:-38}"
RSS_BUDGET_GB="${4:-0}"
LOG="$WD/watchdog.log"
PROBE="$WD/.wd_probe.$$"
HEART="$WD/.wd_heartbeat"
PEAK="$WD/peak_rss.txt"
mkdir -p "$WD"

MEM_AVAIL_MIN_GB=40
SWAP_SO_STREAK_MAX=2
IOWAIT_MAX=25
WRITE_LAT_MULT=10
WRITE_IOWAIT_MIN=10
BASELINE_SAMPLES=5
PREDICT_SAMPLES=3
BREACH_KILL=1
OUR_IO_MIN_MBS=100
SWAP_AVAIL_MAX_GB=60
INTERVAL=5

if [ "$RSS_BUDGET_GB" -eq 0 ]; then
  a=$(free -g | awk 'NR==2{print $7}')
  RSS_BUDGET_GB=$(( a / 2 ))
fi

baseline_ms=""; bl_buf=""; so_streak=0; breach_streak=0; n=0; prev_io_b=""; our_io_mbs=0; our_cores=0
peak_sum=0; peak_max=0; prev_max=0
log(){ echo "$(date '+%F %T') $*" >> "$LOG"; }
log "WATCHDOG v2 START pattern='$PATTERN' LOAD_MAX=$LOAD_MAX RSS_BUDGET=${RSS_BUDGET_GB}G MEM_MIN=${MEM_AVAIL_MIN_GB}G interval=${INTERVAL}s"

while true; do
  n=$((n+1))
  date +%s > "$HEART"

  if ! pgrep -f "$PATTERN" >/dev/null 2>&1; then
    log "target gone after $n samples — peak_sum=${peak_sum}G peak_max=${peak_max}G — watchdog exit clean"
    printf 'peak_sum_gb=%s\npeak_max_gb=%s\nsamples=%s\n' "$peak_sum" "$peak_max" "$n" > "$PEAK"
    rm -f "$PROBE"; exit 0
  fi

  read sum_gb max_gb <<EOF
$(ps -o rss= -p $(pgrep -d, -f "$PATTERN") 2>/dev/null | awk '{s+=$1; if($1>m) m=$1} END{printf "%.1f %.1f", s/1048576, m/1048576}')
EOF
  sum_gb=${sum_gb:-0}; max_gb=${max_gb:-0}
  io_b=$(for p in $(pgrep -f "$PATTERN"); do cat /proc/$p/io 2>/dev/null | awk '/^(read_bytes|write_bytes)/{s+=$2} END{print s+0}'; done | awk '{s+=$1} END{print s+0}')
  if [ -n "$prev_io_b" ]; then our_io_mbs=$(awk -v a="$io_b" -v b="$prev_io_b" -v i="$INTERVAL" 'BEGIN{d=a-b; if(d<0)d=0; printf "%.0f", d/1048576/i}'); fi
  prev_io_b=$io_b
  our_cores=$(ps -o pcpu= -p $(pgrep -d, -f "$PATTERN") 2>/dev/null | awk '{s+=$1} END{printf "%.0f", s/100}')
  awk -v a="$sum_gb" -v b="$peak_sum" 'BEGIN{exit !(a>b)}' && peak_sum=$sum_gb
  awk -v a="$max_gb" -v b="$peak_max" 'BEGIN{exit !(a>b)}' && peak_max=$max_gb

  t0=$(date +%s%N)
  dd if=/dev/zero of="$PROBE" bs=1M count=10 oflag=direct >/dev/null 2>&1
  t1=$(date +%s%N); write_ms=$(( (t1-t0)/1000000 )); rm -f "$PROBE"

  load1=$(cut -d' ' -f1 /proc/loadavg); load1i=${load1%.*}
  avail_gb=$(free -g | awk 'NR==2{print $7}')
  vm=$(vmstat 1 2 | tail -1)
  so=$(echo "$vm" | awk '{print $8}'); wa=$(echo "$vm" | awk '{print $16}')

  if [ -z "$baseline_ms" ]; then
    bl_buf="$bl_buf $write_ms"
    if [ "$(echo $bl_buf | wc -w)" -ge "$BASELINE_SAMPLES" ]; then
      baseline_ms=$(echo $bl_buf | tr ' ' '\n' | sort -n | awk '{a[NR]=$1} END{print a[int((NR+1)/2)]}')
      log "baseline_ms=$baseline_ms (median of $BASELINE_SAMPLES:$bl_buf)"
    fi
  fi
  [ "${so:-0}" -gt 0 ] && so_streak=$((so_streak+1)) || so_streak=0

  breach=0; reason=""
  [ "${avail_gb:-999}" -lt "$MEM_AVAIL_MIN_GB" ] && { breach=1; reason="$reason mem_avail=${avail_gb}G<${MEM_AVAIL_MIN_GB}G"; }
  [ "$so_streak" -ge "$SWAP_SO_STREAK_MAX" ] && [ "${avail_gb:-999}" -lt "$SWAP_AVAIL_MAX_GB" ] && { breach=1; reason="$reason swap_out_streak=$so_streak(avail=${avail_gb}G)"; }
  [ "${our_cores:-0}" -gt "$LOAD_MAX" ]          && { breach=1; reason="$reason our_cores=$our_cores>$LOAD_MAX"; }
  [ "${wa:-0}" -gt "$IOWAIT_MAX" ] && [ "${our_io_mbs:-0}" -ge "$OUR_IO_MIN_MBS" ] && { breach=1; reason="$reason iowait=${wa}%(our_io=${our_io_mbs}MB/s)"; }
  awk -v s="$sum_gb" -v b="$RSS_BUDGET_GB" 'BEGIN{exit !(s>b)}' \
     && { breach=1; reason="$reason worker_rss_sum=${sum_gb}G>${RSS_BUDGET_GB}G"; }

  if awk -v c="$max_gb" -v p="$prev_max" 'BEGIN{exit !(c>p && p>0)}'; then
    proj=$(awk -v c="$max_gb" -v p="$prev_max" -v k="$PREDICT_SAMPLES" 'BEGIN{printf "%.1f", c+(c-p)*k}')
    if awk -v pr="$proj" -v b="$RSS_BUDGET_GB" -v c="$max_gb" 'BEGIN{exit !(pr>b && c>=0.9*b)}'; then
      breach=1; reason="$reason rss_trend_proj=${proj}G>${RSS_BUDGET_GB}G(now=${max_gb}G)"
    fi
  fi
  prev_max=$max_gb

  if [ -n "$baseline_ms" ] && [ "$write_ms" -gt $((baseline_ms * WRITE_LAT_MULT + 500)) ] \
     && [ "${wa:-0}" -ge "$WRITE_IOWAIT_MIN" ] && [ "${our_io_mbs:-0}" -ge "$OUR_IO_MIN_MBS" ]; then
    breach=1; reason="$reason write_ms=$write_ms(base=$baseline_ms,wa=${wa}%)"
  fi

  log "n=$n load=$load1 our_cores=$our_cores our_io=${our_io_mbs}MB/s avail=${avail_gb}G rss_sum=${sum_gb}G rss_max=${max_gb}G so=$so wa=${wa}% write_ms=$write_ms breach=$breach$reason"

  if [ "$breach" -eq 1 ]; then
    breach_streak=$((breach_streak+1))
    if [ "$breach_streak" -ge "$BREACH_KILL" ]; then
      log "*** KILL TRIGGERED (streak=$breach_streak):$reason"
      pkill -9 -f "$PATTERN" 2>>"$LOG"
      sleep 2
      log "*** KILL DONE — survivors=$(pgrep -f "$PATTERN" | wc -l) peak_sum=${peak_sum}G peak_max=${peak_max}G"
      printf 'peak_sum_gb=%s\npeak_max_gb=%s\nkilled=1\nreason=%s\n' "$peak_sum" "$peak_max" "$reason" > "$PEAK"
      rm -f "$PROBE"; exit 9
    fi
  else
    breach_streak=0
  fi

  sleep "$INTERVAL"
done
