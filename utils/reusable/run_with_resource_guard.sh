#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -euo pipefail

ROOT="${ROOT:-${PROJECT_ROOT}}"
JOB_NAME="guarded_job"
MIN_AVAIL_MEM_GB=80
MAX_JOB_MEM_GB=120
MIN_DATA_FREE_GB=500
MIN_TMP_FREE_GB=50
MAX_LOAD_PER_CORE="1.50"
MAX_OTHER_USER_RSS_GB=32
MAX_TOTAL_OTHER_RSS_GB=64
MAX_OTHER_USER_CPU_PCT=400
MAX_TOTAL_OTHER_CPU_PCT=800
THREADS=8
TMPDIR_PATH="$ROOT/work/tmp/resource_guard"
REPORT=""
COMMAND=()

usage() {
  cat <<'EOF'
Usage:
  run_with_resource_guard.sh [options] -- COMMAND [ARGS...]

Options:
  --job-name NAME           Label used in reports.
  --min-avail-mem-gb N      Required available RAM before launch. Default: 80.
  --max-job-mem-gb N        Hard virtual-memory cap applied via ulimit. Default: 120.
  --min-data-free-gb N      Required free space under /data. Default: 500.
  --min-tmp-free-gb N       Required free space under TMPDIR filesystem. Default: 50.
  --max-load-per-core X     Refuse launch if load1 / nproc is above X. Default: 1.50.
  --max-other-user-rss-gb N Refuse if any other user is using >N GiB RSS. Default: 32.
  --max-total-other-rss-gb N
                            Refuse if all other users are using >N GiB RSS. Default: 64.
  --max-other-user-cpu-pct N
                            Refuse if any other user is using >N ps CPU percent. Default: 400.
  --max-total-other-cpu-pct N
                            Refuse if all other users are using >N ps CPU percent. Default: 800.
  --threads N              Export thread env vars and pass as safety metadata. Default: 8.
  --tmpdir PATH             TMPDIR to create/export. Default: ROOT/work/tmp/resource_guard.
  --report PATH             TSV report path. Default: reports/RESOURCE_GUARD.<job>.tsv.
  --check-only              Run checks and write report, but do not launch a command.
EOF
}

CHECK_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --min-avail-mem-gb) MIN_AVAIL_MEM_GB="$2"; shift 2 ;;
    --max-job-mem-gb) MAX_JOB_MEM_GB="$2"; shift 2 ;;
    --min-data-free-gb) MIN_DATA_FREE_GB="$2"; shift 2 ;;
    --min-tmp-free-gb) MIN_TMP_FREE_GB="$2"; shift 2 ;;
    --max-load-per-core) MAX_LOAD_PER_CORE="$2"; shift 2 ;;
    --max-other-user-rss-gb) MAX_OTHER_USER_RSS_GB="$2"; shift 2 ;;
    --max-total-other-rss-gb) MAX_TOTAL_OTHER_RSS_GB="$2"; shift 2 ;;
    --max-other-user-cpu-pct) MAX_OTHER_USER_CPU_PCT="$2"; shift 2 ;;
    --max-total-other-cpu-pct) MAX_TOTAL_OTHER_CPU_PCT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --tmpdir) TMPDIR_PATH="$2"; shift 2 ;;
    --report) REPORT="$2"; shift 2 ;;
    --check-only) CHECK_ONLY=1; shift ;;
    --help|-h) usage; exit 0 ;;
    --) shift; COMMAND=("$@"); break ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ "$CHECK_ONLY" != "1" && "${#COMMAND[@]}" -eq 0 ]]; then
  echo "No command provided. Use --check-only or pass -- COMMAND." >&2
  usage >&2
  exit 2
fi

mkdir -p "$TMPDIR_PATH"
chmod 700 "$TMPDIR_PATH"
REPORT="${REPORT:-$ROOT/reports/RESOURCE_GUARD.${JOB_NAME}.tsv}"
mkdir -p "$(dirname "$REPORT")"

mem_avail_kb="$(awk '/MemAvailable:/ {print $2}' /proc/meminfo)"
mem_total_kb="$(awk '/MemTotal:/ {print $2}' /proc/meminfo)"
mem_avail_gb="$((mem_avail_kb / 1024 / 1024))"
mem_total_gb="$((mem_total_kb / 1024 / 1024))"
nproc_value="$(nproc --all)"
load1="$(awk '{print $1}' /proc/loadavg)"
load_per_core="$(awk -v l="$load1" -v n="$nproc_value" 'BEGIN {printf "%.3f", l/n}')"
data_free_gb="$(df -BG /data | awk 'NR==2 {gsub("G","",$4); print $4}')"
tmp_free_gb="$(df -BG "$TMPDIR_PATH" | awk 'NR==2 {gsub("G","",$4); print $4}')"
slash_tmp_use_pct="$(df -P /tmp | awk 'NR==2 {gsub("%","",$5); print $5}')"
other_users="$(who 2>/dev/null | awk -v me="${USER:-}" '$1 != me {print $1}' | sort -u | paste -sd, - || true)"
other_usage_tsv="$(
  ps -eo user:32=,pcpu=,rss= 2>/dev/null \
    | awk -v me="${USER:-}" '
        $1 != me {
          cpu[$1] += $2
          rss[$1] += $3
        }
        END {
          for (u in rss) {
            printf "%s\t%.1f\t%.3f\n", u, cpu[u], rss[u] / 1024 / 1024
          }
        }
      ' \
    | sort || true
)"
other_usage_summary="$(
  printf '%s\n' "$other_usage_tsv" \
    | awk 'NF == 3 {printf "%s:cpu_pct=%.1f,rss_gb=%.3f;", $1, $2, $3}'
)"
other_usage_totals="$(
  printf '%s\n' "$other_usage_tsv" \
    | awk '
        NF == 3 {
          total_cpu += $2
          total_rss += $3
          if ($2 > max_cpu) {
            max_cpu = $2
            max_cpu_user = $1
          }
          if ($3 > max_rss) {
            max_rss = $3
            max_rss_user = $1
          }
        }
        END {
          if (max_cpu_user == "") {
            max_cpu_user = "NONE"
          }
          if (max_rss_user == "") {
            max_rss_user = "NONE"
          }
          printf "%.1f\t%.3f\t%s\t%.1f\t%s\t%.3f\n", total_cpu, total_rss, max_cpu_user, max_cpu, max_rss_user, max_rss
        }
      '
)"
read -r total_other_cpu_pct total_other_rss_gb max_other_cpu_user max_other_cpu_pct max_other_rss_user max_other_rss_gb <<< "$other_usage_totals"
gpu_summary="NA"
if command -v nvidia-smi >/dev/null 2>&1; then
  gpu_summary="$(nvidia-smi --query-gpu=index,name,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null | tr '\n' ';' || echo NA)"
fi

state="PASS"
detail="ok"
fail() {
  state="FAIL"
  if [[ "$detail" == "ok" ]]; then
    detail="$1"
  else
    detail="$detail; $1"
  fi
}

if (( mem_avail_gb < MIN_AVAIL_MEM_GB )); then
  fail "available_mem_gb=${mem_avail_gb}<${MIN_AVAIL_MEM_GB}"
fi
if (( MAX_JOB_MEM_GB >= mem_total_gb - 20 )); then
  fail "max_job_mem_gb=${MAX_JOB_MEM_GB} leaves less than 20GB OS headroom on ${mem_total_gb}GB host"
fi
if (( data_free_gb < MIN_DATA_FREE_GB )); then
  fail "data_free_gb=${data_free_gb}<${MIN_DATA_FREE_GB}"
fi
if (( tmp_free_gb < MIN_TMP_FREE_GB )); then
  fail "tmpdir_free_gb=${tmp_free_gb}<${MIN_TMP_FREE_GB}"
fi
awk -v got="$load_per_core" -v max="$MAX_LOAD_PER_CORE" 'BEGIN {exit !(got > max)}' && fail "load_per_core=${load_per_core}>${MAX_LOAD_PER_CORE}"
awk -v got="${max_other_rss_gb:-0}" -v max="$MAX_OTHER_USER_RSS_GB" 'BEGIN {exit !(got > max)}' && fail "max_other_user_rss_gb=${max_other_rss_user:-unknown}:${max_other_rss_gb}>${MAX_OTHER_USER_RSS_GB}"
awk -v got="${total_other_rss_gb:-0}" -v max="$MAX_TOTAL_OTHER_RSS_GB" 'BEGIN {exit !(got > max)}' && fail "total_other_rss_gb=${total_other_rss_gb}>${MAX_TOTAL_OTHER_RSS_GB}"
awk -v got="${max_other_cpu_pct:-0}" -v max="$MAX_OTHER_USER_CPU_PCT" 'BEGIN {exit !(got > max)}' && fail "max_other_user_cpu_pct=${max_other_cpu_user:-unknown}:${max_other_cpu_pct}>${MAX_OTHER_USER_CPU_PCT}"
awk -v got="${total_other_cpu_pct:-0}" -v max="$MAX_TOTAL_OTHER_CPU_PCT" 'BEGIN {exit !(got > max)}' && fail "total_other_cpu_pct=${total_other_cpu_pct}>${MAX_TOTAL_OTHER_CPU_PCT}"

tmp_report="${REPORT}.tmp.$$"
{
  printf 'key\tvalue\n'
  printf 'timestamp\t%s\n' "$(date -Is)"
  printf 'job_name\t%s\n' "$JOB_NAME"
  printf 'state\t%s\n' "$state"
  printf 'detail\t%s\n' "$detail"
  printf 'host\t%s\n' "$(hostname)"
  printf 'user\t%s\n' "${USER:-unknown}"
  printf 'mem_total_gb\t%s\n' "$mem_total_gb"
  printf 'mem_available_gb\t%s\n' "$mem_avail_gb"
  printf 'min_available_mem_gb\t%s\n' "$MIN_AVAIL_MEM_GB"
  printf 'max_job_mem_gb\t%s\n' "$MAX_JOB_MEM_GB"
  printf 'nproc\t%s\n' "$nproc_value"
  printf 'threads\t%s\n' "$THREADS"
  printf 'load1\t%s\n' "$load1"
  printf 'load_per_core\t%s\n' "$load_per_core"
  printf 'max_load_per_core\t%s\n' "$MAX_LOAD_PER_CORE"
  printf 'total_other_cpu_pct\t%s\n' "${total_other_cpu_pct:-0.0}"
  printf 'max_other_user_cpu_pct\t%s\n' "${max_other_cpu_user:-NONE}:${max_other_cpu_pct:-0.0}"
  printf 'max_other_user_cpu_pct_limit\t%s\n' "$MAX_OTHER_USER_CPU_PCT"
  printf 'total_other_cpu_pct_limit\t%s\n' "$MAX_TOTAL_OTHER_CPU_PCT"
  printf 'total_other_rss_gb\t%s\n' "${total_other_rss_gb:-0.000}"
  printf 'max_other_user_rss_gb\t%s\n' "${max_other_rss_user:-NONE}:${max_other_rss_gb:-0.000}"
  printf 'max_other_user_rss_gb_limit\t%s\n' "$MAX_OTHER_USER_RSS_GB"
  printf 'total_other_rss_gb_limit\t%s\n' "$MAX_TOTAL_OTHER_RSS_GB"
  printf 'other_user_usage_summary\t%s\n' "${other_usage_summary:-NONE}"
  printf 'data_free_gb\t%s\n' "$data_free_gb"
  printf 'tmpdir\t%s\n' "$TMPDIR_PATH"
  printf 'tmpdir_free_gb\t%s\n' "$tmp_free_gb"
  printf 'slash_tmp_use_pct\t%s\n' "$slash_tmp_use_pct"
  printf 'other_logged_in_users\t%s\n' "${other_users:-NONE}"
  printf 'gpu_summary\t%s\n' "$gpu_summary"
  printf 'command\t%s\n' "${COMMAND[*]:-CHECK_ONLY}"
} > "$tmp_report"
mv "$tmp_report" "$REPORT"
chmod 600 "$REPORT"
cat "$REPORT"

if [[ "$state" != "PASS" ]]; then
  exit 75
fi
if [[ "$CHECK_ONLY" == "1" ]]; then
  exit 0
fi

export TMPDIR="$TMPDIR_PATH"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"
export BLIS_NUM_THREADS="$THREADS"

max_job_mem_kb="$((MAX_JOB_MEM_GB * 1024 * 1024))"
ulimit -Sv "$max_job_mem_kb"
ulimit -Hv "$max_job_mem_kb" 2>/dev/null || true

if command -v /usr/bin/time >/dev/null 2>&1; then
  exec /usr/bin/time -v "${COMMAND[@]}"
fi
exec "${COMMAND[@]}"
