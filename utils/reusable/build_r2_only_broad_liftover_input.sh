#!/usr/bin/env bash
: "${DATA_ROOT:?Set DATA_ROOT}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -euo pipefail

ROOT="${PROJECT_ROOT}"
RAW=""
CHROM="22"
LABEL=""
THREADS="${LIFTOVER_THREADS:-12}"
R2_MIN="0.8"
OUT_DIR=""
BCF_DIR=""
REPORT=""
DRY_RUN="NO"
ALLOW_HEAVY="${ALLOW_HEAVY:-NO}"
SORT_TMP_DIR=""

BCFTOOLS="${BCFTOOLS:-bcftools}"
CROSSMAP="${CROSSMAP:-CrossMap}"
CHAIN="${CHAIN:-${DATA_ROOT}/Tools/liftover/hg19ToHg38.over.chain.gz}"
REF37="${REF37:-}"
REF38="${REF38:-}"

usage() {
  cat <<'USAGE'
Usage: build_r2_only_broad_liftover_input.sh --raw VCF.gz [options]

Build a GRCh38 broad technical input from a GRCh37 source VCF using:
  biallelic variants with INFO/R2 >= threshold, without source INFO/MAF cutoff.

Options:
  --raw VCF.gz       Required source VCF/BCF.
  --chrom CHROM      Contig to process. Default: 22.
  --label LABEL      Output label. Default: chr${chrom}.r2only.r2_${threshold}.
  --threads N        bcftools threads. Default: ${LIFTOVER_THREADS:-12}.
  --r2-min FLOAT     INFO/R2 minimum. Default: 0.8.
  --out-dir DIR      Output VCF directory.
  --bcf-dir DIR      Output BCF directory.
  --report TSV       Output aggregate report TSV.
  --sort-temp-dir DIR
                    Directory for bcftools sort temporary files. Default:
                    OUT_DIR/.sort_tmp on the same /data filesystem.
  --ref37 FASTA      Override GRCh37 FASTA.
  --ref38 FASTA      Override GRCh38 FASTA.
  --chain CHAIN      Override hg19ToHg38 chain.
  --bcftools PATH    Override bcftools.
  --crossmap PATH    Override CrossMap.
  --allow-heavy      Allow heavy preprocessing in this invocation.
  --dry-run          Validate paths and print aggregate plan only.
  -h, --help         Show this help.

Safety:
  - No source INFO/MAF cutoff is applied.
  - No sample/person-level rows are printed.
  - Raw logs go to logs/restricted_raw with mode 700.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw) RAW="$2"; shift 2 ;;
    --chrom) CHROM="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --r2-min) R2_MIN="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --bcf-dir) BCF_DIR="$2"; shift 2 ;;
    --report) REPORT="$2"; shift 2 ;;
    --sort-temp-dir) SORT_TMP_DIR="$2"; shift 2 ;;
    --ref37) REF37="$2"; shift 2 ;;
    --ref38) REF38="$2"; shift 2 ;;
    --chain) CHAIN="$2"; shift 2 ;;
    --bcftools) BCFTOOLS="$2"; shift 2 ;;
    --crossmap) CROSSMAP="$2"; shift 2 ;;
    --allow-heavy) ALLOW_HEAVY="YES"; shift ;;
    --dry-run) DRY_RUN="YES"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$RAW" ]]; then
  echo "missing required --raw" >&2
  usage >&2
  exit 2
fi

if [[ -z "$REF37" ]]; then
  REF37="$(cat "$ROOT/resources/reference/grch37_local_fasta.path")"
fi
if [[ -z "$REF38" ]]; then
  REF38="$(cat "$ROOT/resources/reference/grch38_local_fasta.path")"
fi

chrom_clean="${CHROM#chr}"
r2_label="${R2_MIN//./_}"
if [[ -z "$LABEL" ]]; then
  LABEL="chr${chrom_clean}.r2only.r2_${r2_label}"
fi
if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="$ROOT/work/grch38_liftover/r2_only_broad/${LABEL}"
fi
if [[ -z "$BCF_DIR" ]]; then
  BCF_DIR="$ROOT/work/grch38_liftover/r2_only_broad_bcf/${LABEL}"
fi
if [[ -z "$REPORT" ]]; then
  REPORT="$ROOT/reports/${LABEL}.R2_ONLY_BROAD_LIFTOVER_QC.tsv"
fi
if [[ -z "$SORT_TMP_DIR" ]]; then
  SORT_TMP_DIR="$OUT_DIR/.sort_tmp"
fi

LOG_DIR="$ROOT/logs/restricted_raw"
SAN_LOG_DIR="$ROOT/logs/sanitized"
TMP_PARENT="$OUT_DIR/.tmp"
LOCK="$OUT_DIR/.build.lock"
FINAL_VCF="$OUT_DIR/${LABEL}.grch38.norm.vcf.gz"
FINAL_BCF="$BCF_DIR/${LABEL}.grch38.norm.bcf"

umask 077
mkdir -p "$OUT_DIR" "$BCF_DIR" "$TMP_PARENT" "$SORT_TMP_DIR" "$LOG_DIR" "$SAN_LOG_DIR" "$(dirname "$REPORT")"
chmod 700 "$LOG_DIR" "$SAN_LOG_DIR"

for path in "$RAW" "$RAW.tbi" "$CHAIN" "$REF37" "$REF37.fai" "$REF38" "$REF38.fai" "$BCFTOOLS" "$CROSSMAP"; do
  test -e "$path" || { echo "MISSING_REQUIRED_PATH $path" >&2; exit 3; }
done

if [[ "$DRY_RUN" == "YES" ]]; then
  {
    printf 'metric\tstatus\tvalue\tdetail\n'
    printf 'DRY_RUN\tPASS\tYES\tno heavy command executed\n'
    printf 'SOURCE_POLICY\tPASS\tR2_ONLY_NO_SOURCE_INFO_MAF_CUTOFF\tINFO/R2>=%s;biallelic\n' "$R2_MIN"
    printf 'CHROM\tPASS\t%s\trequested contig\n' "$CHROM"
    printf 'RAW\tPASS\t%s\tpath exists; records not printed\n' "$RAW"
    printf 'CHAIN\tPASS\t%s\tGRCh37_to_GRCh38\n' "$CHAIN"
    printf 'REF37\tPASS\t%s\tpath exists\n' "$REF37"
    printf 'REF38\tPASS\t%s\tpath exists\n' "$REF38"
    printf 'FINAL_VCF\tPLANNED\t%s\twill be bgzip-indexed VCF\n' "$FINAL_VCF"
    printf 'FINAL_BCF\tPLANNED\t%s\twill be indexed BCF\n' "$FINAL_BCF"
    printf 'REPORT\tPLANNED\t%s\taggregate QC only\n' "$REPORT"
  } > "$REPORT"
  chmod 600 "$REPORT"
  printf 'DRY_RUN_REPORT\t%s\n' "$REPORT"
  exit 0
fi

if [[ "$ALLOW_HEAVY" != "YES" ]]; then
  echo "BLOCKED: set ALLOW_HEAVY=YES or pass --allow-heavy after operator approval for broad liftover preprocessing" >&2
  exit 4
fi

if ! mkdir "$LOCK" 2>/dev/null; then
  echo "BLOCKED: another writer appears active: $LOCK" >&2
  exit 9
fi
cleanup() { rmdir "$LOCK" 2>/dev/null || true; }
trap cleanup EXIT

RUN_ID="$(date +%Y%m%d_%H%M%S)_$$"
TMP="$TMP_PARENT/$RUN_ID"
mkdir -p "$TMP"
RAW_LOG="$LOG_DIR/${LABEL}.r2_only_broad.$RUN_ID.raw.log"
SAN_LOG="$SAN_LOG_DIR/${LABEL}.r2_only_broad.$RUN_ID.sanitized.log"
: > "$RAW_LOG"
chmod 600 "$RAW_LOG"

log() {
  printf '%s\t%s\n' "$(date -Is)" "$*" | tee -a "$SAN_LOG" >&2
  printf '%s\t%s\n' "$(date -Is)" "$*" >> "$RAW_LOG"
}

count_records() {
  "$BCFTOOLS" index -n "$1"
}

count_ref_mismatch() {
  local log_file="$1"
  if [[ ! -s "$log_file" ]]; then
    echo 0
  else
    grep -Eic 'REF_MISMATCH|different from the reference|does not match|Reference allele mismatch' "$log_file" || true
  fi
}

duplicate_keys() {
  local vcf="$1"
  "$BCFTOOLS" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$vcf" \
    | LC_ALL=C sort -S 2G \
    | uniq -d \
    | wc -l
}

write_report() {
  local status="$1"
  local reason="$2"
  local tmp_report="$TMP/report.tsv"
  {
    printf 'metric\tstatus\tvalue\tdetail\n'
    printf 'SOURCE_POLICY\tPASS\tR2_ONLY_NO_SOURCE_INFO_MAF_CUTOFF\tINFO/R2>=%s;biallelic\n' "$R2_MIN"
    printf 'SOURCE_BUILD_VERIFIED\tPASS\tGRCh37\tUSER_CONFIRMED_PROVENANCE\n'
    printf 'LIFTOVER_DIRECTION\tPASS\tGRCh37_TO_GRCh38\thg19ToHg38.over.chain.gz\n'
    printf 'CHROM\tPASS\t%s\trequested contig\n' "$CHROM"
    printf 'SOURCE_INPUT_N\t%s\t%s\tR2-only biallelic source count\n' "${SOURCE_INPUT_STATUS:-NOT_REACHED}" "${SOURCE_INPUT_N:-NA}"
    printf 'SOURCE_NORMALIZED_N\t%s\t%s\t%s\n' "${SOURCE_NORMALIZED_STATUS:-NOT_REACHED}" "${SOURCE_NORMALIZED_N:-NA}" "${SOURCE_NORM_PATH:-NA}"
    printf 'SOURCE_REF_MISMATCH_N\t%s\t%s\trestricted_log=%s.source_norm\n' "${SOURCE_REF_STATUS:-NOT_REACHED}" "${SOURCE_REF_MISMATCH_N:-NA}" "$RAW_LOG"
    printf 'SOURCE_DUPLICATE_KEY_N\t%s\t%s\tsource CHROM_POS_REF_ALT duplicate count\n' "${SOURCE_DUP_STATUS:-NOT_REACHED}" "${SOURCE_DUPLICATE_KEY_N:-NA}"
    printf 'LIFTOVER_MAPPED_N\t%s\t%s\t%s\n' "${LIFTOVER_STATUS:-NOT_REACHED}" "${LIFTOVER_MAPPED_N:-NA}" "${LIFTED_PATH:-NA}"
    printf 'LIFTOVER_REJECTED_N\t%s\t%s\tCrossMap unmapped records if emitted\n' "${LIFTOVER_REJECT_STATUS:-NOT_REACHED}" "${REJECT_N:-NA}"
    printf 'TARGET_NORMALIZED_N\t%s\t%s\t%s\n' "${TARGET_NORMALIZED_STATUS:-NOT_REACHED}" "${TARGET_NORMALIZED_N:-NA}" "${TARGET_NORM_PATH:-NA}"
    printf 'TARGET_REF_MISMATCH_N\t%s\t%s\trestricted_log=%s.target_norm\n' "${TARGET_REF_STATUS:-NOT_REACHED}" "${TARGET_REF_MISMATCH_N:-NA}" "$RAW_LOG"
    printf 'DUPLICATE_CANONICAL_KEY_N\t%s\t%s\tGRCh38 CHROM_POS_REF_ALT duplicate count\n' "${TARGET_DUP_STATUS:-NOT_REACHED}" "${DUPLICATE_CANONICAL_KEY_N:-NA}"
    printf 'FINAL_VCF\t%s\t%s\tindexed VCF\n' "${FINAL_STATUS:-NOT_REACHED}" "$FINAL_VCF"
    printf 'FINAL_BCF\t%s\t%s\tindexed BCF\n' "${FINAL_STATUS:-NOT_REACHED}" "$FINAL_BCF"
    printf 'FINAL_COHORT_MAF_GATE\tBLOCKED_PRODUCTION\tNEEDS_FROZEN_COHORT_AF_MAF\tR2-only broad input must be intersected with frozen cohort AF/MAF before production association\n'
    printf 'PERSON_LEVEL_STDOUT\tPASS\tNONE\tsample columns never printed\n'
    printf 'RUN_STATUS\t%s\t%s\t%s\n' "$status" "$reason" "$SAN_LOG"
  } > "$tmp_report"
  mv -f "$tmp_report" "$REPORT"
  chmod 600 "$REPORT"
}

log "START build R2-only broad liftover input label=$LABEL chrom=$CHROM threads=$THREADS"
log "SOURCE_POLICY=R2_ONLY_NO_SOURCE_INFO_MAF_CUTOFF R2_MIN=$R2_MIN"

log "STEP source broad prefilter: biallelic INFO/R2>=$R2_MIN without INFO/MAF cutoff"
"$BCFTOOLS" view \
  --threads "$THREADS" \
  -r "$CHROM" \
  -m2 -M2 \
  -i "INFO/R2>=$R2_MIN" \
  -Oz \
  -o "$TMP/source.r2_only.vcf.gz.tmp" \
  "$RAW" \
  > "$RAW_LOG.bcftools_view" 2>&1
mv -f "$TMP/source.r2_only.vcf.gz.tmp" "$TMP/source.r2_only.vcf.gz"
"$BCFTOOLS" index --threads "$THREADS" -f "$TMP/source.r2_only.vcf.gz"
SOURCE_INPUT_N="$(count_records "$TMP/source.r2_only.vcf.gz")"
SOURCE_INPUT_STATUS="PASS"
log "SOURCE_INPUT_N=$SOURCE_INPUT_N"

log "STEP source GRCh37 normalization"
"$BCFTOOLS" norm \
  --threads "$THREADS" \
  -f "$REF37" \
  -m -any \
  -c w \
  -Oz \
  -o "$TMP/source.grch37.norm.vcf.gz.tmp" \
  "$TMP/source.r2_only.vcf.gz" \
  > "$RAW_LOG.source_norm" 2>&1
mv -f "$TMP/source.grch37.norm.vcf.gz.tmp" "$TMP/source.grch37.norm.vcf.gz"
"$BCFTOOLS" index --threads "$THREADS" -f "$TMP/source.grch37.norm.vcf.gz"
SOURCE_NORMALIZED_N="$(count_records "$TMP/source.grch37.norm.vcf.gz")"
SOURCE_REF_MISMATCH_N="$(count_ref_mismatch "$RAW_LOG.source_norm")"
SOURCE_DUPLICATE_KEY_N="$(duplicate_keys "$TMP/source.grch37.norm.vcf.gz")"
SOURCE_NORM_PATH="$TMP/source.grch37.norm.vcf.gz"
SOURCE_NORMALIZED_STATUS="PASS"
SOURCE_REF_STATUS="$([[ "$SOURCE_REF_MISMATCH_N" == "0" ]] && echo PASS || echo FAIL)"
SOURCE_DUP_STATUS="$([[ "$SOURCE_DUPLICATE_KEY_N" == "0" ]] && echo PASS || echo FAIL)"
log "SOURCE_NORMALIZED_N=$SOURCE_NORMALIZED_N"
log "SOURCE_REF_MISMATCH_N=$SOURCE_REF_MISMATCH_N"
log "SOURCE_DUPLICATE_KEY_N=$SOURCE_DUPLICATE_KEY_N"
if [[ "$SOURCE_REF_STATUS" != "PASS" || "$SOURCE_DUP_STATUS" != "PASS" ]]; then
  write_report "FAIL" "source_ref_or_duplicate_gate"
  exit 5
fi

log "STEP CrossMap GRCh37/hg19 to GRCh38/hg38"
"$CROSSMAP" vcf \
  --chromid s \
  --ref-consistent \
  "$CHAIN" \
  "$TMP/source.grch37.norm.vcf.gz" \
  "$REF38" \
  "$TMP/lifted.crossmap.vcf" \
  > "$RAW_LOG.crossmap" 2>&1

REJECT_N=0
if [[ -s "$TMP/lifted.crossmap.vcf.unmap" ]]; then
  REJECT_N="$(grep -vc '^#' "$TMP/lifted.crossmap.vcf.unmap" || true)"
elif [[ -s "$TMP/lifted.crossmap.vcf.unmapped" ]]; then
  REJECT_N="$(grep -vc '^#' "$TMP/lifted.crossmap.vcf.unmapped" || true)"
fi
LIFTOVER_REJECT_STATUS="PASS"

log "STEP sort lifted VCF"
"$BCFTOOLS" sort \
  -T "$SORT_TMP_DIR/${LABEL}.sort" \
  -Oz \
  -o "$TMP/lifted.sorted.vcf.gz.tmp" \
  "$TMP/lifted.crossmap.vcf" \
  > "$RAW_LOG.sort" 2>&1
mv -f "$TMP/lifted.sorted.vcf.gz.tmp" "$TMP/lifted.sorted.vcf.gz"
"$BCFTOOLS" index --threads "$THREADS" -f "$TMP/lifted.sorted.vcf.gz"
LIFTOVER_MAPPED_N="$(count_records "$TMP/lifted.sorted.vcf.gz")"
LIFTED_PATH="$TMP/lifted.sorted.vcf.gz"
LIFTOVER_STATUS="PASS"
log "LIFTOVER_MAPPED_N=$LIFTOVER_MAPPED_N"
log "LIFTOVER_REJECTED_N=$REJECT_N"

log "STEP target GRCh38 normalization"
"$BCFTOOLS" norm \
  --threads "$THREADS" \
  -f "$REF38" \
  -m -any \
  -c w \
  -Oz \
  -o "$TMP/${LABEL}.grch38.norm.vcf.gz.tmp" \
  "$TMP/lifted.sorted.vcf.gz" \
  > "$RAW_LOG.target_norm" 2>&1
mv -f "$TMP/${LABEL}.grch38.norm.vcf.gz.tmp" "$TMP/${LABEL}.grch38.norm.vcf.gz"
"$BCFTOOLS" index --threads "$THREADS" -f "$TMP/${LABEL}.grch38.norm.vcf.gz"
TARGET_NORMALIZED_N="$(count_records "$TMP/${LABEL}.grch38.norm.vcf.gz")"
TARGET_REF_MISMATCH_N="$(count_ref_mismatch "$RAW_LOG.target_norm")"
DUPLICATE_CANONICAL_KEY_N="$(duplicate_keys "$TMP/${LABEL}.grch38.norm.vcf.gz")"
TARGET_NORM_PATH="$TMP/${LABEL}.grch38.norm.vcf.gz"
TARGET_NORMALIZED_STATUS="PASS"
TARGET_REF_STATUS="$([[ "$TARGET_REF_MISMATCH_N" == "0" ]] && echo PASS || echo FAIL)"
TARGET_DUP_STATUS="$([[ "$DUPLICATE_CANONICAL_KEY_N" == "0" ]] && echo PASS || echo FAIL)"
log "TARGET_NORMALIZED_N=$TARGET_NORMALIZED_N"
log "TARGET_REF_MISMATCH_N=$TARGET_REF_MISMATCH_N"
log "DUPLICATE_CANONICAL_KEY_N=$DUPLICATE_CANONICAL_KEY_N"
if [[ "$TARGET_REF_STATUS" != "PASS" || "$TARGET_DUP_STATUS" != "PASS" ]]; then
  write_report "FAIL" "target_ref_or_duplicate_gate"
  exit 6
fi

log "STEP BCF conversion"
"$BCFTOOLS" view \
  --threads "$THREADS" \
  -Ob \
  -o "$TMP/${LABEL}.grch38.norm.bcf.tmp" \
  "$TMP/${LABEL}.grch38.norm.vcf.gz" \
  > "$RAW_LOG.bcf_view" 2>&1
mv -f "$TMP/${LABEL}.grch38.norm.bcf.tmp" "$TMP/${LABEL}.grch38.norm.bcf"
"$BCFTOOLS" index --threads "$THREADS" -f "$TMP/${LABEL}.grch38.norm.bcf"

mv -f "$TMP/${LABEL}.grch38.norm.vcf.gz" "$FINAL_VCF"
mv -f "$TMP/${LABEL}.grch38.norm.vcf.gz.csi" "$FINAL_VCF.csi"
mv -f "$TMP/${LABEL}.grch38.norm.bcf" "$FINAL_BCF"
mv -f "$TMP/${LABEL}.grch38.norm.bcf.csi" "$FINAL_BCF.csi"
FINAL_STATUS="PASS"
write_report "PASS" "done"
log "DONE final_vcf=$FINAL_VCF final_bcf=$FINAL_BCF report=$REPORT"
