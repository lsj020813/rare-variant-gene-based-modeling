#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
export CODE_ROOT
set -uo pipefail
: "${PRS_CS_ROOT:?Set PRS_CS_ROOT to the official PRS-CS checkout}"
[ -f "${PRS_CS_ROOT}/PRScs.py" ] || { printf "%s\n" "PRS-CS checkout missing PRScs.py: ${PRS_CS_ROOT}" >&2; exit 2; }
N=$1; T=$2; LD=$3; S=$4; NG=$5; W=${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p $W/prscs_$T; O=$W/prscs_$T/$T
[ -s $O.chr$N.done ] && exit 0
$PY "${CODE_ROOT}/baselines/prs_cs/mk_bim.py" $N > /dev/null
$PY "${PRS_CS_ROOT}/PRScs.py" --ref_dir=$LD --bim_prefix=$W/geno_hm3/chr$N.rs --sst_file=$S --n_gwas=$NG --out_dir=$O --chrom=$N --seed=20260926 > $W/prscs_$T/chr$N.log 2>&1
rc=$?; f=$(ls ${O}_pst_eff_a1_b0.5_phiauto_chr$N.txt 2>/dev/null)
[ $rc -eq 0 ] && [ -s "$f" ] && echo "ok $(date -Is) $(wc -l < $f)" > $O.chr$N.done
echo "prscs $T chr$N rc=$rc"
