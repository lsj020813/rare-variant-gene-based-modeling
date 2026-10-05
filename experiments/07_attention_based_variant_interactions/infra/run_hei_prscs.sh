#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export CODE_ROOT
cd ${PROJECT_ROOT}/work/prs; PY="${PYTHON_BIN:-python3}"
[ -s out/hei/bbj_hei_prscs.txt ] || { $PY "${CODE_ROOT}/baselines/prs_cs/prep_bbj_hei.py" > logs/prep_bbj_hei.log 2>&1 || exit 1; }
NG=$($PY -c "import json;print(int(round(json.load(open('out/hei/prscs_prep_summary.json'))['N_est_median'])))")
echo "N_gwas=$NG" >> logs/prep_bbj_hei.log
seq 1 22 | xargs -P 6 -I{} bash "${CODE_ROOT}/baselines/prs_cs/prscs_chr.sh" {} hei ref/ldblk_1kg_eas out/hei/bbj_hei_prscs.txt $NG
n=$(ls prscs_hei/*.done 2>/dev/null | wc -l); echo "prscs hei done $n $(date -Is)"; [ $n -eq 22 ] && echo ok > out/hei/prscs.done
