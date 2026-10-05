#!/usr/bin/env bash
: "${CONTAINER_STATE_DIR:?Set CONTAINER_STATE_DIR}"
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
: "${SAIGE_IMAGE:?Set SAIGE_IMAGE}"
: "${UDOCKER_BIN:?Set UDOCKER_BIN}"
set -uo pipefail

W=${PROJECT_ROOT}/work/run_cond/LDLR2nd_v4
B=${PROJECT_ROOT}/work/run_cond/Barm12
M=${PROJECT_ROOT}/work/run_cond/LDLR19/mini.vcf.gz
R=${PROJECT_ROOT}/work/ref
BCF=bcftools
RAW=${GENOTYPE_DIR}/chr19.vcf.gz
UD=${UDOCKER_BIN}

LEAD1=19:11242307:G:C
LEAD2=19:11343795:A:G
GENES="${GENES:-ENSG00000127616 ENSG00000079805 ENSG00000142453}"

mkdir -p "$W"
cd "$W"
export TMPDIR=${PROJECT_ROOT}/work/tmp UDOCKER_DIR=${CONTAINER_STATE_DIR}
export PROOT_NO_SECCOMP=1

echo "== stage 1: 조건 변이가 mini VCF 안에 있는지 확인 =="
n_mini=$($BCF index -n "$M")
echo "  mini records: $n_mini"
[ "$n_mini" -ge 1 ] || { echo "GATE FAIL: mini empty"; exit 6; }

has1=$($BCF query -r "19:11242307-11242307" -f '%CHROM:%POS:%REF:%ALT\n' "$M" 2>/dev/null | grep -cx "$LEAD1" || true)
has2=$($BCF query -r "19:11343795-11343795" -f '%CHROM:%POS:%REF:%ALT\n' "$M" 2>/dev/null | grep -cx "$LEAD2" || true)
echo "  lead1 in mini: $has1 | lead2 in mini: $has2"
[ "$has1" -ge 1 ] || { echo "GATE FAIL: lead1 absent from mini"; exit 6; }

if [ "$has2" -ge 1 ]; then
  VCF="$M"; IDX="$M.csi"
  echo "  -> lead2 already present. concat/sort 생략 (v3 사망 지점 회피)"
else
  echo "== stage 1b: lead2 추출 후 병합 (메모리 상한 명시) =="
  $BCF view -r "19:11343795-11343795" "$RAW" -Oz -o lead2.vcf.gz --threads 2
  $BCF index -f -c lead2.vcf.gz
  n2=$($BCF index -n lead2.vcf.gz); echo "  lead2 records: $n2"
  [ "$n2" -eq 1 ] || { echo "GATE FAIL: lead2 extract != 1 ($n2)"; exit 6; }
  $BCF concat -a "$M" lead2.vcf.gz -Oz -o mini2.vcf.gz --threads 2
  $BCF sort -m 1G -T "$TMPDIR" mini2.vcf.gz -Oz -o mini2.s.vcf.gz
  $BCF index -f -c mini2.s.vcf.gz
  [ -s mini2.s.vcf.gz.csi ] || { echo "GATE FAIL: index missing"; exit 6; }
  n3=$($BCF index -n mini2.s.vcf.gz); echo "  mini2 records: $n3 (expect $((n_mini+1)))"
  [ "$n3" -eq $((n_mini+1)) ] || { echo "GATE FAIL: mini2 count $n3 != $((n_mini+1))"; exit 6; }
  q=$($BCF query -r "19:11343795-11343795" -f '%CHROM:%POS:%REF:%ALT\n' mini2.s.vcf.gz | grep -cx "$LEAD2" || true)
  [ "$q" -ge 1 ] || { echo "GATE FAIL: lead2 not queryable in mini2"; exit 6; }
  VCF="$W/mini2.s.vcf.gz"; IDX="$W/mini2.s.vcf.gz.csi"
fi
echo "  검정 입력 VCF = $VCF"

echo "== stage 2: 유전자별 합동 조건부 =="
fails=0
for g in $GENES; do
  d="$B/$g"
  [ -s "$d/group.txt" ] || { echo "GATE FAIL: $g group.txt missing"; fails=$((fails+1)); continue; }
  lbl=$(awk '$2=="anno"{print $3; exit}' "$d/group.txt")
  [ "$lbl" = "all" ] || { echo "GATE FAIL: $g anno label='$lbl' != all"; fails=$((fails+1)); continue; }

  $UD run --volume=/data:/data ${SAIGE_IMAGE} \
    step2_SPAtests.R \
      --vcfFile="$VCF" --vcfFileIndex="$IDX" --vcfField=DS \
      --chrom=19 --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 \
      --GMMATmodelFile="$R/saige_step1_v4/tchl_v4.rda" \
      --varianceRatioFile="$R/saige_step1_v4/tchl_v4.varianceRatio.txt" \
      --groupFile="$d/group.txt" --annotation_in_groupTest=all \
      --maxMAF_in_groupTest=0.01 --LOCO=FALSE \
      --condition="$LEAD1,$LEAD2" \
      --SAIGEOutputFile="$W/$g.cond12" > "$W/$g.cond12.log" 2>&1

  if [ ! -s "$W/$g.cond12" ]; then
    echo "FAIL $g: output missing/empty"; fails=$((fails+1)); continue
  fi
  rows=$(( $(wc -l < "$W/$g.cond12") - 1 ))
  hdr=$(head -1 "$W/$g.cond12")
  case "$hdr" in *Pvalue_cond*) ;; *) echo "FAIL $g: no Pvalue_cond column"; fails=$((fails+1)); continue;; esac
  [ "$rows" -ge 1 ] || { echo "FAIL $g: 0 rows"; fails=$((fails+1)); continue; }

  pu=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue") p=i; next}{print $p; exit}' "$d/uncond" 2>/dev/null)
  p1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p; exit}' "$d/cond" 2>/dev/null)
  p2=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Pvalue_cond") p=i; next}{print $p; exit}' "$W/$g.cond12")
  echo "OK $g  uncond=$pu  lead1=$p1  lead1+2=$p2  rows=$rows"
done

echo "fails: $fails"
[ "$fails" -eq 0 ] || { echo "LDLR2ND_V4_GATE_FAIL"; exit 7; }
echo "LDLR2ND_V4_DONE"
