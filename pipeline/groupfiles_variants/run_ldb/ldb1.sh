#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
W=${PROJECT_ROOT}/work; D=$W/run_ldb; mkdir -p $D/out $D/logs; cd $D
BCF=${BCFTOOLS:-bcftools}
P2=${PLINK2:-plink2}
P1=${PROJECT_ROOT}/work/ref/plink19_env/bin/plink
PY=${PYTHON:-python3}
COM=$W/ref/common05/chr19.maf05.vcf.gz; BAND=$W/ref/band_vcf/chr19.band.vcf.gz
TH=8; MEM=20000; SEED=20260914
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
ST=${1:-all}
run(){ [ "$ST" = all ] || [ "$ST" = "$1" ]; }
if run sub; then
  [ -s out/sub10k.txt ] || { $BCF query -l $BAND | $PY -c "import sys,random;random.seed($SEED);ids=[l.strip() for l in sys.stdin];random.shuffle(ids);print('\n'.join(sorted(ids[:10000])))" > out/sub10k.txt; }
  echo "sub n=$(wc -l < out/sub10k.txt)"; touch out/sub.done
fi
if run common; then
  [ -s out/common19.bed ] || $P2 --vcf $COM dosage=DS --keep out/sub10k.keep --double-id --max-alleles 2 --snps-only --maf 0.05 --make-bed --out out/common19 --threads $TH --memory $MEM > logs/common19.plink2.log 2>&1 || { echo "FAIL common rc=$?"; tail -3 logs/common19.plink2.log; exit 2; }
  echo "common variants=$(wc -l < out/common19.bim) samples=$(wc -l < out/common19.fam)"; touch out/common.done
fi
if run blocks; then
  [ -s out/blocks19.blocks.det ] || $P1 --bfile out/common19 --blocks no-pheno-req --blocks-max-kb 500 --out out/blocks19 --threads $TH --memory $MEM > logs/blocks19.plink.log 2>&1 || { echo "FAIL blocks rc=$?"; tail -3 logs/blocks19.plink.log; exit 3; }
  echo "blocks=$(tail -n +2 out/blocks19.blocks.det | wc -l)"; touch out/blocks.done
fi
if run band; then
  [ -s out/band19.pgen ] || $P2 --vcf $BAND dosage=DS --keep out/sub10k.keep --double-id --make-pgen --out out/band19 --threads $TH --memory $MEM > logs/band19.plink2.log 2>&1 || { echo "FAIL band rc=$?"; tail -3 logs/band19.plink2.log; exit 4; }
  echo "band variants=$(tail -n +2 out/band19.pvar | grep -vc '^#') samples=$(tail -n +2 out/band19.psam | wc -l)"; touch out/band.done
fi
if run r2; then
  [ -s out/index5k.txt ] || grep -v '^#' out/band19id.pvar | cut -f3 | $PY -c "import sys,random;random.seed($SEED);v=[l.strip() for l in sys.stdin];random.shuffle(v);print('\n'.join(v[:5000]))" > out/index5k.txt
  [ -s out/r2_19.vcor ] || $P2 --pfile out/band19id --r2-unphased --ld-snp-list out/index5k.txt --ld-window-kb 500 --ld-window-r2 0 --ld-window 999999 --out out/r2_19 --threads $TH --memory $MEM > logs/r2_19.plink2.log 2>&1 || { echo "FAIL r2 rc=$?"; tail -3 logs/r2_19.plink2.log; exit 5; }
  echo "r2 pairs=$(tail -n +2 out/r2_19.vcor | wc -l)"; touch out/r2.done
fi
echo "LDB1_STAGE_${ST}_DONE"
