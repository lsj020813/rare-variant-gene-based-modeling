: "${CONDA_INIT_SCRIPT:?Set CONDA_INIT_SCRIPT}"
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR}"
source ${CONDA_INIT_SCRIPT}; conda activate "${CONDA_ENV:?Set CONDA_ENV}"
R=${GENOTYPE_DIR}/chr22.vcf.gz
ls $R.tbi $R.csi 2>&1 | head -2
echo "TYPED rows in first 3000: $(zcat $R | grep -v '^#' | head -3000 | grep -c TYPED)"
echo "filter TYPED=1: $(zcat $R | head -3400 | bcftools view -i 'TYPED=1' -G - 2>&1 | grep -vc '^#')"
echo "filter INFO/TYPED: $(zcat $R | head -3400 | bcftools view -i 'INFO/TYPED' -G - 2>&1 | grep -vc '^#')"
echo "filter TYPED (bare): $(zcat $R | head -3400 | bcftools view -i 'TYPED' -G - 2>&1 | grep -vc '^#')"
zcat $R | grep -v '^#' | grep -m1 TYPED | cut -f1-8 | cut -c1-160
echo DIAG_DONE
