#!/usr/bin/env bash
: "${CONDA_PREFIX:?Set CONDA_PREFIX}"
: "${GENOTYPE_DIR:?Set GENOTYPE_DIR}"
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set +e
export PATH=${CONDA_PREFIX}/bin:$PATH
export TMPDIR=${PROJECT_ROOT}/work/tmp TMP=${PROJECT_ROOT}/work/tmp TEMP=${PROJECT_ROOT}/work/tmp
ORIG=${GENOTYPE_DIR}/chr22.vcf.gz
W=${PROJECT_ROOT}/work/run_ourfm
echo "start $(date)"
echo "== VCF ID shape (region sample) =="; bcftools query -f '%ID\n' -r 22:20000000-20200000 $ORIG | awk '{a=($0 ~ /^rs/)?"rs":(($0 ~ /^[0-9]+:[0-9]+:[A-Z]+:[A-Z]+$/)?"cpra":(($0==".")?"dot":"other")); c[a]++} END{for(k in c) print k, c[k]}'
echo "== SAIGE MarkerID shape =="; awk 'NR>1{a=($3 ~ /^rs/)?"rs":(($3 ~ /^[0-9]+:[0-9]+:[A-Z]+:[A-Z]+$/)?"cpra":"other"); c[a]++} END{for(k in c) print k, c[k]}' ${PROJECT_ROOT}/work/ref/gwas05/tchl.chr22.txt
echo "== MarkerID == CHR:POS:A1:A2 ? =="; awk 'NR>1{k=$1":"$2":"$4":"$5; if(k==$3)m++; else n++} END{print "match",m,"mismatch",n}' ${PROJECT_ROOT}/work/ref/gwas05/tchl.chr22.txt
echo "== R2 dist chr22 3Mb sample =="; S=$(date +%s); bcftools query -f '%INFO/R2\t%INFO/MAF\n' -r 22:20000000-23000000 $ORIG | awk '{n++; if($1>=0.7)r++; if($1>=0.7 && $2>=0.05)c++; if($1>=0.7 && $2>=0.01 && $2<0.05)l++; if($1>=0.7 && $2>=0.001 && $2<0.01)b++} END{print "n_3Mb",n,"r2ge0.7",r,"common",c,"low",l,"band",b}'; echo "query_secs $(( $(date +%s)-S ))"
echo "== rda sampleID =="; Rscript -e '.required_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) stop(name, " must be set", call. = FALSE)
  value
}
for(t in c("tchl","htn","dm","lip")){e=new.env(); load(sprintf(paste0(.required_env("PROJECT_ROOT"), "/work/ref/saige_step1_v4/%s_v4.rda"),t), envir=e); o=get(ls(e)[1],envir=e); cat(t, ls(e)[1], class(o), length(o$sampleID), sum(duplicated(o$sampleID)), o$traitType, "\n")}' 2>&1 | tail -5
echo "== VCF samples vs rda overlap =="; bcftools query -l $ORIG > $W/tmp/_vcf_samples.txt; Rscript -e '.required_env <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) stop(name, " must be set", call. = FALSE)
  value
}
v=readLines(paste0(.required_env("PROJECT_ROOT"), "/work/run_ourfm/tmp/_vcf_samples.txt")); for(t in c("tchl","htn","dm","lip")){e=new.env(); load(sprintf(paste0(.required_env("PROJECT_ROOT"), "/work/ref/saige_step1_v4/%s_v4.rda"),t), envir=e); o=get(ls(e)[1],envir=e); cat(t, "in_vcf", sum(o$sampleID %in% v), "of", length(o$sampleID), "\n")}' 2>&1 | tail -4
rm -f $W/tmp/_vcf_samples.txt
echo "== plink2 dosage r-unphased help =="; plink2 --help r-unphased 2>&1 | head -40
echo "== susieR funcs =="; Rscript -e '.libPaths("~/R/lib"); library(susieR); print(args(susie_rss))' 2>&1 | tail -8
echo "PROBE_DONE $(date)"
