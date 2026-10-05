#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -u; BCF=bcftools; R=${PROJECT_ROOT}/work/ref15; N=$1; O=$R/annot/t4/chr$N.t4.tsv; mkdir -p $R/annot/t4; [ -s $O.done ] && exit 0
{ printf 'key37\tmaf\tr2\tavg_cs\tis_typed\tis_indel\n'; $BCF query -f '%CHROM:%POS:%REF:%ALT\t%MAF\t%R2\t%AVG_CS\t%TYPED\t%REF\t%ALT\n' $R/band_vcf/chr$N.band.vcf.gz | awk -F'\t' 'BEGIN{OFS="\t"}{k=$1; sub(/^chr/,"",k); t=($5=="."||$5=="")?0:1; ind=(length($6)!=length($7))?1:0; print k,$2,$3,$4,t,ind}'; } > $O.tmp
n=$(($(wc -l < $O.tmp)-1)); [ $n -gt 0 ] || { echo "[chr$N] GATE FAIL 0 rows"; exit 1; }; mv $O.tmp $O && echo "ok $n" > $O.done && echo "[chr$N] T4 $n rows"
