#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
BCF=bcftools
S=${PROJECT_ROOT}/work/ref/orig_index/chr22.vcf.gz
echo "source $(ls -lh $S | awk '{print $5}')"
echo '--- MAF bands (chr22 source)'
$BCF query -f '%INFO/MAF\n' $S | awk '$1!="."{n++; if($1>=0.001&&$1<=0.01)a++; else if($1>0.01&&$1<=0.05)b++; else if($1>0.05)c++} END{printf "total %d | 0.1-1%%: %d | 1-5%%: %d | >5%%: %d\n", n,a,b,c}'
echo '--- DS density 1-5% (every 40th variant)'
$BCF view -i 'INFO/MAF>0.01 && INFO/MAF<=0.05' $S -Ou | $BCF query -f '[%DS\t]\n' | awk 'NR%40==0{for(i=1;i<=NF;i++){t++; if($i>0.05)nz++}} END{printf "cells %d nz %d = %.2f%%\n", t, nz, 100*nz/t}'
echo '--- DS density 0.1-1% band (every 40th)'
$BCF query -f '[%DS\t]\n' ${PROJECT_ROOT}/work/ref/band_vcf/chr22.band.vcf.gz | awk 'NR%40==0{for(i=1;i<=NF;i++){t++; if($i>0.05)nz++}} END{printf "cells %d nz %d = %.2f%%\n", t, nz, 100*nz/t}'
echo BAND15_MEASURE_DONE
