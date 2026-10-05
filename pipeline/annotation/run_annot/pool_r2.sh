#!/usr/bin/env bash
set -u
run1(){ N=$1; bash -c "ulimit -v 4000000; bash ext.sh $N" > ext$N.log 2>&1; }
export -f run1
printf '%s\n' 2 5 6 | xargs -P 3 -I{} bash -c 'run1 {}'
echo POOL_R2_DONE
