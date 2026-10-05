#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT to the working data root}"
[[ "$PROJECT_ROOT" =~ [^[:space:]] ]] || { printf "%s\n" "PROJECT_ROOT must be nonblank" >&2; exit 2; }
set -u
D=${PROJECT_ROOT}/work/prs/ref; cd $D
for P in "7ek4lwwf2b7f749/ldblk_1kg_eas" "mt6var0z96vb6fv/ldblk_1kg_eur"; do
  n=$(basename $P); [ -s $n.done ] && { echo "have $n"; continue; }
  curl -sS -L --retry 5 --retry-delay 20 -C - -o $n.tar.gz.part "https://www.dropbox.com/s/$P.tar.gz?dl=1" && mv $n.tar.gz.part $n.tar.gz || { echo "FAIL $n"; continue; }
  tar -tzf $n.tar.gz >/dev/null 2>&1 || { echo "CORRUPT $n"; continue; }
  tar -xzf $n.tar.gz && echo "ok $(date -Is) $(du -sh $n | cut -f1)" > $n.done && echo "DONE $n"
done
echo "ALL $(date -Is)"
