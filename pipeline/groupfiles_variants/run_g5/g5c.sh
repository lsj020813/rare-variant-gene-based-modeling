#!/usr/bin/env bash
: "${PROJECT_ROOT:?Set PROJECT_ROOT}"
set -uo pipefail
O=${PROJECT_ROOT}/work/ref/deductive
L=${PROJECT_ROOT}/work/ref/loftee
V=${PROJECT_ROOT}/work/ref/vep_env
S=$V/share/ensembl-vep-116.1-0
B=https://personal.broadinstitute.org/konradk/loftee_data/GRCh38
mkdir -p "$L"
exec 9>"$O/.g5c.lock"; flock -n 9 || { echo "another stage C"; exit 0; }

get(){
  local U=$1 F=$2
  [ -e "$F.done" ] && { echo "[$(basename $F)] already"; return 0; }
  local WANT=""
  for t in 1 2 3 4 5; do
    WANT=$(curl -sIL -m 60 "$U" | grep -i '^content-length' | tail -1 | tr -d '\r' | awk '{print $2}')
    case "$WANT" in ''|*[!0-9]*) ;; *) break;; esac
    WANT=$(curl -sL -m 60 -r 0-0 -D - -o /dev/null "$U" | grep -i '^content-range' | tr -d '\r' | sed 's@.*/@@')
    case "$WANT" in ''|*[!0-9]*) ;; *) break;; esac
    sleep 5
  done
  case "$WANT" in ''|*[!0-9]*) echo "FAIL size unknown after 5 probes $U"; return 1;; esac
  echo "[$(basename $F)] expect $WANT bytes"
  for i in 1 2 3 4 5; do
    curl -sSL --retry 8 --retry-delay 10 -m 86400 -o "$F.part" "$U" || true
    local S2=$(stat -c %s "$F.part" 2>/dev/null || echo 0)
    if [ "$S2" -eq "$WANT" ]; then mv "$F.part" "$F"; touch "$F.done"
      echo "[$(basename $F)] OK $S2"; return 0; fi
    echo "[$(basename $F)] try$i $S2 != $WANT — restart"; rm -f "$F.part"
  done
  echo "[$(basename $F)] GIVEUP"; return 1
}
get $B/loftee.sql.gz                                   $L/loftee.sql.gz               || exit 4
get $B/human_ancestor.fa.gz                            $L/human_ancestor.fa.gz        || exit 4
get $B/human_ancestor.fa.gz.fai                        $L/human_ancestor.fa.gz.fai    || exit 4
get $B/human_ancestor.fa.gz.gzi                        $L/human_ancestor.fa.gz.gzi    || exit 4
get $B/gerp_conservation_scores.homo_sapiens.GRCh38.bw $L/gerp.bw                     || exit 4
[ -s "$L/loftee.sql" ] || gunzip -kf "$L/loftee.sql.gz" 2>/dev/null || true
echo "PROCUREMENT_OK $(du -sh $L | cut -f1)"

export PERL5LIB="$S/modules:$S:$V/lib/perl5/site_perl:${PERL5LIB:-}"
export PATH="$V/bin:$PATH"
GTF=$O/gencode.sorted.gtf.gz
FA=${PROJECT_ROOT}/work/tmp/ref/GRCh38_no_alt.fa
$V/bin/perl $V/bin/vep --input_file "$O/vep_input.vcf.gz" --output_file "$O/vep_lof.tsv" \
  --gtf "$GTF" --fasta "$FA" --format vcf --tab --no_stats --force_overwrite \
  --fields "Uploaded_variation,Location,Allele,Gene,Feature,Consequence,BIOTYPE,CANONICAL,LoF,LoF_filter,LoF_flags" \
  --canonical --biotype --fork 6 --buffer_size 5000 \
  --dir_plugins "$S" \
  --plugin LoF,loftee_path:$S,human_ancestor_fa:$L/human_ancestor.fa.gz,conservation_file:$L/loftee.sql,gerp_bigwig:$L/gerp.bw \
  2>&1 | grep -vE "Ignoring '(start_codon|UTR|Selenocysteine)'" | tail -15

NL=$(grep -vc '^#' "$O/vep_lof.tsv" 2>/dev/null || echo 0)
[ "$NL" -ge 1000 ] || { echo "GATE FAIL vep_lof only $NL rows"; exit 5; }
echo "LOFTEE annotated rows: $NL"

awk -F'\t' 'BEGIN{OFS="\t"} /^#/{next}
  $7=="protein_coding" && $8=="YES" &&
  $6 ~ /stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant|start_lost|transcript_ablation/ {
    print $1, $2, $4, $6, $9, $10 }' "$O/vep_lof.tsv" | sort -u > "$O/dz_all.tsv"
awk -F'\t' '$5=="HC"' "$O/dz_all.tsv" | sort -u > "$O/dz_hc.tsv"

echo "== candidates (consequence only): $(wc -l < $O/dz_all.tsv)"
echo "== HC (final deductive zone):     $(wc -l < $O/dz_hc.tsv)"
echo "== LoF verdict breakdown =="
awk -F'\t' '{v=($5==""?"(none)":$5); c[v]++} END{for(k in c) printf "   %-8s %s\n", k, c[k]}' "$O/dz_all.tsv" | sort -k2 -rn
echo "== LC filter reasons =="
awk -F'\t' '$5=="LC"{n=split($6,a,","); for(i=1;i<=n;i++) c[a[i]]++} END{for(k in c) printf "   %-18s %s\n", k, c[k]}' "$O/dz_all.tsv" | sort -k2 -rn | head -8
echo "== genes =="
awk -F'\t' '{g[$3]++} END{n1=0;n2=0; for(k in g){n1++; if(g[k]>=2)n2++}
  printf "   HC genes >=1: %d\n   HC genes >=2: %d\n", n1, n2}' "$O/dz_hc.tsv"
echo STAGE_C_COMPLETE
