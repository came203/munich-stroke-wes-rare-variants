#!/usr/bin/env bash
# =====================================================================
# 04_rebuild_hg19_v2.sh  –  clean hg19 rebuild of the AQ20 joint call for RV-EXCALIBER (VM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · RVGRS_Plan.md §7d/7e (24.11.2026)
#
# Fixes the v1 problems (CrossMap 0.7.2 on multiallelic records → 1,419 malformed GTs, duplicate records such as
# TTK 6:80717709:GTTTT:G, key collisions):
#   1  GRCh38: autosomes, keep FORMAT/GT only, split multiallelic records (norm -m -any), drop ALT=*, drop AC=0
#   2  liftover of each record's REF span with UCSC liftOver (-minMatch=1.0); keep only same-length, '+'-strand
#      mappings to chr1–22 (strand flips / indel-span changes / unmapped are dropped and counted)
#   3  rewrite CHROM/POS, hg19 contig header, sort
#   4  hg19: norm -f hg19.fa --check-ref x (drops REF mismatches) → chr1→1 → ID chr:pos:ref:alt
#      → ALL copies of duplicated variants removed (no arbitrary keep-first) → AC≥1
#   5  PLINK 1.9 bed/bim/fam (same options as v1)
#   6  QC: counts per step, MAC spectrum, duplicate IDs, TTK/ATXN3 spot check, KING relatedness (plink2)
# Input : /mnt/hgfs/D/WES/FOR_REPORT/joint_AQ20.vcf.gz (+ .csi)   – GLnexus 1.4.1 DeepVariantWES, min_AQ1 20
# Output: /mnt/hgfs/D/WES/FOR_REPORT/v2/  joint_AQ20_hg19.norm.vcf.gz(.csi), my_samples_plink_hg19_v2.{bed,bim,fam},
#         04_rebuild.log, 04_rebuild_summary.tsv, liftover_dropped.tsv, dup_ids.txt, *_king.kin0
# Usage : bash 04_rebuild_hg19_v2.sh       (≈30–60 min; temp files in $HOME/tmp_rebuild_v2 – local disk, ~15 GB)
# =====================================================================
set -euo pipefail
export LC_ALL=C
WES=/mnt/hgfs/D/WES
IN=${IN:-$WES/FOR_REPORT/joint_AQ20.vcf.gz}
OUT=${OUT:-$WES/FOR_REPORT/v2}
CHAIN=${CHAIN:-$WES/03_references/hg38ToHg19.over.chain.gz}
HG19=${HG19:-$WES/03_references/hg19.fa}
LIFTOVER=${LIFTOVER:-/mnt/hgfs/D/software/liftOver}
BCF=${BCF:-bcftools}; PLINK=${PLINK:-plink}; PLINK2=${PLINK2:-plink2}
TH=${TH:-8}
TMPD=${TMPD:-$HOME/tmp_rebuild_v2}
DS=my_samples_plink_hg19_v2
mkdir -p "$OUT" "$TMPD"
exec > >(tee "$OUT/04_rebuild.log") 2>&1
SUM=$OUT/04_rebuild_summary.tsv; printf "step\tn\tnote\n" > "$SUM"
step(){ echo; echo "== $(date '+%F %T') $*"; }
rec(){ printf "%s\t%s\t%s\n" "$1" "$2" "${3:-}" >> "$SUM"; printf "   %-34s %12s  %s\n" "$1" "$2" "${3:-}"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
nrec(){ "$BCF" index -n "$1" 2>/dev/null || "$BCF" view -H "$1" | wc -l; }

step "0 checks"
for f in "$IN" "$CHAIN" "$HG19" "$HG19.fai"; do [ -s "$f" ] || die "missing $f"; done
[ -x "$LIFTOVER" ] || die "liftOver not executable: $LIFTOVER"
for t in "$BCF" "$PLINK" "$PLINK2"; do command -v "$t" >/dev/null || die "$t not found"; done
echo "   bcftools: $("$BCF" --version | head -1) | plink: $("$PLINK" --version | head -1) | plink2: $("$PLINK2" --version | head -1)"
[ -e "$IN.csi" ] || [ -e "$IN.tbi" ] || { echo "   indexing $IN"; "$BCF" index --threads "$TH" "$IN"; }
echo "   free space TMPD ($TMPD): $(df -h "$TMPD" | awk 'NR==2{print $4}') | OUT: $(df -h "$OUT" | awk 'NR==2{print $4}')"
ns=$("$BCF" query -l "$IN" | wc -l); [ "$ns" -eq 1538 ] || die "expected 1538 samples, found $ns"
rec input_records "$(nrec "$IN")" "$(basename "$IN"), $ns samples"

step "1 GRCh38: autosomes, GT only, split multiallelic, drop ALT=* and AC=0"
REG=$(seq -s, -f 'chr%g' 1 22)
"$BCF" view --threads "$TH" -r "$REG" -Ou "$IN" \
  | "$BCF" annotate --threads "$TH" -x INFO,FORMAT/RNC,FORMAT/DP,FORMAT/AD,FORMAT/GQ,FORMAT/PL -Ou \
  | "$BCF" norm --threads "$TH" -m -any -Ou \
  | "$BCF" view --threads "$TH" -e 'ALT="*"' -c 1 -Ob -o "$TMPD/s1_b38_split.bcf"
"$BCF" index -f "$TMPD/s1_b38_split.bcf"
N1=$(nrec "$TMPD/s1_b38_split.bcf"); rec b38_autosomal_split_biallelic "$N1"

step "2 liftover of REF spans (UCSC liftOver, minMatch 1.0)"
"$BCF" query -f '%CHROM\t%POS0\t%END\n' "$TMPD/s1_b38_split.bcf" | awk -v OFS='\t' '{print $1,$2,$3,NR,0,"+"}' > "$TMPD/s2_b38.bed"
case "$CHAIN" in *.gz) zcat "$CHAIN" > "$TMPD/chain.txt";; *) cp "$CHAIN" "$TMPD/chain.txt";; esac
"$LIFTOVER" -minMatch=1.0 "$TMPD/s2_b38.bed" "$TMPD/chain.txt" "$TMPD/s2_hg19.bed" "$TMPD/s2_unmapped.bed" 2> "$TMPD/liftover.err" || die "liftOver failed: $(tail -2 "$TMPD/liftover.err")"
awk -F'\t' -v OFS='\t' -v O="$OUT/liftover_dropped.tsv" '
  NR==FNR { len[$4]=$3-$2; next }
  { id=$4; reason="";
    if ($6!="+") reason="strand_flip"; else if ($3-$2!=len[id]) reason="span_changed"; else if ($1 !~ /^chr([1-9]|1[0-9]|2[0-2])$/) reason="non_autosomal_target";
    if (reason=="") print id, $1, $2+1; else { print id, reason > O; n[reason]++ } }
  END { for (r in n) print r, n[r] > "/dev/stderr" }' "$TMPD/s2_b38.bed" "$TMPD/s2_hg19.bed" > "$TMPD/s2_map.tsv" 2> "$TMPD/s2_drop_counts.txt"
NUN=$(grep -vc '^#' "$TMPD/s2_unmapped.bed" || true)
rec liftover_unmapped "$NUN" "liftOver could not map (minMatch 1.0)"
while read -r r n; do rec "liftover_dropped_$r" "$n"; done < "$TMPD/s2_drop_counts.txt"
NMAP=$(wc -l < "$TMPD/s2_map.tsv"); rec liftover_kept "$NMAP"

step "3 rewrite CHROM/POS with hg19 contigs, sort"
"$BCF" view -h "$TMPD/s1_b38_split.bcf" | grep -v '^##contig=' | grep -v '^#CHROM' > "$TMPD/h.txt"
awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {printf "##contig=<ID=%s,length=%s,assembly=hg19>\n", $1, $2}' "$HG19.fai" >> "$TMPD/h.txt"
"$BCF" view -h "$TMPD/s1_b38_split.bcf" | tail -1 >> "$TMPD/h.txt"
mkdir -p "$TMPD/sort"
{ cat "$TMPD/h.txt"
  "$BCF" view -H "$TMPD/s1_b38_split.bcf" | awk -F'\t' -v OFS='\t' -v M="$TMPD/s2_map.tsv" '
    BEGIN { while ((getline l < M) > 0) { split(l, a, "\t"); c[a[1]]=a[2]; p[a[1]]=a[3] } }
    (NR in c) { $1=c[NR]; $2=p[NR]; print }'
} | "$BCF" sort -m 6G -T "$TMPD/sort" -Ob -o "$TMPD/s3_hg19_sorted.bcf"
"$BCF" index -f "$TMPD/s3_hg19_sorted.bcf"
rec hg19_rewritten_sorted "$(nrec "$TMPD/s3_hg19_sorted.bcf")"

step "4 hg19: REF check + left-align, chr1→1, ID, remove all copies of duplicated variants"
seq 1 22 | awk '{print "chr"$1"\t"$1}' > "$TMPD/chrmap.txt"
"$BCF" norm --threads "$TH" -f "$HG19" --check-ref x -Ou "$TMPD/s3_hg19_sorted.bcf" 2> "$TMPD/norm.err" \
  | "$BCF" annotate --rename-chrs "$TMPD/chrmap.txt" -Ou \
  | "$BCF" annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ob -o "$TMPD/s4_hg19_norm.bcf"
"$BCF" index -f "$TMPD/s4_hg19_norm.bcf"
echo "   bcftools norm: $(grep -E 'Lines|REF' "$TMPD/norm.err" | tr '\n' ' ' || true)"
N4=$(nrec "$TMPD/s4_hg19_norm.bcf"); rec hg19_after_ref_check "$N4" "REF mismatches excluded by --check-ref x"
"$BCF" query -f '%ID\n' "$TMPD/s4_hg19_norm.bcf" | sort | uniq -d > "$OUT/dup_ids.txt"
ND=$(wc -l < "$OUT/dup_ids.txt"); rec duplicated_variant_ids "$ND" "all copies removed"
if [ "$ND" -gt 0 ]; then
  "$BCF" view --threads "$TH" -e "ID=@$OUT/dup_ids.txt" -c 1 -Oz -o "$OUT/joint_AQ20_hg19.norm.vcf.gz" "$TMPD/s4_hg19_norm.bcf"
else
  "$BCF" view --threads "$TH" -c 1 -Oz -o "$OUT/joint_AQ20_hg19.norm.vcf.gz" "$TMPD/s4_hg19_norm.bcf"
fi
"$BCF" index -f "$OUT/joint_AQ20_hg19.norm.vcf.gz"
NF=$(nrec "$OUT/joint_AQ20_hg19.norm.vcf.gz"); rec final_hg19_records "$NF"

step "5 PLINK 1.9"
"$PLINK" --vcf "$OUT/joint_AQ20_hg19.norm.vcf.gz" --double-id --keep-allele-order --vcf-half-call missing \
  --allow-extra-chr --make-bed --out "$OUT/$DS" > /dev/null
{ grep -E 'variants|people|genotyping rate' "$OUT/$DS.log" | sed 's/^/   /'; } || true
rec plink_variants "$(wc -l < "$OUT/$DS.bim")"; rec plink_samples "$(wc -l < "$OUT/$DS.fam")"
rec plink_duplicate_ids "$(cut -f2 "$OUT/$DS.bim" | sort | uniq -d | wc -l)" "must be 0"

step "6 QC"
"$PLINK" --bfile "$OUT/$DS" --freq counts --allow-no-sex --out "$TMPD/qc" > /dev/null
awk 'NR>1{m=($5<$6?$5:$6); b=(m==0?"0":(m==1?"1":(m==2?"2":(m<=5?"3-5":(m<=10?"6-10":">10"))))); c[b]++; t++}
     END{for(k in c) printf "   MAC %-5s %8d (%.1f%%)\n", k, c[k], 100*c[k]/t}' "$TMPD/qc.frq.counts" | sort -V
rec plink_MAC1 "$(awk 'NR>1 && ($5==1||$6==1){n++} END{print n+0}' "$TMPD/qc.frq.counts")" "v1 preprocessed: 6 of 115,508"
echo "   spot check TTK 6:80717700-80717720 and ATXN3 14:92537340-92537390 (ID, carriers):"
for r in 6:80717700-80717720 14:92537340-92537390; do
  "$BCF" query -r "$r" -f '%ID[\t%GT]\n' "$OUT/joint_AQ20_hg19.norm.vcf.gz" \
    | awk -F'\t' '{c=0; for(i=2;i<=NF;i++) if($i ~ /1/) c++; printf "     %-60s carriers %d\n", substr($1,1,60), c}' | head -25 || true
done
"$PLINK2" --bfile "$OUT/$DS" --autosome --maf 0.01 --make-king-table --king-table-filter 0.0884 --out "$OUT/${DS}_king" > /dev/null 2>&1 || echo "   KING failed – see $OUT/${DS}_king.log"
if [ -s "$OUT/${DS}_king.kin0" ]; then
  NK=$(($(wc -l < "$OUT/${DS}_king.kin0") - 1)); rec king_pairs_ge_2nd_degree "$NK" "kinship >= 0.0884"
  echo "   pairs (kinship = last column, sorted):"; { awk 'NR>1{print $NF"\t"$0}' "$OUT/${DS}_king.kin0" | sort -k1,1gr | head -25 | cut -f2- | sed 's/^/     /'; } || true
fi
echo; echo "== done $(date '+%F %T'). Summary: $SUM | temp: $TMPD (delete when satisfied)"
