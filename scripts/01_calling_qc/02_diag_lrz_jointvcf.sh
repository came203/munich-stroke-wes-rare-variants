#!/usr/bin/env bash
# =====================================================================
# 02_diag_lrz_jointvcf.sh  –  READ-ONLY diagnostic of the GRCh38 joint VCF (run on LRZ)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · RVGRS_Plan.md §7b (Part A findings 24.11.2026)
#
# Questions:
#   1  Are singletons already missing in joint.vcf.gz (GRCh38)? (raw hg19 PLINK set: 60 MAC=1 of 323,489)
#   2  How many multiallelic records / duplicate positions are there (source of the TTK duplicate,
#      the 1,419 malformed GT records after CrossMap + norm -m -)?
#   3  Which filters/steps are recorded in the header (GLnexus config, bcftools commands)?
#   4  Is there an earlier / unfiltered version (GLnexus BCF, other joint*.vcf.gz) and how many gVCFs?
#   5  Which tools are available on LRZ for a clean rebuild (bcftools plugins incl. liftover, CrossMap, java/picard)?
# Usage (interactive node or login node; bcftools stats takes a few minutes):
#   bash 02_diag_lrz_jointvcf.sh            # override: VCF=/path/joint.vcf.gz BCF=/path/bcftools
# Writes only ./02_diag_lrz_<date>.log and a temp dir (removed on exit).
# =====================================================================
set -u
export LC_ALL=C
BCF=${BCF:-/dss/dsshome1/06/ru36ros/bcftools-1.19/bcftools}
VCF=${VCF:-/dss/dsshome1/06/ru36ros/joint.vcf.gz}
DIR=$(dirname "$VCF")
PREFIX=${1:-./02_diag_lrz_$(date +%Y%m%d_%H%M)}
exec > >(tee "$PREFIX.log") 2>&1
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
sec(){ printf "\n==================== %s ====================\n" "$1"; }
echo "LRZ joint-VCF diagnostic | $(hostname) | $(date '+%F %T')"
echo "bcftools: $BCF ($("$BCF" --version 2>/dev/null | head -1))"
echo "VCF:      $VCF"

sec "1 files"
ls -la "$VCF"* 2>&1
echo "  other joint / BCF files in $DIR (depth 2):"
find "$DIR" -maxdepth 2 \( -iname 'joint*' -o -iname '*.bcf' -o -iname '*glnexus*' \) 2>/dev/null | head -30 | sed 's/^/    /'
echo "  gVCFs under $DIR (depth 4): $(find "$DIR" -maxdepth 4 -name '*.g.vcf.gz' 2>/dev/null | wc -l)"
find "$DIR" -maxdepth 4 -name '*.g.vcf.gz' 2>/dev/null | head -3 | sed 's/^/    e.g. /'

sec "2 header provenance"
"$BCF" view -h "$VCF" | grep -E '^##(source|GLnexus|bcftools|DeepVariant|command|FILTER|reference|fileformat)' | grep -v 'viewCommand' | cut -c1-300 | head -40
echo "  FORMAT fields: $("$BCF" view -h "$VCF" | grep -o '^##FORMAT=<ID=[^,]*' | sed 's/##FORMAT=<ID=//' | paste -sd' ' -)"
echo "  INFO fields:   $("$BCF" view -h "$VCF" | grep -o '^##INFO=<ID=[^,]*' | sed 's/##INFO=<ID=//' | paste -sd' ' -)"
echo "  contig naming: $("$BCF" view -h "$VCF" | grep -m3 '^##contig' | sed 's/,.*//; s/##contig=<ID=//' | paste -sd' ' -)"
echo "  samples: $("$BCF" query -l "$VCF" | wc -l)"

sec "3 bcftools stats (records, multiallelic, singletons, AF spectrum)"
"$BCF" stats "$VCF" > "$TMP/stats.txt" 2> "$TMP/stats.err"
grep '^SN' "$TMP/stats.txt" | cut -f3- | sed 's/^/  /'
echo "  singleton stats (SiS: allele count 1 → SNPs / transitions / transversions / indels):"
grep '^SiS' "$TMP/stats.txt" | cut -f3- | sed 's/^/    /'
echo "  AF spectrum (first 6 bins: AF  n_SNPs  n_indels):"
grep '^AF' "$TMP/stats.txt" | awk -F'\t' '{print "    "$3"  "$4"  "$7}' | head -6
[ -s "$TMP/stats.err" ] && { echo "  stats stderr (first 5 lines):"; head -5 "$TMP/stats.err" | sed 's/^/    /'; }

sec "4 FILTER values, duplicate positions, ALT-count distribution"
"$BCF" query -f '%CHROM\t%POS\t%FILTER\t%ALT\n' "$VCF" > "$TMP/q.txt"
echo "  FILTER: $(cut -f3 "$TMP/q.txt" | sort | uniq -c | sort -k1,1nr | awk '{printf "%s:%s  ",$2,$1}')"
echo "  records sharing CHROM:POS with another record: $(cut -f1,2 "$TMP/q.txt" | uniq -d | wc -l) positions"
echo "  ALT alleles per record: $(awk -F'\t' '{n=split($4,a,","); c[(n>5?">5":n)]++} END{for(k in c) printf "%s:%d  ", k, c[k]}' "$TMP/q.txt")"

sec "5 genotype-based singleton count (records with ALT allele count == 1)"
echo "  records with AC==1 (bcftools view -c1 -C1): $("$BCF" view -H -c 1 -C 1 "$VCF" 2> "$TMP/ac.err" | wc -l)"
[ -s "$TMP/ac.err" ] && head -3 "$TMP/ac.err" | sed 's/^/    stderr: /'

sec "6 tools available for a clean rebuild"
echo "  bcftools plugins (liftover / fill-tags / norm-related):"
"$BCF" plugin -l 2>/dev/null | grep -iE 'liftover|fill-tags|fixref|setGT' | sed 's/^/    /'
echo "  BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS:-<unset>}"
for t in CrossMap CrossMap.py java picard gatk liftOver; do printf "  %-10s %s\n" "$t" "$(command -v $t 2>/dev/null || echo '-')"; done
echo; echo "Log: $PREFIX.log"
