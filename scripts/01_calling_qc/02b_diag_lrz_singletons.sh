#!/usr/bin/env bash
# =====================================================================
# 02b_diag_lrz_singletons.sh  –  READ-ONLY follow-up on LRZ: where did the singletons go, can we get them back?
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · after 02_diag_lrz_jointvcf.sh (24.11.2026):
#   joint.vcf.gz: 318,676 records, only 2 with AC==1; GLnexus config min_allele_copy_number=1 → a later step removed them
# Checks:
#   1  ALL ##*Command header lines of joint.vcf.gz and joint.vep.vcf.gz (02 hid them)
#   2  INFO/AC spectrum (first ALT) of joint.vcf.gz
#   3  GLnexus merge lists: N gVCFs per list, do the gVCFs still exist (sizes, dates)?
#   4  GLnexus.DB / other GLnexus outputs (BCF) still present?
#   5  joint.vep.vcf.gz: record count, AC==1 count (is it a different/unfiltered version?)
#   6  job scripts / shell history mentioning glnexus or bcftools view (exact filter command)
# Usage: bash 02b_diag_lrz_singletons.sh     (≈5 min; overrides: BCF=… H=…)
# =====================================================================
set -u
export LC_ALL=C
BCF=${BCF:-/dss/dsshome1/06/ru36ros/bcftools-1.19/bcftools}
H=${H:-/dss/dsshome1/06/ru36ros}
VCF=$H/joint.vcf.gz; VEP=$H/joint.vep.vcf.gz
PREFIX=${1:-./02b_diag_lrz_$(date +%Y%m%d_%H%M)}
exec > >(tee "$PREFIX.log") 2>&1
sec(){ printf "\n==================== %s ====================\n" "$1"; }
echo "LRZ singleton follow-up | $(hostname) | $(date '+%F %T')"

sec "1 processing commands recorded in headers"
for f in "$VCF" "$VEP"; do
  echo "  --- $(basename "$f")"
  "$BCF" view -h "$f" 2>/dev/null | grep -E '^##[A-Za-z_]*(Command|command|cmdline|CommandLine|VEP=)' | grep -v "$(date '+%a %b')" | cut -c1-400 | sed 's/^/    /'
  "$BCF" view -h "$f" 2>/dev/null | grep -E '^##GLnexusConfig=' | cut -c1-600 | sed 's/^/    /'
done

sec "2 INFO/AC spectrum of joint.vcf.gz (first ALT allele)"
"$BCF" query -f '%INFO/AC\t%INFO/AN\n' "$VCF" | awk -F'\t' '{split($1,a,","); ac=a[1]+0;
   b=(ac<=1?ac:(ac==2?"2":(ac<=5?"3-5":(ac<=10?"6-10":">10")))); c[b]++; if(ac<mn||NR==1)mn=ac}
   END{for(k in c) printf "  AC=%s: %d\n", k, c[k]; printf "  min AC: %d\n", mn}' | sort -V

sec "3 GLnexus merge lists and gVCF availability"
for l in "$H"/mergelist_glnexus*.txt; do
  n=$(grep -c . "$l"); ex=0; miss=0
  while read -r g; do [ -z "$g" ] && continue; if [ -e "$g" ]; then ex=$((ex+1)); else miss=$((miss+1)); fi; done < "$l"
  echo "  $(basename "$l"): $n entries | existing $ex | missing $miss"
  head -2 "$l" | sed 's/^/      e.g. /'
done
first=$(grep -m1 . "$H/mergelist_glnexus_filtered.txt" 2>/dev/null)
if [ -n "$first" ] && [ -e "$first" ]; then
  d=$(dirname "$first"); echo "  gVCF dir: $d | files: $(ls "$d" | grep -c 'g.vcf.gz$') | total size: $(du -sh "$d" 2>/dev/null | cut -f1)"
  ls -la --time-style=long-iso "$first" | sed 's/^/    /'
fi

sec "4 GLnexus database / BCF outputs"
du -sh "$H/GLnexus.DB" 2>/dev/null | sed 's/^/  GLnexus.DB: /'
ls -la --time-style=long-iso "$H"/GLnexus.DB 2>/dev/null | head -8 | sed 's/^/    /'
find "$H" /dss/dssfs03/pn52re/pn52re-dss-0001 -maxdepth 4 \( -name '*.bcf' -o -name '*joint*.vcf*' -o -name '*glnexus*.bcf*' \) 2>/dev/null |
  grep -v '/GLnexus/test/' | head -20 | while read -r f; do echo "  $(ls -la --time-style=long-iso "$f" | awk '{print $5, $6, $7}')  $f"; done

sec "5 joint.vep.vcf.gz – different version?"
if [ -e "$VEP" ]; then
  echo "  records: $("$BCF" view -H "$VEP" | wc -l) | AC==1 records: $("$BCF" view -H -c 1 -C 1 "$VEP" 2>/dev/null | wc -l)"
fi

sec "6 job scripts / history mentioning glnexus or bcftools view"
for h in "$H/.bash_history"; do [ -r "$h" ] && grep -nE 'glnexus|bcftools (view|filter)|joint\.(bcf|vcf)' "$h" | tail -25 | cut -c1-300 | sed 's/^/  /'; done
find "$H" -maxdepth 3 -type f \( -name '*.sh' -o -name '*.slurm' -o -name '*.sbatch' -o -name '*.cmd' \) -newermt 2025-01-01 2>/dev/null |
  xargs -r grep -lE 'glnexus|joint\.vcf' 2>/dev/null | head -10 | while read -r f; do
    echo "  --- $f"; grep -nE 'glnexus|bcftools|joint' "$f" | head -12 | cut -c1-300 | sed 's/^/      /'; done
echo; echo "Log: $PREFIX.log"
