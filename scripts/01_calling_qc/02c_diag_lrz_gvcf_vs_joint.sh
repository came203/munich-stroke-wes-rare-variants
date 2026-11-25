#!/usr/bin/env bash
# =====================================================================
# 02c_diag_lrz_gvcf_vs_joint.sh  –  READ-ONLY on LRZ: did GLnexus drop the singletons?
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · after 02b (24.11.2026):
#   joint.vcf.gz = plain `bcftools view joint.bcf` (no filter); INFO/AC empty; genotype-based MAC=1 only 60/323,489 (hg19 PLINK);
#   GLnexus 1.4.1 DeepVariantWES (min_AQ1 35, min_AQ2 20, revise_genotypes true); all 1,559 gVCFs present (.g.vcf)
# Checks:
#   1  full text of the GLnexus job script (script_nexus2.sh) and the post-processing script (script_post2.sh)
#   2  genotype-based allele-count spectrum of joint.vcf.gz on chr22 (INFO/AC is empty)
#   3  for N samples: PASS variant calls in the sample's own gVCF (chr22) vs that sample's genotypes in joint.vcf.gz:
#      recovered as carrier / site present but sample not carrier / site absent; GQ of the lost calls;
#      cohort allele count (from joint) of the recovered ones
# Usage: bash 02c_diag_lrz_gvcf_vs_joint.sh        (≈5–10 min; N=3 samples by default: N=5 bash …)
# Writes only ./02c_diag_lrz_<date>.log and a temp dir (removed on exit).
# =====================================================================
set -u
export LC_ALL=C
BCF=${BCF:-/dss/dsshome1/06/ru36ros/bcftools-1.19/bcftools}
H=${H:-/dss/dsshome1/06/ru36ros}
VCF=$H/joint.vcf.gz
LIST=$H/mergelist_glnexus_filtered.txt
CHR=${CHR:-chr22}
N=${N:-3}
PREFIX=${1:-./02c_diag_lrz_$(date +%Y%m%d_%H%M)}
exec > >(tee "$PREFIX.log") 2>&1
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
sec(){ printf "\n==================== %s ====================\n" "$1"; }
echo "gVCF vs joint | $(hostname) | $(date '+%F %T') | region $CHR | N=$N"

sec "1 GLnexus job script and post-processing script"
for s in "$H/script_nexus2.sh" "$H/script_post2.sh"; do
  echo "  --- $s"; grep -vE '^[[:space:]]*$' "$s" 2>/dev/null | cut -c1-250 | sed 's/^/    /'
done

sec "2 genotype-based allele-count spectrum, joint.vcf.gz $CHR"
"$BCF" query -r "$CHR" -f '%ALT[\t%GT]\n' "$VCF" | awk -F'\t' '{ac=0; for(i=2;i<=NF;i++){g=$i; gsub(/[|]/,"/",g); n=split(g,a,"/"); for(j=1;j<=n;j++) if(a[j]!="." && a[j]!="0") ac++}
   b=(ac<=2?ac:(ac<=5?"3-5":(ac<=10?"6-10":">10"))); c[b]++; t++}
   END{for(k in c) printf "  AC=%s: %d (%.1f%%)\n", k, c[k], 100*c[k]/t; printf "  records: %d\n", t}' | sort -V

sec "3 per-sample: own gVCF PASS calls vs joint genotypes ($CHR)"
"$BCF" query -r "$CHR" -f '%POS\t%REF\t%ALT[\t%GT]\n' "$VCF" > "$TMP/joint_chr.txt"
"$BCF" query -l "$VCF" > "$TMP/samples.txt"
# cohort AC per position (sum over records at a position)
awk -F'\t' '{ac=0; for(i=4;i<=NF;i++){g=$i; gsub(/[|]/,"/",g); n=split(g,a,"/"); for(j=1;j<=n;j++) if(a[j]!="." && a[j]!="0") ac++} s[$1]+=ac}
     END{for(p in s) print p"\t"s[p]}' "$TMP/joint_chr.txt" | sort -k1,1 > "$TMP/pos_ac.txt"
printf "  %-12s %8s %10s %14s %12s %14s %s\n" sample gvcf_PASS recovered site_not_carrier site_absent "lost_GQ>=30" "AC_in_joint_of_recovered(1|2|3-5|>5)"
head -n "$N" "$LIST" | while read -r g; do
  sname=$(grep -m1 '^#CHROM' "$g" | cut -f10)
  col=$(awk -v s="$sname" '$1==s{print NR; exit}' "$TMP/samples.txt")
  [ -z "$col" ] && { echo "  $sname: not in joint sample list"; continue; }
  # PASS variant calls with a non-ref GT in the sample's gVCF
  awk -F'\t' -v c="$CHR" '$1==c && $7=="PASS" && $5!="<*>" {
      n=split($9,f,":"); split($10,v,":"); for(i=1;i<=n;i++) x[f[i]]=v[i];
      g=x["GT"]; gsub(/[|]/,"/",g); if(g ~ /[1-9]/) print $2"\t"$4"\t"$5"\t"g"\t"x["GQ"] }' "$g" | sort -k1,1 > "$TMP/g.txt"
  # sample genotype in joint (column 3+col)
  awk -F'\t' -v k=$((3 + col)) '{g=$k; gsub(/[|]/,"/",g); car=(g ~ /[1-9]/)?1:0; if(!(($1) in c) || car) c[$1]=car} END{for(p in c) print p"\t"c[p]}' \
      "$TMP/joint_chr.txt" | sort -k1,1 > "$TMP/j.txt"
  join -t $'\t' -a1 -e NA -o 1.1,1.5,2.2 "$TMP/g.txt" "$TMP/j.txt" > "$TMP/gj.txt"
  join -t $'\t' -a1 -e NA -o 1.1,1.2,1.3,2.2 "$TMP/gj.txt" "$TMP/pos_ac.txt" > "$TMP/gja.txt"
  awk -F'\t' -v s="$sname" '{t++; if($3=="1"){r++; a=$4+0; b=(a<=1?"1":(a==2?"2":(a<=5?"3-5":">5"))); h[b]++}
       else if($3=="0"){nc++; if($2+0>=30) lq++} else {ab++; if($2+0>=30) lq++}}
       END{printf "  %-12s %8d %10d %14d %12d %14d %s|%s|%s|%s\n", s, t, r, nc, ab, lq, h["1"]+0, h["2"]+0, h["3-5"]+0, h[">5"]+0}' "$TMP/gja.txt"
  echo "    examples of lost calls (pos GQ joint_carrier cohortAC):"
  awk -F'\t' '$3!="1"' "$TMP/gja.txt" | sort -t$'\t' -k2,2nr | head -4 | sed 's/^/      /'
done
echo
echo "Interpretation: if most PASS calls with high GQ are 'site_absent' and the recovered ones almost never have"
echo "cohort AC=1, GLnexus dropped private alleles (unifier/genotyper settings) – the gVCFs still contain them."
echo; echo "Log: $PREFIX.log"
