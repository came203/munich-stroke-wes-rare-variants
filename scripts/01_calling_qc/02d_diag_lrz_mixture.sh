#!/usr/bin/env bash
# =====================================================================
# 02d_diag_lrz_mixture.sh  –  READ-ONLY on LRZ: sample mixture / contamination / duplicate pairs?
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · after 02c (24.11.2026):
#   per-sample gVCFs contain no private calls with GQ>=30; the few private calls have GQ<=22 (dropped by GLnexus
#   min_AQ1=35); chr22: 0 singletons, 22% doubletons → pattern fits mixed / contaminated samples or duplicate pairs
# Checks:
#   1  header of one gVCF (DeepVariant version / command line)
#   2  per sample (N gVCFs): het VAF distribution and GQ of PASS het calls (clean diploid: VAF peak ~0.5;
#      mixture/contamination: excess at 0.1–0.35), DP
#   3  joint.vcf.gz, all autosomes: records with genotype-based AC==2 carried by two different samples →
#      how often does the SAME pair share them? (random: almost never >2; mixtures/duplicates: tens–hundreds)
#      + per-sample count of AC==2 carriers
# Usage: bash 02d_diag_lrz_mixture.sh        (≈10–15 min; N=5 by default)
# Writes only ./02d_diag_lrz_<date>.log and a temp dir (removed on exit).
# =====================================================================
set -u
export LC_ALL=C
BCF=${BCF:-/dss/dsshome1/06/ru36ros/bcftools-1.19/bcftools}
H=${H:-/dss/dsshome1/06/ru36ros}
VCF=$H/joint.vcf.gz
LIST=$H/mergelist_glnexus_filtered.txt
N=${N:-5}
PREFIX=${1:-./02d_diag_lrz_$(date +%Y%m%d_%H%M)}
exec > >(tee "$PREFIX.log") 2>&1
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
sec(){ printf "\n==================== %s ====================\n" "$1"; }
echo "mixture check | $(hostname) | $(date '+%F %T') | N=$N"

sec "1 gVCF header (first sample)"
g1=$(head -1 "$LIST")
grep -m 40 '^##' "$g1" | grep -iE 'deepvariant|command|source|reference|version' | cut -c1-300 | sed 's/^/  /'

sec "2 per-sample het VAF / GQ / DP (PASS 0/1 calls, all chromosomes)"
printf "  %-10s %7s %7s %8s %8s %8s %8s %9s %9s %7s\n" sample n_het "VAF<.1" ".1-.2" ".2-.35" ".35-.65" ">.65" "GQ<25(%)" "medGQ" medDP
head -n "$N" "$LIST" | while read -r g; do
  awk -F'\t' '/^#CHROM/{s=$10} /^#/{next}
    $7=="PASS" && $5!="<*>" {
      n=split($9,f,":"); split($10,v,":"); for(i=1;i<=n;i++) x[f[i]]=v[i]
      gt=x["GT"]; gsub(/[|]/,"/",gt); if(gt!="0/1") next
      split(x["VAF"],vf,","); vaf=vf[1]+0; t++
      if(vaf<.1)a++; else if(vaf<.2)b++; else if(vaf<.35)c++; else if(vaf<=.65)d++; else e++
      if(x["GQ"]+0<25) lq++; gq[t]=x["GQ"]+0; dp[t]=x["DP"]+0 }
    END{ asort(gq); asort(dp)
         printf "  %-10s %7d %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%% %8.1f%% %9d %7d\n", s, t, 100*a/t, 100*b/t, 100*c/t, 100*d/t, 100*e/t,
                100*lq/t, gq[int(t/2)+1], dp[int(t/2)+1] }' "$g"
done
echo "  (gawk needed for asort; if medians show 0, awk is not gawk – percentages are still valid)"

sec "3 sharing of AC==2 variants between sample pairs (joint.vcf.gz, all records)"
"$BCF" query -l "$VCF" > "$TMP/samples.txt"
"$BCF" query -f '[%GT\t]\n' "$VCF" | awk -F'\t' '
  { ac=0; k=0
    for(i=1;i<NF;i++){ g=$i; if(g=="0/0"||g=="0|0"||g=="./."||g==".|."||g=="./0"||g=="0/.") continue
      gsub(/[|]/,"/",g); n=split(g,a,"/"); c=0; for(j=1;j<=n;j++) if(a[j]!="." && a[j]!="0") c++
      if(c>0){ ac+=c; k++; id[k]=i; if(ac>2) break } }
    if(ac==2 && k==2){ pr[id[1]","id[2]]++; per[id[1]]++; per[id[2]]++; np++ } }
  END{ for(p in pr){ d[pr[p]]++; print "PAIR\t"pr[p]"\t"p }
       for(c in d) print "DIST\t"c"\t"d[c]
       for(s in per) print "PER\t"s"\t"per[s]
       print "NP\t"np }' > "$TMP/pairs.txt"
np=$(awk -F'\t' '$1=="NP"{print $2}' "$TMP/pairs.txt")
echo "  AC==2 records carried by two different samples (one het each): $np"
echo "  how many distinct pairs share k such variants (k: n_pairs):"
awk -F'\t' '$1=="DIST"{print $2"\t"$3}' "$TMP/pairs.txt" | sort -k1,1n | awk -F'\t' '{printf "    k=%s: %s\n",$1,$2}' | head -12
echo "    ... (k>12: $(awk -F'\t' '$1=="DIST" && $2>12{s+=$3} END{print s+0}' "$TMP/pairs.txt") pairs)"
echo "  top 15 pairs (shared AC==2 variants, sample A, sample B):"
awk -F'\t' '$1=="PAIR"{print $2"\t"$3}' "$TMP/pairs.txt" | sort -k1,1nr | head -15 |
  while IFS=$'\t' read -r k p; do a=${p%,*}; b=${p#*,}; echo "    $k  $(sed -n "${a}p" "$TMP/samples.txt")  $(sed -n "${b}p" "$TMP/samples.txt")"; done
echo "  samples whose AC==2 variants are mostly shared with ONE partner (share of their AC==2 variants with top partner):"
awk -F'\t' '$1=="PER"{per[$2]=$3} $1=="PAIR"{split($3,s,","); if($2>top[s[1]]){top[s[1]]=$2}; if($2>top[s[2]]){top[s[2]]=$2}}
  END{ for(x in per){ f=top[x]/per[x]; b=(f>=.5?">=50%":(f>=.2?"20-50%":"<20%")); c[b]++ } for(b in c) printf "    %s: %d samples\n", b, c[b] }' "$TMP/pairs.txt"
echo "  per-sample number of AC==2 variants: $(awk -F'\t' '$1=="PER"{print $3}' "$TMP/pairs.txt" | sort -n | awk '{v[NR]=$1} END{printf "min %d | median %d | p90 %d | max %d (n samples %d)", v[1], v[int(NR/2)+1], v[int(NR*.9)], v[NR], NR}')"
echo
echo "Interpretation: clean cohort → VAF peak 0.35–0.65 (>85%), pairs rarely share >2 AC==2 variants."
echo "Mixed/contaminated samples → excess VAF 0.1–0.35 and recurrent pairs sharing many AC==2 variants."
echo; echo "Log: $PREFIX.log"
