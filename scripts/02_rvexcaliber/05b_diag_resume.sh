#!/usr/bin/env bash
# =====================================================================
# 05b_diag_resume.sh  –  READ-ONLY (VM, seconds): basis for resuming RV-EXCALIBER v2 after ANNOVAR
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · after 05 missense run 24.11.2026:
#   internal variant list empty (93 B), all gnomAD variant lists 0 B → suspected: missense_only.R filters
#   ExonicFunc == "nonsynonymous SNV" (space) while RV-EXCALIBER tables use "nonsynonymous_SNV"
# Checks:
#   1  Func / ExonicFunc values and M-CAP availability in the internal pruned annotation and in a gnomAD repo file
#   2  the filter lines of get_Varlist_rvexcaliber_missense_only.R that compare ExonicFunc / M-CAP
#   3  step-1A shell script from the get_Varlist call to the end (what happens after filtering) + function layout
#   4  sub_generate_RVBurdenMatrix_internal script (burden-matrix generation) and step-2 get_Varlist block
# Usage: bash 05b_diag_resume.sh   → ./05b_diag_<date>.log
# =====================================================================
set -u
export LC_ALL=C
W=${W:-$HOME/rvx_v2_missense}
DS=my_samples_plink_hg19_v2u
RVX=/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main
SH=$RVX/scripts/shell; R=$RVX/scripts/Rscripts
PREFIX=${1:-./05b_diag_$(date +%Y%m%d_%H%M)}
exec > >(tee "$PREFIX.log") 2>&1
sec(){ printf "\n==================== %s ====================\n" "$1"; }
tabcols(){ # file, column names...: value counts of named columns (gz or plain)
  local f=$1; shift
  { case "$f" in *.gz) zcat "$f";; *) cat "$f";; esac; } | awk -F'\t' -v cols="$*" '
    NR==1{n=split(cols,c," "); for(i=1;i<=NF;i++) for(j=1;j<=n;j++) if($i==c[j]) k[j]=i; next}
    {for(j=1;j<=n;j++) if(k[j]) t[c[j]"\t"$(k[j])]++; if(k[n]) {if($(k[n])==".") mm++; else mv++}}
    END{for(x in t) if(t[x]>=20 || x ~ /synonymous|frameshift|stop/) print "    "x"\t"t[x]; printf "    (last column: %d non-missing, %d \".\")\n", mv, mm}' | sort
}

sec "1a internal pruned annotation: Func / ExonicFunc / M-CAP"
PA=$W/${DS}_pruned_annotation_for_R_input_hcc.txt.gz
echo "  header: $(zcat "$PA" | head -1)"
tabcols "$PA" Func.refGene ExonicFunc.refGene M-CAP_score | grep -v 'M-CAP_score'
sec "1b gnomAD repo pruned annotation (chr22, hcc)"
G=$(ls "$RVX"/gnomAD_211_exomes_hg19/22_gnomAD_pruned_annotation_for_R_input_hcc.txt.gz 2>/dev/null | head -1)
echo "  file: $G"; echo "  header: $(zcat "$G" | head -1)"
tabcols "$G" Func.refGene ExonicFunc.refGene M-CAP_score | grep -v 'M-CAP_score'

sec "2 comparisons in get_Varlist_rvexcaliber_missense_only.R"
grep -nE 'ExonicFunc|Func\.refGene|M-CAP|MCAP|AF_int|gnomAD_check' "$R/get_Varlist_rvexcaliber_missense_only.R" | cut -c1-200 | sed 's/^/  /'

sec "3 step-1A script: functions and everything from the get_Varlist block to the end"
F=$SH/get_RVBurdenMatrix_internal_rvexcaliber.clean.sh
echo "  function definitions / standalone calls:"
grep -nE '^[[:space:]]*function [A-Za-z_]+|^[[:space:]]*[A-Za-z_]+[[:space:]]*$' "$F" | sed 's/^/    /'
s=$(grep -n 'function get_Varlist' "$F" | head -1 | cut -d: -f1)
echo "  --- lines $s..EOF (comments/blank removed):"
awk -v s="$s" 'NR>=s && $0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/ {printf "%5d  %s\n", NR, substr($0,1,230)}' "$F"

sec "4 burden-matrix generation and step-2 gnomAD get_Varlist block"
S=$SH/sub_generate_RVBurdenMatrix_internal_rvexcaliber.clean.sh
echo "  --- $S (non-comment lines):"
awk '$0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/ {printf "%5d  %s\n", NR, substr($0,1,230)}' "$S" | head -80
echo "  --- step 2 lines 380-480:"
awk 'NR>=380 && NR<=480 && $0 !~ /^[[:space:]]*#/ && $0 !~ /^[[:space:]]*$/ {printf "%5d  %s\n", NR, substr($0,1,230)}' "$SH/get_RVBurdenMatrix_gnomAD_rvexcaliber.clean.sh"
echo; echo "Log: $PREFIX.log"
