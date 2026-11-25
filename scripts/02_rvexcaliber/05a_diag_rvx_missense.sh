#!/usr/bin/env bash
# =====================================================================
# 05a_diag_rvx_missense.sh  –  READ-ONLY (VM): what does get_Varlist_rvexcaliber_missense_only.R change,
#                                and how must the RV-EXCALIBER v2 rerun be called?
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · decision Rainer 24.11.2026: primary = missense-only (M-CAP > 0.025)
# Checks:
#   1  diff get_Varlist_rvexcaliber_original.R → get_Varlist_rvexcaliber_missense_only.R (the R scripts actually used
#      live in /home/rainer/WES_analysis/04_software/RV-EXCALIBER-main)
#   2  where get_Varlist is called (step 1A internal, step 2 gnomAD) and with which arguments
#   3  documented arguments of steps 1A / 2 / 3 (header comment blocks), incl. job_interval
#   4  hard-coded paths in shell + R scripts (path_to_rvexcaliber, plink, bedtools, annovar, scripts dir)
#   5  size of the repo data folders (copy vs. symlink for a separate missense-only installation)
#   6  v2 PLINK set: missingness of the 6 related samples (choose which to drop)
# Usage: bash 05a_diag_rvx_missense.sh     → ./05a_diag_<date>.log
# =====================================================================
set -u
export LC_ALL=C
RVXH=${RVXH:-/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main}
RVXM=${RVXM:-/mnt/hgfs/D/WES/04_software/RV-EXCALIBER-main}
V2=${V2:-/mnt/hgfs/D/WES/FOR_REPORT/v2/my_samples_plink_hg19_v2}
PREFIX=${1:-./05a_diag_$(date +%Y%m%d_%H%M)}
exec > >(tee "$PREFIX.log") 2>&1
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
sec(){ printf "\n==================== %s ====================\n" "$1"; }
R=$RVXH/scripts/Rscripts; S=$RVXM/scripts/shell
echo "RV-EXCALIBER missense-only diagnostic | $(date '+%F %T')"

sec "1 diff original -> missense_only (get_Varlist)"
ls -l --time-style=long-iso "$R"/get_Varlist_rvexcaliber*.R 2>&1 | sed 's/^/  /'
diff "$R/get_Varlist_rvexcaliber_original.R" "$R/get_Varlist_rvexcaliber_missense_only.R" | cut -c1-220 | head -150

sec "2 calls of get_Varlist in shell scripts (both copies)"
for d in "$S" "$RVXH/scripts/shell"; do
  echo "  --- $d"; grep -nE 'get_Varlist|Rscript' "$d"/*.sh 2>/dev/null | grep -v ':[0-9]*:[[:space:]]*#' | cut -c1-260 | sed "s#$d/##; s/^/    /"
done
echo "  shell scripts identical between copies: $(diff -rq "$S" "$RVXH/scripts/shell" >/dev/null 2>&1 && echo yes || echo NO)"
diff -rq "$S" "$RVXH/scripts/shell" 2>/dev/null | head -10 | sed 's/^/    /'

sec "3 documented arguments (header comment blocks)"
for f in get_RVBurdenMatrix_internal_rvexcaliber.clean.sh get_RVBurdenMatrix_gnomAD_rvexcaliber.clean.sh get_SummaryAssociations_iCF_gCF_rvexcaliber.clean.sh; do
  echo "  --- $f"
  awk 'NR<=200 && /^#/ && /(^# *[A-Za-z_]+ *$|Description|Notes|Options|[Ee]xample|Usage)/' "$S/$f" | cut -c1-230 | sed 's/^/    /' | head -60
  echo "    positional assignments:"; grep -nE '^[a-zA-Z_]+=\$[0-9]+' "$S/$f" | sed 's/^/      /'
done

sec "4 hard-coded paths"
grep -nHE '^[[:space:]]*(path_to_rvexcaliber|annovar|plink|bedtools|scripts|gnomAD|ranking)[a-zA-Z_]*=' "$S"/*.sh 2>/dev/null | cut -c1-220 | sed "s#$S/##; s/^/  /"
echo "  paths inside R scripts:"
grep -nHE '/home/|/mnt/' "$R"/*.R 2>/dev/null | cut -c1-200 | sed "s#$R/##; s/^/    /" | head -20

sec "5 repo data folder sizes (copy vs symlink)"
du -sh "$RVXH"/* 2>/dev/null | sed 's/^/  /'
echo "  free space in $(dirname "$RVXH"): $(df -h "$(dirname "$RVXH")" | awk 'NR==2{print $4}')"

sec "6 related samples in v2 set: missingness"
if command -v plink >/dev/null && [ -e "$V2.bed" ]; then
  cp "${REL_PAIRS_FILE:?set REL_PAIRS_FILE (related pairs, IID IID per line)}" "$TMP/rel.txt"
  plink --bfile "$V2" --keep "$TMP/rel.txt" --missing --allow-no-sex --out "$TMP/rel" > /dev/null 2>&1
  awk 'NR==1 || 1' "$TMP/rel.imiss" | sed 's/^/  /'
else echo "  plink or $V2.bed missing"; fi
echo; echo "Log: $PREFIX.log"
