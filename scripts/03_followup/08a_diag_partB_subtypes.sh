#!/usr/bin/env bash
# =====================================================================
# 08a_diag_partB_subtypes.sh  –  read-only diagnostic for Part B (TOAST subtype runs on v2 missense) (VM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · RVGRS_Plan.md §3 (Part B) · 24.11.2026
# Question: which files does RV-EXCALIBER step 3 read, and can the two burden matrices simply be
#           column-subset per subtype (CE / LAS / SVD) and step 3 re-run on them?
# Checks:
#   1 files in the v2 missense work dir (names, sizes)
#   2 step-3 shell script + the R scripts it calls: input-file patterns, arguments
#   3 case + gnomAD burden matrix: separator, dimensions, first/last column names, example rows,
#     how gnomAD column names map to case IIDs (raw / "X"+IID / make.names)
#   4 subtype files (stroke == 2): N, overlap with .fam and with both matrices' columns
#   5 any other per-sample file in the work dir that step 3 could depend on
# Nothing is written except the log: /mnt/hgfs/D/WES/FOR_REPORT/partB_v2/08a_log.txt
# Usage: bash 08a_diag_partB_subtypes.sh
# =====================================================================
set -uo pipefail
export LC_ALL=C
WES=/mnt/hgfs/D/WES
RVX=/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main
SHDIR=$RVX/scripts/shell; RDIR=$RVX/scripts/Rscripts
WORK=${WORK:-$HOME/rvx_v2_missense}
DS=my_samples_plink_hg19_v2u
TAG=0.001_0.025_nfe_hcc
FAM=$WES/FOR_REPORT/v2/plink_in/$DS.fam
SUBDIR=/mnt/hgfs/D/MUNICH_GENOTYPES/NEW_IMP/final
CM=$WORK/${DS}_RVBurdenMatrix_${TAG}.txt.gz
GM=$WORK/gnomAD_RVBurdenMatrix_${DS}_${TAG}.txt.gz
S3=$SHDIR/get_SummaryAssociations_iCF_gCF_rvexcaliber.clean.sh
OUT=$WES/FOR_REPORT/partB_v2; mkdir -p "$OUT"
exec > >(tee "$OUT/08a_log.txt") 2>&1
hdr(){ echo; echo "==== $* ===="; }
echo "08a Part B diagnostic | $(date '+%F %T') | WORK=$WORK"

hdr "1 work dir"
ls -la "$WORK" | awk '{printf "  %10s  %s\n", $5, $NF}'
[ -d "$WORK/shell_fixed" ] && { echo "  shell_fixed/:"; ls "$WORK/shell_fixed" | sed 's/^/    /'; }

hdr "2a step-3 shell script ($S3)"
if [ -s "$S3" ]; then nl -ba "$S3" | sed 's/^/  /'; else echo "  NOT FOUND"; ls "$SHDIR" | sed 's/^/    /'; fi

hdr "2b R scripts called by step 3: input/argument lines"
for r in $(grep -o 'Rscripts/[A-Za-z0-9_.-]*\.R' "$S3" 2>/dev/null | sort -u); do
  f=$RVX/scripts/$r; echo "  --- $r ($( [ -s "$f" ] && wc -l < "$f" || echo missing ) lines)"
  [ -s "$f" ] && grep -n -E 'commandArgs|args\[|fread|read\.(table|delim|csv)|readRDS|readLines|paste0?\(.*(txt|gz|fam|bim)|file\.path|colnames|names\(|grep\(|X' "$f" | head -80 | sed 's/^/    /'
done

hdr "3 burden matrices"
for M in "$CM" "$GM"; do
  echo "  --- $(basename "$M")"
  [ -s "$M" ] || { echo "    MISSING"; continue; }
  { zcat "$M" | head -3 > "$OUT/.m3.tmp"; } || true
  H=$(head -1 "$OUT/.m3.tmp")
  TABS=$(printf '%s' "$H" | tr -cd '\t' | wc -c); SPCS=$(printf '%s' "$H" | tr -cd ' ' | wc -c)
  echo "    tabs in header: $TABS | spaces in header: $SPCS"
  SEP=$'\t'; [ "$TABS" -eq 0 ] && SEP=' '
  NC=$(printf '%s\n' "$H" | awk -F"$SEP" '{print NF}')
  echo "    rows (incl header): $({ zcat "$M" || true; } | wc -l) | columns: $NC"
  echo "    first 5 colnames: $(printf '%s\n' "$H" | awk -F"$SEP" '{for(i=1;i<=5&&i<=NF;i++) printf "%s ", $i}')"
  echo "    last 3 colnames:  $(printf '%s\n' "$H" | awk -F"$SEP" '{for(i=NF-2;i<=NF;i++) if(i>0) printf "%s ", $i}')"
  echo "    rows 2-3, first 6 fields:"; awk -F"$SEP" 'NR>1{for(i=1;i<=6&&i<=NF;i++) printf "%s ", $i; print ""}' "$OUT/.m3.tmp" | sed 's/^/      /'
  printf '%s\n' "$H" | awk -F"$SEP" '{for(i=1;i<=NF;i++) print $i}' > "$OUT/.cols_$(basename "$M" .txt.gz).tmp"
done
rm -f "$OUT/.m3.tmp"

hdr "4 subtype files vs fam vs matrix columns (R)"
Rscript - "$FAM" "$SUBDIR" "$OUT/.cols_${DS}_RVBurdenMatrix_${TAG}.tmp" "$OUT/.cols_gnomAD_RVBurdenMatrix_${DS}_${TAG}.tmp" <<'EOF'
suppressMessages(library(data.table))
a <- commandArgs(TRUE); fam <- fread(a[1], header = FALSE, colClasses = "character")
cc <- if (file.exists(a[3])) readLines(a[3]) else character(); gc <- if (file.exists(a[4])) readLines(a[4]) else character()
read_sub <- function(f) { b <- readBin(f, "raw", 2L)
  enc <- if (identical(b, as.raw(c(0xff, 0xfe))) || identical(b, as.raw(c(0xfe, 0xff)))) "UTF-16" else ""
  con <- file(f, encoding = enc); x <- readLines(con, warn = FALSE); close(con)
  d <- fread(text = sub("^﻿", "", sub("\r$", "", x)), header = TRUE, colClasses = "character")
  list(n_rows = nrow(d), cols = paste(names(d), collapse = ","), stroke_tab = paste(names(table(d$stroke)), table(d$stroke), sep = ":", collapse = " "),
       ids = d[IID != "IID" & stroke == "2", unique(IID)]) }
cat(sprintf("  fam: %d samples (IID examples: %s)\n", nrow(fam), paste(head(fam$V2, 3), collapse = ", ")))
cat(sprintf("  case matrix columns: %d | first non-Gene: %s\n", length(cc), paste(head(cc[-1], 3), collapse = ", ")))
cat(sprintf("  gnomAD matrix columns: %d | first non-Gene: %s\n", length(gc), paste(head(gc[-1], 3), collapse = ", ")))
cat(sprintf("  case cols == fam IIDs: %d of %d fam | case cols not in fam: %d\n", sum(fam$V2 %in% cc), nrow(fam), sum(!cc[-1] %in% fam$V2)))
for (nm in c("raw", "X+IID", "make.names")) {
  m <- switch(nm, raw = fam$V2, "X+IID" = paste0("X", fam$V2), make.names = make.names(fam$V2))
  cat(sprintf("  gnomAD cols matched by %-10s: %d of %d fam | unique mapping: %s\n", nm, sum(m %in% gc), nrow(fam), !anyDuplicated(m))) }
for (s in c("ce", "las", "svd")) { f <- file.path(a[2], paste0(s, ".pheno.txt"))
  if (!file.exists(f)) { cat("  ", s, ": file missing\n"); next }
  r <- read_sub(f); ids <- r$ids
  cat(sprintf("  %-3s rows %d | cols %s | stroke %s | cases %d | in fam %d | in case matrix %d | in gnomAD matrix (make.names) %d\n",
      toupper(s), r$n_rows, r$cols, r$stroke_tab, length(ids), sum(ids %in% fam$V2), sum(ids %in% cc), sum(make.names(ids) %in% gc))) }
EOF
rm -f "$OUT"/.cols_*.tmp

hdr "5 other per-sample / DS-named files in work dir (could be step-3 inputs)"
ls "$WORK" | grep -E "$DS|fam|sample|iCF|gCF|CF" | sed 's/^/  /'
echo; echo "08a done | log: $OUT/08a_log.txt"
