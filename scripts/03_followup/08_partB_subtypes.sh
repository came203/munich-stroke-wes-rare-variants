#!/usr/bin/env bash
# =====================================================================
# 08_partB_subtypes.sh  –  Part B: TOAST subtype RV-EXCALIBER runs on the v2 missense result (VM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · RVGRS_Plan.md §3 (B1–B6) · 24.11.2026
# Approach (RVGRS_Plan 7b #17: step 3 reads only the two burden matrices):
#   per subtype (CE, LAS, SVD; stroke == 2 in {ce,las,svd}.pheno.txt, ∩ the 1,535 v2 samples)
#   1 column-subset the case matrix and the gnomAD matrix of the v2 missense run to the subtype's samples
#     (genes with 0 case alleles in the subtype are dropped from the case matrix, as step 1A would do)
#   2 write them under a new dataset name (my_samples_plink_hg19_v2u_<s>) into a separate work dir
#   3 run RV-EXCALIBER step 3 only (same arguments as the main run: 0.001 0.025 nfe hcc fulladjust, MIGen ranking)
#   4 summary: N, genes tested, λ, p<0.05, Bonferroni, FDR, allele-filtered run, top 3 genes; ERCC2/PLG + all-IS top 10 per subtype
# Note: qualifying variants = those of the all-IS run (internal MAF ≤ 0.01 in N = 1,535); re-running step 1A per subtype
#       would apply the internal MAF filter within e.g. 43 LAS cases (1 allele = 1.2% > 1%) and drop every carrier.
# Inputs only read. Outputs:
#   work dirs  $HOME/rvx_v2_missense_sub_{ce,las,svd}
#   report     /mnt/hgfs/D/WES/FOR_REPORT/partB_v2/  (B_summary_by_subtype.tsv, B_top_genes_<S>.tsv,
#              B_allIS_topgenes_by_subtype.tsv, B_summary.txt, 08_log.txt, results/ = copied step-3 outputs)
# Usage: nohup bash 08_partB_subtypes.sh > 08_partB.out 2>&1 &      (subset only: ONLY_SUBSET=1 bash 08_partB_subtypes.sh)
# =====================================================================
set -euo pipefail
export LC_ALL=C
WES=/mnt/hgfs/D/WES
RVX=/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main
SHDIR=$RVX/scripts/shell                                  # same copy as 05c step 3
WORK=${WORK:-$HOME/rvx_v2_missense}
DS=my_samples_plink_hg19_v2u
GMAF=0.001; MCAP=0.025; ETH=nfe; COV=hcc; ADJ=fulladjust
TAG=${GMAF}_${MCAP}_${ETH}_${COV}
FAM=$WES/FOR_REPORT/v2/plink_in/$DS.fam
SUBDIR=/mnt/hgfs/D/MUNICH_GENOTYPES/NEW_IMP/final
CM=$WORK/${DS}_RVBurdenMatrix_${TAG}.txt.gz
GM=$WORK/gnomAD_RVBurdenMatrix_${DS}_${TAG}.txt.gz
MAIN=$WORK/rvexcaliber_testing_${DS}_SummaryAssociations_${TAG}_iCFgCFadjust_rvexcaliber_base.txt
S3=get_SummaryAssociations_iCF_gCF_rvexcaliber.clean.sh
OUT=$WES/FOR_REPORT/partB_v2; mkdir -p "$OUT/results"
exec > >(tee -a "$OUT/08_log.txt") 2>&1
step(){ echo; echo "== $(date '+%F %T') $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ [ -s "$1" ] || die "missing/empty: $1"; }

step "0 checks"
for f in "$CM" "$GM" "$MAIN" "$FAM" "$SHDIR/$S3" "$SUBDIR/ce.pheno.txt" "$SUBDIR/las.pheno.txt" "$SUBDIR/svd.pheno.txt"; do need "$f"; done
[ "$(wc -l < "$FAM")" -eq 1535 ] || die "expected 1535 samples in $FAM"
Rscript -e 'for(p in c("data.table","R.utils")) if(!requireNamespace(p, quietly=TRUE)) stop("R package missing: ", p)' || die "R packages"
echo "   inputs ok | step-3 script: $SHDIR/$S3"

step "1 subset burden matrices per subtype"
cat > "$OUT/08_subset.R" <<'EOF'
suppressMessages(library(data.table))
a <- commandArgs(TRUE); cm_f <- a[1]; gm_f <- a[2]; fam_f <- a[3]; subdir <- a[4]; ds <- a[5]; tag <- a[6]; home <- a[7]; out <- a[8]
sep_of <- function(f) { h <- readLines(gzcon(file(f, "rb")), n = 1); if (grepl("\t", h)) "\t" else " " }
rd <- function(f) fread(cmd = paste("zcat", shQuote(f)), sep = sep_of(f), header = TRUE, check.names = FALSE, colClasses = list(character = 1))
read_sub <- function(f) { b <- readBin(f, "raw", 2L)
  enc <- if (identical(b, as.raw(c(0xff, 0xfe))) || identical(b, as.raw(c(0xfe, 0xff)))) "UTF-16" else ""
  con <- file(f, encoding = enc); x <- readLines(con, warn = FALSE); close(con)
  d <- fread(text = sub("^﻿", "", sub("\r$", "", x)), header = TRUE, colClasses = "character"); d[IID != "IID" & stroke == "2", unique(IID)] }
fam <- fread(fam_f, header = FALSE, colClasses = "character")$V2
cm <- rd(cm_f); gm <- rd(gm_f); scm <- sep_of(cm_f); sgm <- sep_of(gm_f)
cat(sprintf("   case matrix %d genes x %d cols (first col '%s') | gnomAD matrix %d genes x %d cols (first col '%s')\n",
            nrow(cm), ncol(cm) - 1, names(cm)[1], nrow(gm), ncol(gm) - 1, names(gm)[1]))
## column-name scheme per matrix (raw IID, "X"+IID, make.names(IID))
scheme <- function(cols) { sc <- list(raw = identity, X = function(i) paste0("X", i), make.names = make.names)
  hits <- sapply(sc, function(f) sum(f(fam) %in% cols)); best <- names(which.max(hits))
  cat(sprintf("     matches of the 1,535 fam IIDs: %s -> using '%s'\n", paste(names(hits), hits, sep = ":", collapse = " "), best))
  if (hits[best] != length(fam)) stop("not all fam IIDs found among matrix columns (", hits[best], " of ", length(fam), ")")
  sc[[best]] }
cat("   case matrix:\n");   fc <- scheme(names(cm))
cat("   gnomAD matrix:\n"); fg <- scheme(names(gm))
if (ncol(cm) - 1 != length(fam) || ncol(gm) - 1 != length(fam)) cat("   NOTE: matrices have columns beyond the fam samples – these are dropped in the subsets\n")
num_c <- names(cm)[-1]; cm[, (num_c) := lapply(.SD, as.numeric), .SDcols = num_c]
ns <- list()
for (s in c("ce", "las", "svd")) {
  ids <- intersect(read_sub(file.path(subdir, paste0(s, ".pheno.txt"))), fam)
  wd <- file.path(home, paste0("rvx_v2_missense_sub_", s)); dir.create(wd, showWarnings = FALSE)
  dss <- paste0(ds, "_", s)
  c1 <- cm[, c(names(cm)[1], fc(ids)), with = FALSE]
  rs <- rowSums(c1[, -1]); c1 <- c1[rs > 0]
  g1 <- gm[, c(names(gm)[1], fg(ids)), with = FALSE]
  fwrite(c1, file.path(wd, paste0(dss, "_RVBurdenMatrix_", tag, ".txt.gz")), sep = scm, quote = FALSE)
  fwrite(g1, file.path(wd, paste0("gnomAD_RVBurdenMatrix_", dss, "_", tag, ".txt.gz")), sep = sgm, quote = FALSE)
  writeLines(ids, file.path(out, paste0("ids_", s, ".txt")))
  ns[[s]] <- data.table(subtype = toupper(s), N = length(ids), genes_case_matrix = nrow(c1), case_alleles = sum(rs),
                        genes_gnomAD_matrix = nrow(g1), workdir = wd)
  cat(sprintf("   %-3s N %d | case matrix: %d genes with >=1 allele (%d alleles) | gnomAD matrix %d genes x %d cols\n",
              toupper(s), length(ids), nrow(c1), sum(rs), nrow(g1), ncol(g1) - 1))
}
fwrite(rbindlist(ns), file.path(out, "B_subset_info.tsv"), sep = "\t")
EOF
Rscript "$OUT/08_subset.R" "$CM" "$GM" "$FAM" "$SUBDIR" "$DS" "$TAG" "$HOME" "$OUT" || die "subsetting failed"
[ "${ONLY_SUBSET:-0}" = "1" ] && { echo "ONLY_SUBSET=1 -> stop after subsetting"; exit 0; }

step "2 step 3 per subtype"
declare -A ST
for s in ce las svd; do
  WS=$HOME/rvx_v2_missense_sub_$s; DSS=${DS}_$s
  RES=$WS/rvexcaliber_testing_${DSS}_SummaryAssociations_${TAG}_iCFgCFadjust_rvexcaliber_base.txt
  echo "-- ${s^^} | $WS"
  if [ -s "$RES" ]; then echo "   result present – skipped"; ST[$s]=ok; continue; fi
  echo "   \$ bash $S3 $WS $DSS NA $GMAF $MCAP $ETH $COV $ADJ"
  if ( cd "$SHDIR" && bash "$S3" "$WS" "$DSS" NA "$GMAF" "$MCAP" "$ETH" "$COV" "$ADJ" ) > "$WS/08_step3.log" 2>&1 && [ -s "$RES" ]; then
    ST[$s]=ok; echo "   done"
  else
    ST[$s]=failed; echo "   FAILED – last lines of $WS/08_step3.log:"; { tail -15 "$WS/08_step3.log" | sed 's/^/     /'; } || true
  fi
  cp "$WS"/rvexcaliber_testing_${DSS}_* "$OUT/results/" 2>/dev/null || true
done
echo "   status: CE ${ST[ce]} | LAS ${ST[las]} | SVD ${ST[svd]}"

step "3 summary"
cat > "$OUT/08_summary.R" <<'EOF'
suppressMessages(library(data.table)); options(width = 200)
a <- commandArgs(TRUE); home <- a[1]; ds <- a[2]; tag <- a[3]; main_f <- a[4]; out <- a[5]
rd <- function(f) { d <- fread(f); setnames(d, c("Gene","test_allele_count","test_gnomAD_allele_count_iCFgCFadjust","test_P_rvexcaliber_base_iCFgCFadjust"),
                                         c("Gene","OAC","EXP","P"), skip_absent = TRUE); d[, FDR := p.adjust(P, "BH")][, OR := OAC / EXP][] }
info <- fread(file.path(out, "B_subset_info.tsv")); main <- rd(main_f); TOP <- main[order(P)][1:10, Gene]
rows <- list(); cross <- list()
for (s in c("ce", "las", "svd")) {
  wd <- file.path(home, paste0("rvx_v2_missense_sub_", s)); dss <- paste0(ds, "_", s)
  f <- file.path(wd, paste0("rvexcaliber_testing_", dss, "_SummaryAssociations_", tag, "_iCFgCFadjust_rvexcaliber_base.txt"))
  ff <- sub("SummaryAssociations_", "SummaryAssociations_allele_filter_", f)
  if (!file.exists(f)) { rows[[s]] <- data.table(subtype = toupper(s), N = info[subtype == toupper(s), N], status = "step 3 failed"); next }
  r <- rd(f); n <- nrow(r); lam <- median(qchisq(r$P, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
  t3 <- r[order(P)][1:min(3, n)]
  rf <- if (file.exists(ff)) rd(ff) else NULL
  rows[[s]] <- data.table(subtype = toupper(s), N = info[subtype == toupper(s), N], status = "ok", genes_tested = n, lambda = round(lam, 3),
    n_p05 = r[P < 0.05, .N], exp_p05 = round(0.05 * n), bonf_thr = signif(0.05 / n, 3), n_bonf = r[P < 0.05 / n, .N], n_fdr05 = r[FDR < 0.05, .N],
    genes_allele_filter = if (is.null(rf)) NA_integer_ else nrow(rf),
    n_bonf_allele_filter = if (is.null(rf)) NA_integer_ else rf[P < 0.05 / nrow(rf), .N],
    top3 = paste(sprintf("%s %d/%.1f p=%.2g", t3$Gene, t3$OAC, t3$EXP, t3$P), collapse = "; "))
  fwrite(r[order(P)][1:min(50, n)][, `:=`(EXP = round(EXP, 2), OR = round(OR, 2), in_allele_filter_run = if (is.null(rf)) NA else Gene %in% rf$Gene)],
         file.path(out, paste0("B_top_genes_", toupper(s), ".tsv")), sep = "\t")
  cross[[s]] <- merge(data.table(Gene = TOP), r[, .(Gene, OAC, EXP = round(EXP, 2), OR = round(OR, 2), P = signif(P, 3))], all.x = TRUE)[, subtype := toupper(s)]
}
sm <- rbindlist(rows, fill = TRUE); cr <- rbindlist(cross)
cr <- merge(cr, main[, .(Gene, OAC_allIS = OAC, EXP_allIS = round(EXP, 2), P_allIS = signif(P, 3))], by = "Gene")
cr[, share_of_allIS_alleles := round(OAC / OAC_allIS, 2)]
info_n <- setNames(info$N, info$subtype); cr[, share_of_N := round(info_n[subtype] / 1535, 3)]
fwrite(sm, file.path(out, "B_summary_by_subtype.tsv"), sep = "\t")
fwrite(cr[order(P_allIS, subtype)], file.path(out, "B_allIS_topgenes_by_subtype.tsv"), sep = "\t")
txt <- c(sprintf("Part B v2 (missense, step 3 on column-subset matrices) | %s", format(Sys.time())), "",
         capture.output(print(sm)), "", "all-IS top 10 genes within subtypes (share_of_allIS_alleles vs share_of_N = subtype's share of samples):",
         capture.output(print(dcast(cr, Gene + OAC_allIS + P_allIS ~ subtype, value.var = c("OAC", "P"))[order(P_allIS)])))
writeLines(txt); writeLines(txt, file.path(out, "B_summary.txt"))
EOF
Rscript "$OUT/08_summary.R" "$HOME" "$DS" "$TAG" "$MAIN" "$OUT" || die "summary failed"
echo; echo "08 done | $(date '+%F %T') | report files: $OUT"
