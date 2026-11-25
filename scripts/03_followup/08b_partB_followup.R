#!/usr/bin/env Rscript
# =====================================================================
# 08b_partB_followup.R  –  Part B follow-up: why so few genes are tested per subtype + subtype O/E for the all-IS top genes (VM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · RVGRS_Plan.md §3 (B4–B6) · 24.11.2026
# Background (08 log): genes tested CE 1,039 / LAS 56 / SVD 222 although 4,242 / 1,573 / 2,539 genes carry ≥1 allele;
#   ERCC2, PLG, ERBB3, GCNT2, DNMT3A, PTCHD4 untested in every subtype → step 3 prunes genes by some count rule.
# 1 rule check: print the gene-pruning lines of get_iCF_association_test_rvexcaliber.R; per run (all-IS + 3 subtypes)
#   compare tested vs untested genes on raw gnomAD expected count (row sum of gnomAD matrix), observed alleles, carriers
#   → empirical threshold (max of untested vs min of tested)
# 2 descriptive subtype O/E for the all-IS top 10 (B6): O = subtype case alleles; E = all-IS iCF/gCF-adjusted EXP scaled by the
#   subtype's share of the raw gnomAD expectation; Poisson p (excess) and binomial p (subtype vs remaining samples, two-sided)
# Inputs only read. Output: /mnt/hgfs/D/WES/FOR_REPORT/partB_v2/ (B_rule_check.tsv, B_topgenes_subtype_OE.tsv, 08b_log.txt)
# Usage: Rscript 08b_partB_followup.R
# =====================================================================
suppressMessages(library(data.table)); options(width = 220)
WES  <- "/mnt/hgfs/D/WES"; HOME <- Sys.getenv("HOME")
WORK <- file.path(HOME, "rvx_v2_missense"); DS <- "my_samples_plink_hg19_v2u"; TAG <- "0.001_0.025_nfe_hcc"
RVX  <- "/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main"
OUT  <- file.path(WES, "FOR_REPORT/partB_v2")
zz <- file(file.path(OUT, "08b_log.txt"), open = "wt"); sink(zz, split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = ""); hdr <- function(x) cat("\n==== ", x, " ====\n", sep = "")
rdm <- function(f) fread(cmd = paste("zcat", shQuote(f)), header = TRUE, check.names = FALSE, colClasses = list(character = 1))
rdr <- function(f) { d <- fread(f); setnames(d, c("Gene","test_allele_count","test_gnomAD_allele_count_iCFgCFadjust","test_P_rvexcaliber_base_iCFgCFadjust"),
                                         c("Gene","OAC","EXP","P"), skip_absent = TRUE); d[] }
say("08b Part B follow-up | %s", format(Sys.time()))

## 1a pruning code
hdr("1a gene-pruning lines in get_iCF_association_test_rvexcaliber.R")
src <- readLines(file.path(RVX, "scripts/Rscripts/get_iCF_association_test_rvexcaliber.R"))
hit <- which(grepl("prun|filter|thresh|min_|>=|<=| > | < |subset\\(|rowSums|rowMeans|colSums", src) & !grepl("^\\s*#", src))
for (i in hit[hit > 95 & hit < 250]) say("  %4d: %s", i, trimws(src[i]))

## 1b empirical rule
hdr("1b tested vs untested genes per run")
runs <- data.table(run = c("ALL", "CE", "LAS", "SVD"),
  dir = c(WORK, file.path(HOME, paste0("rvx_v2_missense_sub_", c("ce", "las", "svd")))),
  ds = c(DS, paste0(DS, "_", c("ce", "las", "svd"))))
mats <- list(); rule <- list()
for (k in seq_len(nrow(runs))) { r <- runs[k]
  cm <- rdm(file.path(r$dir, paste0(r$ds, "_RVBurdenMatrix_", TAG, ".txt.gz")))
  gm <- rdm(file.path(r$dir, paste0("gnomAD_RVBurdenMatrix_", r$ds, "_", TAG, ".txt.gz")))
  cs <- names(cm)[-1]; gs <- names(gm)[-1]
  cm[, (cs) := lapply(.SD, as.numeric), .SDcols = cs]; gm[, (gs) := lapply(.SD, as.numeric), .SDcols = gs]
  g <- data.table(Gene = cm[[1]], O = rowSums(cm[, ..cs]), carriers = rowSums(cm[, ..cs] > 0))
  g <- merge(g, data.table(Gene = gm[[1]], E_raw = rowSums(gm[, ..gs])), by = "Gene", all.x = TRUE)
  res <- rdr(file.path(r$dir, paste0("rvexcaliber_testing_", r$ds, "_SummaryAssociations_", TAG, "_iCFgCFadjust_rvexcaliber_base.txt")))
  g[, tested := Gene %in% res$Gene]; g <- g[O > 0]
  mats[[r$run]] <- list(g = g, N = length(cs), res = res)
  q <- function(x) paste(signif(quantile(x, c(0, .05, .5, .95, 1), na.rm = TRUE), 3), collapse = "/")
  rule[[r$run]] <- data.table(run = r$run, N = length(cs), genes_O_gt0 = nrow(g), genes_tested = g[tested == TRUE, .N],
    genes_no_gnomAD_row = g[is.na(E_raw), .N],
    Eraw_tested_min = signif(g[tested == TRUE, min(E_raw, na.rm = TRUE)], 3), Eraw_untested_max = signif(g[tested == FALSE, max(E_raw, na.rm = TRUE)], 3),
    O_tested_min = g[tested == TRUE, min(O)], O_untested_max = g[tested == FALSE, max(O)],
    carriers_tested_min = g[tested == TRUE, min(carriers)],
    untested_with_Eraw_ge_1 = g[tested == FALSE & E_raw >= 1, .N], tested_with_Eraw_lt_1 = g[tested == TRUE & E_raw < 1, .N])
  say("  %-3s N %4d | genes O>0 %5d | tested %5d | E_raw (min/5%%/50%%/95%%/max) tested %s | untested %s",
      r$run, length(cs), nrow(g), g[tested == TRUE, .N], q(g[tested == TRUE, E_raw]), q(g[tested == FALSE, E_raw]))
}
rule <- rbindlist(rule); print(rule); fwrite(rule, file.path(OUT, "B_rule_check.tsv"), sep = "\t")
sep_ok <- rule[, Eraw_untested_max < Eraw_tested_min | (untested_with_Eraw_ge_1 == 0 & tested_with_Eraw_lt_1 == 0)]
say("raw-expected rule separates tested/untested in: %s", paste(rule$run[sep_ok], collapse = ", "))

## 2 subtype O/E for the all-IS top genes
hdr("2 all-IS top 10: descriptive O/E within subtypes")
allg <- mats$ALL$g; main <- mats$ALL$res; TOP <- main[order(P)][1:10, Gene]
oe <- rbindlist(lapply(c("CE", "LAS", "SVD"), function(s) { sg <- mats[[s]]$g
  rbindlist(lapply(TOP, function(gn) {
    O_all <- allg[Gene == gn, O]; Er_all <- allg[Gene == gn, E_raw]; EXP_all <- main[Gene == gn, EXP]
    O_s <- sg[Gene == gn, O]; if (!length(O_s)) O_s <- 0
    Er_s <- mats[[s]]$N / mats$ALL$N * Er_all                        # gnomAD matrix values are identical across samples per gene
    E_s <- EXP_all * Er_s / Er_all; E_rest <- EXP_all - E_s; O_rest <- O_all - O_s
    data.table(Gene = gn, subtype = s, N_sub = mats[[s]]$N, O_sub = O_s, E_sub = round(E_s, 2), OE_sub = round(O_s / E_s, 2),
               O_rest = O_rest, OE_rest = round(O_rest / E_rest, 2), OAC_all = O_all, EXP_all = round(EXP_all, 2),
               p_excess_sub = signif(ppois(O_s - 1, E_s, lower.tail = FALSE), 3),
               p_sub_vs_rest = signif(binom.test(O_s, O_all, E_s / EXP_all)$p.value, 3),
               tested_in_subtype_run = gn %in% mats[[s]]$res$Gene) })) }))
fwrite(oe, file.path(OUT, "B_topgenes_subtype_OE.tsv"), sep = "\t")
print(dcast(oe, Gene ~ subtype, value.var = c("O_sub", "E_sub", "p_sub_vs_rest"))[match(TOP, Gene)])
say("smallest subtype-vs-rest p: %s (30 tests; Bonferroni 0.0017)", oe[which.min(p_sub_vs_rest), sprintf("%s %s O %d vs E %.2f, p %.3g", Gene, subtype, O_sub, E_sub, p_sub_vs_rest)])
say("\n08b done | %s", format(Sys.time())); sink(); close(zz)
