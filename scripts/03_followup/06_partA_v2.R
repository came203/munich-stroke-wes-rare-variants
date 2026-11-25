#!/usr/bin/env Rscript
# =====================================================================
# 06_partA_v2.R  –  Part A on the v2 missense result: variant breakdown of the top genes (VM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · RVGRS_Plan.md §2 (A1, A2) on v2 (24.11.2026)
#   v2 = AQ20 joint call → clean hg19 rebuild → N 1,535 (relatives removed) → RV-EXCALIBER missense-only (M-CAP > 0.025)
# Genes: Bonferroni/FDR genes (ERCC2, PLG) + next nominal top genes (TOPN in total)
# Per variant: AC, het/hom, missing, cohort AF vs gnomAD NFE AF (fold), crude Poisson p, share of gene alleles,
#              M-CAP, AAChange, exon; per gene: n variants, carriers, max single-variant share (A2 flag)
# Reference: single-variant concentration across all tested genes; carriers by TOAST subtype and sample-ID class
# Inputs only read. Output: /mnt/hgfs/D/WES/FOR_REPORT/partA_v2/
#   A2_topgene_summary.tsv, A2_topgene_variants.tsv, A2_allgenes_concentration.tsv, A2_concentration_by_OAC.tsv,
#   A2_carriers_by_subtype.tsv, A2_carriers_by_idclass.tsv, A2_carriers_INDIVIDUAL_LEVEL.tsv (local only),
#   A2_summary.txt, A2_log.txt
# Usage: Rscript 06_partA_v2.R            (TOPN via env: TOPN=12 Rscript 06_partA_v2.R)
# =====================================================================
suppressMessages(library(data.table))
options(width = 230)
WES   <- "/mnt/hgfs/D/WES"
FINAL <- file.path(WES, "05_rvexcaliber_output_v2_missense")
WORK  <- file.path(Sys.getenv("HOME"), "rvx_v2_missense")
DS    <- "my_samples_plink_hg19_v2u"
TAG   <- "0.001_0.025_nfe_hcc"
PLK   <- file.path(WES, "FOR_REPORT/v2/plink_in", DS)                 # genotypes (IDs chr:pos:ref:alt)
pick  <- function(...) { for (f in c(...)) if (file.exists(f)) return(f); stop("none found: ", paste(c(...), collapse = " | ")) }
MAIN  <- pick(file.path(FINAL, paste0("rvexcaliber_testing_", DS, "_SummaryAssociations_", TAG, "_iCFgCFadjust_rvexcaliber_base.txt")),
              file.path(WORK,  paste0("rvexcaliber_testing_", DS, "_SummaryAssociations_", TAG, "_iCFgCFadjust_rvexcaliber_base.txt")))
MAINF <- sub("SummaryAssociations_", "SummaryAssociations_allele_filter_", MAIN)
VLF   <- pick(file.path(FINAL, paste0(DS, "_variant_list_extracted_", TAG, ".txt.gz")), file.path(WORK, paste0(DS, "_variant_list_extracted_", TAG, ".txt.gz")))
MAF   <- pick(file.path(WORK, paste0(DS, ".hg19_multianno.txt.gz")), file.path(WORK, paste0(DS, ".hg19_multianno.txt")),
              file.path(FINAL, paste0(DS, ".hg19_multianno.txt")))
SUBDIR <- "/mnt/hgfs/D/MUNICH_GENOTYPES/NEW_IMP/final"
OUT   <- file.path(WES, "FOR_REPORT/partA_v2"); TMPD <- file.path(OUT, "tmp")
PLINK <- Sys.which("plink")
TOPN  <- as.integer(Sys.getenv("TOPN", "10"))
SHARE_FLAG <- 0.5
dir.create(TMPD, recursive = TRUE, showWarnings = FALSE)
zz <- file(file.path(OUT, "A2_log.txt"), open = "wt"); sink(zz, split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
hdr <- function(x) cat("\n==== ", x, " ====\n", sep = "")
zread <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = paste("zcat", shQuote(f)), ...) else fread(f, ...)
run_plink <- function(args, tag) { so <- file.path(TMPD, paste0(tag, ".stdout"))
  if (system2(PLINK, args, stdout = so, stderr = so) != 0) stop("plink failed (", tag, ") – see ", so) }
say("Part A v2 | %s | main: %s", format(Sys.time()), basename(MAIN))

## 1 main result + gene selection
hdr("1 main result (v2 missense)")
main <- fread(MAIN)
setnames(main, c("Gene", "test_allele_count", "test_gnomAD_allele_count_iCFgCFadjust", "test_P_rvexcaliber_base_iCFgCFadjust"), c("Gene", "OAC", "EXP", "P"))
main[, `:=`(FDR = p.adjust(P, "BH"), OR = OAC / EXP, in_allele_filter_run = Gene %in% fread(MAINF)$Gene)]
TOP <- main[order(P)][1:TOPN, Gene]
say("genes %d | Bonferroni %d | FDR<0.05 %d | top %d: %s", nrow(main), main[P < 0.05 / .N, .N], main[FDR < 0.05, .N], TOPN, paste(TOP, collapse = ", "))

## 2 qualifying variants + cohort counts
hdr("2 qualifying variants and cohort allele counts")
vl <- zread(VLF, header = FALSE, colClasses = list(character = c(1, 2, 4)))
setnames(vl, c("ID", "A_cnt", "V3", "Gene", "AF_gnomAD_nfe", "gene_cumAF"))
vl[, AF_gnomAD_nfe := as.numeric(AF_gnomAD_nfe)]
vl[, c("chr", "pos", "ref", "alt") := tstrsplit(ID, ":", fixed = TRUE)][, pos := as.integer(pos)]
bim <- fread(paste0(PLK, ".bim"), header = FALSE, select = 2, col.names = "ID")
vl[, in_bim := ID %in% bim$ID]
say("variant list rows %d | genes %d | in bim %d | not in bim %d | gnomAD-absent %d",
    nrow(vl), uniqueN(vl$Gene), vl[in_bim == TRUE, .N], vl[in_bim == FALSE, .N], vl[AF_gnomAD_nfe == 0, .N])
writeLines(unique(vl[in_bim == TRUE, ID]), file.path(TMPD, "qv.ids"))
run_plink(c("--bfile", PLK, "--extract", file.path(TMPD, "qv.ids"), "--freq", "counts", "--allow-no-sex", "--out", file.path(TMPD, "qv")), "freq_qv")
frq <- fread(file.path(TMPD, "qv.frq.counts"), colClasses = list(character = c("SNP", "A1", "A2")))[!duplicated(SNP)]
vl <- merge(vl, frq[, .(ID = SNP, A1, A2, C1, C2, G0)], by = "ID", all.x = TRUE)
vl[, AC := fifelse(A1 == A_cnt, C1, fifelse(A2 == A_cnt, C2, NA_integer_))][, AN := C1 + C2][, AF_cohort := AC / AN]
t1 <- table(cut(vl$AC, c(-Inf, 0, 1, 2, 5, 10, Inf), labels = c("0", "1", "2", "3-5", "6-10", ">10")))
say("AC spectrum qualifying: %s | max cohort AF %.4f", paste(names(t1), t1, sep = ":", collapse = " "), vl[, max(AF_cohort, na.rm = TRUE)])

## 3 concentration, all genes
hdr("3 single-variant concentration, all genes")
g <- vl[!is.na(AC) & AC > 0, .(n_var_obs = .N, OAC_recomp = sum(AC), max_AC = max(AC), top_var = ID[which.max(AC)],
                              AC_gnomAD_absent = sum(AC[AF_gnomAD_nfe == 0])), by = Gene]
g[, `:=`(max_share = max_AC / OAC_recomp, share_gnomAD_absent = AC_gnomAD_absent / OAC_recomp)]
g <- merge(main[, .(Gene, OAC, EXP, OR, P, FDR, in_allele_filter_run)], g, by = "Gene", all.x = TRUE)
say("recomputed OAC == main OAC: %d of %d", g[OAC == OAC_recomp, .N], nrow(g))
conc <- rbindlist(lapply(c(5, 10, 20, 40), function(k) g[OAC >= k & !is.na(max_share), .(min_OAC = k, n_genes = .N,
  pct_share_gt50 = round(100 * mean(max_share > SHARE_FLAG), 1), median_max_share = round(median(max_share), 2),
  median_share_gnomAD_absent = round(median(share_gnomAD_absent), 2))]))
print(conc)

## 4 top genes per variant
hdr("4 top genes – per variant")
tv <- vl[Gene %in% TOP & in_bim == TRUE]
fwrite(unique(tv[, .(ID, A_cnt)]), file.path(TMPD, "top.ra"), sep = "\t", col.names = FALSE)
writeLines(unique(tv$ID), file.path(TMPD, "top.ids"))
run_plink(c("--bfile", PLK, "--extract", file.path(TMPD, "top.ids"), "--recode", "A", "--recode-allele", file.path(TMPD, "top.ra"),
            "--allow-no-sex", "--out", file.path(TMPD, "top")), "recode_top")
raw <- fread(file.path(TMPD, "top.raw"), colClasses = list(character = c("FID", "IID")))
long <- melt(raw, id.vars = "IID", measure.vars = setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")),
             variable.name = "col", value.name = "dos")
long[, ID := sub("_[^_]+$", "", as.character(col))]
long <- merge(long, unique(tv[, .(ID, Gene)]), by = "ID", allow.cartesian = TRUE)
vstat <- long[, .(n_het = sum(dos == 1, na.rm = TRUE), n_hom = sum(dos == 2, na.rm = TRUE), n_miss = sum(is.na(dos))), by = .(Gene, ID)]
ma <- zread(MAF, select = c("Chr", "Start", "Ref", "Alt", "ExonicFunc.refGene", "AAChange.refGene", "M-CAP_score"),
            colClasses = list(character = c("Chr", "Ref", "Alt", "M-CAP_score")))
ma[, ID := paste(Chr, Start, Ref, Alt, sep = ":")]
ma <- ma[ID %in% tv$ID][!duplicated(ID)]
ma[, AAChange_first := sub(",.*$", "", AAChange.refGene)][, exon := sub("^[^:]*:[^:]*:(exon[0-9]+).*$", "\\1", AAChange_first)]
tvar <- merge(tv[, .(Gene, ID, AC, AN, AF_cohort, AF_gnomAD_nfe)], vstat, by = c("Gene", "ID"), all.x = TRUE)
tvar <- merge(tvar, ma[, .(ID, AAChange_first, exon, MCAP = `M-CAP_score`)], by = "ID", all.x = TRUE)
tvar <- merge(tvar, g[, .(Gene, OAC_gene = OAC_recomp)], by = "Gene", all.x = TRUE)
tvar[, `:=`(share_of_gene = AC / OAC_gene, exp_AC_gnomAD = AN * AF_gnomAD_nfe, gnomAD_absent = AF_gnomAD_nfe == 0,
            fold_vs_gnomAD = fifelse(AF_gnomAD_nfe > 0, AF_cohort / AF_gnomAD_nfe, Inf))]
tvar[, p_poisson_crude := fifelse(exp_AC_gnomAD > 0 & AC > 0, ppois(AC - 1, exp_AC_gnomAD, lower.tail = FALSE), NA_real_)]
tvar <- tvar[AC > 0]
setorderv(tvar, c("Gene", "AC"), c(1L, -1L), na.last = TRUE)
for (gn in TOP[TOP %in% tvar$Gene][1:min(4, length(TOP))]) {
  cat(sprintf("\n  -- %s (OAC %d, EXP %.2f, p %.2e)\n", gn, main[Gene == gn, OAC], main[Gene == gn, EXP], main[Gene == gn, P]))
  print(tvar[Gene == gn, .(ID, AAChange = substr(AAChange_first, 1, 48), MCAP, AC, share = round(share_of_gene, 2), n_hom, n_miss,
                           AF_cohort = signif(AF_cohort, 2), AF_gnomAD = signif(AF_gnomAD_nfe, 2), fold = round(fold_vs_gnomAD))])
  ex <- tvar[Gene == gn, .(AC = sum(AC), n_var = .N), by = exon][order(-AC)]
  cat("     alleles by exon: ", paste(sprintf("%s:%d(%d var)", ex$exon, ex$AC, ex$n_var), collapse = "  "), "\n")
}
carr <- long[dos > 0, .(n_carriers = uniqueN(IID)), by = Gene]
topv <- tvar[, .SD[1], by = Gene]
tsum <- merge(g[Gene %in% TOP], carr, by = "Gene", all.x = TRUE)
tsum <- merge(tsum, topv[, .(Gene, top_AAChange = AAChange_first, top_AF_cohort = AF_cohort, top_AF_gnomAD = AF_gnomAD_nfe, top_gnomAD_absent = gnomAD_absent)], by = "Gene", all.x = TRUE)
tsum[, flag_single_variant_gt50 := max_share > SHARE_FLAG]
setorder(tsum, P)
hdr("4b top-gene summary")
print(tsum[, .(Gene, OAC, EXP = round(EXP, 2), OR = round(OR, 1), P = signif(P, 3), in_allele_filter_run, n_var_obs, n_carriers, max_AC,
               max_share = round(max_share, 2), flag_single_variant_gt50, share_gnomAD_absent = round(share_gnomAD_absent, 2))])

## 5 carriers by subtype / ID class
hdr("5 carriers by TOAST subtype and sample-ID class")
read_sub <- function(f) { b <- readBin(f, "raw", 2L)
  enc <- if (identical(b, as.raw(c(0xff, 0xfe))) || identical(b, as.raw(c(0xfe, 0xff)))) "UTF-16" else ""
  con <- file(f, encoding = enc); x <- readLines(con, warn = FALSE); close(con)
  d <- fread(text = sub("^﻿", "", sub("\r$", "", x)), header = TRUE, colClasses = "character"); d[IID != "IID" & stroke == "2", unique(IID)] }
fam <- fread(paste0(PLK, ".fam"), header = FALSE, colClasses = "character")
samp <- data.table(IID = fam$V2, subtype = "none")
for (s in c("ce", "las", "svd")) { ids <- read_sub(file.path(SUBDIR, paste0(s, ".pheno.txt"))); samp[IID %in% ids, subtype := toupper(s)] }
samp[, id_class := sub("^[0-9]+_(.*)_[^_]+$", "\\1", IID)][id_class == IID, id_class := "other"][grepl("^[0-9]+$", IID), id_class := "numeric"]
say("samples by subtype: %s | by ID class: %s", paste(samp[, .N, by = subtype][, paste(subtype, N, sep = ":")], collapse = " "),
    paste(samp[, .N, by = id_class][, paste(id_class, N, sep = ":")], collapse = " "))
cs <- merge(unique(long[dos > 0, .(Gene, IID)]), samp, by = "IID")
tab_by <- function(var) rbindlist(lapply(TOP, function(gn) {
  tot <- samp[, .N, by = var]; setnames(tot, var, "level"); cg <- cs[Gene == gn, .N, by = var]; setnames(cg, var, "level")
  x <- merge(tot, cg, by = "level", all.x = TRUE, suffixes = c("_total", "_carriers")); x[is.na(N_carriers), N_carriers := 0L]
  nc <- sum(x$N_carriers); nt <- sum(x$N_total); pc <- x$N_carriers; pt <- x$N_total
  x[, exp_carriers := round(nc * N_total / nt, 2)]
  x[, p_fisher := sapply(seq_along(pc), function(i) fisher.test(matrix(c(pc[i], pt[i] - pc[i], nc - pc[i], nt - pt[i] - (nc - pc[i])), 2))$p.value)]
  x[, Gene := gn][] }))
tab_sub <- tab_by("subtype"); tab_idc <- tab_by("id_class")
print(dcast(tab_sub, Gene ~ level, value.var = "N_carriers")); print(dcast(tab_idc, Gene ~ level, value.var = "N_carriers"))
say("smallest Fisher p (subtype): %s", tab_sub[which.min(p_fisher), sprintf("%s %s p=%.3g", Gene, level, p_fisher)])

## outputs + summary
fwrite(tsum, file.path(OUT, "A2_topgene_summary.tsv"), sep = "\t"); fwrite(tvar, file.path(OUT, "A2_topgene_variants.tsv"), sep = "\t")
fwrite(g[order(P)], file.path(OUT, "A2_allgenes_concentration.tsv"), sep = "\t"); fwrite(conc, file.path(OUT, "A2_concentration_by_OAC.tsv"), sep = "\t")
fwrite(tab_sub, file.path(OUT, "A2_carriers_by_subtype.tsv"), sep = "\t"); fwrite(tab_idc, file.path(OUT, "A2_carriers_by_idclass.tsv"), sep = "\t")
fwrite(merge(long[dos > 0, .(Gene, ID, IID, dos)], samp, by = "IID")[order(Gene, ID, IID)], file.path(OUT, "A2_carriers_INDIVIDUAL_LEVEL.tsv"), sep = "\t")
sig <- main[P < 0.05 / nrow(main), Gene]
summ <- c(sprintf("v2 missense (N=1,535, %d genes): exome-wide genes %s.", nrow(main), paste(sig, collapse = ", ")),
  unlist(lapply(sig, function(gn) { x <- tsum[Gene == gn]
    sprintf("%s: %d alleles vs %.1f expected (OR %.1f, p %.1e) from %d variants in %d carriers; largest single-variant share %.0f%% (%s); %.0f%% of alleles from gnomAD-absent variants; allele-filter run: %s.",
            gn, x$OAC, x$EXP, x$OR, x$P, x$n_var_obs, x$n_carriers, 100 * x$max_share, x$top_AAChange, 100 * x$share_gnomAD_absent, x$in_allele_filter_run) })),
  sprintf("Reference: among genes with OAC>=20 (n=%d), %.1f%% have a single-variant share >50%%.", conc[min_OAC == 20, n_genes], conc[min_OAC == 20, pct_share_gt50]))
writeLines(summ, file.path(OUT, "A2_summary.txt")); hdr("SUMMARY"); cat(summ, sep = "\n")
sink(); close(zz)
