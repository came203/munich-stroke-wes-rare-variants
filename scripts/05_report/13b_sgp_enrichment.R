#!/usr/bin/env Rscript
## =============================================================================
## 13b_sgp_enrichment.R – monogenic stroke gene panels SGP1 / SGP2 (Lund; Ilinca et al. 2019/2023) on the v2 result
## Per panel (SGP1, SGP2, union): genes in panel / matched to refGene / tested in v2 (8,790 genes)
##   O/E = sum OAC / sum EXP over tested panel genes (vs gnomAD expectation, as in the main test)
##   permutation p: 10,000 random gene sets of the same size, matched on EXP decile (one-sided, O/E >= observed)
##   nominal genes (p < 0.05 & OR > 1) in panel vs expectation (hypergeometric)
##   same O/E with UKB as reference (11b, exp_from_ukb) – robustness to the gnomAD calibration issue
## Run: Rscript /mnt/hgfs/D/WES/13b_sgp_enrichment.R [SGP1.xlsx] [SGP2.xlsx]
##      (without arguments: env SGP1/SGP2, else searched under /mnt/hgfs/D/WES and ~ up to depth 4)
## Needs R package readxl
## Output: FOR_REPORT/partC/ C_sgp_enrichment.tsv, C_sgp_genes.tsv, 13b_sgp.log
## =============================================================================
suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("readxl", quietly = TRUE)) stop("install.packages('readxl') first")
set.seed(2611)
OUT <- "/mnt/hgfs/D/WES/FOR_REPORT/partC"
REFGENE <- "/mnt/hgfs/D/WES/04_software/annovar/humandb/hg19_refGene.txt"
B <- 10000
sink(file.path(OUT, "13b_sgp.log"), split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
say("13b_sgp_enrichment.R  %s", format(Sys.time()))

find_panel <- function(i, tag) {
  a <- commandArgs(TRUE)
  if (length(a) >= i && file.exists(a[i])) return(a[i])
  e <- Sys.getenv(tag); if (nzchar(e) && file.exists(e)) return(e)
  hit <- system(sprintf("find /mnt/hgfs/D/WES %s -maxdepth 4 -iname '*%s*.xlsx' 2>/dev/null", shQuote(path.expand("~")), tag), intern = TRUE)
  if (!length(hit)) stop(sprintf("%s xlsx not found – pass the path as argument %d", tag, i))
  if (length(hit) > 1) say("several %s files, using the first: %s", tag, paste(hit, collapse = " | "))
  hit[1]
}
sym_ref <- unique(fread(REFGENE, header = FALSE, select = 13)$V13)

## gene column = the column (any sheet) with the largest overlap with refGene symbols
read_panel <- function(f) {
  best <- NULL; bo <- 0
  for (s in readxl::excel_sheets(f)) {
    x <- as.data.table(readxl::read_excel(f, sheet = s, col_names = FALSE, col_types = "text"))
    for (j in seq_along(x)) {
      v <- unique(trimws(na.omit(x[[j]])))
      o <- sum(v %chin% sym_ref)
      if (o > bo) { bo <- o; best <- list(sheet = s, col = j, v = v) }
    }
  }
  v <- best$v[grepl("^[A-Za-z0-9][A-Za-z0-9.-]*$", best$v)]
  v <- v[!(toupper(v) %chin% c("GENE", "GENES", "SYMBOL", "HGNC"))]
  say("%s: sheet '%s', column %d -> %d symbols (%d in refGene)", basename(f), best$sheet, best$col, length(v), sum(v %chin% sym_ref))
  if (any(!v %chin% sym_ref)) say("  not in refGene (alias?): %s", paste(v[!v %chin% sym_ref], collapse = ", "))
  v
}
p1 <- read_panel(find_panel(1, "SGP1")); p2 <- read_panel(find_panel(2, "SGP2"))
panels <- list(SGP1 = p1, SGP2 = p2, SGP_union = union(p1, p2))

w  <- fread(file.path(OUT, "C_weights_v2.tsv"))[, .(Gene, OAC, EXP, OR, P, P_BH)]
u  <- fread(file.path(OUT, "C_11b_gene_ukb_vs_munich.tsv"))[, .(Gene, exp_from_ukb, or_vs_ukb, p_vs_ukb_gt)]
w[, bin := cut(EXP, unique(quantile(EXP, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]
u <- merge(u, w[, .(Gene, OAC)], by = "Gene")
u[, bin := cut(exp_from_ukb, unique(quantile(exp_from_ukb, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)]

perm_oe <- function(univ, genes, o_col, e_col) {
  sel <- univ[Gene %chin% genes]
  obs <- sum(sel[[o_col]]) / sum(sel[[e_col]])
  need <- sel[, .N, by = bin]
  idx <- split(seq_len(nrow(univ)), univ$bin)
  O <- univ[[o_col]]; E <- univ[[e_col]]
  null <- vapply(seq_len(B), function(b) {
    k <- unlist(lapply(seq_len(nrow(need)), function(r) sample(idx[[as.character(need$bin[r])]], need$N[r])))
    sum(O[k]) / sum(E[k])
  }, 0)
  c(n = nrow(sel), O = sum(sel[[o_col]]), E = sum(sel[[e_col]]), OE = obs,
    perm_p = (1 + sum(null >= obs)) / (B + 1), null_median = median(null))
}
nom_all <- w[, sum(P < 0.05 & OR > 1)]
res <- rbindlist(lapply(names(panels), function(pn) {
  gs <- panels[[pn]]
  a <- perm_oe(w, gs, "OAC", "EXP"); b <- perm_oe(u, gs, "OAC", "exp_from_ukb")
  k <- w[Gene %chin% gs, sum(P < 0.05 & OR > 1)]; n <- a[["n"]]
  data.table(panel = pn, genes_in_panel = length(gs), in_refGene = sum(gs %chin% sym_ref), tested = n,
             O = a[["O"]], E_gnomad = a[["E"]], OE_gnomad = a[["OE"]], perm_p_gnomad = a[["perm_p"]],
             E_ukb = b[["E"]], OE_ukb = b[["OE"]], perm_p_ukb = b[["perm_p"]],
             nominal_in_panel = k, nominal_expected = n * nom_all / nrow(w),
             hyper_p = phyper(k - 1, nom_all, nrow(w) - nom_all, n, lower.tail = FALSE),
             nominal_genes = paste(w[Gene %chin% gs & P < 0.05 & OR > 1][order(P), sprintf("%s (p %.2g)", Gene, P)], collapse = ", "))
}))
say("\nresult (O/E = observed / expected qualifying alleles over tested panel genes; perm p one-sided, EXP-decile matched, B = %d):", B)
print(res[, !"nominal_genes"], digits = 3)
for (i in seq_len(nrow(res))) say("%s nominal genes: %s", res$panel[i], if (nzchar(res$nominal_genes[i])) res$nominal_genes[i] else "none")
fwrite(res, file.path(OUT, "C_sgp_enrichment.tsv"), sep = "\t")

gl <- data.table(Gene = panels$SGP_union)[, `:=`(SGP1 = Gene %chin% p1, SGP2 = Gene %chin% p2, in_refGene = Gene %chin% sym_ref)]
gl <- merge(gl, w[, .(Gene, OAC, EXP, OR, P, P_BH)], by = "Gene", all.x = TRUE)
gl <- merge(gl, u[, .(Gene, exp_from_ukb, or_vs_ukb, p_vs_ukb_gt)], by = "Gene", all.x = TRUE)
gl[, tested := !is.na(P)]
setorder(gl, -tested, P)
fwrite(gl, file.path(OUT, "C_sgp_genes.tsv"), sep = "\t")
nom <- gl[tested == TRUE & P < 0.05 & OR > 1, .(Gene, SGP1, SGP2, OAC, EXP = round(EXP, 2), OR = round(OR, 2), P, P_BH,
                                              OR_vs_ukb = round(or_vs_ukb, 2), p_vs_ukb = p_vs_ukb_gt)][order(P)]
say("\nnominal panel genes (p < 0.05, OR > 1 vs gnomAD) with UKB-referenced OR:"); print(nom, digits = 3)
fwrite(nom, file.path(OUT, "C_sgp_nominal.tsv"), sep = "\t")
nomu <- gl[tested == TRUE & !is.na(p_vs_ukb_gt) & p_vs_ukb_gt < 0.05 & or_vs_ukb > 1, .(Gene, OAC, exp_from_ukb = round(exp_from_ukb, 2), OR_vs_ukb = round(or_vs_ukb, 2), p_vs_ukb = p_vs_ukb_gt, OR_gnomad = round(OR, 2), P)][order(p_vs_ukb)]
say("panel genes nominal vs UKB reference (descriptive):"); print(nomu, digits = 3)
say("\nkey genes:"); print(gl[Gene %chin% c("NOTCH3", "HTRA1", "COL4A1", "COL4A2", "TREX1", "PLG", "ERCC2")], digits = 3)
say("untested panel genes (no RV-EXCALIBER test, e.g. expected count < 1): %d", gl[tested == FALSE, .N])
say("wrote C_sgp_enrichment.tsv, C_sgp_genes.tsv")
sink()
