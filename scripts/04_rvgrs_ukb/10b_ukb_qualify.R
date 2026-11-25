#!/usr/bin/env Rscript
## =============================================================================
## 10b_ukb_qualify.R – Part C step 2b: apply the D4 rule (partC_qualify_fun.R, C_filter_rule.tsv,
##   C_class_set.txt) to the UKB ANNOVAR output of 10_ukb_qualify.sh -> qualifying variants
## Run on the VM after "10_ukb_qualify.sh all":
##   Rscript /mnt/hgfs/D/WES/10b_ukb_qualify.R            (all chr)   or   ... 10b_ukb_qualify.R 22
## Writes: ~/partC_ukb/C_ukb_qualifying_variants.tsv.gz (variant level, local; input to 11)
##         FOR_REPORT/partC/C_ukb_gene_counts.tsv, C_ukb_funnel.tsv, 10b_ukb_qualify.log
## Note: internal AF in UKB = ALT_CTS/OBS_CT in the analysis set (alt = minor after pre-filter),
##       counted allele = ALT (as in Munich D4: counted allele == minor allele for all 53,489)
## =============================================================================
suppressPackageStartupMessages(library(data.table))
setDTthreads(8)

BASE    <- "/mnt/hgfs/D/WES"
OUT     <- file.path(BASE, "FOR_REPORT/partC")
W       <- path.expand("~/partC_ukb")
HCC_BED <- "/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main/gnomAD_211_exomes_hg19_coverage/ALL_CHROM_gnomad_refGene_exons_highcov20X.mod.bed"
F_VL    <- path.expand("~/rvx_v2_missense/my_samples_plink_hg19_v2u_variant_list_extracted_0.001_0.025_nfe_hcc.txt.gz")
N_MUN   <- 1535L
TOP2    <- c("ERCC2", "PLG")

a  <- commandArgs(FALSE); sf <- sub("^--file=", "", a[grep("^--file=", a)])
SCRIPT_DIR <- if (length(sf) == 1) dirname(normalizePath(sf)) else getwd()
source(file.path(SCRIPT_DIR, "partC_qualify_fun.R"))
CHRS <- as.integer(commandArgs(TRUE)); if (!length(CHRS)) CHRS <- 1:22

sink(file.path(OUT, "10b_ukb_qualify.log"), split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
say("10b_ukb_qualify.R  %s  chr: %s", format(Sys.time()), paste(CHRS, collapse = ","))

## ---------- rule ----------
rule <- fread(file.path(OUT, "C_filter_rule.tsv"), colClasses = "character")
rv <- setNames(rule$value, rule$param)
if (rv[["D4_verdict"]] != "PASS") stop("D4 not PASS")
class_set <- readLines(file.path(OUT, "C_class_set.txt"))
hcc_col <- if (rv[["hcc_mode"]] == "rvx") "in_hcc_rvx" else "in_hcc"
say("rule: M-CAP > %s | gnomAD %s %s | internal %s %s %s | classes %d | hcc %s | require_mac %s | gene %s",
    rv[["mcap_thr"]], rv[["gnomad_mode"]], rv[["gnomad_thr"]], rv[["int_col"]], rv[["int_cmp"]], rv[["int_thr"]],
    length(class_set), hcc_col, rv[["require_mac"]], rv[["gene_mode"]])
if (rv[["gene_mode"]] != "asis") stop("gene_mode ", rv[["gene_mode"]], " not implemented here")

hcc <- read_hcc_bed(HCC_BED)
w   <- fread(file.path(OUT, "C_weights_v2.tsv"))
n_ukb <- nrow(fread(file.path(W, "keep_analysis.txt")))

## ---------- per chromosome ----------
res <- lapply(CHRS, function(c) {
  D <- file.path(W, paste0("chr", c))
  if (!file.exists(file.path(D, "DONE"))) stop("chr", c, " not finished (no DONE) – run 10_ukb_qualify.sh")
  ma <- read_multianno(file.path(D, "B.hg19_multianno.txt.gz"))
  mp <- fread(file.path(D, "map.tsv"), colClasses = list(character = c("ID38", "chr", "ref", "alt")))
  mp[, key := make_key(chr, pos, ref, alt)]
  ma[, key := make_key(chr, pos, ref, alt)]
  ma <- merge(ma, mp[, .(key, ID38, ALT_CTS, OBS_CT)], by = "key", all.x = TRUE)
  n_nomap <- ma[is.na(ID38), .N]
  ma[, `:=`(mac = as.numeric(ALT_CTS), maf_obs = as.numeric(ALT_CTS) / as.numeric(OBS_CT))]
  ma[, `:=`(in_hcc = flag_in_hcc(chr, pos, hcc, "strict"), in_hcc_rvx = flag_in_hcc(chr, pos, hcc, "rvx"))]
  ma[, class := make_class(func, exfunc)]
  ma[, reason := qualify_reason(ma, mcap_thr = as.numeric(rv[["mcap_thr"]]), gnomad_thr = as.numeric(rv[["gnomad_thr"]]),
                                int_thr = as.numeric(rv[["int_thr"]]), gnomad_mode = rv[["gnomad_mode"]],
                                int_col = rv[["int_col"]], int_cmp = rv[["int_cmp"]], class_keep = class_set,
                                hcc_col = hcc_col, require_mac = as.logical(rv[["require_mac"]]))]
  q <- ma[reason == "qualifying", .(ID38, ALT = alt, Gene = gene, key_hg19 = key, ALT_CTS, OBS_CT,
                                     af_alt = maf_obs, af_nfe, mcap, class)]
  say("chr%-2d annotated %7d | no map %d | qualifying %6d variants, %7d alleles", c, nrow(ma), n_nomap, nrow(q), q[, sum(ALT_CTS)])
  list(funnel = ma[, .N, by = reason][, chr := c], q = q)
})
q      <- rbindlist(lapply(res, `[[`, "q"))
funnel <- rbindlist(lapply(res, `[[`, "funnel"))
if (anyDuplicated(q$ID38)) stop("duplicated ID38 in qualifying list")

## ---------- outputs ----------
q[, in_weights := Gene %chin% w$Gene]
fwrite(q, file.path(W, "C_ukb_qualifying_variants.tsv.gz"), sep = "\t")
ft <- dcast(funnel, reason ~ chr, value.var = "N", fill = 0)
ft[, total := rowSums(.SD), .SDcols = setdiff(names(ft), "reason")]
fwrite(ft, file.path(OUT, "C_ukb_funnel.tsv"), sep = "\t")
print(ft[, .(reason, total)])

g <- q[, .(n_var = .N, AC = sum(ALT_CTS)), by = Gene]
g <- merge(g, w[, .(Gene, OAC_munich = OAC, EXP_munich = EXP, beta_raw, set_allelefilt)], by = "Gene", all = TRUE)
g[is.na(n_var), `:=`(n_var = 0L, AC = 0)]
g[, `:=`(AC_per_1000 = 1000 * AC / n_ukb, in_weights = !is.na(beta_raw))]
setorder(g, -AC)
fwrite(g, file.path(OUT, "C_ukb_gene_counts.tsv"), sep = "\t")

## ---------- summary / sanity ----------
mun_oac <- fread(file.path(OUT, "C_D4_gene_oac_compare.tsv"))[, sum(OAC_bm)]
say("\n==== SUMMARY ====")
say("UKB analysis set n = %d | qualifying variants %d in %d genes | alleles %d", n_ukb, nrow(q), uniqueN(q$Gene), q[, sum(ALT_CTS)])
say("weighted genes (8,790) with >= 1 UKB qualifying variant: %d (%.1f%%) | UKB alleles in weighted genes: %.1f%%",
    g[in_weights & n_var > 0, .N], 100 * g[in_weights & n_var > 0, .N] / g[in_weights == TRUE, .N],
    100 * q[in_weights == TRUE, sum(ALT_CTS)] / q[, sum(ALT_CTS)])
say("qualifying alleles per person: UKB %.2f | Munich cases %.2f (all genes; Munich 53,489 variants)",
    q[, sum(ALT_CTS)] / n_ukb, mun_oac / N_MUN)
gg <- g[in_weights == TRUE & n_var > 0]
say("gene-level Spearman rho, UKB AC/n vs Munich gnomAD-expected EXP/1535 (weighted genes): %.3f",
    cor(gg$AC / n_ukb, gg$EXP_munich / N_MUN, method = "spearman"))
vl <- fread(cmd = paste("zcat", shQuote(F_VL)), header = FALSE, colClasses = "character")
vl[, key := parse_id_key(V1)$vkey]
for (gn in TOP2) {
  mk <- vl[V4 == gn, key]; uq <- q[Gene == gn]
  say("%-5s UKB: %d variants, %d alleles (expected carriers ~%.0f) | Munich variants %d, of these qualifying in UKB: %d",
      gn, nrow(uq), uq[, sum(ALT_CTS)], uq[, sum(ALT_CTS)], length(mk), sum(mk %chin% uq$key_hg19))
}
say("\nwrote %s", file.path(W, "C_ukb_qualifying_variants.tsv.gz"))
say("wrote %s", file.path(OUT, c("C_ukb_gene_counts.tsv", "C_ukb_funnel.tsv", "10b_ukb_qualify.log")))
sink()
