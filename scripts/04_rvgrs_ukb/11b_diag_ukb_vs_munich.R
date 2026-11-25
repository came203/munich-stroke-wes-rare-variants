#!/usr/bin/env Rscript
## =============================================================================
## 11b_diag_ukb_vs_munich.R – READ-ONLY diagnostic before 12 (25.11.2026)
## Question: is the gnomAD-based expectation too low for ERCC2/PLG (UKB population carries
##   3.8x / 7.2x the gnomAD-calibrated rate)? -> technical explanation of the discovery signal?
## A gene level (8,790 weighted genes): UKB rate / gnomAD-expected rate (r_ukb_exp); where do ERCC2/PLG
##   rank; does the Munich OR (vs gnomAD) track r_ukb_exp across genes; Munich vs UKB as external
##   control (calibrated by the global Munich/UKB allele ratio c), Poisson test + lambda
## B variant level ERCC2/PLG: Munich case AC vs UKB AC per variant, gnomAD presence (all 5 pops missing
##   = absent), share of alleles from gnomAD-absent variants; O/E vs gnomAD on gnomAD-present variants only
## C genome-wide share of qualifying alleles from gnomAD-absent variants (Munich vs UKB)
## Run: Rscript /mnt/hgfs/D/WES/11b_diag_ukb_vs_munich.R
## Writes FOR_REPORT/partC/: 11b_diag.log, C_11b_gene_ukb_vs_munich.tsv, C_11b_topgenes.tsv,
##   C_11b_ERCC2_PLG_variants.tsv (variant-level allele counts, no individual data – keep local)
## Note: UKB (IDT capture, UKB pipeline) vs Munich (Twist, DeepVariant/GLnexus) differ technically –
##   UKB-as-control comparisons are descriptive, not a formal replication
## =============================================================================
suppressPackageStartupMessages(library(data.table))
setDTthreads(8)
BASE <- "/mnt/hgfs/D/WES"; OUT <- file.path(BASE, "FOR_REPORT/partC"); W <- path.expand("~/partC_ukb")
WORK <- path.expand("~/rvx_v2_missense"); DS <- "my_samples_plink_hg19_v2u"; TAG <- "0.001_0.025_nfe_hcc"
N_MUN <- 1535L; TOP2 <- c("ERCC2", "PLG")
a <- commandArgs(FALSE); sf <- sub("^--file=", "", a[grep("^--file=", a)])
source(file.path(if (length(sf) == 1) dirname(normalizePath(sf)) else getwd(), "partC_qualify_fun.R"))
sink(file.path(OUT, "11b_diag.log"), split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
hdr <- function(x) cat("\n================ ", x, " ================\n", sep = "")
say("11b_diag_ukb_vs_munich.R  %s", format(Sys.time()))
n_ukb <- nrow(fread(file.path(W, "keep_analysis.txt")))
lam <- function(p) { p <- p[!is.na(p)]; qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1) }

## =============================================================================
hdr("A  gene level: UKB population vs gnomAD expectation vs Munich cases")
g <- fread(file.path(OUT, "C_ukb_gene_counts.tsv"))[in_weights == TRUE]
w <- fread(file.path(OUT, "C_weights_v2.tsv"))[, .(Gene, OR, P, P_BH)]
g <- merge(g, w, by = "Gene")
g[, `:=`(rate_mun = OAC_munich / N_MUN, rate_ukb = AC / n_ukb, rate_exp = EXP_munich / N_MUN)]
c_glob <- g[, sum(rate_mun) / sum(rate_ukb)]
say("weighted genes %d (UKB AC > 0: %d) | qualifying alleles/person: Munich %.2f, UKB %.2f, gnomAD-expected %.2f",
    nrow(g), g[AC > 0, .N], g[, sum(rate_mun)], g[, sum(rate_ukb)], g[, sum(rate_exp)])
say("global ratios: Munich/UKB c = %.3f | Munich/gnomAD-exp = %.3f | UKB/gnomAD-exp = %.3f",
    c_glob, g[, sum(rate_mun) / sum(rate_exp)], g[, sum(rate_ukb) / sum(rate_exp)])

g <- g[AC > 0]
g[, `:=`(r_ukb_exp = rate_ukb / rate_exp, exp_from_ukb = rate_ukb * N_MUN * c_glob)]
g[, or_vs_ukb := OAC_munich / exp_from_ukb]
g[, `:=`(p_vs_ukb_2s = mapply(function(o, e) poisson.test(o, T = e)$p.value, OAC_munich, exp_from_ukb),
         p_vs_ukb_gt = mapply(function(o, e) poisson.test(o, T = e, alternative = "greater")$p.value, OAC_munich, exp_from_ukb))]
q_r <- quantile(g$r_ukb_exp, c(.05, .25, .5, .75, .95))
say("r_ukb_exp (UKB rate / gnomAD-expected rate) quantiles 5/25/50/75/95%%: %s", paste(sprintf("%.2f", q_r), collapse = " / "))
for (gn in TOP2) say("  %-5s r_ukb_exp %.2f = percentile %.1f%% | OR vs gnomAD %.1f -> OR vs UKB (calibrated) %.2f, p(>) %.2g",
                      gn, g[Gene == gn, r_ukb_exp], 100 * mean(g$r_ukb_exp <= g[Gene == gn, r_ukb_exp]),
                      g[Gene == gn, OR], g[Gene == gn, or_vs_ukb], g[Gene == gn, p_vs_ukb_gt])
gs <- g[EXP_munich >= 3]
say("genes with EXP >= 3 (n %d): Spearman log OR(gnomAD) vs log r_ukb_exp = %.3f | vs log OR(UKB) = %.3f",
    nrow(gs), cor(log(gs$OR + 1e-3), log(gs$r_ukb_exp), method = "spearman"),
    cor(log(gs$OR + 1e-3), log(gs$or_vs_ukb + 1e-3), method = "spearman"))
say("Munich vs UKB-as-control (two-sided Poisson, genes with expected >= 1, n %d): lambda %.3f | p<0.05 %d | Bonferroni %d | FDR<0.05 %d",
    g[exp_from_ukb >= 1, .N], lam(g[exp_from_ukb >= 1, p_vs_ukb_2s]), g[exp_from_ukb >= 1 & p_vs_ukb_2s < 0.05, .N],
    g[exp_from_ukb >= 1 & p_vs_ukb_gt < 0.05 / .N, .N], g[exp_from_ukb >= 1, sum(p.adjust(p_vs_ukb_gt, "BH") < 0.05)])
say("same test vs gnomAD (RV-EXCALIBER P), same genes: lambda %.3f", lam(g[exp_from_ukb >= 1, P]))

top <- unique(rbind(g[order(P)][1:20], g[Gene %chin% TOP2]))
top <- top[order(P), .(Gene, OAC_munich, EXP_munich, OR_gnomad = OR, P_gnomad = P, UKB_AC = AC,
                       r_ukb_exp, exp_from_ukb, OR_vs_ukb = or_vs_ukb, p_vs_ukb_gt, p_vs_ukb_2s)]
say("\ntop 20 genes (gnomAD-based P) – Munich vs gnomAD and vs UKB:"); print(top, digits = 3)
top2 <- g[order(p_vs_ukb_gt)][1:15, .(Gene, OAC_munich, exp_from_ukb, OR_vs_ukb = or_vs_ukb, p_vs_ukb_gt, OR_gnomad = OR, P_gnomad = P)]
say("\ntop 15 genes by Munich vs UKB-as-control:"); print(top2, digits = 3)
fwrite(g, file.path(OUT, "C_11b_gene_ukb_vs_munich.tsv"), sep = "\t")
fwrite(top, file.path(OUT, "C_11b_topgenes.tsv"), sep = "\t")

## =============================================================================
hdr("B  variant level ERCC2 / PLG")
q <- fread(cmd = paste("zcat", shQuote(file.path(W, "C_ukb_qualifying_variants.tsv.gz"))))
q[, absent := is.na(af_nfe)]                       # all-pop rule: qualifying with NFE missing => all 5 pops missing
vl <- fread(cmd = paste("zcat", shQuote(file.path(WORK, sprintf("%s_variant_list_extracted_%s.txt.gz", DS, TAG)))),
            header = FALSE, colClasses = "character")
vl[, key := parse_id_key(V1)$vkey]
fq <- fread(path.expand("~/partC_work/munich_v2u_freq.acount")); setnames(fq, sub("^#", "", names(fq)))
fq[, `:=`(key = parse_id_key(ID)$vkey, mac_mun = pmin(ALT_CTS, OBS_CT - ALT_CTS))]
ma <- read_multianno(file.path(WORK, sprintf("%s.hg19_multianno.txt.gz", DS)))
ma[, key := make_key(chr, pos, ref, alt)]
ma[, absent_mun := Reduce(`&`, lapply(POPS, function(p) is.na(ma[[p]])))]
mun <- merge(vl[, .(key, Gene = V4)], fq[, .(key, mac_mun)], by = "key", all.x = TRUE)
mun <- merge(mun, ma[, .(key, af_nfe_mun = af_nfe, mcap_mun = mcap, absent_mun, class_mun = make_class(func, exfunc))], by = "key", all.x = TRUE)

vt <- rbindlist(lapply(TOP2, function(gn) {
  m <- mun[Gene == gn]; u <- q[Gene == gn]
  x <- merge(m[, .(key, AC_munich = mac_mun, af_nfe_mun, absent_mun, mcap_mun, class_mun)],
             u[, .(key = key_hg19, AC_ukb = ALT_CTS, AF_ukb = af_alt, af_nfe, absent, mcap, class)], by = "key", all = TRUE)
  x[, `:=`(Gene = gn, in_munich = !is.na(AC_munich), in_ukb = !is.na(AC_ukb))]
  x[is.na(AC_munich), AC_munich := 0]; x[is.na(AC_ukb), AC_ukb := 0]
  x[, gnomad_absent := fifelse(in_ukb, absent, absent_mun)]
  x[, gnomad_nfe := fifelse(in_ukb, af_nfe, af_nfe_mun)]
  x[]
}))
for (gn in TOP2) {
  x <- vt[Gene == gn]
  say("\n%s: variants Munich %d / UKB %d / both %d | Munich alleles %d (from gnomAD-absent variants: %d = %.0f%%) | UKB alleles %d (absent: %.0f%%)",
      gn, x[in_munich == TRUE, .N], x[in_ukb == TRUE, .N], x[in_munich & in_ukb, .N], x[, sum(AC_munich)],
      x[gnomad_absent == TRUE, sum(AC_munich)], 100 * x[gnomad_absent == TRUE, sum(AC_munich)] / x[, sum(AC_munich)],
      x[, sum(AC_ukb)], 100 * x[gnomad_absent == TRUE, sum(AC_ukb)] / x[, sum(AC_ukb)])
  e <- g[Gene == gn, EXP_munich]; op <- x[gnomad_absent == FALSE, sum(AC_munich)]
  say("  O/E vs gnomAD: all variants %d / %.2f = %.1f | gnomAD-present variants only %d / %.2f = %.1f (EXP counts gnomAD-present alleles only)",
      x[, sum(AC_munich)], e, x[, sum(AC_munich)] / e, op, e, op / e)
  say("  Munich variants sorted by Munich AC (with UKB AC / AF, gnomAD NFE AF; NA = absent from gnomAD):")
  print(x[in_munich == TRUE][order(-AC_munich), .(key, AC_munich, AC_ukb, AF_ukb = signif(AF_ukb, 2), gnomad_nfe, gnomad_absent, mcap = fifelse(is.na(mcap), mcap_mun, mcap))], digits = 3)
  say("  top 10 UKB variants by UKB AC:")
  print(x[in_ukb == TRUE][order(-AC_ukb)][1:min(.N, 10), .(key, AC_ukb, AF_ukb = signif(AF_ukb, 2), AC_munich, gnomad_nfe, gnomad_absent, mcap, class)], digits = 3)
}
fwrite(vt, file.path(OUT, "C_11b_ERCC2_PLG_variants.tsv"), sep = "\t")

## =============================================================================
hdr("C  genome-wide share of qualifying alleles from gnomAD-absent variants")
say("Munich cases: %.1f%% of %d alleles | UKB: %.1f%% of %d alleles",
    100 * mun[absent_mun == TRUE, sum(mac_mun, na.rm = TRUE)] / mun[, sum(mac_mun, na.rm = TRUE)], mun[, sum(mac_mun, na.rm = TRUE)],
    100 * q[absent == TRUE, sum(ALT_CTS)] / q[, sum(ALT_CTS)], q[, sum(ALT_CTS)])
ga <- mun[, .(absent_share_mun = sum(mac_mun[absent_mun %in% TRUE], na.rm = TRUE) / sum(mac_mun, na.rm = TRUE)), by = Gene]
gs2 <- merge(gs, ga, by = "Gene")
say("genes with EXP >= 3: Spearman log OR(gnomAD) vs Munich absent-allele share = %.3f", cor(log(gs2$OR + 1e-3), gs2$absent_share_mun, method = "spearman"))
say("\nwrote C_11b_gene_ukb_vs_munich.tsv, C_11b_topgenes.tsv, C_11b_ERCC2_PLG_variants.tsv (variant-level, local), 11b_diag.log")
sink()
