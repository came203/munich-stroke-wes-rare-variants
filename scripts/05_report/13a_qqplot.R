#!/usr/bin/env Rscript
## =============================================================================
## 13a_qqplot.R – QQ plot of the main v2 RV-EXCALIBER run (all ischemic stroke, missense, unfiltered)
## RV-EXCALIBER itself wrote only the allele-filtered QQ plot. This script draws
##   (1) report figure: unfiltered run (8,790 genes), 95% CI band, lambda, Bonferroni line, top genes labelled
##   (2) companion: unfiltered vs allele-filtered overlay
## Run: Rscript /mnt/hgfs/D/WES/13a_qqplot.R
## Output: /mnt/hgfs/D/WES/FOR_REPORT/figures/QQ_v2_main.{png,pdf}, QQ_v2_main_vs_allelefilter.{png,pdf}
## Report style: no en/em dashes in figure text
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
RES <- "/mnt/hgfs/D/WES/05_rvexcaliber_output_v2_missense"
OUT <- "/mnt/hgfs/D/WES/FOR_REPORT/figures"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
f_main <- file.path(RES, "rvexcaliber_testing_my_samples_plink_hg19_v2u_SummaryAssociations_0.001_0.025_nfe_hcc_iCFgCFadjust_rvexcaliber_base.txt")
f_af   <- sub("SummaryAssociations_", "SummaryAssociations_allele_filter_", f_main, fixed = TRUE)
COL1 <- "#257de0"; COL2 <- "#f47230"; INK <- "#191919"; INK2 <- "#666666"; BAND <- "#e6e6e6"   # dataviz reference palette
N_LABEL <- 2

rd <- function(f, set) {
  x <- fread(f)
  data.table(Gene = x$Gene, P = as.numeric(x$test_P_rvexcaliber_base_iCFgCFadjust), set = set)
}
qq <- function(x) {
  x <- x[order(P)]; n <- nrow(x); i <- seq_len(n)
  x[, `:=`(exp = -log10(i / (n + 1)), obs = -log10(P),
           lo = -log10(qbeta(0.975, i, n - i + 1)), hi = -log10(qbeta(0.025, i, n - i + 1)))]
  x[]
}
lam <- function(p) qchisq(median(p), 1, lower.tail = FALSE) / qchisq(0.5, 1)
m  <- qq(rd(f_main, "All genes (unfiltered)"))
af <- qq(rd(f_af, "Allele-filtered"))
n_m <- nrow(m); bonf <- -log10(0.05 / n_m)
cat(sprintf("unfiltered: %d genes, lambda %.3f | allele-filtered: %d genes, lambda %.3f\n", n_m, lam(m$P), nrow(af), lam(af$P)))
print(m[1:5, .(Gene, P)])

th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#eeeeee", linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 12), plot.subtitle = element_text(colour = INK2),
        legend.position = "top", legend.title = element_blank(), legend.text = element_text(colour = INK))
lim <- max(c(m$obs, m$exp, bonf)) * 1.05
save2 <- function(g, name, w = 5.5, h = 5.5) {
  ggsave(file.path(OUT, paste0(name, ".png")), g, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".pdf")), g, width = w, height = h)
  cat("wrote", file.path(OUT, paste0(name, c(".png", ".pdf"))), "\n")
}

## (1) report figure
g1 <- ggplot(m, aes(exp, obs)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = BAND) +
  geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.4) +
  geom_hline(yintercept = bonf, colour = INK2, linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = 0.05, y = bonf, label = sprintf("Bonferroni (0.05 / %s genes)", format(n_m, big.mark = ",")),
           hjust = 0, vjust = -0.5, size = 3, colour = INK2) +
  geom_point(colour = COL1, size = 1.3) +
  geom_text(data = m[1:N_LABEL], aes(label = Gene), hjust = 1.2, vjust = 0.4, size = 3.4, fontface = "italic", colour = INK) +
  annotate("text", x = 0.05, y = lim * 0.97, hjust = 0, size = 3.6, colour = INK, parse = TRUE,
           label = sprintf("lambda == %.2f", lam(m$P))) +
  coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) +
  labs(x = expression(Expected ~ -log[10](italic(P))), y = expression(Observed ~ -log[10](italic(P))),
       title = "Gene-based rare missense variant burden",
       subtitle = sprintf("1,535 ischemic stroke cases vs gnomAD NFE (RV-EXCALIBER), %s genes", format(n_m, big.mark = ","))) + th
save2(g1, "QQ_v2_main")

## (2) unfiltered vs allele-filtered
both <- rbind(m, af)
both[, set := factor(set, levels = c("All genes (unfiltered)", "Allele-filtered"))]
g2 <- ggplot(both, aes(exp, obs, colour = set)) +
  geom_abline(slope = 1, intercept = 0, colour = INK2, linewidth = 0.4) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c(COL1, COL2),
                      labels = c(sprintf("All genes (n = %s, lambda %.2f)", format(n_m, big.mark = ","), lam(m$P)),
                                 sprintf("Allele-filtered (n = %s, lambda %.2f)", format(nrow(af), big.mark = ","), lam(af$P)))) +
  coord_cartesian(xlim = c(0, lim), ylim = c(0, lim)) +
  labs(x = expression(Expected ~ -log[10](italic(P))), y = expression(Observed ~ -log[10](italic(P))),
       title = "Gene-based burden test, unfiltered vs allele-filtered") + th
save2(g2, "QQ_v2_main_vs_allelefilter")
