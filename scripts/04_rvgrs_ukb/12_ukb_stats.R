#!/usr/bin/env Rscript
## =============================================================================
## 12_ukb_stats.R – Part C step 4: RVGRS association with ischemic stroke in UKB (RVGRS_Plan.md 4b)
## Input : ~/partC_ukb/C_ukb_scores.tsv.gz (11), IS.pheno, IS.EUR (meta-PRS), incident.is.txt
## Models: M1 logistic IS ~ RVGRS_SD + total qualifying count (D2) + age + age^2 + sex + PC1-10
##         M2 M1 + meta-PRS_SD          M3 top 1% vs rest (M2 covariates)
##         M4 Cox, incident only (prevalent cases excluded; baseline age)      M5 dC (DeLong), NRI, IDI
##         M6 ERCC2 / PLG qualifying-variant carriers (M2 covariates without RVGRS; Firth if <10 carrier cases)
##         Sens.: no count covariate; beta_shrunk; nominal genes; allele-filtered genes; incident-only logistic; + chip
## Run   : Rscript /mnt/hgfs/D/WES/12_ukb_stats.R        (R: data.table, survival, pROC, ggplot2; logistf optional)
## Output: FOR_REPORT/partC/ C_results.tsv, C_m5_discrimination.tsv, C_m6_carriers.tsv, C_deciles.tsv,
##         C_figure_forest.{png,pdf}, C_figure_deciles.{png,pdf}, 12_ukb_stats.log (summary for the report)
## =============================================================================
suppressPackageStartupMessages({ library(data.table); library(survival); library(pROC); library(ggplot2) })
setDTthreads(8)
OUT <- "/mnt/hgfs/D/WES/FOR_REPORT/partC"; W <- path.expand("~/partC_ukb")
F_PHENO <- "/cluster2/dichganslab/rmalik2/REGENIE/FORCCR2/IS.pheno"
F_EUR   <- "/mnt/hgfs/D/data/UKB_IS_RISK/IS.EUR"
F_INC   <- "/mnt/hgfs/D/data/UKB_IS_RISK/incident.is.txt"
COL <- "#257de0"; INK <- "#191919"; INK2 <- "#666666"   # dataviz reference palette (series-1, text)
sink(file.path(OUT, "12_ukb_stats.log"), split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
hdr <- function(x) cat("\n================ ", x, " ================\n", sep = "")
say("12_ukb_stats.R  %s", format(Sys.time()))

## ---------- data ----------
sc  <- fread(cmd = paste("zcat", shQuote(file.path(W, "C_ukb_scores.tsv.gz"))), colClasses = list(character = "IID"))
ph  <- fread(F_PHENO, colClasses = list(character = c("FID", "IID")))
eur <- fread(F_EUR, colClasses = list(character = c("FID", "IID")))[, .(IID, PRS)]
inc <- fread(F_INC, colClasses = list(character = "eid"))[, .(IID = eid, event, days, blage)]
d <- merge(merge(sc, ph, by = "IID"), eur, by = "IID")
say("analysis set n = %d, cases = %d (expected 438,694 / 9,308)", nrow(d), d[ispheno == 1, .N])
z <- function(x) (x - mean(x)) / sd(x)
d[, `:=`(rv = z(w_raw), rv_shr = z(w_shrunk), rv_nom = z(w_nom), rv_af = z(w_af),
         cnt = z(n_w), cnt_nom = z(n_nom), cnt_af = z(n_af), prs = z(PRS),
         age = isage, age2 = (isage - mean(isage))^2 / 100, sexf = factor(sex), chipf = factor(chip),
         c_ercc2 = as.integer(ERCC2 > 0), c_plg = as.integer(PLG > 0), c_any = as.integer(ERCC2 + PLG > 0))]
inc_ev <- inc[event == 1, IID]
d[, prevalent := ispheno == 1 & !(IID %chin% inc_ev)]
PCS <- paste0("PC", 1:10)
COV <- c("age", "age2", "sexf", PCS)
say("cor(RVGRS raw score, total count n_w) = %.3f | cor(RVGRS, meta-PRS) = %.4f | prevalent cases (no incident event) = %d",
    cor(d$w_raw, d$n_w), cor(d$rv, d$prs), d[prevalent == TRUE, .N])
say("RVGRS (w_raw) mean cases %.3f vs non-cases %.3f", d[ispheno == 1, mean(w_raw)], d[ispheno == 0, mean(w_raw)])

## ---------- helpers ----------
res <- list()
lr <- function(model, term, rhs, data = d, note = "") {
  f <- as.formula(paste("ispheno ~", paste(rhs, collapse = " + ")))
  m <- glm(f, data = data, family = binomial)
  cf <- summary(m)$coefficients
  tn <- grep(paste0("^", term), rownames(cf), value = TRUE)[1]
  b <- cf[tn, 1]; se <- cf[tn, 2]
  r <- data.table(model = model, term = term, type = "OR", n = nobs(m), cases = sum(data$ispheno[!is.na(fitted(m))] == 1),
                  est = exp(b), lci = exp(b - 1.96 * se), uci = exp(b + 1.96 * se), p = cf[tn, 4], note = note)
  res[[length(res) + 1]] <<- r
  say("%-34s %s %.3f (%.3f-%.3f) p %.3g  [n %d, cases %d]", model, r$type, r$est, r$lci, r$uci, r$p, r$n, r$cases)
  invisible(m)
}

## =============================================================================
hdr("M1-M3 logistic (prevalent + incident IS)")
m1 <- lr("M1 RVGRS (primary)", "rv", c("rv", "cnt", COV))
m2 <- lr("M2 RVGRS + meta-PRS", "rv", c("rv", "cnt", "prs", COV))
lr("M2 meta-PRS (same model)", "prs", c("rv", "cnt", "prs", COV))
d[, top1 := as.integer(rv >= quantile(rv, 0.99))]
lr("M3 RVGRS top 1% vs rest", "top1", c("top1", "cnt", "prs", COV), note = sprintf("%d in top 1%%", d[, sum(top1)]))

hdr("Sensitivity")
lr("S1 no count covariate", "rv", c("rv", COV))
lr("S2 beta_shrunk", "rv_shr", c("rv_shr", "cnt", COV))
lr("S3 nominal genes (239)", "rv_nom", c("rv_nom", "cnt_nom", COV))
lr("S4 allele-filtered genes (5,105)", "rv_af", c("rv_af", "cnt_af", COV))
lr("S5 incident-only logistic", "rv", c("rv", "cnt", COV), data = d[prevalent == FALSE])
lr("S6 M2 + genotyping chip", "rv", c("rv", "cnt", "prs", "chipf", COV))

## =============================================================================
hdr("M4 Cox, incident IS (prevalent cases excluded)")
dc <- merge(d[prevalent == FALSE], inc[, .(IID, event, days, blage)], by = "IID")
dc <- dc[!is.na(days) & days > 0]
dc[, `:=`(ev = as.integer(event == 1), blage2 = (blage - mean(blage))^2 / 100)]
fc <- as.formula(paste("Surv(days, ev) ~ rv + cnt + prs + blage + blage2 + sexf +", paste(PCS, collapse = " + ")))
m4 <- coxph(fc, data = dc)
s4 <- summary(m4)$coefficients["rv", ]
r4 <- data.table(model = "M4 Cox incident (HR per SD)", term = "rv", type = "HR", n = m4$n, cases = m4$nevent,
                 est = exp(s4[["coef"]]), lci = exp(s4[["coef"]] - 1.96 * s4[["se(coef)"]]),
                 uci = exp(s4[["coef"]] + 1.96 * s4[["se(coef)"]]), p = s4[["Pr(>|z|)"]], note = "baseline age")
res[[length(res) + 1]] <- r4
say("%-34s HR %.3f (%.3f-%.3f) p %.3g  [n %d, events %d]", r4$model, r4$est, r4$lci, r4$uci, r4$p, r4$n, r4$cases)

## =============================================================================
hdr("M5 discrimination: base (age, age2, sex, PCs, count, meta-PRS) vs + RVGRS")
fb <- glm(as.formula(paste("ispheno ~", paste(c("cnt", "prs", COV), collapse = " + "))), data = d, family = binomial)
fn <- glm(as.formula(paste("ispheno ~", paste(c("rv", "cnt", "prs", COV), collapse = " + "))), data = d, family = binomial)
pb <- fitted(fb); pn <- fitted(fn); y <- d$ispheno
rb <- roc(y, pb, quiet = TRUE, direction = "<"); rn <- roc(y, pn, quiet = TRUE, direction = "<")
rt <- roc.test(rb, rn, method = "delong", paired = TRUE)
cib <- ci.auc(rb); cin <- ci.auc(rn)
dci <- if (!is.null(rt$conf.int)) as.numeric(rt$conf.int) else c(NA, NA)
ev <- y == 1; up <- pn > pb; dn <- pn < pb
nri_e <- mean(up[ev]) - mean(dn[ev]); nri_n <- mean(dn[!ev]) - mean(up[!ev])
se_e <- sqrt((mean(up[ev]) + mean(dn[ev]) - nri_e^2) / sum(ev)); se_n <- sqrt((mean(up[!ev]) + mean(dn[!ev]) - nri_n^2) / sum(!ev))
nri <- nri_e + nri_n; se_nri <- sqrt(se_e^2 + se_n^2)
idi <- (mean(pn[ev]) - mean(pb[ev])) - (mean(pn[!ev]) - mean(pb[!ev]))
se_idi <- sqrt(var(pn[ev] - pb[ev]) / sum(ev) + var(pn[!ev] - pb[!ev]) / sum(!ev))
m5 <- data.table(metric = c("C base", "C base + RVGRS", "delta C (new - base)", "continuous NRI", "NRI events", "NRI non-events", "IDI"),
                 est = c(auc(rb), auc(rn), auc(rn) - auc(rb), nri, nri_e, nri_n, idi),
                 lci = c(cib[1], cin[1], -dci[2], nri - 1.96 * se_nri, nri_e - 1.96 * se_e, nri_n - 1.96 * se_n, idi - 1.96 * se_idi),
                 uci = c(cib[3], cin[3], -dci[1], nri + 1.96 * se_nri, nri_e + 1.96 * se_e, nri_n + 1.96 * se_n, idi + 1.96 * se_idi),
                 p   = c(NA, NA, rt$p.value, 2 * pnorm(-abs(nri / se_nri)), NA, NA, 2 * pnorm(-abs(idi / se_idi))))
## pROC reports the CI of (AUC1 - AUC2) = base - new; sign flipped above to new - base
print(m5, digits = 4)
say("note: dC/NRI/IDI in-sample (1 added parameter); NRI continuous (any up/down movement), own code (Pencina 2011)")
fwrite(m5, file.path(OUT, "C_m5_discrimination.tsv"), sep = "\t")

## =============================================================================
hdr("M6 ERCC2 / PLG carrier replication (within-UKB case-control; independent of gnomAD)")
m6 <- rbindlist(lapply(c("c_ercc2", "c_plg", "c_any"), function(v) {
  tab <- d[, .(carriers = sum(get(v)), carrier_cases = sum(get(v) * ispheno))]
  f <- as.formula(paste("ispheno ~", paste(c(v, "prs", COV), collapse = " + ")))
  if (tab$carrier_cases < 10 && requireNamespace("logistf", quietly = TRUE)) {
    m <- logistf::logistf(f, data = d); i <- which(names(coef(m)) == v)
    out <- data.table(est = exp(coef(m)[i]), lci = exp(m$ci.lower[i]), uci = exp(m$ci.upper[i]), p = m$prob[i], method = "Firth")
  } else {
    m <- glm(f, data = d, family = binomial); cf <- summary(m)$coefficients[v, ]
    out <- data.table(est = exp(cf[1]), lci = exp(cf[1] - 1.96 * cf[2]), uci = exp(cf[1] + 1.96 * cf[2]), p = cf[4], method = "logistic")
  }
  cbind(data.table(carrier = sub("c_", "", v), tab,
                   freq_cases = tab$carrier_cases / d[ispheno == 1, .N], freq_controls = (tab$carriers - tab$carrier_cases) / d[ispheno == 0, .N]), out)
}))
print(m6, digits = 3)
say("discovery (Munich vs gnomAD): ERCC2 OR 7.2, PLG OR 13.7 (vs UKB as reference, 11b: 1.68 / 1.70)")
fwrite(m6, file.path(OUT, "C_m6_carriers.tsv"), sep = "\t")

## =============================================================================
hdr("Deciles of RVGRS (residualised on total count), M2 covariates, reference = decile 1")
d[, rv_res := resid(lm(w_raw ~ n_w, data = d))]
d[, dec := cut(rv_res, quantile(rv_res, 0:10 / 10), include.lowest = TRUE, labels = 1:10)]
md <- glm(as.formula(paste("ispheno ~ dec + prs +", paste(COV, collapse = " + "))), data = d, family = binomial)
cf <- summary(md)$coefficients
dd <- data.table(decile = 1:10, n = d[, .N, keyby = dec]$N, cases = d[, sum(ispheno), keyby = dec]$V1, est = 1, lci = 1, uci = 1, p = NA_real_)
for (k in 2:10) { r <- cf[paste0("dec", k), ]; dd[decile == k, `:=`(est = exp(r[1]), lci = exp(r[1] - 1.96 * r[2]), uci = exp(r[1] + 1.96 * r[2]), p = r[4])] }
print(dd, digits = 3)
fwrite(dd, file.path(OUT, "C_deciles.tsv"), sep = "\t")

## =============================================================================
rt_all <- rbindlist(res)
fwrite(rt_all, file.path(OUT, "C_results.tsv"), sep = "\t")

## ---------- figures ----------
th <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e6e6", linewidth = 0.3),
        axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK),
        plot.title = element_text(colour = INK, face = "bold", size = 12), plot.subtitle = element_text(colour = INK2))
fp <- rt_all[term %chin% c("rv", "rv_shr", "rv_nom", "rv_af")][!grepl("^M2 meta-PRS", model)]   # M4 (Cox) is already in rt_all (term rv)
fp[, lab := sprintf("%s  %.3f (%.3f-%.3f)", sub("^(M[0-9]|S[0-9]) ", "\\1  ", model), est, lci, uci)]
fp[, lab := factor(lab, levels = rev(lab))]
g1 <- ggplot(fp, aes(x = est, y = lab)) +
  geom_vline(xintercept = 1, colour = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_errorbar(aes(xmin = lci, xmax = uci), width = 0, orientation = "y", colour = COL, linewidth = 0.7) +
  geom_point(colour = COL, size = 2.6) +
  scale_x_log10() + labs(x = "OR (HR for M4) per SD of RVGRS, 95% CI", y = NULL,
                         title = "Rare-variant risk score and ischemic stroke in UK Biobank",
                         subtitle = sprintf("n = %s, %s cases; weights from Munich WES (1,535 cases) vs gnomAD",
                                            format(nrow(d), big.mark = ","), format(d[ispheno == 1, .N], big.mark = ","))) + th
ggsave(file.path(OUT, "C_figure_forest.png"), g1, width = 8, height = 4.2, dpi = 300, bg = "white")
ggsave(file.path(OUT, "C_figure_forest.pdf"), g1, width = 8, height = 4.2)
g2 <- ggplot(dd, aes(x = decile, y = est)) +
  geom_hline(yintercept = 1, colour = INK2, linewidth = 0.4, linetype = "dashed") +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, colour = COL, linewidth = 0.7) +
  geom_point(colour = COL, size = 2.6) +
  scale_x_continuous(breaks = 1:10) + scale_y_log10() +
  labs(x = "RVGRS decile (residualised on total qualifying allele count)", y = "OR vs decile 1 (95% CI)",
       title = "Ischemic stroke by RVGRS decile", subtitle = "adjusted for meta-PRS, age, age^2, sex, PC1-10") + th
ggsave(file.path(OUT, "C_figure_deciles.png"), g2, width = 6.5, height = 4, dpi = 300, bg = "white")
ggsave(file.path(OUT, "C_figure_deciles.pdf"), g2, width = 6.5, height = 4)

## ---------- report summary ----------
hdr("SUMMARY FOR THE REPORT (Aim 5)")
g <- function(m) rt_all[model == m]
say("UKB WES analysis set: n = %s, %s ischemic stroke cases (incident %d events in Cox set)",
    format(nrow(d), big.mark = ","), format(d[ispheno == 1, .N], big.mark = ","), r4$cases)
for (m in c("M1 RVGRS (primary)", "M2 RVGRS + meta-PRS", "M3 RVGRS top 1% vs rest")) say("%-28s OR %.3f (%.3f-%.3f), p %.2g", m, g(m)$est, g(m)$lci, g(m)$uci, g(m)$p)
say("%-28s HR %.3f (%.3f-%.3f), p %.2g", "M4 Cox incident", r4$est, r4$lci, r4$uci, r4$p)
say("M5 C-statistic %.4f -> %.4f (dC %.5f, DeLong p %.2g); continuous NRI %.3f (p %.2g); IDI %.2e (p %.2g)",
    m5$est[1], m5$est[2], m5$est[3], m5$p[3], m5$est[4], m5$p[4], m5$est[7], m5$p[7])
for (i in seq_len(nrow(m6))) say("M6 %-5s carriers %d (cases %d): OR %.2f (%.2f-%.2f), p %.2g [%s]",
                                 m6$carrier[i], m6$carriers[i], m6$carrier_cases[i], m6$est[i], m6$lci[i], m6$uci[i], m6$p[i], m6$method[i])
say("\nwrote C_results.tsv, C_m5_discrimination.tsv, C_m6_carriers.tsv, C_deciles.tsv, C_figure_forest/deciles .png/.pdf, 12_ukb_stats.log")
sink()
