#!/usr/bin/env Rscript
## =============================================================================
## 09_weights.R – Part C step 1 (RVGRS_Plan.md 4b; D1–D6 confirmed 25.11.2026)
##  (1) D1: gene weight table from the v2 missense main result (all 8,790 genes)
##      OR = OAC/EXP; beta_raw = ln(OR) (primary); beta_shrunk = ln((OAC+0.5)/(EXP+0.5));
##      set flags: all / nominal (p<0.05 & OR>1) / allele-filtered / FDR / ERCC2+PLG
##  (2) D4: re-implement the qualifying rule (partC_qualify_fun.R) on the Munich ANNOVAR
##      output and check it reproduces RV-EXCALIBER's qualifying variants (53,489) and the
##      per-gene OAC (burden matrix + main result) exactly. Tests 24 rule variants
##      (nonsyn yes/no x gnomAD mode 3 x internal-MAF denominator obs/2N x '<='/'<')
##      and 3 gene modes; writes the winning rule to C_filter_rule.tsv for 10_ukb_qualify.
## v3 (25.11.2026): grid = class set (pruned/nonsyn/any) x hcc mode x MAC>=1; gnomAD all-pop <, internal <
## v2 (25.11.2026): key fix (data.table key=), ordered ID keys, getBinary-type modes,
##      AF_int cross-check vs pruned annotation, prints the get_Varlist filter lines.
##
## Run on the VM after 09a (partC_qualify_fun.R must sit next to this script):
##   Rscript /mnt/hgfs/D/WES/09_weights.R 2>&1 | tee ~/partC_work/09_weights.stdout
##   optional: PLINK2=/path/to/plink2 Rscript ...
## Reads only. Writes FOR_REPORT/partC/: C_weights_v2.tsv (Zenodo), C_D4_combos.tsv,
##   C_D4_reason_funnel.tsv, C_D4_gene_oac_compare.tsv, C_D4_mismatch_variants.tsv (local),
##   C_filter_rule.tsv, 09_weights.log; scratch ~/partC_work/munich_v2u_freq.*
## =============================================================================
suppressPackageStartupMessages(library(data.table))
setDTthreads(8)

## ---------- config (RVGRS_Plan H3) ----------
BASE    <- "/mnt/hgfs/D/WES"
RES_DIR <- file.path(BASE, "05_rvexcaliber_output_v2_missense")
WORK    <- path.expand("~/rvx_v2_missense")
DS      <- "my_samples_plink_hg19_v2u"
TAG     <- "0.001_0.025_nfe_hcc"
V2U     <- file.path(BASE, "FOR_REPORT/v2/plink_in", DS)
HCC_BED <- "/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main/gnomAD_211_exomes_hg19_coverage/ALL_CHROM_gnomad_refGene_exons_highcov20X.mod.bed"
OUT     <- file.path(BASE, "FOR_REPORT/partC")
SCR     <- path.expand("~/partC_work")
TOP2    <- c("ERCC2", "PLG")
THR     <- list(mcap = 0.025, gnomad_maf = 0.001, int_maf = 0.01)
EXPECT  <- list(n_samples = 1535L, n_genes = 8790L, n_qual = 53489L, n_bm_genes = 11585L, lambda = 0.768)

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(SCR, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "09_weights.log"), split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
hdr <- function(x) cat("\n================ ", x, " ================\n", sep = "")
chk <- function(ok, what) say("%-6s %s", if (isTRUE(ok)) "OK" else "CHECK", what)
say("09_weights.R  %s", format(Sys.time()))

## ---------- helpers / shared rule ----------
a  <- commandArgs(FALSE); sf <- sub("^--file=", "", a[grep("^--file=", a)])
SCRIPT_DIR <- if (length(sf) == 1) dirname(normalizePath(sf)) else getwd()
FUN <- file.path(SCRIPT_DIR, "partC_qualify_fun.R")
if (!file.exists(FUN)) stop("partC_qualify_fun.R not found next to this script: ", FUN)
source(FUN)

rd <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = paste("zcat", shQuote(f)), ...) else fread(f, ...)
one_file <- function(dir, pattern) {
  x <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(x) > 1 && any(grepl("\\.gz$", x))) x <- x[grepl("\\.gz$", x)]
  if (length(x) != 1) stop(sprintf("expected 1 file matching '%s' in %s, found %d: %s",
                                   pattern, dir, length(x), paste(basename(x), collapse = ", ")))
  x
}
cand <- c(Sys.getenv("PLINK2"), Sys.which("plink2"),
          Sys.glob("/home/rainer/WES_analysis/04_software/plink2"),
          Sys.glob("/home/rainer/WES_analysis/04_software/*/plink2"))
cand <- cand[nzchar(cand) & file.exists(cand)]
if (!length(cand)) stop("plink2 not found – run with PLINK2=/path/to/plink2")
PLINK2 <- unname(cand[1]); say("plink2: %s", PLINK2)

## ---------- inputs ----------
f_main <- file.path(RES_DIR, sprintf("rvexcaliber_testing_%s_SummaryAssociations_%s_iCFgCFadjust_rvexcaliber_base.txt", DS, TAG))
f_af   <- one_file(RES_DIR, "SummaryAssociations_allele_filter_.*rvexcaliber_base\\.txt$")
f_vl   <- file.path(WORK, sprintf("%s_variant_list_extracted_%s.txt.gz", DS, TAG))
f_bm   <- file.path(WORK, sprintf("%s_RVBurdenMatrix_%s.txt.gz", DS, TAG))
f_ma   <- one_file(WORK, "hg19_multianno\\.txt(\\.gz)?$")
for (f in c(f_main, f_af, f_vl, f_bm, f_ma, paste0(V2U, c(".bed", ".bim", ".fam")), HCC_BED))
  if (!file.exists(f)) stop("missing input: ", f) else say("input  %s", f)

## =============================================================================
hdr("PART 1  Gene weights (D1)")
res  <- rd(f_main)
need <- c("Gene", "test_allele_count", "test_gnomAD_allele_count_iCFgCFadjust", "test_P_rvexcaliber_base_iCFgCFadjust")
if (length(miss <- setdiff(need, names(res)))) stop("main result lacks columns: ", paste(miss, collapse = ", "))
w <- res[, .(Gene, OAC = as.numeric(test_allele_count),
             EXP = as.numeric(test_gnomAD_allele_count_iCFgCFadjust),
             P   = as.numeric(test_P_rvexcaliber_base_iCFgCFadjust))]
if (anyDuplicated(w$Gene)) stop("duplicated genes in main result")
af <- rd(f_af)

w[, OR := OAC / EXP]
w[, beta_raw    := fifelse(OAC > 0 & EXP > 0, log(OAC / EXP), NA_real_)]
w[, beta_shrunk := log((OAC + 0.5) / (EXP + 0.5))]
w[, P_BH := p.adjust(P, "BH")]
w[, `:=`(set_all        = !is.na(beta_raw),
         set_nominal    = P < 0.05 & OR > 1,
         set_allelefilt = Gene %chin% af$Gene,
         set_fdr        = P_BH < 0.05,
         set_top2       = Gene %chin% TOP2)]
lam <- qchisq(median(w$P), 1, lower.tail = FALSE) / qchisq(0.5, 1)

chk(nrow(w) == EXPECT$n_genes, sprintf("genes in main result: %d (expected %d)", nrow(w), EXPECT$n_genes))
chk(abs(lam - EXPECT$lambda) < 0.005, sprintf("lambda(median p): %.3f (expected %.3f)", lam, EXPECT$lambda))
chk(w[OAC == 0, .N] == 0, sprintf("genes with OAC = 0 (beta_raw undefined, weight 0 in primary score): %d", w[OAC == 0, .N]))
say("genes EXP < 1: %d | EXP < 3.79 (allele filter): %d | allele-filtered set: %d (file %d rows)",
    w[EXP < 1, .N], w[EXP < 3.79, .N], w[set_allelefilt == TRUE, .N], nrow(af))
say("p<0.05: %d | nominal (p<0.05 & OR>1): %d | FDR<0.05: %d | Bonferroni: %d",
    w[P < 0.05, .N], w[set_nominal == TRUE, .N], w[set_fdr == TRUE, .N], w[P < 0.05 / .N, .N])
say("beta_raw: median %.3f, mean %.3f, IQR %.3f to %.3f, range %.2f to %.2f, share < 0: %.1f%%",
    median(w$beta_raw, na.rm = TRUE), mean(w$beta_raw, na.rm = TRUE),
    quantile(w$beta_raw, .25, na.rm = TRUE), quantile(w$beta_raw, .75, na.rm = TRUE),
    min(w$beta_raw, na.rm = TRUE), max(w$beta_raw, na.rm = TRUE), 100 * mean(w$beta_raw < 0, na.rm = TRUE))
say("beta_shrunk: median %.3f, range %.2f to %.2f", median(w$beta_shrunk), min(w$beta_shrunk), max(w$beta_shrunk))
say("(D2: global shift of beta is absorbed by the total-qualifying-count covariate in 12)")
print(w[Gene %chin% TOP2])

## =============================================================================
hdr("PART 2  D4 – reproduce RV-EXCALIBER qualifying variants on Munich")
N <- nrow(fread(paste0(V2U, ".fam"), header = FALSE))
chk(N == EXPECT$n_samples, sprintf("v2u samples: %d", N))

## internal allele counts (plink2)
fq_pre <- file.path(SCR, "munich_v2u_freq")
rc <- system(sprintf("%s --bfile %s --freq counts cols=+pos --threads 8 --memory 6000 --out %s > %s.stdout 2>&1",
                     shQuote(PLINK2), shQuote(V2U), shQuote(fq_pre), shQuote(fq_pre)))
if (rc != 0) stop("plink2 --freq failed, see ", fq_pre, ".log")
fq <- fread(paste0(fq_pre, ".acount"))
setnames(fq, sub("^#", "", names(fq)))
if (length(miss <- setdiff(c("CHROM", "POS", "REF", "ALT", "ALT_CTS", "OBS_CT"), names(fq))))
  stop("acount lacks: ", paste(miss, collapse = ", "))
fq <- fq[!grepl(",", ALT)]
fq[, `:=`(ALT_CTS = as.numeric(ALT_CTS), OBS_CT = as.numeric(OBS_CT))]
fq[, `:=`(mac   = pmin(ALT_CTS, OBS_CT - ALT_CTS),
          minor = fifelse(ALT_CTS <= OBS_CT - ALT_CTS, ALT, REF))]
fq[, `:=`(maf_obs = fifelse(OBS_CT > 0, mac / OBS_CT, NA_real_), maf_2N = mac / (2 * N))]
fq[, key := parse_id_key(ID)$vkey]                       # key = bim ID (chr:pos:ref:alt)
say("plink freq: %d variants | duplicated ID keys: %d | plink REF/ALT swapped vs ID: %d",
    nrow(fq), fq[duplicated(key), .N], fq[key != make_key(CHROM, POS, REF, ALT), .N])
fq <- unique(fq, by = "key")

## ANNOVAR annotation
ma <- read_multianno(f_ma)
say("multianno rows: %d | exonic: %d | nonsynonymous_SNV: %d | M-CAP non-missing: %d",
    nrow(ma), ma[grepl("exonic", func), .N], ma[exfunc %chin% NONSYN, .N], ma[!is.na(mcap), .N])
ma[, key := make_key(chr, pos, ref, alt)]
n_sw <- ma[!(key %chin% fq$key) & make_key(chr, pos, alt, ref) %chin% fq$key, .N]
ma[!(key %chin% fq$key) & make_key(chr, pos, alt, ref) %chin% fq$key, key := make_key(chr, pos, alt, ref)]
ma[, poskey := paste(chr, pos, sep = ":")]
say("multianno keys matched to bim IDs after allele swap: %d | duplicated keys (dropped): %d (SNVs: %d)",
    n_sw, ma[!is.na(key) & duplicated(key), .N], ma[!is.na(key) & duplicated(key) & nchar(ref) == 1 & nchar(alt) == 1, .N])
ma <- unique(ma, by = "key")
ma <- merge(ma, fq[, .(key, mac, minor, maf_obs, maf_2N)], by = "key", all.x = TRUE, sort = FALSE)
say("multianno SNVs without plink freq match: %d",
    ma[nchar(ref) == 1 & nchar(alt) == 1 & is.na(mac), .N])

## RV-EXCALIBER's own internal AF (pruned annotation, column AF_int) vs ours
f_pr <- file.path(WORK, sprintf("%s_pruned_annotation_for_R_input_hcc.txt.gz", DS))
if (!file.exists(f_pr)) stop("pruned annotation needed for the class set: ", f_pr)
{
  pr <- rd(f_pr, sep = "\t", colClasses = "character", na.strings = c(".", ""))
  pr[, key := make_key(Chr, Pos, Ref, Alt)]
  pr[!(key %chin% ma$key), key := make_key(Chr, Pos, Alt, Ref)]
  pc <- merge(pr[, .(key, AF_int = as.numeric(AF_int))], ma[, .(key, maf_obs, maf_2N)], by = "key")
  rel <- function(a, b) abs(a - b) / pmax(a, 1e-12)
  say("pruned annotation: %d rows, %d matched | AF_int ~ maf_2N (rel diff < 1e-3): %.1f%% | AF_int ~ maf_obs: %.1f%% | AF_int > 0.5: %d",
      nrow(pr), nrow(pc), 100 * pc[AF_int > 0, mean(rel(AF_int, maf_2N) < 1e-3, na.rm = TRUE)],
      100 * pc[AF_int > 0, mean(rel(AF_int, maf_obs) < 1e-3, na.rm = TRUE)], pc[AF_int > 0.5, .N])
}

## variant classes (Func|ExonicFunc) that RV-EXCALIBER keeps = classes present in the pruned annotation
pr[, class := make_class(Func.refGene, ExonicFunc.refGene)]
ma[, class := make_class(func, exfunc)]
ma[, in_pruned := key %chin% pr$key]
CLASS_PRUNED <- sort(unique(pr$class))
say("classes in pruned annotation (n = %d):", length(CLASS_PRUNED)); print(pr[, .N, by = class][order(-N)])
say("multianno class x in_pruned (top 40):")
print(dcast(ma[, .N, by = .(class, in_pruned)], class ~ in_pruned, value.var = "N", fill = 0)[order(-`TRUE`)][1:min(.N, 40)])

## filter lines of the RV-EXCALIBER variant rule actually used
f_gv <- file.path(RES_DIR, "get_Varlist_rvexcaliber_missense_fixed.R")
if (file.exists(f_gv)) {
  gv <- readLines(f_gv); hit <- grep("getBinary\\(|AF_int|internal_MAF|M.CAP|MCAP|nonsynonymous|exonic|eth", gv)
  hit <- hit[hit > 520]
  say("get_Varlist filter lines (line > 520, max 80):")
  for (i in head(hit, 80)) say("  %4d: %s", i, trimws(gv[i]))
}

hcc <- read_hcc_bed(HCC_BED)
ma[, in_hcc       := flag_in_bed(chr, pos, hcc)]
ma[, in_hcc_shift := flag_in_bed(chr, pos + 1L, hcc)]   # diagnostic: off-by-one BED handling
ma[, in_hcc_rvx   := flag_in_hcc(chr, pos, hcc, "rvx")]
say("in hcc 'rvx' mode (s <= pos <= e): %d", ma[in_hcc_rvx == TRUE, .N])
say("hcc intervals: %d | multianno variants in hcc: %d (with +1 shift: %d)",
    nrow(hcc), ma[in_hcc == TRUE, .N], ma[in_hcc_shift == TRUE, .N])

## RV-EXCALIBER qualifying list (no header: ID, counted allele, 1, Gene, gnomAD NFE AF, gene cumAF)
vl <- rd(f_vl, header = FALSE, colClasses = "character")
nm <- c("ID", "counted", "one", "Gene", "gnomAD_AF", "gene_cumAF")
k  <- min(6L, ncol(vl)); setnames(vl, seq_len(k), nm[seq_len(k)])
pk <- parse_id_key(vl$ID)
vl[, `:=`(key = pk$vkey, poskey = pk$poskey)]
match_col <- if (mean(is.na(vl$key)) < 0.001) "key" else "poskey"
ma[, mkey := get(match_col)]; vl[, mkey := get(match_col)]
th <- unique(na.omit(vl$mkey))
chk(length(th) == EXPECT$n_qual,
    sprintf("RV-EXCALIBER list: %d rows, %d unique variants (expected %d), %d genes; matching on %s",
            nrow(vl), length(th), EXPECT$n_qual, uniqueN(vl$Gene), match_col))

## burden matrix row sums = per-gene OAC (all genes with >= 1 allele)
bm <- rd(f_bm, header = TRUE)
oac_bm <- data.table(Gene = bm[[1]], OAC_bm = rowSums(as.matrix(bm[, -1, with = FALSE]), na.rm = TRUE))
chk(ncol(bm) - 1 == N, sprintf("burden matrix: %d genes x %d samples", nrow(bm), ncol(bm) - 1))
rm(bm); invisible(gc())

## --- rule variants ---
## v3: gnomAD all-pop '<' and internal '<' fixed (v2 grid: best; '<=' vs '<' and obs vs 2N identical)
CLASS_SETS <- list(pruned = CLASS_PRUNED,
                   nonsyn = grep("\\|nonsynonymous_SNV$", unique(ma$class), value = TRUE),
                   any    = NULL)
combos <- CJ(class_mode = c("pruned", "nonsyn", "any"), hcc_col = c("in_hcc_rvx", "in_hcc"),
             require_mac = c(FALSE, TRUE), sorted = FALSE)
combos[, `:=`(gnomad_mode = "allpop_af_lt", int_col = "maf_obs", int_cmp = "lt")]
run_rule <- function(cb) qualify_reason(ma, mcap_thr = THR$mcap, gnomad_thr = THR$gnomad_maf, int_thr = THR$int_maf,
                                        gnomad_mode = cb$gnomad_mode, int_col = cb$int_col, int_cmp = cb$int_cmp,
                                        class_keep = CLASS_SETS[[cb$class_mode]], hcc_col = cb$hcc_col,
                                        require_mac = cb$require_mac)
cres <- rbindlist(lapply(seq_len(nrow(combos)), function(i) {
  r <- run_rule(combos[i])
  ours <- unique(na.omit(ma$mkey[r == "qualifying"]))
  data.table(combos[i], n_ours = length(ours), n_rvx = length(th), n_both = length(intersect(ours, th)),
             only_ours = length(setdiff(ours, th)), only_rvx = length(setdiff(th, ours)))
}))
cres[, n_mismatch := only_ours + only_rvx]
print(cres[order(n_mismatch)])
best <- cres[order(n_mismatch)][1]                       # ties: first in grid order
say("chosen rule: class %s | hcc %s | require_mac %s | gnomAD %s | internal %s %s (mismatches %d)",
    best$class_mode, best$hcc_col, best$require_mac, best$gnomad_mode, best$int_col, best$int_cmp, best$n_mismatch)
## extras of the class-free rule, by class (shows what the pruning removes)
r_any <- run_rule(modifyList(as.list(best), list(class_mode = "any")))
say("variants qualifying without class filter but not in RV-EXCALIBER list, by class:")
print(ma[r_any == "qualifying" & !(mkey %chin% th), .N, by = class][order(-N)])

ma[, reason := run_rule(best)]
funnel <- ma[, .N, by = reason][order(reason)]
funnel[, remaining := nrow(ma) - cumsum(N) + N]
print(funnel)

## --- gene assignment mode ---
q <- ma[reason == "qualifying"]
gene_oac <- function(mode) {
  ag <- assign_genes(q$gene, mode)
  ag[, `:=`(mkey = q$mkey[i], mac = q$mac[i])]
  ag
}
gm <- rbindlist(lapply(c("split", "first", "asis"), function(m) {
  g <- gene_oac(m)[, .(OAC_ours = sum(mac)), by = Gene]
  cmp <- merge(g, oac_bm, by = "Gene", all = TRUE)
  cmp[is.na(OAC_ours), OAC_ours := 0]; cmp[is.na(OAC_bm), OAC_bm := 0]
  data.table(gene_mode = m, n_genes_ours = nrow(g), n_genes_bm = nrow(oac_bm),
             n_genes_diff = cmp[OAC_ours != OAC_bm, .N], sum_abs_diff = cmp[, sum(abs(OAC_ours - OAC_bm))])
}))
print(gm)
best_mode <- gm[order(n_genes_diff, sum_abs_diff)][1, gene_mode]
ag <- gene_oac(best_mode)
say("chosen gene mode: %s", best_mode)

## --- per-gene comparison (ours vs burden matrix vs main result) ---
gc <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE),
             list(ag[, .(OAC_ours = sum(mac), n_var_ours = uniqueN(mkey)), by = Gene],
                  oac_bm,
                  vl[, .(n_var_rvx = uniqueN(mkey)), by = Gene],
                  w[, .(Gene, OAC_main = OAC)]))
for (v in c("OAC_ours", "OAC_bm", "n_var_ours", "n_var_rvx")) set(gc, which(is.na(gc[[v]])), v, 0)
gc[, `:=`(match_bm = OAC_ours == OAC_bm, match_main = is.na(OAC_main) | OAC_ours == OAC_main)]

## --- further consistency checks ---
pairs_ours <- unique(paste(ag$mkey, ag$Gene)); pairs_rvx <- unique(paste(vl$mkey, vl$Gene))
cm <- merge(unique(vl[, .(mkey, counted, gnomAD_AF)], by = "mkey"),
            unique(ma[, .(mkey, minor, af_nfe)], by = "mkey"), by = "mkey")
cm[, af0 := fifelse(is.na(af_nfe), 0, af_nfe)]
n_cnt_mm <- cm[counted != minor, .N]
af_diff  <- cm[, max(abs(suppressWarnings(as.numeric(gnomAD_AF)) - af0), na.rm = TRUE)]
exf_rvx  <- ma[mkey %chin% th, .N, by = class][order(-N)]

hdr("D4 RESULT")
chk(best$n_mismatch == 0, sprintf("variant set: ours %d vs RV-EXCALIBER %d | only ours %d | only RV-EXCALIBER %d",
                                  best$n_ours, best$n_rvx, best$only_ours, best$only_rvx))
chk(length(setdiff(pairs_ours, pairs_rvx)) + length(setdiff(pairs_rvx, pairs_ours)) == 0,
    sprintf("variant-gene pairs: only ours %d | only RV-EXCALIBER %d",
            length(setdiff(pairs_ours, pairs_rvx)), length(setdiff(pairs_rvx, pairs_ours))))
chk(gc[match_bm == FALSE, .N] == 0, sprintf("per-gene OAC vs burden matrix: %d of %d genes differ", gc[match_bm == FALSE, .N], nrow(gc)))
chk(gc[match_main == FALSE, .N] == 0, sprintf("per-gene OAC vs main result (tested genes): %d of %d differ",
                                              gc[match_main == FALSE, .N], gc[!is.na(OAC_main), .N]))
chk(n_cnt_mm == 0, sprintf("counted allele == cohort minor allele: %d mismatches of %d", n_cnt_mm, nrow(cm)))
chk(af_diff < 1e-6, sprintf("gnomAD NFE AF (list vs multianno) max abs diff: %.2g", af_diff))
say("Func|ExonicFunc of RV-EXCALIBER qualifying variants:"); print(exf_rvx)
print(gc[Gene %chin% TOP2])

d4_pass <- best$n_mismatch == 0 && gc[match_bm == FALSE, .N] == 0 && gc[match_main == FALSE, .N] == 0
say("\nD4 VERDICT: %s", if (d4_pass) "PASS – rule reproduces RV-EXCALIBER exactly; use it for UKB (10)"
                       else "FAIL – inspect C_D4_mismatch_variants.tsv (reason column) before writing 10")

## =============================================================================
hdr("OUTPUTS")
w <- merge(w, gc[, .(Gene, OAC_ours, n_var_ours, n_var_rvx, oac_match = match_main)], by = "Gene", all.x = TRUE)
setcolorder(w, c("Gene", "OAC", "EXP", "OR", "P", "P_BH", "beta_raw", "beta_shrunk", "n_var_rvx", "n_var_ours",
                 "OAC_ours", "oac_match", "set_all", "set_nominal", "set_allelefilt", "set_fdr", "set_top2"))
setorder(w, P)
fwrite(w, file.path(OUT, "C_weights_v2.tsv"), sep = "\t")
fwrite(cres, file.path(OUT, "C_D4_combos.tsv"), sep = "\t")
fwrite(funnel, file.path(OUT, "C_D4_reason_funnel.tsv"), sep = "\t")
fwrite(gc[order(match_bm, -abs(OAC_ours - OAC_bm))], file.path(OUT, "C_D4_gene_oac_compare.tsv"), sep = "\t")

ours_best <- unique(na.omit(q$mkey))
only_o <- setdiff(ours_best, th); only_r <- setdiff(th, ours_best)
mcols <- c("mkey", "chr", "pos", "ref", "alt", "gene", "exfunc", POPS, "mcap", "mac", "maf_obs", "maf_2N",
           "class", "in_pruned", "in_hcc", "in_hcc_shift", "in_hcc_rvx", "reason")
mm <- rbind(ma[mkey %chin% only_o, c(list(side = "only_ours"), .SD), .SDcols = mcols],
            ma[mkey %chin% only_r, c(list(side = "only_rvx"), .SD), .SDcols = mcols],
            data.table(side = "only_rvx", mkey = setdiff(only_r, ma$mkey), reason = "not_in_multianno"),
            fill = TRUE)
fwrite(mm, file.path(OUT, "C_D4_mismatch_variants.tsv"), sep = "\t")
if (nrow(mm)) { say("mismatch reasons:"); print(mm[, .N, by = .(side, reason)]) }

rule <- data.table(
  param = c("mcap_thr", "gnomad_thr", "gnomad_mode", "int_thr", "int_col", "int_cmp", "class_mode", "hcc_mode", "require_mac",
            "gene_mode", "match_on", "hcc_bed", "annovar_dbs", "D4_verdict", "n_qualifying_munich", "date"),
  value = c(THR$mcap, THR$gnomad_maf, best$gnomad_mode, THR$int_maf, best$int_col, best$int_cmp,
            best$class_mode, sub("in_hcc_?", "", best$hcc_col), best$require_mac, best_mode, match_col, HCC_BED,
            "refGene, gnomad211_exome (AF_nfe/afr/sas/eas/amr), dbnsfp35c (M-CAP_score)",
            if (d4_pass) "PASS" else "FAIL", length(ours_best), format(Sys.Date())),
  note  = c("M-CAP_score > thr (strict)", "", "see partC_qualify_fun.R header", "internal AF of cohort",
            "maf_obs = MAC/OBS_CT; maf_2N = MAC/(2N)", "le: <= thr; lt: < thr",
            "class set -> C_class_set.txt", "rvx: s <= pos <= e; '' = strict", "FALSE: MAC-0 variants kept (no effect on OAC)", "asis = Gene.refGene string as is",
            "", "", "", "", "", ""))
fwrite(rule, file.path(OUT, "C_filter_rule.tsv"), sep = "\t")
if (!is.null(CLASS_SETS[[best$class_mode]])) writeLines(CLASS_SETS[[best$class_mode]], file.path(OUT, "C_class_set.txt"))
for (f in c("C_weights_v2.tsv", "C_D4_combos.tsv", "C_D4_reason_funnel.tsv", "C_D4_gene_oac_compare.tsv",
            "C_D4_mismatch_variants.tsv", "C_filter_rule.tsv", "C_class_set.txt", "09_weights.log")) say("wrote %s", file.path(OUT, f))
say("NOTE: C_D4_mismatch_variants.tsv is variant-level – keep local, not for Zenodo")
sink()
