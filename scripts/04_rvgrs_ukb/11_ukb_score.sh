#!/usr/bin/env bash
# =============================================================================
# 11_ukb_score.sh – Part C step 3: RVGRS scoring in UKB (RVGRS_Plan.md 4b, D1/D2/D5)
# Run on the VM after 10 all + 10b:   bash /mnt/hgfs/D/WES/11_ukb_score.sh
#
# 1 R: variant-level score file from C_ukb_qualifying_variants.tsv.gz x C_weights_v2.tsv
#      (scored allele = ALT; genes without a Munich weight get weight 0)
#      columns  w_raw    beta_raw, all 8,790 genes                       (D1 primary)
#               w_shrunk beta_shrunk, all genes                          (sens.)
#               w_nom    beta_raw, nominal genes (p<0.05 & OR>1)         (sens.)
#               w_af     beta_raw, allele-filtered genes (5,105)         (sens.)
#               n_w      1 per qualifying allele in the 8,790 genes      (D2 covariate for w_raw/w_shrunk)
#               n_all    1 per qualifying allele, all genes              (descriptive)
#               n_nom    1 per allele in nominal genes                   (covariate for w_nom)
#               n_af     1 per allele in allele-filtered genes           (covariate for w_af)
#               ERCC2, PLG  1 per allele in that gene                    (M6 carrier replication)
# 2 plink2 --score per chr, no-mean-imputation (missing = 0, D5), sums only
# 3 R: sum over chr -> ~/partC_ukb/C_ukb_scores.tsv.gz (individual level, local only)
#      + FOR_REPORT/partC/C_ukb_score_summary.tsv, 11_ukb_score.log
# =============================================================================
set -u
export LC_ALL=C
BASE=/mnt/hgfs/D/WES
OUT=$BASE/FOR_REPORT/partC
W=$HOME/partC_ukb
S=$W/score; mkdir -p "$S"
UKB=/cluster2/dichganslab/rmalik2/UKB_FINAL
KEEP=$W/keep_analysis.txt
PLINK2=${PLINK2:-$(command -v plink2)}
LOG=$OUT/11_ukb_score.log
exec > >(tee "$LOG") 2>&1
die(){ echo "ERROR: $*"; exit 1; }
echo "11_ukb_score.sh  $(date)"
for f in "$W/C_ukb_qualifying_variants.tsv.gz" "$OUT/C_weights_v2.tsv" "$KEEP"; do [ -s "$f" ] || die "missing $f"; done

# ---- 1 score file ----------------------------------------------------------
W="$W" OUT="$OUT" S="$S" Rscript - <<'EOF' || exit 1
suppressPackageStartupMessages(library(data.table))
q <- fread(cmd = paste("zcat", shQuote(file.path(Sys.getenv("W"), "C_ukb_qualifying_variants.tsv.gz"))))
w <- fread(file.path(Sys.getenv("OUT"), "C_weights_v2.tsv"))
x <- merge(q[, .(ID38, ALT, Gene)], w[, .(Gene, beta_raw, beta_shrunk, set_nominal, set_allelefilt)], by = "Gene", all.x = TRUE)
inw <- !is.na(x$beta_raw)
z <- function(v) fifelse(is.na(v), 0, v)
x[, `:=`(w_raw    = z(beta_raw),
         w_shrunk = z(beta_shrunk),
         w_nom    = fifelse(inw & set_nominal, beta_raw, 0),
         w_af     = fifelse(inw & set_allelefilt, beta_raw, 0),
         n_w      = as.numeric(inw),
         n_all    = 1,
         n_nom    = as.numeric(inw & set_nominal),
         n_af     = as.numeric(inw & set_allelefilt),
         ERCC2    = as.numeric(Gene == "ERCC2"),
         PLG      = as.numeric(Gene == "PLG"))]
cols <- c("ID38", "ALT", "w_raw", "w_shrunk", "w_nom", "w_af", "n_w", "n_all", "n_nom", "n_af", "ERCC2", "PLG")
x[, chr := sub(":.*", "", ID38)]
cat(sprintf("score file: %d variants | in weighted genes %d | nominal %d | allele-filtered %d | ERCC2 %d | PLG %d\n",
            nrow(x), sum(x$n_w), sum(x$n_nom), sum(x$n_af), sum(x$ERCC2), sum(x$PLG)))
for (c in unique(x$chr)) fwrite(x[chr == c, ..cols], file.path(Sys.getenv("S"), paste0("score_chr", c, ".txt")), sep = "\t")
fwrite(x[, ..cols][, chr := x$chr], file.path(Sys.getenv("S"), "score_all.txt.gz"), sep = "\t")
EOF

# ---- 2 plink2 --score per chr ----------------------------------------------
PF(){ [ -e "$UKB/chr$1.pvar" ] && echo "" || echo "vzs"; }
for c in $(seq 1 22); do
  sf=$S/score_chr$c.txt
  [ -s "$sf" ] || { echo "chr$c: no qualifying variants – skip"; continue; }
  t0=$(date +%s)
  "$PLINK2" --pfile "$UKB/chr$c" $(PF $c) --keep "$KEEP" \
            --score "$sf" 1 2 header-read no-mean-imputation cols=maybefid,nallele,scoresums \
            --score-col-nums 3-12 --threads 16 --memory 10000 --out "$S/chr$c" > "$S/chr$c.stdout" 2>&1 \
    || { tail -20 "$S/chr$c.log"; die "plink2 --score chr$c"; }
  nv=$(grep -m1 -oE '[0-9]+ variants processed' "$S/chr$c.log" || echo "? variants processed")
  echo "chr$c: $nv, $(( $(date +%s) - t0 )) s $(grep -i warning "$S/chr$c.log" | head -2 | tr '\n' ' ')"
done

# ---- 3 combine ---------------------------------------------------------------
W="$W" OUT="$OUT" S="$S" Rscript - <<'EOF' || exit 1
suppressPackageStartupMessages(library(data.table))
S <- Sys.getenv("S"); OUT <- Sys.getenv("OUT")
fs <- list.files(S, pattern = "^chr[0-9]+\\.sscore$", full.names = TRUE)
cat("sscore files:", length(fs), "\n")
tot <- NULL
for (f in fs) {
  d <- fread(f); setnames(d, sub("^#", "", names(d)))
  sc <- grep("_SUM$", names(d), value = TRUE)
  d <- d[, c("IID", sc), with = FALSE]
  if (is.null(tot)) tot <- d else {
    if (!identical(tot$IID, d$IID)) { d <- d[match(tot$IID, IID)]; if (anyNA(d$IID)) stop("IID mismatch in ", f) }
    for (s in sc) set(tot, j = s, value = tot[[s]] + d[[s]])
  }
}
setnames(tot, sub("_SUM$", "", names(tot)))
fwrite(tot, file.path(Sys.getenv("W"), "C_ukb_scores.tsv.gz"), sep = "\t")
sc <- setdiff(names(tot), "IID")
sm <- rbindlist(lapply(sc, function(s) data.table(score = s, mean = mean(tot[[s]]), sd = sd(tot[[s]]),
                                                  min = min(tot[[s]]), max = max(tot[[s]]),
                                                  pct_nonzero = 100 * mean(tot[[s]] != 0))))
cat(sprintf("individuals: %d\n", nrow(tot))); print(sm, digits = 4)
cat(sprintf("cor(w_raw, n_w) = %.3f | cor(w_raw, n_all) = %.3f | ERCC2 carriers %d | PLG carriers %d\n",
            cor(tot$w_raw, tot$n_w), cor(tot$w_raw, tot$n_all), sum(tot$ERCC2 > 0), sum(tot$PLG > 0)))
fwrite(sm, file.path(OUT, "C_ukb_score_summary.tsv"), sep = "\t")
cat("wrote", file.path(Sys.getenv("W"), "C_ukb_scores.tsv.gz"), "(individual level – local only) and C_ukb_score_summary.tsv\n")
EOF
echo "11 done $(date) – next: 12_ukb_stats.R"
