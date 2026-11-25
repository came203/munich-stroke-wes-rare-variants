#!/usr/bin/env bash
# =====================================================================
# 05_run_rvexcaliber_v2.sh  –  RV-EXCALIBER rerun on the rebuilt v2 set (VM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · RVGRS_Plan.md §7f (24.11.2026)
#
# CLASS=missense (primary, decision Rainer 24.11.2026): shell scripts of the /home/rainer/WES_analysis copy, which call
#                get_Varlist_rvexcaliber_missense_only.R (nonsynonymous SNV, M-CAP >= 0.025) in step 1A and step 2
# CLASS=default  (sensitivity): shell scripts of the /mnt/hgfs copy → get_Varlist_rvexcaliber.R (original classes:
#                missense M-CAP + splicing/stopgain/stoploss/frameshift)
# Both copies use the same R scripts, gnomAD files and MIGen ranking set (path_to_rvexcaliber = /home copy).
#
# Steps: 0 input = v2 PLINK set minus 1 sample per related pair (IDs from REMOVE_FILE; higher missingness) → N 1,535
#        1A internal burden matrix | 2 gnomAD burden matrix | 3 iCF/gCF summary associations (default MIGen ranking)
#        4 summary (λ, genes, Bonferroni/FDR, top genes, v1 top genes) + exact commands for the report
# Parameters (as main run): gnomAD MAF 0.001, internal MAF 0.01, M-CAP 0.025, nfe, hcc, fulladjust, ranking NA
# Usage:  nohup bash 05_run_rvexcaliber_v2.sh missense > 05_missense.out 2>&1 &
#         nohup bash 05_run_rvexcaliber_v2.sh default  > 05_default.out  2>&1 &   (after the first has finished)
# Work dir: $HOME/rvx_v2_<CLASS> (local disk; gzip on the VMware shared folder fails with "Value too large for defined data type")
# Output: copied to /mnt/hgfs/D/WES/05_rvexcaliber_output_v2_<CLASS>/  (results, matrices, variant lists, 05_summary_<CLASS>.txt)
# =====================================================================
set -euo pipefail
export LC_ALL=C
CLASS=${1:-missense}
WES=/mnt/hgfs/D/WES
V2=$WES/FOR_REPORT/v2/my_samples_plink_hg19_v2
INDIR=$WES/FOR_REPORT/v2/plink_in
DS=my_samples_plink_hg19_v2u
OUTDIR=${WORK:-$HOME/rvx_v2_${CLASS}}                 # RV-EXCALIBER works on local VM disk (gzip on /mnt/hgfs fails: EOVERFLOW)
FINAL=$WES/05_rvexcaliber_output_v2_${CLASS}          # results are copied here at the end
V1MAIN=$WES/05_rvexcaliber_output/rvexcaliber_testing_my_samples_plink_hg19_SummaryAssociations_0.001_0.025_nfe_hcc_iCFgCFadjust_rvexcaliber_base.txt
case "$CLASS" in
  missense) SHDIR=/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main/scripts/shell; VARLIST=get_Varlist_rvexcaliber_missense_only.R ;;
  default)  SHDIR=$WES/04_software/RV-EXCALIBER-main/scripts/shell;                     VARLIST=get_Varlist_rvexcaliber.R ;;
  *) echo "CLASS must be missense or default"; exit 1 ;;
esac
GMAF=0.001; IMAF=0.01; MCAP=0.025; ETH=nfe; COV=hcc; JOBS=${JOBS:-8}; ADJ=fulladjust
TAG=${GMAF}_${MCAP}_${ETH}_${COV}
REMOVE="$(tr '\n' ' ' < "${REMOVE_FILE:?set REMOVE_FILE (one IID per line, one sample per related pair)}")"
mkdir -p "$INDIR" "$OUTDIR" "$FINAL"
LOG=$OUTDIR/05_run_${CLASS}.log; exec > >(tee -a "$LOG") 2>&1
CMDS=$OUTDIR/05_commands_${CLASS}.txt; : > "$CMDS"
step(){ echo; echo "== $(date '+%F %T') $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ [ -s "$1" ] || die "expected output missing: $1"; echo "   ok: $(basename "$1") ($(du -h "$1" | cut -f1))"; }
run(){ echo "   \$ $*" | tee -a "$CMDS"; ( cd "$SHDIR" && "$@" ) || die "step failed (exit $?): $2 – see messages above and $OUTDIR"; }

step "0 checks + input (CLASS=$CLASS)"
for f in get_RVBurdenMatrix_internal_rvexcaliber.clean.sh get_RVBurdenMatrix_gnomAD_rvexcaliber.clean.sh get_SummaryAssociations_iCF_gCF_rvexcaliber.clean.sh; do
  [ -s "$SHDIR/$f" ] || die "missing $SHDIR/$f"
done
n_vl=$( { grep -c "Rscripts/$VARLIST " "$SHDIR/get_RVBurdenMatrix_internal_rvexcaliber.clean.sh" "$SHDIR/get_RVBurdenMatrix_gnomAD_rvexcaliber.clean.sh" || true; } | awk -F: '{s+=$NF} END{print s+0}')
[ "$n_vl" -ge 2 ] || die "$SHDIR scripts do not call $VARLIST in both step 1A and step 2"
echo "   shell dir: $SHDIR | variant list script: $VARLIST (calls found: $n_vl)"
echo "   shell-script differences between the two copies (provenance):"
diff "$WES/04_software/RV-EXCALIBER-main/scripts/shell/get_RVBurdenMatrix_internal_rvexcaliber.clean.sh" \
     "/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main/scripts/shell/get_RVBurdenMatrix_internal_rvexcaliber.clean.sh" | head -20 | sed 's/^/     /' || true
Rscript -e 'for(p in c("data.table","dplyr","ggplot2","R.utils")) if(!requireNamespace(p, quietly=TRUE)) stop("R package missing: ", p, " – install.packages(\"", p, "\")")' || die "R packages"
if [ ! -s "$INDIR/$DS.bed" ]; then
  printf '%s\n' $REMOVE | awk '{print $1, $1}' > "$INDIR/related_removed.txt"
  plink --bfile "$V2" --remove "$INDIR/related_removed.txt" --keep-allele-order --allow-no-sex --make-bed --out "$INDIR/$DS" > /dev/null
fi
echo "   input $INDIR/$DS: $(wc -l < "$INDIR/$DS.fam") samples, $(wc -l < "$INDIR/$DS.bim") variants (removed: $REMOVE)"
[ "$(wc -l < "$INDIR/$DS.fam")" -eq 1535 ] || die "expected 1535 samples"

step "1A internal burden matrix (ANNOVAR annotation, filters)"
if [ ! -s "$OUTDIR/${DS}_RVBurdenMatrix_${TAG}.txt.gz" ]; then
  run bash get_RVBurdenMatrix_internal_rvexcaliber.clean.sh "$INDIR" "$OUTDIR" "$DS" "$GMAF" "$IMAF" "$MCAP" "$ETH" "$COV" "$JOBS"
fi
need "$OUTDIR/${DS}_RVBurdenMatrix_${TAG}.txt.gz"

step "2 gnomAD burden matrix"
if [ ! -s "$OUTDIR/gnomAD_RVBurdenMatrix_${DS}_${TAG}.txt.gz" ]; then
  run bash get_RVBurdenMatrix_gnomAD_rvexcaliber.clean.sh "$OUTDIR" "$DS" "$GMAF" "$MCAP" "$ETH" "$COV"
fi
need "$OUTDIR/gnomAD_RVBurdenMatrix_${DS}_${TAG}.txt.gz"

step "3 summary associations (iCF/gCF, default MIGen ranking)"
RES=$OUTDIR/rvexcaliber_testing_${DS}_SummaryAssociations_${TAG}_iCFgCFadjust_rvexcaliber_base.txt
if [ ! -s "$RES" ]; then
  run bash get_SummaryAssociations_iCF_gCF_rvexcaliber.clean.sh "$OUTDIR" "$DS" NA "$GMAF" "$MCAP" "$ETH" "$COV" "$ADJ"
fi
need "$RES"
RESF=${RES/SummaryAssociations_/SummaryAssociations_allele_filter_}
LAMF=$(ls "$OUTDIR"/rvexcaliber_testing_${DS}_GenomicInflationFactor_allele_filter_${TAG}*iCFgCFadjust* 2>/dev/null | head -1 || true)

step "4 summary"
cat > "$OUTDIR/05_summary.R" <<'EOF'
suppressMessages(library(data.table))
a <- commandArgs(TRUE); res <- a[1]; resf <- a[2]; lamf <- a[3]; v1 <- a[4]; out <- a[5]; cls <- a[6]
rd <- function(f) { d <- fread(f); setnames(d, c("Gene","test_allele_count","test_gnomAD_allele_count_iCFgCFadjust","test_P_rvexcaliber_base_iCFgCFadjust"),
                                         c("Gene","OAC","EXP","P"), skip_absent = TRUE); d[, FDR := p.adjust(P, "BH")][, OR := OAC / EXP][] }
r <- rd(res); n <- nrow(r)
lam_med <- median(qchisq(r$P, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
s <- c(sprintf("RV-EXCALIBER v2 | class = %s | %s", cls, format(Sys.time())),
       sprintf("unfiltered: genes %d | lambda(median p) %.3f | p<0.05: %d (exp %.0f) | p>0.95: %d | Bonferroni (p<%.2e): %d | FDR<0.05: %d",
               n, lam_med, r[P < 0.05, .N], 0.05 * n, r[P > 0.95, .N], 0.05 / n, r[P < 0.05 / n, .N], r[FDR < 0.05, .N]))
if (file.exists(resf)) { f <- rd(resf); s <- c(s, sprintf("allele-filtered: genes %d | Bonferroni (p<%.2e): %d | FDR<0.05: %d",
                                                         nrow(f), 0.05 / nrow(f), f[P < 0.05 / nrow(f), .N], f[FDR < 0.05, .N])) }
if (nzchar(lamf) && file.exists(lamf)) s <- c(s, paste("GenomicInflationFactor file:", paste(readLines(lamf), collapse = " | ")))
top <- r[order(P)][1:20]
if (file.exists(resf)) top[, in_allele_filter_run := Gene %in% f$Gene]
fwrite(top, out, sep = "\t")
v1g <- c("TTK","ATXN3","ZNF266","ERCC2","C10orf82","DENND4B","PLG","PTPRN2")
cmp <- merge(data.table(Gene = v1g), r[, .(Gene, OAC_v2 = OAC, EXP_v2 = round(EXP, 2), P_v2 = signif(P, 3), FDR_v2 = signif(FDR, 3))], all.x = TRUE)
if (file.exists(v1)) cmp <- merge(cmp, rd(v1)[, .(Gene, OAC_v1 = OAC, EXP_v1 = round(EXP, 2), P_v1 = signif(P, 3))], by = "Gene", all.x = TRUE)
writeLines(s); cat("\ntop 20 genes (v2):\n"); print(top[, .(Gene, OAC, EXP = round(EXP, 2), OR = round(OR, 1), P = signif(P, 3), FDR = signif(FDR, 3))])
cat("\nv1 top genes in v2:\n"); print(cmp)
writeLines(c(s, "", capture.output(print(top)), "", capture.output(print(cmp))), sub("_top_genes_", "_summary_", sub("\\.tsv$", ".txt", out)))
EOF
Rscript "$OUTDIR/05_summary.R" "$RES" "$RESF" "${LAMF:-}" "$V1MAIN" "$OUTDIR/05_top_genes_${CLASS}.tsv" "$CLASS"
step "5 copy results to $FINAL"
cp -f "$OUTDIR"/rvexcaliber_testing_* "$OUTDIR"/rvexcaliber_gnomAD_* "$OUTDIR"/*RVBurdenMatrix_* "$OUTDIR"/*variant_list* \
      "$OUTDIR"/*pruned_annotation* "$OUTDIR"/*.hg19_multianno.txt "$OUTDIR"/05_* "$FINAL"/ 2>/dev/null || true
cp -f "$LOG" "$FINAL"/ 2>/dev/null || true
ls -la "$FINAL" | sed 's/^/   /'
echo; echo "== done $(date '+%F %T') | commands: $CMDS | summary: $FINAL/05_summary_${CLASS}.txt | work dir: $OUTDIR"
