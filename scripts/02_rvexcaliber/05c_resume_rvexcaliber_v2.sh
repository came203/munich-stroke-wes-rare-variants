#!/usr/bin/env bash
# =====================================================================
# 05c_resume_rvexcaliber_v2.sh  –  resume the RV-EXCALIBER v2 missense run WITHOUT re-running ANNOVAR (VM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · after 05b (24.11.2026)
#
# Cause of the failed run: get_Varlist_rvexcaliber_missense_only.R filters ExonicFunc == "nonsynonymous SNV" (space),
#   RV-EXCALIBER tables (internal + gnomAD) use "nonsynonymous_SNV" → 0 variants on both sides.
# v2 of 05c: missense_only.R also changed the variant-list format (4 instead of 6 columns) → step-2 crash; the fixed
#   script is now derived from the ORIGINAL get_Varlist (correct format) with the LoF-or-M-CAP clause reduced to M-CAP.
# Fix (original scripts untouched):
#   R     : new get_Varlist_rvexcaliber_missense_fixed.R = missense_only.R + at the start of its getVarFilter():
#           ExonicFunc '_' → ' ', gnomAD AF '.'/NA → 0 (absent from gnomAD qualifies, as in v1), AF_int numeric,
#           filter counts printed; M-CAP '>= MCAP' → '> MCAP' (as original RV-EXCALIBER / plan: M-CAP > 0.025)
#   shell : patched copies in $WORK/shell_fixed: step 1A + step 2 call the fixed R script; in step 1A the
#           annovar_annotation call is skipped when ${dataset}.hg19_multianno.txt.gz exists (reused from the 05 run).
#           All other step-1A functions run again (they rebuild the ID-renamed preprocessed PLINK set that step 1A
#           deletes at its end and that the burden-matrix generation needs).
# Then: step 2 (gnomAD) and step 3 (summary associations) as in 05, summary, copy to /mnt/hgfs.
# Usage : nohup bash 05c_resume_rvexcaliber_v2.sh > 05c_missense.out 2>&1 &
# =====================================================================
set -euo pipefail
export LC_ALL=C
WES=/mnt/hgfs/D/WES
RVX=/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main
SHDIR=$RVX/scripts/shell; RDIR=$RVX/scripts/Rscripts
INDIR=$WES/FOR_REPORT/v2/plink_in
DS=my_samples_plink_hg19_v2u
WORK=${WORK:-$HOME/rvx_v2_missense}
FINAL=$WES/05_rvexcaliber_output_v2_missense
V1MAIN=$WES/05_rvexcaliber_output/rvexcaliber_testing_my_samples_plink_hg19_SummaryAssociations_0.001_0.025_nfe_hcc_iCFgCFadjust_rvexcaliber_base.txt
GMAF=0.001; IMAF=0.01; MCAP=0.025; ETH=nfe; COV=hcc; JOBS=${JOBS:-8}; ADJ=fulladjust
TAG=${GMAF}_${MCAP}_${ETH}_${COV}
RSRC=$RDIR/get_Varlist_rvexcaliber_original.R          # = get_Varlist_rvexcaliber.R of the v1 main run (md5 identical)
RFIX=$RDIR/get_Varlist_rvexcaliber_missense_fixed.R
SHFIX=$WORK/shell_fixed
mkdir -p "$FINAL" "$SHFIX"
LOG=$WORK/05c_resume.log; exec > >(tee -a "$LOG") 2>&1
CMDS=$WORK/05_commands_missense.txt
step(){ echo; echo "== $(date '+%F %T') $*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ [ -s "$1" ] || die "expected output missing/empty: $1"; echo "   ok: $(basename "$1") ($(du -h "$1" | cut -f1))"; }
run(){ echo "   \$ $*" | tee -a "$CMDS"; ( cd "$SHDIR" && "$@" ) || die "step failed (exit $?): $2"; }
zrows(){ { zcat -f "$1" 2>/dev/null || true; } | wc -l; }

step "0 checks"
need "$WORK/$DS.hg19_multianno.txt.gz"
need "$INDIR/$DS.bed"
[ "$(wc -l < "$INDIR/$DS.fam")" -eq 1535 ] || die "expected 1535 samples in $INDIR/$DS.fam"
Rscript -e 'for(p in c("data.table","dplyr","ggplot2","R.utils")) if(!requireNamespace(p, quietly=TRUE)) stop("R package missing: ", p)' || die "R packages"

step "1 fixed R script: $RFIX  (original get_Varlist, qualifying clause reduced to M-CAP > threshold)"
OLDC='df$Func.refGene %in% disruptive | df$ExonicFunc.refGene %in% disruptive | df$M.CAP_score>y'
n_old=$(grep -cF "$OLDC" "$RSRC" || true)
echo "   occurrences of the LoF-or-M-CAP clause in $(basename "$RSRC"): $n_old | lines mentioning 'disruptive':"
grep -n 'disruptive' "$RSRC" | cut -c1-160 | sed 's/^/     /'
[ "$n_old" -eq 1 ] || die "expected exactly one qualifying clause in $RSRC"
awk -v old="$OLDC" '{ i = index($0, old); if (i) { $0 = substr($0, 1, i-1) "df$M.CAP_score>y" substr($0, i+length(old)) "  # FIX 05c: missense-only (M-CAP > threshold)" } print }' "$RSRC" > "$RFIX.tmp"
[ "$(grep -c 'disruptive |' "$RFIX.tmp" || true)" -eq 0 ] || die "patch incomplete"
mv "$RFIX.tmp" "$RFIX"
Rscript -e 'invisible(parse(file = commandArgs(TRUE)[1])); cat("   R syntax check passed\n")' "$RFIX" || die "syntax error in $RFIX"
diff "$RSRC" "$RFIX" | sed 's/^/     /' || true

step "2 patched shell copies in $SHFIX"
F1=get_RVBurdenMatrix_internal_rvexcaliber.clean.sh; F2=get_RVBurdenMatrix_gnomAD_rvexcaliber.clean.sh
cp -f "$SHDIR/$F1" "$SHDIR/$F2" "$SHFIX/"
sed -i 's#Rscripts/get_Varlist_rvexcaliber_missense_only\.R #Rscripts/get_Varlist_rvexcaliber_missense_fixed.R #' "$SHFIX/$F1" "$SHFIX/$F2"
sed -i 's#^\([[:space:]]*\)annovar_annotation[[:space:]]*$#\1if [[ -s ${outdir}/${internal_dataset}.hg19_multianno.txt.gz ]]; then echo "Reusing existing ANNOVAR output (05c)"; [[ -s ${outdir}/${internal_dataset}.hg19_multianno.txt ]] || zcat ${outdir}/${internal_dataset}.hg19_multianno.txt.gz > ${outdir}/${internal_dataset}.hg19_multianno.txt; else annovar_annotation; fi#' "$SHFIX/$F1"
c1=$(grep -c 'get_Varlist_rvexcaliber_missense_fixed.R' "$SHFIX/$F1" || true); c2=$(grep -c 'get_Varlist_rvexcaliber_missense_fixed.R' "$SHFIX/$F2" || true)
c3=$(grep -c 'Reusing existing ANNOVAR output' "$SHFIX/$F1" || true)
echo "   step 1A calls fixed R: $c1 | step 2 calls fixed R: $c2 | ANNOVAR reuse patch: $c3"
[ "$c1" -eq 1 ] && [ "$c2" -eq 1 ] && [ "$c3" -eq 1 ] || die "patching failed"
diff "$SHDIR/$F1" "$SHFIX/$F1" | sed 's/^/     /' || true

VL0=$WORK/${DS}_variant_list_extracted_${TAG}.txt.gz; BM0=$WORK/${DS}_RVBurdenMatrix_${TAG}.txt.gz
DONE1A=no
if [ -s "$VL0" ] && [ -s "$BM0" ] && [ "$(zrows "$VL0")" -gt 1000 ] && [ "$( { zcat "$VL0" | awk '{print NF}' | sort -u | paste -sd, -; } || true)" = "6" ]; then DONE1A=yes; fi
echo "   valid step-1A output already present: $DONE1A"
if [ "$DONE1A" = no ]; then
step "3 clean stale outputs of the failed run (annotation kept)"
rm -f "$WORK"/*_gnomAD_variant_list_extracted_* "$WORK"/${DS}_variant_list_extracted_* "$WORK"/${DS}_RVBurdenMatrix_* \
      "$WORK"/rvexcaliber_weights_* "$WORK"/gnomAD_RVBurdenMatrix_* "$WORK"/gnomAD_variant_list_extracted_*
ls -la "$WORK" | sed 's/^/   /'

step "4 step 1A (all functions except ANNOVAR) with fixed variant filter"
run bash "$SHFIX/$F1" "$INDIR" "$WORK" "$DS" "$GMAF" "$IMAF" "$MCAP" "$ETH" "$COV" "$JOBS"
else
  step "3+4 skipped (valid step-1A output from previous 05c run)"
fi
need "$WORK/${DS}_RVBurdenMatrix_${TAG}.txt.gz"
VLR=$(zrows "$WORK/${DS}_variant_list_extracted_${TAG}.txt.gz")
GN=$(( $(zrows "$WORK/${DS}_RVBurdenMatrix_${TAG}.txt.gz") - 1 ))
echo "   qualifying case variants: $VLR | genes in case matrix: $GN"
[ "$VLR" -gt 1000 ] || die "variant list still (almost) empty"
NFV=$(zcat "$WORK/${DS}_variant_list_extracted_${TAG}.txt.gz" | awk '{print NF}' | sort -u | paste -sd, -)
echo "   variant list columns: $NFV (expected 6: ID allele 1 gene gnomAD_AF gene_cumAF)"; [ "$NFV" = "6" ] || die "unexpected variant-list format"
echo "   example rows:"; { zcat "$WORK/${DS}_variant_list_extracted_${TAG}.txt.gz" | head -3 | sed 's/^/     /'; } || true
echo "   case variant classes (from pruned annotation, qualifying IDs):"
{ zcat "$WORK/${DS}_variant_list_extracted_${TAG}.txt.gz" | awk '{print $1}' | sort -u > "$WORK/05c_q.ids"; } || true
{ zcat "$WORK/${DS}_pruned_annotation_for_R_input_${COV}.txt.gz" | awk -F'\t' 'NR==FNR{q[$1]=1; next} FNR>1{id=$1":"$2":"$3":"$4; if(id in q) c[$8]++} END{for(k in c) printf "     %s %d\n", k, c[k]}' "$WORK/05c_q.ids" - | sort -k2,2nr; } || true

step "5 step 2 gnomAD burden matrix (fixed variant filter)"
[ -s "$WORK/gnomAD_RVBurdenMatrix_${DS}_${TAG}.txt.gz" ] && echo "   already present – skipped" || run bash "$SHFIX/$F2" "$WORK" "$DS" "$GMAF" "$MCAP" "$ETH" "$COV"
need "$WORK/gnomAD_RVBurdenMatrix_${DS}_${TAG}.txt.gz"
echo "   gnomAD qualifying variants: $(zrows "$WORK/gnomAD_variant_list_extracted_${TAG}.txt.gz") | genes in gnomAD matrix: $(( $(zrows "$WORK/gnomAD_RVBurdenMatrix_${DS}_${TAG}.txt.gz") - 1 ))"

step "6 step 3 summary associations (iCF/gCF, default MIGen ranking)"
RES=$WORK/rvexcaliber_testing_${DS}_SummaryAssociations_${TAG}_iCFgCFadjust_rvexcaliber_base.txt
[ -s "$RES" ] && echo "   already present – skipped" || run bash get_SummaryAssociations_iCF_gCF_rvexcaliber.clean.sh "$WORK" "$DS" NA "$GMAF" "$MCAP" "$ETH" "$COV" "$ADJ"
need "$RES"
RESF=${RES/SummaryAssociations_/SummaryAssociations_allele_filter_}
LAMF=$(ls "$WORK"/rvexcaliber_testing_${DS}_GenomicInflationFactor_allele_filter_${TAG}*iCFgCFadjust* 2>/dev/null | head -1 || true)

step "7 summary"
cat > "$WORK/05_summary.R" <<'EOF'
suppressMessages(library(data.table))
a <- commandArgs(TRUE); res <- a[1]; resf <- a[2]; lamf <- a[3]; v1 <- a[4]; out <- a[5]; cls <- a[6]
rd <- function(f) { d <- fread(f); setnames(d, c("Gene","test_allele_count","test_gnomAD_allele_count_iCFgCFadjust","test_P_rvexcaliber_base_iCFgCFadjust"),
                                         c("Gene","OAC","EXP","P"), skip_absent = TRUE); d[, FDR := p.adjust(P, "BH")][, OR := OAC / EXP][] }
r <- rd(res); n <- nrow(r)
lam_med <- median(qchisq(r$P, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
s <- c(sprintf("RV-EXCALIBER v2 | class = %s | %s", cls, format(Sys.time())),
       sprintf("unfiltered: genes %d | lambda(median p) %.3f | p<0.05: %d (exp %.0f) | p>0.95: %d | Bonferroni (p<%.2e): %d | FDR<0.05: %d",
               n, lam_med, r[P < 0.05, .N], 0.05 * n, r[P > 0.95, .N], 0.05 / n, r[P < 0.05 / n, .N], r[FDR < 0.05, .N]))
f <- NULL
if (file.exists(resf)) { f <- rd(resf); s <- c(s, sprintf("allele-filtered: genes %d | Bonferroni (p<%.2e): %d | FDR<0.05: %d",
                                                         nrow(f), 0.05 / nrow(f), f[P < 0.05 / nrow(f), .N], f[FDR < 0.05, .N])) }
if (nzchar(lamf) && file.exists(lamf)) s <- c(s, paste("GenomicInflationFactor file:", paste(readLines(lamf), collapse = " | ")))
top <- r[order(P)][1:min(20, n)]
if (!is.null(f)) top[, in_allele_filter_run := Gene %in% f$Gene]
fwrite(top, out, sep = "\t")
v1g <- c("TTK","ATXN3","ZNF266","ERCC2","C10orf82","DENND4B","PLG","PTPRN2")
cmp <- merge(data.table(Gene = v1g), r[, .(Gene, OAC_v2 = OAC, EXP_v2 = round(EXP, 2), P_v2 = signif(P, 3), FDR_v2 = signif(FDR, 3))], all.x = TRUE)
if (file.exists(v1)) cmp <- merge(cmp, rd(v1)[, .(Gene, OAC_v1 = OAC, EXP_v1 = round(EXP, 2), P_v1 = signif(P, 3))], by = "Gene", all.x = TRUE)
writeLines(s); cat("\ntop 20 genes (v2):\n"); print(top[, .(Gene, OAC, EXP = round(EXP, 2), OR = round(OR, 1), P = signif(P, 3), FDR = signif(FDR, 3))])
cat("\nv1 top genes in v2:\n"); print(cmp)
writeLines(c(s, "", capture.output(print(top)), "", capture.output(print(cmp))), sub("_top_genes_", "_summary_", sub("\\.tsv$", ".txt", out)))
EOF
Rscript "$WORK/05_summary.R" "$RES" "$RESF" "${LAMF:-}" "$V1MAIN" "$WORK/05_top_genes_missense.tsv" "missense (fixed filter)"

step "8 copy results to $FINAL"
cp -f "$WORK"/rvexcaliber_testing_* "$WORK"/rvexcaliber_gnomAD_* "$WORK"/*RVBurdenMatrix_* "$WORK"/*variant_list* \
      "$WORK"/*pruned_annotation* "$WORK"/05_* "$RFIX" "$FINAL"/ 2>/dev/null || true
cp -f "$LOG" "$CMDS" "$FINAL"/ 2>/dev/null || true
ls -la "$FINAL" | sed 's/^/   /'
echo; echo "== done $(date '+%F %T') | summary: $FINAL/05_summary_missense.txt"
