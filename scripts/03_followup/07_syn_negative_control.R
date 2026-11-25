#!/usr/bin/env Rscript
# =====================================================================
# 07_syn_negative_control.R  –  synonymous negative control for the v2 top genes (VM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · after 06_partA_v2.R (24.11.2026)
# Question: is the missense excess in ERCC2 / PLG (and other top genes) specific to missense variants, or does it
#           also appear for rare SYNONYMOUS variants (→ gene-level technical difference vs gnomAD: capture, coverage,
#           mapping, e.g. PLG/LPA paralogy)?
# Method (raw, i.e. without RV-EXCALIBER iCF/gCF calibration; calibrated by control genes instead):
#   per gene g and class c (missense M-CAP > 0.025 | synonymous):
#     O = cohort alleles (gnomAD NFE AF ≤ 0.001, cohort AF ≤ 0.01, hcc exons, PASS-extracted set)
#     E = 2 × 1,535 × Σ gnomAD NFE AF over gnomAD variants of class c in g (AF ≤ 0.001, hcc exons)
#   O/E per gene relative to the median O/E of 60 control genes (random, EXP 2–30 in the v2 result, seed 2611)
#   → rel_missense vs rel_synonymous; Poisson p for the synonymous excess
# Inputs (read only): v2 result + variant lists, cohort multianno (ANNOVAR), ANNOVAR humandb (refGene, gnomad211_exome),
#   RV-EXCALIBER hcc BED, v2 PLINK set. ANNOVAR geneanno is run on a small gnomAD subset (temp dir).
# Output: /mnt/hgfs/D/WES/FOR_REPORT/syn_control/  S_gene_table.tsv, S_summary.txt, S_log.txt
# Usage: Rscript 07_syn_negative_control.R        (≈5–15 min; most time = scanning hg19_gnomad211_exome.txt)
# =====================================================================
suppressMessages(library(data.table))
options(width = 220)
WES   <- "/mnt/hgfs/D/WES"; HOME <- Sys.getenv("HOME")
DS    <- "my_samples_plink_hg19_v2u"; TAG <- "0.001_0.025_nfe_hcc"
FINAL <- file.path(WES, "05_rvexcaliber_output_v2_missense"); WORK <- file.path(HOME, "rvx_v2_missense")
pick  <- function(...) { for (f in c(...)) if (file.exists(f)) return(f); stop("none found: ", paste(c(...), collapse = " | ")) }
MAIN  <- pick(file.path(FINAL, paste0("rvexcaliber_testing_", DS, "_SummaryAssociations_", TAG, "_iCFgCFadjust_rvexcaliber_base.txt")))
VLF   <- pick(file.path(FINAL, paste0(DS, "_variant_list_extracted_", TAG, ".txt.gz")), file.path(WORK, paste0(DS, "_variant_list_extracted_", TAG, ".txt.gz")))
GVLF  <- pick(file.path(FINAL, paste0("gnomAD_variant_list_extracted_", TAG, ".txt.gz")), file.path(WORK, paste0("gnomAD_variant_list_extracted_", TAG, ".txt.gz")))
MAF   <- pick(file.path(WORK, paste0(DS, ".hg19_multianno.txt.gz")), file.path(WORK, paste0(DS, ".hg19_multianno.txt")))
PLK   <- file.path(WES, "FOR_REPORT/v2/plink_in", DS)
ANNO  <- file.path(WES, "04_software/annovar"); HDB <- file.path(ANNO, "humandb")
GNOM  <- file.path(HDB, "hg19_gnomad211_exome.txt"); REFG <- file.path(HDB, "hg19_refGene.txt")
HCC   <- "/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main/gnomAD_211_exomes_hg19_coverage/ALL_CHROM_gnomad_refGene_exons_highcov20X.mod.bed"
OUT   <- file.path(WES, "FOR_REPORT/syn_control"); TMPD <- file.path(OUT, "tmp")
N2    <- 2 * 1535; GMAF <- 0.001; IMAF <- 0.01; MCAP <- 0.025
NTOP  <- 10; NCTRL <- 60
dir.create(TMPD, recursive = TRUE, showWarnings = FALSE)
zz <- file(file.path(OUT, "S_log.txt"), open = "wt"); sink(zz, split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
zread <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = paste("zcat", shQuote(f)), ...) else fread(f, ...)
sys <- function(cmd) { r <- system(cmd); if (r != 0) stop("command failed: ", cmd) }
for (f in c(GNOM, REFG, HCC, paste0(PLK, ".bim"), file.path(ANNO, "annotate_variation.pl"))) if (!file.exists(f)) stop("missing ", f)
say("synonymous negative control | %s", format(Sys.time()))

## 1 genes: top + controls
main <- fread(MAIN)
setnames(main, c("Gene", "test_allele_count", "test_gnomAD_allele_count_iCFgCFadjust", "test_P_rvexcaliber_base_iCFgCFadjust"), c("Gene", "OAC", "EXP", "P"))
top <- main[order(P)][1:NTOP, Gene]
set.seed(2611); ctrl <- sample(main[EXP >= 2 & EXP <= 30 & !(Gene %in% top), Gene], NCTRL)
genes <- data.table(Gene = c(top, ctrl), set = c(rep("top", NTOP), rep("control", NCTRL)))
say("top genes: %s | %d control genes (EXP 2-30)", paste(top, collapse = ", "), NCTRL)

## 2 gene regions (refGene, primary contigs)
rg <- fread(REFG, header = FALSE, select = c(3, 5, 6, 13), col.names = c("chr", "s", "e", "Gene"))
rg <- rg[Gene %in% genes$Gene & !grepl("_", chr)][, chr := sub("^chr", "", chr)]
reg <- rg[, .(s = min(s) - 50L, e = max(e) + 50L), by = .(Gene, chr)][, .SD[which.max(e - s)], by = Gene]
fwrite(reg[, .(chr, s, e, Gene)], file.path(TMPD, "regions.tsv"), sep = "\t", col.names = FALSE)
say("regions for %d of %d genes", nrow(reg), nrow(genes))

## 3 gnomAD variants in these regions (AF_nfe), ANNOVAR geneanno, hcc filter
gsub_f <- file.path(TMPD, "gnomad_sub.avinput")
if (!file.exists(gsub_f)) {
  say("scanning %s ...", GNOM)
  awk <- paste0("awk -F'\\t' -v OFS='\\t' -v R=", shQuote(file.path(TMPD, "regions.tsv")),
    " 'BEGIN{while((getline l<R)>0){split(l,a,\"\\t\"); n[a[1]]++; S[a[1],n[a[1]]]=a[2]; E[a[1],n[a[1]]]=a[3]; G[a[1],n[a[1]]]=a[4]}}",
    " NR==1{for(i=1;i<=NF;i++) if($i==\"AF_nfe\") k=i; next}",
    " {c=$1; sub(/^chr/,\"\",c); if(!(c in n)) next; for(j=1;j<=n[c];j++) if($2>=S[c,j] && $2<=E[c,j]){print c,$2,$3,$4,$5,$k,G[c,j]; break}}' ",
    shQuote(GNOM), " > ", shQuote(gsub_f))
  sys(awk)
}
say("gnomAD variants in regions: %d", length(readLines(gsub_f)))
sys(paste("perl", shQuote(file.path(ANNO, "annotate_variation.pl")), "-geneanno -buildver hg19 -dbtype refGene -outfile",
          shQuote(file.path(TMPD, "gnomad_sub")), shQuote(gsub_f), shQuote(HDB), "> /dev/null 2>&1"))
ev <- fread(file.path(TMPD, "gnomad_sub.exonic_variant_function"), header = FALSE, sep = "\t")
# columns: line#, type, AAChange, chr, start, end, ref, alt, AF_nfe, region_gene
setnames(ev, c("line", "type", "aa", "chr", "start", "end", "ref", "alt", "AF_nfe", "Gene_region"))
ev[, AF_nfe := suppressWarnings(as.numeric(AF_nfe))]
ev[, gene_anno := sub(":.*$", "", aa)]
ev <- ev[gene_anno == Gene_region & type %in% c("synonymous SNV", "nonsynonymous SNV") & !is.na(AF_nfe) & AF_nfe > 0 & AF_nfe <= GMAF]
fwrite(ev[, .(chr, start - 1L, end, line)], file.path(TMPD, "ev.bed"), sep = "\t", col.names = FALSE)
hcc_f <- file.path(TMPD, "hcc_nochr.bed")
if (!file.exists(hcc_f)) sys(paste("sed 's/^chr//'", shQuote(HCC), ">", shQuote(hcc_f)))
sys(paste("bedtools intersect -u -a", shQuote(file.path(TMPD, "ev.bed")), "-b", shQuote(hcc_f), ">", shQuote(file.path(TMPD, "ev_hcc.bed"))))
keep <- fread(file.path(TMPD, "ev_hcc.bed"), header = FALSE)$V4
gsyn <- ev[line %in% keep & type == "synonymous SNV", .(E_syn = N2 * sum(AF_nfe), n_gnomad_syn = .N), by = .(Gene = Gene_region)]
say("gnomAD rare synonymous variants (hcc): %d in %d genes", ev[line %in% keep & type == "synonymous SNV", .N], nrow(gsyn))

## 4 gnomAD missense expected (raw) from the RV-EXCALIBER gnomAD variant list (same filter as the main run)
gvl <- zread(GVLF, header = FALSE, select = c(4, 5), col.names = c("Gene", "AF"))
gmis <- gvl[Gene %in% genes$Gene, .(E_mis = N2 * sum(as.numeric(AF), na.rm = TRUE), n_gnomad_mis = .N), by = Gene]

## 5 cohort: synonymous from multianno, missense from v2 variant list → allele counts (plink)
ma <- zread(MAF, select = c("Chr", "Start", "Ref", "Alt", "Gene.refGene", "ExonicFunc.refGene", "AF_nfe"),
            colClasses = list(character = c("Chr", "Ref", "Alt", "AF_nfe")))
ma[, AFn := suppressWarnings(as.numeric(AF_nfe))][is.na(AFn), AFn := 0]
csyn <- ma[ExonicFunc.refGene == "synonymous SNV" & Gene.refGene %in% genes$Gene & AFn <= GMAF, .(ID = paste(Chr, Start, Ref, Alt, sep = ":"), Gene = Gene.refGene)]
vl <- zread(VLF, header = FALSE, select = c(1, 2, 4), col.names = c("ID", "A_cnt", "Gene"))[Gene %in% genes$Gene]
ids <- unique(c(csyn$ID, vl$ID)); writeLines(ids, file.path(TMPD, "cohort.ids"))
sys(paste("plink --bfile", shQuote(PLK), "--extract", shQuote(file.path(TMPD, "cohort.ids")), "--freq counts --allow-no-sex --out",
          shQuote(file.path(TMPD, "cohort")), "> /dev/null 2>&1"))
fq <- fread(file.path(TMPD, "cohort.frq.counts"), colClasses = list(character = c("SNP", "A1", "A2")))[!duplicated(SNP)]
fq[, alt := sub("^.*:", "", SNP)][, AC := fifelse(A1 == alt, C1, C2)][, AF := AC / (C1 + C2)]
csyn <- merge(csyn, fq[, .(ID = SNP, AC, AF)], by = "ID")[AF <= IMAF & AC > 0]
osyn <- csyn[, .(O_syn = sum(AC), n_cohort_syn = .N), by = Gene]
omis <- merge(vl, fq[, .(ID = SNP, AC)], by = "ID")[, .(O_mis = sum(AC)), by = Gene]

## 6 table + calibration by control genes
tab <- Reduce(function(a, b) merge(a, b, by = "Gene", all.x = TRUE), list(genes, main[, .(Gene, OAC, EXP, P)], omis, gmis, osyn, gsyn))
for (v in c("O_mis", "O_syn", "n_cohort_syn")) set(tab, which(is.na(tab[[v]])), v, 0L)
tab[, `:=`(OE_mis = O_mis / E_mis, OE_syn = O_syn / E_syn)]
cal_mis <- tab[set == "control" & is.finite(OE_mis), median(OE_mis)]
cal_syn <- tab[set == "control" & is.finite(OE_syn), median(OE_syn)]
tab[, `:=`(rel_mis = OE_mis / cal_mis, rel_syn = OE_syn / cal_syn)]
tab[, p_syn_excess := fifelse(E_syn > 0, ppois(O_syn - 1, E_syn * cal_syn, lower.tail = FALSE), NA_real_)]
setorder(tab, set, P)
say("control-gene calibration: median O/E missense %.2f | synonymous %.2f (n controls with data: %d / %d)",
    cal_mis, cal_syn, tab[set == "control" & is.finite(OE_mis), .N], tab[set == "control" & is.finite(OE_syn), .N])
say("control genes: IQR rel_missense %.2f–%.2f | rel_synonymous %.2f–%.2f",
    tab[set == "control", quantile(rel_mis, .25, na.rm = TRUE)], tab[set == "control", quantile(rel_mis, .75, na.rm = TRUE)],
    tab[set == "control", quantile(rel_syn, .25, na.rm = TRUE)], tab[set == "control", quantile(rel_syn, .75, na.rm = TRUE)])
cat("\n"); print(tab[set == "top", .(Gene, P = signif(P, 3), OAC, EXP = round(EXP, 1), O_mis, E_mis = round(E_mis, 1), rel_mis = round(rel_mis, 2),
                                     O_syn, E_syn = round(E_syn, 1), rel_syn = round(rel_syn, 2), p_syn_excess = signif(p_syn_excess, 2))])
fwrite(tab, file.path(OUT, "S_gene_table.tsv"), sep = "\t")
summ <- unlist(lapply(c("ERCC2", "PLG"), function(gn) { x <- tab[Gene == gn]
  sprintf("%s: missense O/E %.1f× control median (O %d, E %.1f) vs synonymous O/E %.2f× (O %d, E %.1f; Poisson p for excess %.2g) → %s",
          gn, x$rel_mis, x$O_mis, x$E_mis, x$rel_syn, x$O_syn, x$E_syn, x$p_syn_excess,
          ifelse(is.na(x$rel_syn), "no synonymous data",
                 ifelse(x$rel_syn < 1.5 & x$p_syn_excess > 0.05, "excess specific to missense (no technical excess detected)",
                        "synonymous excess present – technical explanation possible"))) }))
summ <- c(summ, sprintf("Calibration: median raw O/E in %d control genes = %.2f (missense), %.2f (synonymous).", NCTRL, cal_mis, cal_syn))
writeLines(summ, file.path(OUT, "S_summary.txt")); cat("\n"); cat(summ, sep = "\n")
sink(); close(zz)
