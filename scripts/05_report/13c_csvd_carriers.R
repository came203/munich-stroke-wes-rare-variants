#!/usr/bin/env Rscript
## =============================================================================
## 13c_csvd_carriers.R – rare variants and carriers in NOTCH3, HTRA1, COL4A1, COL4A2, TREX1 (Munich v2, N 1,535)
## Independent of the hcc restriction used by RV-EXCALIBER: all bim variants inside the gene bounds (hg19 refGene)
## are annotated with the same ANNOVAR setup (refGene, gnomad211_exome, dbnsfp35c).
## Rare protein-altering = Func exonic / exonic;splicing / splicing, not synonymous/unknown;
##   gnomAD 2.1.1 exome AF < 0.001 in all of NFE/AFR/SAS/EAS/AMR (absent = rare); cohort alt AF < 0.01, AC >= 1
## Flags: qualifies in the main RV-EXCALIBER test (variant in the 53,489 list), in hcc, NOTCH3 cysteine-altering
##   missense in exons 2-24 (EGFr domains, CADASIL-typical)
## Run: Rscript /mnt/hgfs/D/WES/13c_csvd_carriers.R   (partC_qualify_fun.R next to the script; plink2, perl/ANNOVAR)
## Output: FOR_REPORT/partC/ C_csvd_variants.tsv (rare protein-altering variants, counts only),
##         C_csvd_all_coding.tsv (all coding/splice variants incl. common), C_csvd_summary.tsv, 13c_csvd.log
##         ~/partC_work/13c/ carriers_INDIVIDUAL_LEVEL.tsv (local only – consent)
## =============================================================================
suppressPackageStartupMessages(library(data.table))
GENES <- c("NOTCH3", "HTRA1", "COL4A1", "COL4A2", "TREX1")
CANON <- c(NOTCH3 = "NM_000435", HTRA1 = "NM_002775", COL4A1 = "NM_001845", COL4A2 = "NM_001846", TREX1 = "NM_033629")
BASE <- "/mnt/hgfs/D/WES"; OUT <- file.path(BASE, "FOR_REPORT/partC")
V2U <- file.path(BASE, "FOR_REPORT/v2/plink_in/my_samples_plink_hg19_v2u")
ANNOVAR <- file.path(BASE, "04_software/annovar")
HCC_BED <- "/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main/gnomAD_211_exomes_hg19_coverage/ALL_CHROM_gnomad_refGene_exons_highcov20X.mod.bed"
F_VL <- path.expand("~/rvx_v2_missense/my_samples_plink_hg19_v2u_variant_list_extracted_0.001_0.025_nfe_hcc.txt.gz")
SCR <- path.expand("~/partC_work/13c"); dir.create(SCR, recursive = TRUE, showWarnings = FALSE)
a <- commandArgs(FALSE); sf <- sub("^--file=", "", a[grep("^--file=", a)])
source(file.path(if (length(sf) == 1) dirname(normalizePath(sf)) else getwd(), "partC_qualify_fun.R"))
PLINK2 <- Sys.which("plink2"); if (!nzchar(PLINK2)) stop("plink2 not found")
sink(file.path(OUT, "13c_csvd.log"), split = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
say("13c_csvd_carriers.R  %s  genes: %s", format(Sys.time()), paste(GENES, collapse = ", "))

## ---------- gene bounds (hg19 refGene, all transcripts) ----------
rg <- fread(file.path(ANNOVAR, "humandb/hg19_refGene.txt"), header = FALSE, select = c(3, 5, 6, 13),
            col.names = c("chrom", "txStart", "txEnd", "Gene"))
rg <- rg[Gene %chin% GENES & !grepl("_", chrom)]
gb <- rg[, .(chr = sub("^chr", "", chrom[1]), start = min(txStart) + 1L, end = max(txEnd)), by = Gene]
print(gb)

## ---------- bim variants in bounds -> ANNOVAR avinput ----------
bim <- fread(paste0(V2U, ".bim"), header = FALSE, col.names = c("chr", "ID", "cm", "pos", "a1", "a2"), colClasses = c("character", "character", "numeric", "integer", "character", "character"))
p <- tstrsplit(bim$ID, ":", fixed = TRUE); bim[, `:=`(ref = p[[3]], alt = p[[4]])]
sel <- rbindlist(lapply(seq_len(nrow(gb)), function(i) bim[chr == gb$chr[i] & pos >= gb$start[i] & pos <= gb$end[i]][, region := gb$Gene[i]]))
to_av <- function(pos, ref, alt) {
  n <- 0L; m <- min(nchar(ref), nchar(alt))
  while (n < m && substr(ref, n + 1, n + 1) == substr(alt, n + 1, n + 1)) n <- n + 1L
  r <- substring(ref, n + 1); al <- substring(alt, n + 1)
  if (r == "")       list(s = pos + n - 1L, e = pos + n - 1L, r = "-", a = al)
  else if (al == "") list(s = pos + n, e = pos + n + nchar(r) - 1L, r = r, a = "-")
  else               list(s = pos + n, e = pos + n + nchar(r) - 1L, r = r, a = al)
}
av <- rbindlist(Map(to_av, sel$pos, sel$ref, sel$alt))
sel[, `:=`(av_start = av$s, av_end = av$e, av_ref = av$r, av_alt = av$a)]
sel[, akey := make_key(chr, av_start, av_ref, av_alt)]
say("bim variants inside gene bounds: %d (SNVs %d, indels %d)", nrow(sel), sel[nchar(ref) == 1 & nchar(alt) == 1, .N], sel[nchar(ref) != 1 | nchar(alt) != 1, .N])
fwrite(sel[, .(chr, av_start, av_end, av_ref, av_alt)], file.path(SCR, "csvd.avinput"), sep = "\t", col.names = FALSE)

## ---------- ANNOVAR (same protocols as the main run / script 10) ----------
pre <- file.path(SCR, "csvd")
rc <- system(sprintf("perl %s %s %s -buildver hg19 -out %s -remove -protocol refGene,gnomad211_exome,dbnsfp35c -operation g,f,f -nastring . > %s.annovar.log 2>&1",
                     shQuote(file.path(ANNOVAR, "table_annovar.pl")), shQuote(file.path(SCR, "csvd.avinput")),
                     shQuote(file.path(ANNOVAR, "humandb")), shQuote(pre), shQuote(pre)))
if (rc != 0 || !file.exists(paste0(pre, ".hg19_multianno.txt"))) stop("ANNOVAR failed, see ", pre, ".annovar.log")
ma <- fread(paste0(pre, ".hg19_multianno.txt"), sep = "\t", quote = "", colClasses = "character", na.strings = c(".", ""), check.names = FALSE)
ma[, akey := make_key(Chr, Start, Ref, Alt)]
ma <- merge(ma, sel[, .(akey, ID, region)], by = "akey")
for (pcol in c("AF_nfe", "AF_afr", "AF_sas", "AF_eas", "AF_amr", "AF_popmax", "M-CAP_score")) set(ma, j = pcol, value = suppressWarnings(as.numeric(ma[[pcol]])))

## ---------- genotype counts + carriers (plink2) ----------
writeLines(ma$ID, file.path(SCR, "ids.txt"))
rc <- system(sprintf("%s --bfile %s --extract %s --geno-counts --export A --out %s > %s.plink.stdout 2>&1",
                     shQuote(PLINK2), shQuote(V2U), shQuote(file.path(SCR, "ids.txt")), shQuote(file.path(SCR, "g")), shQuote(file.path(SCR, "g"))))
if (rc != 0) stop("plink2 failed, see ", file.path(SCR, "g.log"))
gc <- fread(file.path(SCR, "g.gcount")); setnames(gc, sub("^#", "", names(gc)))
hetc <- grep("^HET", names(gc), value = TRUE)[1]; homc <- grep("^TWO_ALT", names(gc), value = TRUE)[1]
gc <- gc[, .(ID, plink_alt = ALT, HOM_REF = HOM_REF_CT, HET = get(hetc), HOM_ALT = get(homc), MISSING = MISSING_CT)]
ma <- merge(ma, gc, by = "ID")
ma[, `:=`(AC = HET + 2 * HOM_ALT, AN = 2 * (HOM_REF + HET + HOM_ALT))]
ma[, AF_cohort := AC / AN]

## ---------- classification ----------
hcc <- read_hcc_bed(HCC_BED)
vl <- fread(cmd = paste("zcat", shQuote(F_VL)), header = FALSE, colClasses = "character")
ma[, `:=`(in_hcc = flag_in_hcc(Chr, as.integer(Start), hcc, "rvx"), qualifies_main_test = ID %chin% vl$V1)]
pops <- c("AF_nfe", "AF_afr", "AF_sas", "AF_eas", "AF_amr")
ma[, gnomad_max5 := do.call(pmax, c(lapply(pops, function(p) ma[[p]]), na.rm = TRUE))]
ma[, gnomad_absent := is.na(gnomad_max5)]
ma[, coding := `Func.refGene` %chin% c("exonic", "exonic;splicing", "splicing") &
              !(`ExonicFunc.refGene` %chin% c("synonymous SNV", "unknown"))]
ma[, rare := (gnomad_absent | gnomad_max5 < 0.001) & AF_cohort < 0.01 & AC >= 1]
## gene + canonical transcript annotation
pick_tx <- function(gene_str, aach, gdet) {
  g <- intersect(strsplit(gene_str, ";", fixed = TRUE)[[1]], GENES)[1]
  if (is.na(g)) return(list(Gene = NA_character_, tx = NA_character_, exon = NA_character_, cDNA = NA_character_, protein = NA_character_))
  ent <- c(if (!is.na(aach)) strsplit(aach, ",", fixed = TRUE)[[1]], if (!is.na(gdet)) strsplit(gdet, "[,;]")[[1]])
  all_ent <- ent
  ent <- ent[grepl(CANON[[g]], ent, fixed = TRUE)]
  if (!length(ent)) ent <- all_ent
  if (!length(ent) || is.na(ent[1])) return(list(Gene = g, tx = NA_character_, exon = NA_character_, cDNA = NA_character_, protein = NA_character_))
  f <- strsplit(ent[1], ":", fixed = TRUE)[[1]]
  list(Gene = g, tx = f[grepl("^NM_", f)][1], exon = sub("exon", "", f[grepl("^exon", f)][1]),
       cDNA = f[grepl("^c\\.", f)][1], protein = f[grepl("^p\\.", f)][1])
}
ann <- rbindlist(Map(pick_tx, ma$`Gene.refGene`, ma$`AAChange.refGene`, ma$`GeneDetail.refGene`), use.names = TRUE)
ma <- cbind(ma, ann)
ma <- ma[!is.na(Gene)]
mm <- regmatches(ma$protein, regexec("^p\\.([A-Z])([0-9]+)([A-Z])$", ma$protein))
ma[, `:=`(aa_ref = vapply(mm, function(x) if (length(x) == 4) x[2] else NA_character_, ""),
          aa_alt = vapply(mm, function(x) if (length(x) == 4) x[4] else NA_character_, ""))]
ma[, notch3_cys_egfr := Gene == "NOTCH3" & `ExonicFunc.refGene` %chin% "nonsynonymous SNV" &
                        (aa_ref %chin% "C" | aa_alt %chin% "C") & suppressWarnings(as.integer(exon)) %in% 2:24]

cols <- c("Gene", "ID", "tx", "exon", "cDNA", "protein", "Func.refGene", "ExonicFunc.refGene", "AF_nfe", "gnomad_max5",
          "gnomad_absent", "AF_popmax", "M-CAP_score", "HET", "HOM_ALT", "MISSING", "AC", "AF_cohort",
          "qualifies_main_test", "in_hcc", "notch3_cys_egfr")
allc <- ma[coding == TRUE][order(Gene, as.integer(sub(".*:([0-9]+):.*", "\\1", ID))), ..cols]
rv <- ma[coding == TRUE & rare == TRUE][order(Gene, -AC), ..cols]
fwrite(allc, file.path(OUT, "C_csvd_all_coding.tsv"), sep = "\t")
fwrite(rv, file.path(OUT, "C_csvd_variants.tsv"), sep = "\t")

## ---------- carriers per gene (distinct individuals) ----------
raw <- fread(file.path(SCR, "g.raw"))
vcols <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
## plink2 --export A may count REF: carrier = dosage of the counted allele > 0 if counted == ALT (from ID), else < 2
idmap <- data.table(col = vcols, ID = sub("_[^_]+$", "", vcols), counted = sub("^.*_", "", vcols))
idmap[, alt := tstrsplit(ID, ":", fixed = TRUE)[[4]]]
idmap[, flip := counted != alt]
say("export A counted allele: ALT for %d, REF for %d variants", idmap[flip == FALSE, .N], idmap[flip == TRUE, .N])
carrier_matrix <- function(cols) {
  M <- as.matrix(raw[, cols, with = FALSE])
  fl <- idmap[match(cols, col), flip]
  C <- M > 0
  if (any(fl)) C[, fl] <- M[, fl, drop = FALSE] < 2
  C[is.na(C)] <- FALSE
  C
}
carr <- function(ids) {
  cc <- idmap[ID %chin% ids, col]
  if (!length(cc)) return(0L)
  sum(rowSums(carrier_matrix(cc)) > 0)
}
sm <- rbindlist(lapply(GENES, function(g) {
  x <- rv[Gene == g]
  data.table(Gene = g, rare_coding_variants = nrow(x), carriers = carr(x$ID),
             missense = x[`ExonicFunc.refGene` == "nonsynonymous SNV", .N],
             lof = x[`ExonicFunc.refGene` %chin% c("stopgain", "frameshift deletion", "frameshift insertion", "startloss", "stoploss") |
                     `Func.refGene` %chin% c("splicing", "exonic;splicing"), .N],
             inframe = x[`ExonicFunc.refGene` %chin% c("nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution"), .N],
             qualifying_variants = x[qualifies_main_test == TRUE, .N], carriers_qualifying = carr(x[qualifies_main_test == TRUE, ID]),
             outside_hcc = x[in_hcc == FALSE, .N],
             notch3_cys_variants = x[notch3_cys_egfr == TRUE, .N], notch3_cys_carriers = carr(x[notch3_cys_egfr == TRUE, ID]))
}))
say("\nper gene (rare = gnomAD AF < 0.001 in all 5 populations or absent, cohort AF < 0.01; protein-altering incl. splice +-2 bp):")
print(sm)
fwrite(sm, file.path(OUT, "C_csvd_summary.tsv"), sep = "\t")
say("\nNOTCH3 cysteine-altering EGFr variants (exons 2-24):")
print(rv[notch3_cys_egfr == TRUE, .(ID, exon, cDNA, protein, het = HET, hom = HOM_ALT, gnomAD_NFE = AF_nfe, MCAP = `M-CAP_score`, main = qualifies_main_test, hcc = in_hcc)])
say("\nrare protein-altering variants:")
options(datatable.print.nrows = 500)
print(rv[, .(Gene, ID, exon, cDNA, protein, cls = fifelse(is.na(`ExonicFunc.refGene`), `Func.refGene`, `ExonicFunc.refGene`),
             gnomAD_NFE = AF_nfe, het = HET, hom = HOM_ALT, MCAP = `M-CAP_score`, main = qualifies_main_test, hcc = in_hcc, cys = notch3_cys_egfr)])

## individual-level carrier file (local only)
cc <- idmap[ID %chin% rv$ID, col]
Cm <- carrier_matrix(cc)
w_ <- which(Cm, arr.ind = TRUE)
long <- data.table(IID = raw$IID[w_[, 1]], col = cc[w_[, 2]])
long <- merge(long, idmap[, .(col, ID)], by = "col")
long <- merge(long, rv[, .(ID, Gene, protein, cDNA)], by = "ID")
fwrite(long[, .(IID, Gene, ID, cDNA, protein)], file.path(SCR, "carriers_INDIVIDUAL_LEVEL.tsv"), sep = "\t")
say("carriers with >1 rare variant in the same gene: %d", long[, .N, by = .(IID, Gene)][N > 1, .N])
say("\nwrote C_csvd_variants.tsv, C_csvd_all_coding.tsv, C_csvd_summary.tsv, 13c_csvd.log; %s (individual level, local only)",
    file.path(SCR, "carriers_INDIVIDUAL_LEVEL.tsv"))
say("note: no ClinVar annotation in the local ANNOVAR humandb – pathogenicity classification to be done manually")
sink()
