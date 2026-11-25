## =============================================================================
## partC_qualify_fun.R – shared qualifying-variant rule for Part C (RVGRS_Plan 4b, D3/D4)  v3
## v3 (25.11.2026): class set (Func|ExonicFunc) instead of nonsyn-only, hcc "rvx" mode (s <= pos <= e),
##   MAC >= 1 optional (RV-EXCALIBER keeps AF_int = 0 variants); defaults = D4 rule
## Sourced by 09_weights.R (Munich reproduction) and 10_ukb_qualify (UKB), so that both
## cohorts are filtered by identical code. Requires data.table.
## v2 (25.11.2026): ordered keys (chr:pos:ref:alt, as bim/pvar IDs); column 'key' renamed
##   (data.table() reserves key=); gnomAD/internal comparison modes to match RV-EXCALIBER
##   getBinary() (strict '<', optionally all 5 gnomAD populations)
##
## Criteria, checked in this order (first failing criterion = "reason"):
##   1 autosome (1–22)        4 optional ExonicFunc.refGene == nonsynonymous_SNV
##   2 SNV (A/C/G/T)          5 M-CAP_score > 0.025 (strict)
##   3 inside hcc BED (hg19)  6 gnomAD AF (gnomad_mode)   7 MAC >= 1   8 internal MAF (int_col/int_cmp)
## gnomad_mode: nfe_maf_le   min(AF_nfe,1-AF_nfe) <= thr, missing = 0      (plan D3)
##              nfe_af_lt    AF_nfe < thr, missing = 0
##              allpop_af_lt all of AF_nfe/afr/sas/eas/amr < thr; all missing = pass,
##                           partly missing = fail                          (getBinary logic)
## =============================================================================

NONSYN <- c("nonsynonymous_SNV", "nonsynonymous SNV")
POPS   <- c("af_nfe", "af_afr", "af_sas", "af_eas", "af_amr")

## ordered key "chr:pos:ref:alt" (chr without "chr"), same scheme as bim/pvar IDs
make_key <- function(chr, pos, ref, alt) {
  chr <- sub("^chr", "", as.character(chr))
  pos <- suppressWarnings(as.integer(pos))
  out <- paste(chr, pos, ref, alt, sep = ":")
  out[is.na(chr) | is.na(pos) | is.na(ref) | is.na(alt)] <- NA_character_
  out
}

## IDs like "6:80717709:G:A" (optional "chr") -> normalised key + position key
parse_id_key <- function(id) {
  p <- tstrsplit(sub("^chr", "", id), ":", fixed = TRUE, fill = NA_character_)
  while (length(p) < 4) p[[length(p) + 1L]] <- rep(NA_character_, length(id))
  list(vkey   = make_key(p[[1]], p[[2]], p[[3]], p[[4]]),
       poskey = paste(p[[1]], suppressWarnings(as.integer(p[[2]])), sep = ":"))
}

## BED (0-based half-open) -> keyed 1-based closed intervals, chr without "chr"
read_hcc_bed <- function(f) {
  b <- fread(f, header = FALSE, select = 1:3, colClasses = "character", col.names = c("chr", "s", "e"))
  b <- b[!grepl("^(track|browser|#)", chr)]
  b[, `:=`(chr = sub("^chr", "", chr),
           start = suppressWarnings(as.integer(s)) + 1L,
           end   = suppressWarnings(as.integer(e)))]
  b <- b[!is.na(start) & !is.na(end) & start <= end, .(chr, start, end)]
  setkey(b, chr, start, end)
  b
}

## TRUE if 1-based position lies inside any interval
flag_in_bed <- function(chr, pos, bed) {
  q <- data.table(xid = seq_along(pos), chr = sub("^chr", "", as.character(chr)),
                  start = as.integer(pos), end = as.integer(pos))
  q <- q[!is.na(start)]
  ov <- foverlaps(q, bed, by.x = c("chr", "start", "end"), type = "within", nomatch = 0L, which = TRUE)
  out <- logical(length(pos))
  out[q$xid[unique(ov$xid)]] <- TRUE
  out
}

## numeric conversion; multi-valued cells ("0.1;0.3") -> max
num1 <- function(x, what = "") {
  y <- suppressWarnings(as.numeric(x))
  bad <- !is.na(x) & is.na(y)
  if (any(bad)) {
    y[bad] <- vapply(strsplit(x[bad], "[;,|]"),
                     function(z) suppressWarnings(max(as.numeric(z), na.rm = TRUE)), 0)
    y[is.infinite(y)] <- NA_real_
    message(sprintf("num1(%s): %d multi-valued/non-numeric cells -> max", what, sum(bad)))
  }
  y
}

## needed columns from ANNOVAR hg19_multianno.txt(.gz) (refGene, gnomad211_exome, dbnsfp35c)
read_multianno <- function(f) {
  con <- gzfile(f); h <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]; close(con)
  pick <- function(pat) {
    i <- grep(pat, h)
    if (length(i) != 1) stop(sprintf("multianno: %d columns match '%s'", length(i), pat))
    i
  }
  idx <- c(chr = pick("^Chr$"), start = pick("^Start$"), ref = pick("^Ref$"), alt = pick("^Alt$"),
           func = pick("^Func\\.refGene$"), gene = pick("^Gene\\.refGene$"),
           exfunc = pick("^ExonicFunc\\.refGene$"),
           af_nfe = pick("^AF_nfe$"), af_afr = pick("^AF_afr$"), af_sas = pick("^AF_sas$"),
           af_eas = pick("^AF_eas$"), af_amr = pick("^AF_amr$"), mcap = pick("^M.CAP_score$"))
  o <- order(idx)                                   # cut -f returns fields in file order
  catc <- if (grepl("\\.gz$", f)) "zcat" else "cat"
  dt <- fread(cmd = sprintf("%s %s | cut -f %s", catc, shQuote(f), paste(idx[o], collapse = ",")),
              sep = "\t", quote = "", header = TRUE, colClasses = "character", na.strings = c(".", ""))
  setnames(dt, names(idx)[o])
  dt[, `:=`(chr = sub("^chr", "", chr), pos = as.integer(start), mcap = num1(mcap, "M-CAP"))]
  for (p in POPS) set(dt, j = p, value = num1(dt[[p]], p))
  dt[, start := NULL]
  dt[]
}

## hcc membership. "rvx": BED start s <= pos <= e (RV-EXCALIBER/bedtools point overlap, found in D4)
##                 "strict": s+1 <= pos <= e (correct 1-based)
flag_in_hcc <- function(chr, pos, bed, mode = c("rvx", "strict")) {
  mode <- match.arg(mode)
  x <- flag_in_bed(chr, pos, bed)
  if (mode == "rvx") x <- x | flag_in_bed(chr, as.integer(pos) + 1L, bed)
  x
}

## variant class "Func.refGene|ExonicFunc.refGene" with spaces -> "_" and missing -> "."
make_class <- function(func, exfunc) {
  f <- fifelse(is.na(func), ".", func)
  e <- fifelse(is.na(exfunc), ".", gsub(" ", "_", exfunc, fixed = TRUE))
  paste(f, e, sep = "|")
}

## first failing criterion per variant ("qualifying" if none)
## needs: chr, ref, alt, <hcc_col>, class (if class_keep), mcap, POPS, mac, <int_col>
## v3 defaults = rule found in D4 (09_weights.R, 25.11.2026)
qualify_reason <- function(dt, mcap_thr = 0.025, gnomad_thr = 0.001, int_thr = 0.01,
                           gnomad_mode = c("allpop_af_lt", "nfe_maf_le", "nfe_af_lt"),
                           int_col = "maf_obs", int_cmp = c("lt", "le"),
                           class_keep = NULL, hcc_col = "in_hcc_rvx", require_mac = FALSE) {
  gnomad_mode <- match.arg(gnomad_mode); int_cmp <- match.arg(int_cmp)
  class_fail <- if (is.null(class_keep)) rep(FALSE, nrow(dt)) else !(dt$class %chin% class_keep)
  mac_fail   <- if (require_mac) is.na(dt$mac) | dt$mac < 1 else rep(FALSE, nrow(dt))
  acgt <- c("A", "C", "G", "T")
  snv  <- nchar(dt$ref) == 1L & nchar(dt$alt) == 1L & dt$ref %chin% acgt & dt$alt %chin% acgt
  nfe0 <- fifelse(is.na(dt$af_nfe), 0, dt$af_nfe)
  g_fail <- switch(gnomad_mode,
    nfe_maf_le   = pmin(nfe0, 1 - nfe0) > gnomad_thr,
    nfe_af_lt    = nfe0 >= gnomad_thr,
    allpop_af_lt = {
      nna <- Reduce(`+`, lapply(POPS, function(p) is.na(dt[[p]])))
      mx  <- do.call(pmax, c(lapply(POPS, function(p) dt[[p]]), na.rm = TRUE))
      (nna > 0 & nna < length(POPS)) | (nna == 0 & mx >= gnomad_thr)
    })
  imaf   <- dt[[int_col]]
  i_fail <- if (int_cmp == "le") imaf > int_thr else imaf >= int_thr
  fcase(!(dt$chr %chin% as.character(1:22)),            "1_non_autosomal",
        !snv,                                           "2_not_SNV",
        !dt[[hcc_col]],                                 "3_not_hcc",
        class_fail,                                     "4_class",
        is.na(dt$mcap) | dt$mcap <= mcap_thr,           "5_mcap",
        g_fail,                                         "6_gnomad",
        is.na(imaf) | mac_fail,                         "7_not_in_cohort",
        i_fail,                                         "8_internal_maf",
        default =                                       "qualifying")
}

## variant -> gene rows; mode "split" (one row per ';'-gene), "first", "asis"
assign_genes <- function(gene, mode = c("split", "first", "asis")) {
  mode <- match.arg(mode)
  if (mode == "asis") return(data.table(i = seq_along(gene), Gene = gene))
  s <- strsplit(gene, ";", fixed = TRUE)
  if (mode == "first") return(data.table(i = seq_along(gene), Gene = vapply(s, `[`, "", 1L)))
  data.table(i = rep(seq_along(s), lengths(s)), Gene = unlist(s))
}
