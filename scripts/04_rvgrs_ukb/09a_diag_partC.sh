#!/usr/bin/env bash
# =============================================================================
# 09a_diag_partC.sh  –  Part C (RVGRS in UKB), READ-ONLY diagnostics
# RVGRS_Plan.md 4b, script 09a (decisions D1–D6 confirmed 25.11.2026)
#
# Run on the VM:   bash /mnt/hgfs/D/WES/09a_diag_partC.sh
#   runtime ~5–15 min (chr22 plink2 --freq on the UKB analysis set is the slow part)
#
# Inputs are only read. Writes:
#   log      -> /mnt/hgfs/D/WES/FOR_REPORT/partC/09a_diag.log
#   scratch  -> ~/partC_diag/ (UKB keep file, chr22 freq, liftOver/REF/hcc test)
#
# Sections
#   0 tools + R packages        4 UKB pgen/pvar/psam format
#   1 Munich inputs for 09      5 UKB sample overlap -> keep file (n 438,694 / 9,308?)
#   2 hcc BED                   6 chr22: freq timing, rare SNVs, liftOver, REF check, hcc
#   3 ANNOVAR                   7 resources
# =============================================================================
set -u
export LC_ALL=C

BASE=/mnt/hgfs/D/WES
OUT=$BASE/FOR_REPORT/partC
SCR=$HOME/partC_diag
mkdir -p "$OUT" "$SCR"
LOG=$OUT/09a_diag.log
exec > >(tee "$LOG") 2>&1

# ---- paths (RVGRS_Plan H3) ---------------------------------------------------
UKB=/cluster2/dichganslab/rmalik2/UKB_FINAL
ISPHENO=/cluster2/dichganslab/rmalik2/REGENIE/FORCCR2/IS.pheno
ISEUR=/mnt/hgfs/D/data/UKB_IS_RISK/IS.EUR
INCID=/mnt/hgfs/D/data/UKB_IS_RISK/incident.is.txt
RVX=/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main
RVX2=$BASE/04_software/RV-EXCALIBER-main
HCC=$RVX/gnomAD_211_exomes_hg19_coverage/ALL_CHROM_gnomad_refGene_exons_highcov20X.mod.bed
ANNOVAR=$BASE/04_software/annovar
CHAIN=$BASE/03_references/hg38ToHg19.over.chain.gz
HG19=$BASE/03_references/hg19.fa
LIFTOVER=/mnt/hgfs/D/software/liftOver
WORK=$HOME/rvx_v2_missense
RES=$BASE/05_rvexcaliber_output_v2_missense
DS=my_samples_plink_hg19_v2u
TAG=0.001_0.025_nfe_hcc
V2U=$BASE/FOR_REPORT/v2/plink_in/$DS

sec(){ echo; echo "================ $* ================"; }
chk(){ if [ -e "$1" ]; then echo "OK      $1 ($(du -h "$1" 2>/dev/null | cut -f1))"; else echo "MISSING $1"; fi; }
find_tool(){
  local p; p=$(command -v "$1" 2>/dev/null)
  [ -z "$p" ] && p=$(find /home/rainer/WES_analysis/04_software -maxdepth 3 -type f -name "$1" -perm -u+x 2>/dev/null | head -1)
  echo "$p"
}

echo "09a_diag_partC.sh  $(date)  host $(hostname)"

# =============================================================================
sec "0. Tools"
PLINK2=${PLINK2:-$(find_tool plink2)}
BEDTOOLS=$(find_tool bedtools)
echo "plink2:   ${PLINK2:-NOT FOUND}";  [ -n "$PLINK2" ] && "$PLINK2" --version 2>&1 | head -1
echo "bedtools: ${BEDTOOLS:-NOT FOUND}"; [ -n "$BEDTOOLS" ] && "$BEDTOOLS" --version 2>&1 | head -1
echo "liftOver: $LIFTOVER  executable=$([ -x "$LIFTOVER" ] && echo yes || echo NO)"
echo "zstdcat:  $(command -v zstdcat || echo 'not found')"
echo "perl:     $(perl -e 'print $^V' 2>/dev/null)"
echo "R:        $(Rscript --version 2>&1 | head -1)"
Rscript -e 'for (p in c("data.table","R.utils","survival","pROC","nricens","logistf","ggplot2")) cat(sprintf("  R pkg %-10s %s\n", p, requireNamespace(p, quietly = TRUE)))'

# =============================================================================
sec "1. Munich inputs for 09_weights.R"
echo "--- work dir $WORK"
ls -la "$WORK" 2>&1 | head -60
echo "--- result dir $RES"
ls -la "$RES" 2>&1 | head -30

MAIN=$RES/rvexcaliber_testing_${DS}_SummaryAssociations_${TAG}_iCFgCFadjust_rvexcaliber_base.txt
chk "$MAIN"
[ -e "$MAIN" ] && { echo "rows (incl header): $(wc -l < "$MAIN")"; head -3 "$MAIN"; }
echo "--- allele-filter / lambda files:"
ls "$RES" 2>/dev/null | grep -E 'allele_filter|GenomicInflation' || echo "  none found"

echo "--- multianno"
MA=$(ls "$WORK"/*hg19_multianno.txt* 2>/dev/null)
echo "candidates: $MA"
MA1=$(echo "$MA" | head -1)
if [ -n "$MA1" ]; then
  echo "header (index:name):"
  zcat -f "$MA1" 2>/dev/null | head -1 | tr '\t' '\n' | nl -ba | head -100
  zcat -f "$MA1" 2>/dev/null | awk -F'\t' 'NR==1{print "header fields:",NF} NR==2{print "row2 fields:",NF; exit}'
  echo "rows 2-3, cols 1-10,20,52:"
  zcat -f "$MA1" 2>/dev/null | sed -n '2,3p' | cut -f1-10,20,52
  echo "data rows: $(zcat -f "$MA1" 2>/dev/null | tail -n +2 | wc -l)"
  echo "ExonicFunc.refGene (col 9) table:"
  zcat -f "$MA1" 2>/dev/null | tail -n +2 | cut -f9 | sort | uniq -c | sort -rn | head -15
  echo "Gene.refGene containing ';' among exonic: $(zcat -f "$MA1" 2>/dev/null | awk -F'\t' 'NR>1 && $6 ~ /exonic/ && $7 ~ /;/' | wc -l)"
fi

VL=$WORK/${DS}_variant_list_extracted_${TAG}.txt.gz
chk "$VL"
if [ -e "$VL" ]; then
  zcat "$VL" 2>/dev/null | head -5
  echo "rows: $(zcat "$VL" | wc -l)   unique col1 IDs: $(zcat "$VL" | awk '{print $1}' | sort -u | wc -l)   unique genes: $(zcat "$VL" | awk '{print $4}' | sort -u | wc -l)"
  echo "ID ':'-part counts:"; zcat "$VL" | awk '{n=split($1,a,":"); c[n]++} END{for(k in c) print "  "k" parts: "c[k]}'
  echo "field separator tab? $(zcat "$VL" | head -1 | grep -c $'\t')"
fi

BM=$WORK/${DS}_RVBurdenMatrix_${TAG}.txt.gz
chk "$BM"
[ -e "$BM" ] && { zcat "$BM" 2>/dev/null | head -1 | awk '{print "fields:",NF,"  first:",$1,$2,$3}'; echo "rows: $(zcat "$BM" | wc -l)"; }

PR=$(ls "$WORK"/*pruned_annotation_for_R_input_hcc.txt* 2>/dev/null | head -1)
echo "--- pruned annotation: ${PR:-not found}"
[ -n "$PR" ] && { zcat -f "$PR" 2>/dev/null | head -3 | cut -c1-300; echo "rows: $(zcat -f "$PR" | wc -l)"; }

echo "--- v2u PLINK set"
for e in bed bim fam; do chk "$V2U.$e"; done
[ -e "$V2U.bim" ] && { echo "bim rows: $(wc -l < "$V2U.bim")  fam rows: $(wc -l < "$V2U.fam")"; head -2 "$V2U.bim"; }

echo "--- RV-EXCALIBER variant rule (fixed get_Varlist)"
GV=$RVX/scripts/Rscripts/get_Varlist_rvexcaliber_missense_fixed.R
chk "$GV"; [ -e "$GV" ] && head -300 "$GV"
echo "--- exact step-1A ANNOVAR call"
[ -e "$WORK/05_commands_missense.txt" ] && cat "$WORK/05_commands_missense.txt" || echo "05_commands_missense.txt not found"
echo "--- hcc / internal-MAF / freq logic in shell scripts"
for d in "$WORK/shell_fixed" "$RVX2/scripts/shell" "$RVX/scripts/shell"; do
  [ -d "$d" ] || continue
  for f in "$d"/*.sh; do
    [ -e "$f" ] || continue
    hits=$(grep -nEi 'hcc|highcov|bedtools|intersect|internal|maf|--freq|--max-maf' "$f" | head -25)
    [ -n "$hits" ] && { echo "## $f"; echo "$hits"; }
  done
done

# =============================================================================
sec "2. hcc BED (hg19)"
chk "$HCC"
if [ -e "$HCC" ]; then
  head -3 "$HCC"
  awk -F'\t' '{n++; bp+=$3-$2; c[$1]++; if($1=="chr22"||$1=="22") bp22+=$3-$2}
              END{print "intervals:",n,"  total bp:",bp,"  chr22 bp:",bp22,"  contigs:",length(c)}' "$HCC"
  echo "contig names: $(cut -f1 "$HCC" | sort -u | head -30 | tr '\n' ' ')"
fi

# =============================================================================
sec "3. ANNOVAR"
chk "$ANNOVAR/table_annovar.pl"
chk "$ANNOVAR/annotate_variation.pl"
for db in hg19_refGene.txt hg19_refGeneMrna.fa hg19_gnomad211_exome.txt hg19_gnomad211_exome.txt.idx hg19_dbnsfp35c.txt hg19_dbnsfp35c.txt.idx; do chk "$ANNOVAR/humandb/$db"; done
echo "gnomad211_exome header: $(head -1 "$ANNOVAR/humandb/hg19_gnomad211_exome.txt" 2>/dev/null | cut -c1-300)"
echo "chain: "; chk "$CHAIN"; chk "$BASE/03_references/hg19ToHg38.over.chain.gz"
echo "hg19 ref:"; chk "$HG19"; chk "$HG19.fai"

# =============================================================================
sec "4. UKB pgen/pvar/psam"
ls -l "$UKB"/chr22.* 2>&1
miss=0
for c in $(seq 1 22); do
  for e in pgen psam; do [ -e "$UKB/chr$c.$e" ] || { echo "MISSING chr$c.$e"; miss=1; }; done
  [ -e "$UKB/chr$c.pvar" ] || [ -e "$UKB/chr$c.pvar.zst" ] || { echo "MISSING chr$c.pvar"; miss=1; }
done
[ $miss -eq 0 ] && echo "chr1-22 pgen/pvar/psam all present"
echo "psam md5 groups across chr1-22 (must be 1 group for summing scores):"
md5sum "$UKB"/chr{1..22}.psam 2>/dev/null | awk '{print $1}' | sort | uniq -c

PV=$UKB/chr22.pvar; PVCAT=cat; PFLAG=""
[ -e "$PV" ] || { PV=$UKB/chr22.pvar.zst; PVCAT=zstdcat; PFLAG=vzs; }
echo "pvar: $PV"
echo "## header lines: $($PVCAT "$PV" | grep -c '^##')"
$PVCAT "$PV" | grep -m1 '^#CHROM'
$PVCAT "$PV" | grep -v '^#' | head -5 | cut -c1-200
$PVCAT "$PV" | awk -F'\t' '!/^#/{
    n++; if($5 ~ /,/) multi++; if(length($4)==1 && length($5)==1) snv++
    id=$1":"$2":"$4":"$5; if($3==id) idok++; else if(bad<3){print "  ID != CHROM:POS:REF:ALT e.g.:",$1,$2,$3,$4,$5; bad++}
    k=$1":"$2; if(k in s) samepos++; s[k]=1; if(d[$3]++) dupid++ }
  END{printf "chr22 records %d | ALT with comma %d | SNV %d | ID==CHROM:POS:REF:ALT %d | records sharing a position %d | duplicate IDs %d\n", n,multi,snv,idok,samepos,dupid}'
echo "psam head:"; head -3 "$UKB/chr22.psam"; echo "psam rows: $(wc -l < "$UKB/chr22.psam")"

# =============================================================================
sec "5. UKB sample overlap (IS.pheno ∩ psam ∩ IS.EUR) -> keep file"
export PSAM="$UKB/chr22.psam" ISPHENO ISEUR INCID SCR
Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table))
ps  <- fread(Sys.getenv("PSAM")); setnames(ps, sub("^#", "", names(ps)))
ph  <- fread(Sys.getenv("ISPHENO")); eur <- fread(Sys.getenv("ISEUR")); inc <- fread(Sys.getenv("INCID"))
cat("psam cols:", names(ps), "| N", nrow(ps), "| FID==IID:", if ("FID" %in% names(ps)) all(ps$FID == ps$IID) else "no FID", "\n")
cat("IS.pheno cols:", names(ph), "| N", nrow(ph), "| ispheno:", paste(names(table(ph$ispheno)), table(ph$ispheno), sep = "=", collapse = " "), "\n")
cat("IS.EUR cols:", names(eur), "| N", nrow(eur), "| PRS NA:", sum(is.na(eur$PRS)), "\n")
cat("incident cols:", names(inc), "| N", nrow(inc), "\n")
for (d in list(ph = ph$IID, ps = ps$IID, eur = eur$IID)) if (anyDuplicated(d)) cat("WARNING duplicated IIDs\n")
ids <- Reduce(intersect, list(as.character(ph$IID), as.character(ps$IID), as.character(eur$IID)))
a <- ph[as.character(IID) %in% ids]
cat(sprintf("ANALYSIS SET n = %d, cases = %d (plan: 438,694 / 9,308)\n", nrow(a), sum(a$ispheno == 1)))
inc_ev <- as.character(inc[event == 1, eid])
cat(sprintf("  cases with incident event: %d | analysis set in incident file: %d | non-cases with incident event (should be 0): %d\n",
            sum(a$ispheno == 1 & as.character(a$IID) %in% inc_ev), sum(as.character(a$IID) %in% as.character(inc$eid)),
            sum(a$ispheno == 0 & as.character(a$IID) %in% inc_ev)))
k <- ps[as.character(IID) %in% ids]
cols <- intersect(c("FID", "IID"), names(k))
f <- file.path(Sys.getenv("SCR"), "keep_analysis.txt")
writeLines(paste0("#", paste(cols, collapse = "\t")), f)
fwrite(k[, ..cols], f, sep = "\t", append = TRUE, col.names = FALSE)
cat("keep file:", f, "rows", nrow(k), "\n")
EOF

# =============================================================================
sec "6. chr22: plink2 --freq timing, rare SNVs, liftOver, REF check, hcc"
if [ -n "$PLINK2" ] && [ -s "$SCR/keep_analysis.txt" ]; then
  t0=$(date +%s)
  "$PLINK2" --pfile "$UKB/chr22" $PFLAG --keep "$SCR/keep_analysis.txt" --snps-only just-acgt \
            --freq counts cols=+pos --threads 16 --memory 10000 --out "$SCR/ukb_chr22" > "$SCR/ukb_chr22.stdout" 2>&1
  rc=$?; echo "plink2 exit $rc, $(( $(date +%s) - t0 )) s"
  grep -Ei 'samples|variants|error|warning' "$SCR/ukb_chr22.log" | head -15
  echo "acount header: $(head -1 "$SCR/ukb_chr22.acount" 2>/dev/null)"
  awk -F'\t' -v bed="$SCR/c22_rare.hg38.bed" 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
     { a=$h["ALT_CTS"]; o=$h["OBS_CT"]; if(o==0) next; m=(a<o-a)?a:o-a; tot++
       if(m>=1){ poly++; if(m/o<=0.01){ rare++; c=$h["#CHROM"]; if(c!~/^chr/) c="chr"c
         print c"\t"($h["POS"]-1)"\t"$h["POS"]"\t"$h["ID"]"\t0\t+" > bed } } }
     END{print "chr22 SNVs:",tot,"| polymorphic in analysis set:",poly,"| MAF<=0.01 & MAC>=1:",rare}' "$SCR/ukb_chr22.acount"
else
  echo "SKIPPED (no plink2 or no keep file)"
fi

if [ -s "$SCR/c22_rare.hg38.bed" ] && [ -x "$LIFTOVER" ]; then
  t0=$(date +%s)
  "$LIFTOVER" "$SCR/c22_rare.hg38.bed" "$CHAIN" "$SCR/c22_rare.hg19.bed" "$SCR/c22_rare.unmapped" 2> "$SCR/liftover.stderr"
  echo "liftOver exit $?, $(( $(date +%s) - t0 )) s"
  echo "mapped $(wc -l < "$SCR/c22_rare.hg19.bed") | unmapped $(grep -vc '^#' "$SCR/c22_rare.unmapped") | minus strand $(awk '$6=="-"' "$SCR/c22_rare.hg19.bed" | wc -l) | other contig $(awk '$1!="chr22"' "$SCR/c22_rare.hg19.bed" | wc -l)"

  if [ -n "$BEDTOOLS" ] && [ -e "$HG19.fai" ]; then
    awk 'BEGIN{OFS="\t"} $6=="+" && $3-$2==1 {print $1,$2,$3,$4}' "$SCR/c22_rare.hg19.bed" | head -5000 > "$SCR/c22_sub.bed"
    "$BEDTOOLS" getfasta -fi "$HG19" -bed "$SCR/c22_sub.bed" -name -tab 2>/dev/null |
      awk -F'\t' '{split($1,n,"::"); split(n[1],p,":"); s=toupper($2); if(s==p[3]) r++; else if(s==p[4]) sw++; else nn++}
                  END{print "REF check (first 5000 + strand): hg19 base == REF",r+0,"| == ALT",sw+0,"| neither",nn+0}'
  fi
  if [ -n "$BEDTOOLS" ] && [ -e "$HCC" ]; then
    awk 'BEGIN{OFS="\t"} {if($1!~/^chr/) $1="chr"$1; print $1,$2,$3}' "$HCC" > "$SCR/hcc_chr.bed"
    nh=$("$BEDTOOLS" intersect -u -a "$SCR/c22_rare.hg19.bed" -b "$SCR/hcc_chr.bed" | wc -l)
    fac=$(awk -F'\t' '{bp+=$3-$2; if($1=="chr22"||$1=="22") b22+=$3-$2} END{printf "%.1f", bp/b22}' "$HCC")
    echo "chr22 rare SNVs in hcc exons: $nh | extrapolation factor (total/chr22 hcc bp): $fac | est. genome-wide ANNOVAR input: $(awk -v a="$nh" -v f="$fac" 'BEGIN{printf "%d", a*f}')"
  fi
else
  echo "liftOver test SKIPPED"
fi

# =============================================================================
sec "7. Resources"
df -h "$HOME" /mnt/hgfs/D /cluster2 /tmp 2>/dev/null
free -g | head -2
echo "cores: $(nproc)"

sec "DONE – send back $LOG"
