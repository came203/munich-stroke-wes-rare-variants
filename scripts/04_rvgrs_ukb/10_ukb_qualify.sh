#!/usr/bin/env bash
# =============================================================================
# 10_ukb_qualify.sh – Part C step 2: UKB variants -> hg19 -> ANNOVAR (identical to Munich)
# RVGRS_Plan.md 4b, D3/D4 (rule reproduced exactly on Munich by 09_weights.R, 25.11.2026)
#
# Usage (VM, scripts in /mnt/hgfs/D/WES/):
#   bash 10_ukb_qualify.sh validate        # 1) ANNOVAR call reproduces Munich chr22 annotation -> must PASS
#   bash 10_ukb_qualify.sh chr 22          # 2) test one chromosome, note runtime
#   nohup bash 10_ukb_qualify.sh all 3 > ~/partC_ukb/10_all.out 2>&1 &   # 3) chr1-22, 3 in parallel
#   Rscript 10b_ukb_qualify.R              # 4) apply the D4 rule -> qualifying variants
#   (finished chromosomes are skipped on rerun; FORCE=1 skips the validate gate;
#    ANNOVAR_EXTRA="..." adds table_annovar arguments if validate shows a difference)
#
# Per chromosome (work dir ~/partC_ukb/chrN, local disk – no gzip on hgfs):
#   1 plink2 --freq on the analysis set (n 438,694; keep file from 09a), SNVs only
#   2 AC >= 1 and alt AF <= 0.01 (pre-filter; exact rule '<' applied in 10b)
#   3 UCSC liftOver hg38 -> hg19; keep + strand, span 1, chr1-22 (as Munich script 04)
#   4 hcc pre-filter in 'rvx' mode (BED s <= pos <= e) via interval [pos-1, pos+1)
#   5 hg19 REF check (drop mismatches); 6 drop all copies of duplicated hg19 keys
#   7 ANNOVAR A: refGene only -> keep classes in C_class_set.txt
#   8 ANNOVAR B: refGene,gnomad211_exome,dbnsfp35c (same layout as Munich multianno)
# Outputs: chrN/map.tsv (ID38 -> hg19, ALT_CTS, OBS_CT), chrN/B.hg19_multianno.txt.gz,
#          chrN/summary.tsv, chrN/DONE; logs copied to FOR_REPORT/partC/
# =============================================================================
set -u
export LC_ALL=C

BASE=/mnt/hgfs/D/WES
OUT=$BASE/FOR_REPORT/partC
W=$HOME/partC_ukb
SELF=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SELF")
FUN=$SCRIPT_DIR/partC_qualify_fun.R

UKB=/cluster2/dichganslab/rmalik2/UKB_FINAL
KEEP_SRC=$HOME/partC_diag/keep_analysis.txt
KEEP=$W/keep_analysis.txt
RVX=/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main
HCC=$RVX/gnomAD_211_exomes_hg19_coverage/ALL_CHROM_gnomad_refGene_exons_highcov20X.mod.bed
HCC_CHR=$W/hcc_chr.bed
ANNOVAR=$BASE/04_software/annovar
ANNOVAR_EXTRA=${ANNOVAR_EXTRA:-}
CHAIN=$BASE/03_references/hg38ToHg19.over.chain.gz
HG19=$BASE/03_references/hg19.fa
LIFTOVER=/mnt/hgfs/D/software/liftOver
PLINK2=${PLINK2:-$(command -v plink2)}
BEDTOOLS=$(command -v bedtools)
THREADS=${THREADS:-4}
MEM=${MEM:-3000}
INT_THR=0.01
CLASS_SET=$OUT/C_class_set.txt
RULE=$OUT/C_filter_rule.tsv
MUNICH_AVI=$HOME/rvx_v2_missense/my_samples_plink_hg19_v2u.annovarInput
MUNICH_MA=$HOME/rvx_v2_missense/my_samples_plink_hg19_v2u.hg19_multianno.txt.gz

die(){ echo "ERROR: $*"; exit 1; }
nl_(){ [ -e "$1" ] && awk 'END{print NR}' "$1" || echo 0; }

usage(){ sed -n '6,14p' "$SELF"; }

setup(){
  mkdir -p "$W"
  for t in "$PLINK2" "$BEDTOOLS" "$LIFTOVER"; do [ -x "$t" ] || die "tool missing: $t"; done
  for f in "$CHAIN" "$HG19" "$HG19.fai" "$HCC" "$CLASS_SET" "$RULE" "$FUN" "$ANNOVAR/table_annovar.pl"; do
    [ -e "$f" ] || die "missing: $f"; done
  [ "$(awk -F'\t' '$1=="D4_verdict"{print $2}' "$RULE")" = "PASS" ] || die "C_filter_rule.tsv: D4_verdict is not PASS"
  if [ ! -s "$KEEP" ]; then [ -s "$KEEP_SRC" ] || die "keep file missing (run 09a)"; cp "$KEEP_SRC" "$KEEP"; fi
  [ -s "$HCC_CHR" ] || awk 'BEGIN{OFS="\t"} {if($1!~/^chr/) $1="chr"$1; print $1,$2,$3}' "$HCC" > "$HCC_CHR"
}

# ANNOVAR: $1 avinput, $2 out prefix, $3 protocols, $4 operations
annovar(){
  perl "$ANNOVAR/table_annovar.pl" "$1" "$ANNOVAR/humandb" -buildver hg19 -out "$2" -remove \
       -protocol "$3" -operation "$4" -nastring . $ANNOVAR_EXTRA > "$2.annovar.log" 2>&1
  [ -s "$2.hg19_multianno.txt" ] || { tail -20 "$2.annovar.log"; die "ANNOVAR failed: $2"; }
}

# keep avinput rows whose class "Func.refGene|ExonicFunc.refGene" (spaces -> _, empty -> .) is in the class set
class_filter(){
  awk -F'\t' -v OFS='\t' 'NR==FNR{keep[$1]=1; next} FNR==1{next}
     { f=$6; e=$9; gsub(/ /,"_",e); if(f=="") f="."; if(e=="") e="."; if((f"|"e) in keep) print $1,$2,$3,$4,$5 }' \
     "$CLASS_SET" "$1" > "$2"
}

annotate_two_stage(){   # $1 avinput, $2 dir
  annovar "$1" "$2/A" refGene g
  class_filter "$2/A.hg19_multianno.txt" "$2/B.avinput"
  annovar "$2/B.avinput" "$2/B" refGene,gnomad211_exome,dbnsfp35c g,f,f
}

# -----------------------------------------------------------------------------
validate(){
  V=$W/validate; mkdir -p "$V"
  exec > >(tee "$V/validate.log") 2>&1
  echo "10 validate  $(date)"
  echo "--- ANNOVAR calls inside RV-EXCALIBER (for reference):"
  grep -hn "table_annovar\|annotate_variation" "$RVX"/scripts/shell/*.sh "$HOME"/rvx_v2_missense/shell_fixed/*.sh 2>/dev/null | head -10
  echo "--- this script: table_annovar.pl <in> humandb -buildver hg19 -remove -protocol refGene[,gnomad211_exome,dbnsfp35c] -operation g[,f,f] -nastring . $ANNOVAR_EXTRA"
  [ -s "$MUNICH_AVI" ] || die "Munich annovarInput missing: $MUNICH_AVI"
  echo "Munich annovarInput head:"; head -2 "$MUNICH_AVI"
  awk -v OFS='\t' '($1=="22"||$1=="chr22"){print $1,$2,$3,$4,$5}' "$MUNICH_AVI" > "$V/m22.avinput"
  echo "Munich chr22 avinput rows: $(nl_ "$V/m22.avinput")"
  t0=$(date +%s); annotate_two_stage "$V/m22.avinput" "$V"; echo "ANNOVAR two-stage: $(( $(date +%s)-t0 )) s"
  FUN="$FUN" MA="$MUNICH_MA" B="$V/B.hg19_multianno.txt" CS="$CLASS_SET" VD="$V/VERDICT" Rscript - <<'EOF'
suppressPackageStartupMessages(library(data.table)); source(Sys.getenv("FUN"))
m <- read_multianno(Sys.getenv("MA"))[chr == "22"]
b <- read_multianno(Sys.getenv("B"))
cs <- readLines(Sys.getenv("CS"))
m[, `:=`(key = make_key(chr, pos, ref, alt), class = make_class(func, exfunc))]
b[, key := make_key(chr, pos, ref, alt)]
m_in <- m[class %chin% cs]
x <- merge(m_in, b, by = "key", suffixes = c(".m", ".u"))
cat(sprintf("Munich chr22 rows %d | in class set %d | UKB-path stage B rows %d | matched %d | missing in B %d | extra in B %d\n",
            nrow(m), nrow(m_in), nrow(b), nrow(x), nrow(m_in) - nrow(x), nrow(b) - nrow(x)))
eq <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
cols <- c("func", "gene", "exfunc", POPS, "mcap")
d <- sapply(cols, function(cl) sum(!eq(x[[paste0(cl, ".m")]], x[[paste0(cl, ".u")]])))
cat("differing values per column:\n"); print(d)
if (any(d > 0)) { cl <- names(d)[d > 0][1]; i <- which(!eq(x[[paste0(cl, ".m")]], x[[paste0(cl, ".u")]]))[1:5]
  print(x[i, c("key", paste0(cl, c(".m", ".u"))), with = FALSE]) }
ok <- all(d == 0) && nrow(x) == nrow(m_in) && nrow(b) == nrow(m_in)
writeLines(if (ok) "PASS" else "FAIL", Sys.getenv("VD"))
cat("VALIDATE:", if (ok) "PASS – ANNOVAR path reproduces Munich annotation" else "FAIL", "\n")
EOF
  cp "$V/validate.log" "$OUT/10_validate.log"
}

need_pass(){
  [ "${FORCE:-0}" = 1 ] && return 0
  [ "$(cat "$W/validate/VERDICT" 2>/dev/null)" = "PASS" ] || die "run 'bash 10_ukb_qualify.sh validate' first (must PASS)"
}

# -----------------------------------------------------------------------------
do_chr(){
  c=$1; D=$W/chr$c; mkdir -p "$D"
  [ -e "$D/DONE" ] && { echo "chr$c already done – skip"; return 0; }
  exec > >(tee "$D/10_chr$c.log") 2>&1
  echo "chr$c start $(date)"; t0=$(date +%s)
  PF=""; [ -e "$UKB/chr$c.pvar" ] || PF=vzs

  # 1 freq
  "$PLINK2" --pfile "$UKB/chr$c" $PF --keep "$KEEP" --snps-only just-acgt --freq counts cols=+pos \
            --threads "$THREADS" --memory "$MEM" --out "$D/freq" > "$D/freq.stdout" 2>&1 || die "plink2 chr$c"

  # 2 rare pre-filter -> hg38 BED (+ strand) and counts
  awk -F'\t' -v OFS='\t' -v thr="$INT_THR" -v bed="$D/rare.hg38.bed" -v cts="$D/rare.counts" '
      NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
      { a=$h["ALT_CTS"]; o=$h["OBS_CT"]
        if(o>0 && a>=1 && a/o<=thr){ c=$h["#CHROM"]; if(c!~/^chr/) c="chr"c
          print c,$h["POS"]-1,$h["POS"],$h["ID"],0,"+" > bed; print $h["ID"],a,o > cts } }' "$D/freq.acount"

  # 3 liftOver; 4 + strand, span 1, autosomes; hcc rvx pre-filter with [pos-1, pos+1)
  "$LIFTOVER" "$D/rare.hg38.bed" "$CHAIN" "$D/rare.hg19.bed" "$D/rare.unmapped" 2> "$D/liftover.stderr" || die "liftOver chr$c"
  awk -v OFS='\t' '$6=="+" && $3-$2==1 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print $1,$2,$3+1,$4}' \
      "$D/rare.hg19.bed" > "$D/lifted.wide.bed"
  "$BEDTOOLS" intersect -u -a "$D/lifted.wide.bed" -b "$HCC_CHR" > "$D/hcc.wide.bed"

  # 5 REF check on hg19
  awk -v OFS='\t' '{print $1,$2,$3-1,$4}' "$D/hcc.wide.bed" > "$D/hcc.bed"
  "$BEDTOOLS" getfasta -fi "$HG19" -bed "$D/hcc.bed" -name -tab > "$D/hcc.fa.tab" 2> "$D/getfasta.stderr"
  awk -F'\t' -v OFS='\t' 'NR==FNR{split($1,n,"::"); base[n[1]]=toupper($2); next}
      { split($4,p,":"); c=$1; sub(/^chr/,"",c)
        if(base[$4]==p[3]) print $4,c,$3,p[3],p[4]; else if(base[$4]==p[4]) sw++; else nn++ }
      END{print "REF check: REF==ALT(swap) " sw+0 ", neither " nn+0 > "/dev/stderr"}' \
      "$D/hcc.fa.tab" "$D/hcc.bed" > "$D/refok.tsv"

  # 6 drop all copies of duplicated hg19 keys; join counts -> map.tsv
  awk -F'\t' -v OFS='\t' 'NR==FNR{n[$2":"$3":"$4":"$5]++; next} n[$2":"$3":"$4":"$5]==1' \
      "$D/refok.tsv" "$D/refok.tsv" > "$D/uniq.tsv"
  awk -F'\t' -v OFS='\t' 'BEGIN{print "ID38","chr","pos","ref","alt","ALT_CTS","OBS_CT"}
      NR==FNR{row[$1]=$0; next} ($1 in row){print row[$1],$2,$3}' "$D/uniq.tsv" "$D/rare.counts" > "$D/map.tsv"

  # 7/8 ANNOVAR two-stage
  awk -v OFS='\t' 'NR>1{print $2,$3,$3,$4,$5}' "$D/map.tsv" > "$D/all.avinput"
  annotate_two_stage "$D/all.avinput" "$D"

  printf "chr\tsnvs\trare\tlifted\thcc\trefok\tuniq\tA_rows\tB_rows\tseconds\n%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" \
    "$c" $(( $(nl_ "$D/freq.acount")-1 )) "$(nl_ "$D/rare.counts")" "$(nl_ "$D/lifted.wide.bed")" "$(nl_ "$D/hcc.bed")" \
    "$(nl_ "$D/refok.tsv")" $(( $(nl_ "$D/map.tsv")-1 )) $(( $(nl_ "$D/A.hg19_multianno.txt")-1 )) \
    $(( $(nl_ "$D/B.hg19_multianno.txt")-1 )) $(( $(date +%s)-t0 )) > "$D/summary.tsv"
  cat "$D/summary.tsv"

  gzip -f "$D/B.hg19_multianno.txt"
  rm -f "$D/A.hg19_multianno.txt" "$D/rare.hg38.bed" "$D/rare.hg19.bed" "$D/lifted.wide.bed" "$D/hcc.wide.bed" \
        "$D/hcc.fa.tab" "$D/freq.acount" "$D/all.avinput" "$D/B.avinput"
  touch "$D/DONE"
  echo "chr$c done $(date)"
}

run_all(){
  PAR=$1; export THREADS=$(( 16 / PAR > 1 ? 16 / PAR : 1 ))
  order=$(ls -S "$UKB"/chr*.pvar* | sed -E 's#.*/chr([0-9]+)\.pvar.*#\1#' | awk '$1>=1 && $1<=22' | awk '!s[$1]++')
  echo "10 all: PAR=$PAR THREADS=$THREADS order: $(echo $order)  $(date)"
  printf '%s\n' $order | xargs -P "$PAR" -I{} bash "$SELF" chr {}
  echo "=== summary $(date)"
  ok=1; first=1
  for c in $(seq 1 22); do
    [ -e "$W/chr$c/DONE" ] || { echo "chr$c NOT DONE"; ok=0; continue; }
    if [ $first = 1 ]; then cat "$W/chr$c/summary.tsv"; first=0; else tail -n +2 "$W/chr$c/summary.tsv"; fi
  done | tee "$OUT/C_ukb_prefilter_summary.tsv"
  cat "$W"/chr*/10_chr*.log > "$OUT/10_ukb_qualify_chr.log" 2>/dev/null
  echo "next: Rscript $SCRIPT_DIR/10b_ukb_qualify.R"
}

# -----------------------------------------------------------------------------
case "${1:-}" in
  validate) setup; validate ;;
  chr)      [ -n "${2:-}" ] || { usage; exit 1; }; setup; need_pass; do_chr "$2" ;;
  all)      setup; need_pass; run_all "${2:-3}" ;;
  *)        usage; exit 1 ;;
esac
