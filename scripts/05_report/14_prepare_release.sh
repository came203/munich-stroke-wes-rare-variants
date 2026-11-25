#!/usr/bin/env bash
# =============================================================================
# 14_prepare_release.sh – build the two public release packages (DFG final report, 4.2 outputs)  v2
#   A) /mnt/hgfs/D/WES/UPLOAD/zenodo_dataset_v2/  gene-level results + RVGRS weights + UKB evaluation (summary level only)
#   B) /mnt/hgfs/D/WES/UPLOAD/github_repo/        pipeline FASTQ -> gVCF -> joint call -> RV-EXCALIBER -> RVGRS (scripts only), NOT pushed
#   FASTQ->gVCF and GLnexus scripts come from the LRZ bundle (14a_collect_lrz_pipeline.sh), extracted to
#   /mnt/hgfs/D/WES/UPLOAD/lrz_pipeline_bundle/ – skipped with a warning if absent
# Includes a privacy scan (sample IDs from the .fam files must not appear anywhere in the release).
# Run on the VM:  bash /mnt/hgfs/D/WES/14_prepare_release.sh
# Nothing is uploaded. Upload/publish steps are printed at the end (GitHub push + Zenodo GitHub integration;
#   Zenodo dataset via web UI or 14b_zenodo_draft.sh).
# =============================================================================
set -u
export LC_ALL=C
BASE=/mnt/hgfs/D/WES
P=$BASE/FOR_REPORT/partC
RES=$BASE/05_rvexcaliber_output_v2_missense
RVX=/home/rainer/WES_analysis/04_software/RV-EXCALIBER-main
REL=$BASE/UPLOAD
LRZB=$REL/lrz_pipeline_bundle
Z=$REL/zenodo_dataset_v2
G=$REL/github_repo
VERSION=1.0.0
DATE=$(date +%Y-%m-%d)
rm -rf "$Z" "$G"   # generated folders: rebuilt from scratch on every run (clean git history)
mkdir -p "$Z" "$G"
LOG=$REL/14_prepare_release.log
exec > >(tee "$LOG") 2>&1
echo "14_prepare_release.sh  $(date)"

# =============================================================================
echo "=== A) Zenodo dataset"
P="$P" RES="$RES" Z="$Z" Rscript - <<'EOF' || exit 1
suppressPackageStartupMessages(library(data.table))
P <- Sys.getenv("P"); Z <- Sys.getenv("Z")
w <- fread(file.path(P, "C_weights_v2.tsv"))
gl <- w[, .(gene = Gene, case_alleles = OAC, expected_alleles = round(EXP, 4), OR = round(OR, 4), P, P_BH,
            in_allele_filtered_set = set_allelefilt, n_qualifying_variants = n_var_rvx)][order(P)]
fwrite(gl, file.path(Z, "gene_level_results_IS_v2.tsv"), sep = "\t")
wt <- w[, .(gene = Gene, beta_raw = round(beta_raw, 6), beta_shrunk = round(beta_shrunk, 6),
            set_all = set_all, set_nominal = set_nominal, set_allele_filtered = set_allelefilt,
            set_fdr = set_fdr, set_ERCC2_PLG = set_top2)][order(gene)]
fwrite(wt, file.path(Z, "rvgrs_gene_weights_v2.tsv"), sep = "\t")
r <- fread(file.path(P, "C_results.tsv"))
fwrite(r[, .(model, term, type, n, cases, estimate = signif(est, 5), lci = signif(lci, 5), uci = signif(uci, 5), P = signif(p, 4), note)],
       file.path(Z, "rvgrs_ukb_association.tsv"), sep = "\t")
fwrite(fread(file.path(P, "C_m5_discrimination.tsv")), file.path(Z, "rvgrs_ukb_discrimination.tsv"), sep = "\t")
fwrite(fread(file.path(P, "C_m6_carriers.tsv")), file.path(Z, "ERCC2_PLG_ukb_carrier_test.tsv"), sep = "\t")
fwrite(fread(file.path(P, "C_deciles.tsv")), file.path(Z, "rvgrs_ukb_deciles.tsv"), sep = "\t")
fwrite(fread(file.path(P, "C_filter_rule.tsv"))[!param %chin% c("hcc_bed", "date", "D4_verdict", "match_on")],
       file.path(Z, "qualifying_variant_rule.tsv"), sep = "\t")
file.copy(file.path(P, "C_class_set.txt"), file.path(Z, "qualifying_variant_classes.txt"), overwrite = TRUE)
cat(sprintf("gene-level results %d genes | weights %d genes | UKB tables written\n", nrow(gl), nrow(wt)))
EOF
for f in rvexcaliber_testing_my_samples_plink_hg19_v2u_GenomicInflationFactor_allele_filter_0.001_0.025_nfe_hcc_iCFgCFadjust_rvexcaliber_base.txt; do
  [ -e "$RES/$f" ] && cp "$RES/$f" "$Z/genomic_inflation_allele_filter_v2.txt"
done
[ -e "$BASE/FOR_REPORT/figures/QQ_v2_main.png" ] && cp "$BASE/FOR_REPORT/figures/QQ_v2_main.png" "$Z/QQ_gene_level_v2.png"

cat > "$Z/README.md" <<EOF
# Rare-variant burden in PRS-selected ischemic stroke exomes (Munich) – gene-level results and rare-variant risk score weights

Version $VERSION ($DATE). Summary-level data only. No individual-level data are included (consent).

## Study
Whole-exome sequencing (Twist Exome 2.0, NCCT Tübingen) of ischemic stroke patients selected for low common-variant polygenic risk,
Institute for Stroke and Dementia Research (ISD), LMU University Hospital Munich. Analysis set: 1,535 unrelated patients of European ancestry.

## Methods (short)
| Step | Setting |
|---|---|
| Calling | DeepVariant 1.8.0 (WES model) per sample, GLnexus 1.4.1 joint calling (DeepVariantWES configuration with min_AQ1 = 20), GRCh38 |
| Liftover | multiallelic records split on GRCh38, UCSC liftOver to GRCh37, reference-allele check, duplicated records removed |
| Gene-based test | RV-EXCALIBER (Lali et al., Nat Commun 2021, https://github.com/GMELab/RV-EXCALIBER), base model, iCF and gCF adjustment, gnomAD v2.1.1 non-Finnish European exomes as reference, bundled MIGen ranking cohort |
| Qualifying variants | M-CAP > 0.025 (dbNSFP 3.5c via ANNOVAR 2019-10-24, refGene), gnomAD 2.1.1 exome AF < 0.001 in NFE, AFR, SAS, EAS and AMR, cohort AF < 0.01, high-coverage coding exons, autosomes. Exact rule: qualifying_variant_rule.tsv and qualifying_variant_classes.txt |
| Risk score (RVGRS) | weight per gene = ln(case alleles / expected alleles). Evaluated in UK Biobank WES (n = 438,694, 9,308 ischemic stroke cases), adjusted for total qualifying allele count, age, age squared, sex and 10 PCs (+ 2019 meta-PRS where stated) |

## Files
| File | Content |
|---|---|
| gene_level_results_IS_v2.tsv | 8,790 genes: case_alleles (observed qualifying alleles in 1,535 cases), expected_alleles (iCF/gCF-adjusted expectation from gnomAD), OR = case/expected, P (RV-EXCALIBER), P_BH (Benjamini-Hochberg), in_allele_filtered_set (gene passes the expected-allele filter, 3.79), n_qualifying_variants |
| rvgrs_gene_weights_v2.tsv | weights per gene: beta_raw = ln(OR) (primary), beta_shrunk = ln((case+0.5)/(expected+0.5)); gene-set flags (all, nominal P < 0.05 and OR > 1, allele-filtered, FDR < 0.05, ERCC2+PLG) |
| rvgrs_ukb_association.tsv | UK Biobank association models (OR or HR per SD, 95% CI, P) |
| rvgrs_ukb_discrimination.tsv | C-statistic with and without RVGRS (DeLong), continuous NRI, IDI |
| rvgrs_ukb_deciles.tsv | OR by RVGRS decile (reference decile 1) |
| ERCC2_PLG_ukb_carrier_test.tsv | within-UK Biobank case-control test of ERCC2 and PLG qualifying-variant carriers |
| qualifying_variant_rule.tsv, qualifying_variant_classes.txt | parameters of the qualifying-variant rule (reproduced exactly against RV-EXCALIBER) |
| genomic_inflation_allele_filter_v2.txt | RV-EXCALIBER genomic inflation output (allele-filtered run) |
| QQ_gene_level_v2.png | QQ plot of the gene-level test (lambda 0.77) |

## Key results
Gene-level test: genomic inflation 0.77 (conservative), exome-wide significant: ERCC2, PLG. The rare-variant risk score was not associated with ischemic stroke in UK Biobank (OR 0.993 per SD, 95% CI 0.972 to 1.014).

## Code
Analysis pipeline: GitHub repository (see related identifiers), archived on Zenodo.

## Funding
Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project number 497256604 (MA 7973/2-1, DI 722/20-1).
UK Biobank analyses under application 2532.

## License
CC BY 4.0
EOF
( cd "$Z" && md5sum $(ls | grep -v md5sums.txt) > md5sums.txt )
echo "Zenodo dataset files:"; ls -la "$Z"

# =============================================================================
echo "=== B) GitHub repository"
mkdir -p "$G"/{scripts/00_fastq_to_gvcf,scripts/01_calling_qc,scripts/02_rvexcaliber,scripts/03_followup,scripts/04_rvgrs_ukb,scripts/05_report,rvexcaliber_patches/shell}
cp_if(){ for f in "$@"; do if [ -e "$f" ]; then cp "$f" "$DEST/"; echo "  + $(basename "$f")"; else echo "  - missing: $f"; fi; done; }
for t in "$LRZB"/*.tar.gz; do [ -e "$t" ] && { tar -xzf "$t" -C "$LRZB" && echo "  extracted $(basename "$t")"; }; done
echo "  bundle content:"; ls -la "$LRZB" 2>/dev/null | sed "s#^#    #"
if [ -d "$LRZB" ] && [ -n "$(ls -A "$LRZB" 2>/dev/null)" ]; then
  echo "  LRZ pipeline bundle -> scripts/00_fastq_to_gvcf"
  ( cd "$LRZB" && find . -type f ! -name '*.txt' ! -name '*.tsv' ! -name '*.csv' ! -name '*.list' ! -name '*.tar*' -size -1024k ) |
  while read -r f; do
    case "$(basename "$f")" in 02*|03_*) d=01_calling_qc ;; *) d=00_fastq_to_gvcf ;; esac
    mkdir -p "$G/scripts/$d/$(dirname "$f")"; cp "$LRZB/$f" "$G/scripts/$d/$f"; echo "    + $d/${f#./}"
  done
else
  echo "  !! no LRZ bundle in $LRZB – run 14a_collect_lrz_pipeline.sh on LRZ, copy + extract the tar there, rerun 14"
fi
DEST=$G/scripts/01_calling_qc; cp_if "$BASE"/04_rebuild_hg19_v2.sh
for f in 02_diag_lrz_jointvcf.sh 03_glnexus_rerun_AQ20.sh; do [ -e "$DEST/$f" ] || echo "  !! $f not found (expected from the LRZ bundle)"; done
DEST=$G/scripts/02_rvexcaliber; cp_if "$BASE"/05_run_rvexcaliber_v2.sh "$BASE"/05a_diag_rvx_missense.sh "$BASE"/05b_diag_resume.sh "$BASE"/05c_resume_rvexcaliber_v2.sh
DEST=$G/scripts/03_followup; cp_if "$BASE"/06_partA_v2.R "$BASE"/07_syn_negative_control.R "$BASE"/08a_diag_partB_subtypes.sh "$BASE"/08_partB_subtypes.sh "$BASE"/08b_partB_followup.R
DEST=$G/scripts/04_rvgrs_ukb; cp_if "$BASE"/09a_diag_partC.sh "$BASE"/partC_qualify_fun.R "$BASE"/09_weights.R "$BASE"/10_ukb_qualify.sh "$BASE"/10b_ukb_qualify.R "$BASE"/11_ukb_score.sh "$BASE"/11b_diag_ukb_vs_munich.R "$BASE"/12_ukb_stats.R
DEST=$G/scripts/05_report; cp_if "$BASE"/13a_qqplot.R "$BASE"/13b_sgp_enrichment.R "$BASE"/13c_csvd_carriers.R "$BASE"/14_prepare_release.sh
DEST=$G/rvexcaliber_patches; cp_if "$RVX"/scripts/Rscripts/get_Varlist_rvexcaliber_missense_fixed.R
DEST=$G/rvexcaliber_patches/shell; cp_if "$HOME"/rvx_v2_missense/shell_fixed/*.sh

cat > "$G/README.md" <<EOF
# Munich stroke WES – rare-variant analysis pipeline (DFG project 497256604)

Analysis code for the gene-based rare-variant analysis of 1,535 PRS-selected ischemic stroke exomes and the evaluation of a
rare-variant risk score (RVGRS) in UK Biobank. Version $VERSION ($DATE). Code only, no data. Paths are those of the analysis
environment (LRZ cluster, ISD VM) and must be adapted.

| Folder | Step | Scripts |
|---|---|---|
| scripts/00_fastq_to_gvcf | LRZ (SLURM), run in this order: (1) script_exis_minimal_lrz.sh = FASTQ to analysis-ready BAM (ngs-bits read QC, trimming and merging of 4 lanes, BWA-MEM2, samblaster duplicate marking, samtools sort/index, ABRA2 2.19 indel realignment on Twist Exome 2.0 targets), (2) script_deepvariant.sh = DeepVariant 1.8.0 WES model, per-sample gVCF, (3) joint calling: scripts/01_calling_qc/03_glnexus_rerun_AQ20.sh (final call set). script_nexus2.sh = original GLnexus call (DeepVariantWES preset, min_AQ1 35), superseded because it removed all singletons | script_exis_minimal_lrz.sh, script_deepvariant.sh, script_nexus2.sh |
| scripts/01_calling_qc | joint-calling diagnostics, GLnexus rerun (min_AQ1 20), GRCh38 split, liftOver to GRCh37, REF check, deduplication, KING relatedness | 02*, 03, 04 |
| scripts/02_rvexcaliber | RV-EXCALIBER run (missense M-CAP > 0.025, gnomAD NFE, hcc, iCF/gCF) | 05* |
| scripts/03_followup | variant breakdown of top genes, synonymous negative control, TOAST subtype analyses | 06, 07, 08* |
| scripts/04_rvgrs_ukb | gene weights and exact reproduction of the RV-EXCALIBER qualifying rule (09), UK Biobank variant qualification with identical ANNOVAR annotation (10, 10b), scoring (11), UKB vs gnomAD reference check (11b), association models (12). partC_qualify_fun.R = shared rule | 09a to 12 |
| scripts/05_report | QQ plot, monogenic stroke gene panel enrichment, cSVD gene carrier lookup, release packaging | 13a to 14 |
| rvexcaliber_patches | modified RV-EXCALIBER files used in this analysis (variant filter reduced to M-CAP > 0.025, fixed shell copies) | – |

Dependencies: ngs-bits, BWA-MEM2, samblaster, samtools, ABRA2 2.19, RV-EXCALIBER (https://github.com/GMELab/RV-EXCALIBER), ANNOVAR 2019-10-24 with refGene, gnomad211_exome, dbnsfp35c (hg19),
GLnexus 1.4.1, DeepVariant 1.8.0, bcftools 1.19, PLINK 1.9 and 2.0, bedtools 2.31, UCSC liftOver, R 4.6 (data.table, survival, pROC, ggplot2, readxl).

Results: Zenodo dataset (see related identifiers).

## Funding
Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project number 497256604.

## License
GPL-3.0 (the repository contains modified RV-EXCALIBER files, which are GPL-3.0).
EOF
cat > "$G/LICENSE" <<'EOF'
This repository is licensed under the GNU General Public License v3.0.
Full text: https://www.gnu.org/licenses/gpl-3.0.txt
(Replace this file with the full GPL-3.0 text before publishing, e.g. via GitHub "Add license" -> GPL-3.0.)
EOF
cat > "$G/CITATION.cff" <<EOF
cff-version: 1.2.0
message: "If you use this code, please cite it as below."
title: "Munich stroke WES rare-variant analysis pipeline"
version: "$VERSION"
date-released: "$DATE"
authors:
  - family-names: "Malik"
    given-names: "Rainer"
    affiliation: "Institute for Stroke and Dementia Research (ISD), LMU University Hospital, LMU Munich"
  - family-names: "Dichgans"
    given-names: "Martin"
    affiliation: "Institute for Stroke and Dementia Research (ISD), LMU University Hospital, LMU Munich"
license: GPL-3.0
EOF
cat > "$G/.zenodo.json" <<EOF
{
  "title": "Munich stroke WES rare-variant analysis pipeline",
  "upload_type": "software",
  "description": "Analysis code for the gene-based rare-variant analysis (RV-EXCALIBER) of 1,535 PRS-selected ischemic stroke exomes and the evaluation of a rare-variant risk score in UK Biobank. Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project number 497256604.",
  "creators": [
    {"name": "Malik, Rainer", "affiliation": "Institute for Stroke and Dementia Research (ISD), LMU University Hospital, LMU Munich"},
    {"name": "Dichgans, Martin", "affiliation": "Institute for Stroke and Dementia Research (ISD), LMU University Hospital, LMU Munich"}
  ],
  "license": "GPL-3.0-only",
  "keywords": ["ischemic stroke", "whole-exome sequencing", "rare variants", "RV-EXCALIBER", "UK Biobank", "polygenic risk"],
  "access_right": "open"
}
EOF
printf '*.tsv\n*.tsv.gz\n*.raw\n*.fam\n*.bim\n*.bed\n*.pgen\n*.vcf*\n*.log\n' > "$G/.gitignore"

echo "=== B2) remove hard-coded sample IDs from the repository copies (originals unchanged)"
FAMS="$BASE/FOR_REPORT/v2/plink_in/my_samples_plink_hg19_v2u.fam $BASE/FOR_REPORT/v2/my_samples_plink_hg19_v2.fam"
cat $FAMS 2>/dev/null | awk '{print $2}' | sort -u | awk 'length($0)>=5' > "$REL/.ids_tmp"
cat > "$REL/.sanitise.pl" <<'PL'
s/^REMOVE="[0-9 ]+"/REMOVE="\$(tr '\\n' ' ' < "\${REMOVE_FILE:?set REMOVE_FILE (one IID per line, one sample per related pair)}")"/;
s/\((?:[0-9]+, )+[0-9]+; higher missingness/(IDs from REMOVE_FILE; higher missingness/;
s/printf '(?:[0-9]+ [0-9]+\\n)+' > "\$TMP\/rel.txt"/cp "\${REL_PAIRS_FILE:?set REL_PAIRS_FILE (related pairs, IID IID per line)}" "\$TMP\/rel.txt"/;
PL
for f in "$G/scripts/02_rvexcaliber/05_run_rvexcaliber_v2.sh" "$G/scripts/02_rvexcaliber/05a_diag_rvx_missense.sh"; do
  [ -e "$f" ] && perl -pi "$REL/.sanitise.pl" "$f" && echo "  sanitised $(basename "$f")"
done
grep -n 'REMOVE=\|REMOVE_FILE\|REL_PAIRS_FILE' "$G/scripts/02_rvexcaliber/"05*.sh | sed 's#^#    #'
rm -f "$REL/.sanitise.pl"
echo "  generic: e-mail addresses -> <your.email>, standalone sample-ID tokens -> <SAMPLE_ID> (all files under scripts/ and rvexcaliber_patches/)"
find "$G/scripts" "$G/rvexcaliber_patches" -type f ! -name 14_prepare_release.sh | while read -r f; do
  IDS="$REL/.ids_tmp" perl -pi -e 'BEGIN{ open(my $F, "<", $ENV{IDS}) or die; while(<$F>){ chomp; $id{$_}=1 } }
    s/[A-Za-z0-9._%+-]+\@[A-Za-z0-9.-]+\.[A-Za-z]{2,}/<your.email>/g;
    BEGIN{ $re = join("|", map { quotemeta } sort { length($b) <=> length($a) } keys %id) }
    s/(?<![A-Za-z0-9_])(?<!\d\.)($re)(?![A-Za-z0-9_])(?!\.\d)/<SAMPLE_ID>/g;' "$f"
done
grep -rn '<SAMPLE_ID>\|<your.email>' "$G/scripts" "$G/rvexcaliber_patches" | grep -v 14_prepare_release.sh | sed 's#^#    #'

# =============================================================================
echo "=== C) privacy scan (sample IDs, individual-level files)"
echo "sample IDs checked: $(wc -l < "$REL/.ids_tmp")"
hits=$(find "$Z" "$G" -type f ! -path '*/.git/*' ! -name '*.png' ! -name '*.pdf' | while read -r f; do
  IDS="$REL/.ids_tmp" F="$f" perl -ne 'BEGIN{ open(my $H, "<", $ENV{IDS}) or die; while(<$H>){ chomp; $id{$_}=1 }
      $re = join("|", map { quotemeta } sort { length($b) <=> length($a) } keys %id) }
    print "$ENV{F}:$.: $1\n" if /(?<![A-Za-z0-9_])(?<!\d\.)($re)(?![A-Za-z0-9_])(?!\.\d)/;' "$f"
done | head -50)
rm -f "$REL/.ids_tmp"
if [ -n "$hits" ]; then
  echo "!!! POSSIBLE SAMPLE IDs IN RELEASE (review/remove before publishing; numbers may also be false positives):"; echo "$hits"
else
  echo "OK: no sample ID from the .fam files found in the release"
fi
mails=$(grep -rnoE '[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}' "$Z" "$G" 2>/dev/null | grep -v '/\.git/' | grep -v 'example.org' | grep -viE 'noreply|git@github.com')
[ -n "$mails" ] && { echo "!!! e-mail addresses in release (e.g. SBATCH --mail-user) – remove before publishing:"; echo "$mails"; } || echo "OK: no e-mail addresses"
secr=$(grep -rniE '(password|passwd|token|secret|api[_-]?key)[[:space:]]*[=:]' "$G" 2>/dev/null | grep -v '/\.git/' | grep -v 'ZENODO_TOKEN')
[ -n "$secr" ] && { echo "!!! possible credentials:"; echo "$secr"; } || echo "OK: no credentials pattern"
big=$(find "$Z" "$G" -type f -size +5M 2>/dev/null); [ -n "$big" ] && echo "large files (check): $big"
indiv=$(find "$Z" "$G" -type f \( -iname '*INDIVIDUAL*' -o -name '*.raw' -o -name '*scores*' \) 2>/dev/null)
[ -n "$indiv" ] && echo "!!! individual-level files present: $indiv" || echo "OK: no individual-level files"

# =============================================================================
echo "note: git on the VMware share (hgfs) can fail (chmod/locks) – if so, copy $G to ~ and run git there"
if command -v git >/dev/null; then
  gitc(){ git -c safe.directory="$G" "$@"; }
  ( cd "$G" && { [ -d .git ] || gitc init -q -b main; } && gitc add -A && { gitc config user.email >/dev/null || gitc config user.email "set-me@example.org"; } && { gitc config user.name >/dev/null || gitc config user.name "Rainer Malik"; } && gitc commit -q -m "Pipeline v$VERSION (DFG final report)" && echo "git: committed locally ($(gitc rev-parse --short HEAD))" ) || echo "git commit skipped (nothing new or error) – alternative: git config --global --add safe.directory $G"
fi

cat <<EOF

=== NEXT STEPS (manual; nothing has been uploaded)
GitHub + code DOI
  1 zenodo.org -> GitHub settings -> switch ON the new repository (do this BEFORE the release)
  2 create empty repo on github.com (e.g. munich-stroke-wes-rare-variants), replace LICENSE with full GPL-3.0 text
  3 cd $G && git remote add origin git@github.com:<user>/munich-stroke-wes-rare-variants.git && git push -u origin main
  4 GitHub -> Releases -> "v$VERSION" -> Zenodo archives it and mints the code DOI (.zenodo.json supplies the metadata)
Zenodo dataset
  5 zenodo.org -> New upload -> files from $Z (or: ZENODO_TOKEN=... bash 14b_zenodo_draft.sh creates a draft, not published)
  6 type Dataset, license CC BY 4.0, related identifier = code DOI / GitHub URL, publish -> dataset DOI
  7 enter both DOIs in the report (Section 4.2) and in claude/Abschlussbericht_Plan.md
Check the git author e-mail of the local commit (git log -1) before pushing; fix with git commit --amend --reset-author if needed.
log: $LOG
EOF
