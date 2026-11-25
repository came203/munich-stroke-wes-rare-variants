# Munich stroke WES – rare-variant analysis pipeline (DFG project 497256604)

Analysis code for the gene-based rare-variant analysis of 1,535 PRS-selected ischemic stroke exomes and the evaluation of a
rare-variant risk score (RVGRS) in UK Biobank. Version 1.0.0 (2026-11-25). Code only, no data. Paths are those of the analysis
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
