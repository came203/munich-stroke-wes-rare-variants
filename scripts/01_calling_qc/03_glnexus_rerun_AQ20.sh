#!/bin/bash
# =====================================================================
# 03_glnexus_rerun_AQ20.sh  –  GLnexus joint calling rerun with singleton-preserving config (LRZ, SLURM)
# Round 2, DFG MA 7973/2-1 / DI 722/20-1 · after 02d (24.11.2026):
#   DeepVariant 1.8.0 gVCFs: all het calls GQ < 25 (median 22) → DeepVariantWES preset (min_AQ1 35 for alleles seen in
#   one sample, min_AQ2 20 for ≥2 samples) dropped every singleton. Fix: identical config, only min_AQ1 35 → 20.
# Inputs (unchanged from original run, script_nexus2.sh): mergelist_glnexus_filtered.txt (1,538 gVCFs),
#   twistCustomExomeV2_2021_12_14.bed. Original outputs (joint.bcf / joint.vcf.gz / GLnexus.DB) are NOT touched.
# Outputs (DSS project area – home is full, scratch not visible from login/Globus): $OUT/joint_AQ20.bcf, $OUT/joint_AQ20.vcf.gz(+.csi), $OUT/DeepVariantWES_AQ20.yml,
#   $OUT/joint_AQ20_QC.txt   (DB in $OUT/GLnexus_AQ20.DB – delete after success)
# Submit:  sbatch 03_glnexus_rerun_AQ20.sh      (cm4_tiny, 32 cores, 96G; original run on teramem: ~1.5 h)
# Status:  squeue -M cm4 -u $USER
# =====================================================================
#SBATCH -J join_AQ20
#SBATCH --output=/dss/dsshome1/06/ru36ros/LOGS/R-%x.%j.out
#SBATCH --error=/dss/dsshome1/06/ru36ros/LOGS/R-%x.%j.err
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<your.email>
#SBATCH --export=NONE
#SBATCH --time=5:00:00
set -uo pipefail
module load charliecloud/0.30 2>/dev/null || true

HOMEDIR=/dss/dsshome1/06/ru36ros
OUT=${OUT:-/dss/dssfs03/pn52re/pn52re-dss-0001/EXOME/v2}   # DSS project area (Rainer, 24.11.2026)
BCF=${HOMEDIR}/bcftools-1.19/bcftools
GLX=${HOMEDIR}/glnexus_cli
LIST=${HOMEDIR}/mergelist_glnexus_filtered.txt
BED=${HOMEDIR}/twistCustomExomeV2_2021_12_14.bed
ORIG=${HOMEDIR}/joint.vcf.gz
THREADS=${SLURM_CPUS_PER_TASK:-32}
CFG=${OUT}/DeepVariantWES_AQ20.yml
DB=${OUT}/GLnexus_AQ20.DB
QC=${OUT}/joint_AQ20_QC.txt
die(){ echo "ERROR: $*" >&2; exit 1; }

echo "== start $(date) | host $(hostname) | threads ${THREADS}"
for f in "$BCF" "$GLX"; do [ -x "$f" ] || die "$f not executable"; done
for f in "$LIST" "$BED" "$ORIG"; do [ -s "$f" ] || die "$f missing"; done
mkdir -p "$OUT" || die "cannot create $OUT"
echo "== free space at $(dirname "$OUT"):"; df -h "$(dirname "$OUT")" | tail -1   # needs ~15 GB (DB ~10 GB, BCF + VCF ~4 GB)
[ -e "$DB" ] && die "$DB exists – remove it or set OUT to a new directory (GLnexus needs a fresh DB dir)"
[ "$(grep -c . "$LIST")" -eq 1538 ] || die "merge list does not have 1538 entries"

## 1 config: exact preset of the original run (from its VCF header), min_AQ1 35 -> 20
"$BCF" view -h "$ORIG" | grep -m1 '^##GLnexusConfig=' | sed 's/^##GLnexusConfig=//' > "$CFG"
grep -q 'min_AQ1: 35' "$CFG" || die "min_AQ1: 35 not found in original config"
sed -i 's/min_AQ1: 35/min_AQ1: 20/' "$CFG"
echo "== config ${CFG}:"; grep -o 'unifier_config: {[^}]*}' "$CFG"

## 2 pre-check: GQ distribution of PASS het calls in the first gVCF (documents the GQ scale)
g1=$(head -1 "$LIST")
echo "== GQ of PASS 0/1 calls in $(basename "$g1"):"
awk -F'\t' '/^#/{next} $7=="PASS" && $5!="<*>" {n=split($9,f,":"); split($10,v,":"); for(i=1;i<=n;i++) x[f[i]]=v[i];
     g=x["GT"]; gsub(/[|]/,"/",g); if(g=="0/1"){q=x["GQ"]+0; b=(q<10?"<10":(q<20?"10-19":(q<25?"20-24":(q<35?"25-34":">=35")))); c[b]++; t++}}
     END{for(k in c) printf "   GQ %-6s %6d (%.1f%%)\n", k, c[k], 100*c[k]/t}' "$g1" | sort

## 3 GLnexus
echo "== GLnexus start $(date)"
"$GLX" --config "$CFG" --bed "$BED" --dir "$DB" --threads "$THREADS" --list "$LIST" > "${OUT}/joint_AQ20.bcf" \
  || die "GLnexus failed"
[ -s "${OUT}/joint_AQ20.bcf" ] || die "joint_AQ20.bcf empty"
echo "== GLnexus done $(date)"

## 4 BCF -> VCF.GZ + index
"$BCF" view --threads 4 -Oz -o "${OUT}/joint_AQ20.vcf.gz" "${OUT}/joint_AQ20.bcf" || die "bcftools view failed"
"$BCF" index --threads 4 -f "${OUT}/joint_AQ20.vcf.gz" || die "index failed"

## 5 QC vs original (genotype-based; INFO/AC is empty in GLnexus output)
{
  echo "QC joint_AQ20 vs joint (original) | $(date)"
  for v in "$ORIG" "${OUT}/joint_AQ20.vcf.gz"; do
    echo "--- $(basename "$v"): records $("$BCF" index -n "$v") | samples $("$BCF" query -l "$v" | wc -l)"
    echo "    chr22 genotype-based AC spectrum:"
    "$BCF" query -r chr22 -f '[%GT\t]\n' "$v" | awk -F'\t' '{ac=0; for(i=1;i<NF;i++){g=$i; gsub(/[|]/,"/",g); n=split(g,a,"/");
        for(j=1;j<=n;j++) if(a[j]!="." && a[j]!="0") ac++} b=(ac<=2?ac:(ac<=5?"3-5":(ac<=10?"6-10":">10"))); c[b]++; t++}
        END{for(k in c) printf "      AC=%-5s %7d (%.1f%%)\n", k, c[k], 100*c[k]/t; printf "      records %d\n", t}' | sort -V
  done
} | tee "$QC"
echo "== done $(date). After checking ${QC}: rm -rf ${DB}  (and optionally ${OUT}/joint_AQ20.bcf)"
