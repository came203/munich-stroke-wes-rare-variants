#!/bin/bash
#set -e  # uncomment to exit on any error

#SBATCH -J join
#SBATCH --output=/dss/dsshome1/06/ru36ros/LOGS/R-%x.%j.out
#SBATCH --error=/dss/dsshome1/06/ru36ros/LOGS/R-%x.%j.err
#SBATCH -D ./
#SBATCH --get-user-env

#SBATCH --clusters=inter
#SBATCH --partition=teramem_inter
#SBATCH --mem=20G # Consider if 20G is enough for ~1500 samples, might need more
#SBATCH --mail-type=FAIL,END # Added END notification
#SBATCH --mail-user=<your.email>
#SBATCH --export=NONE
#SBATCH --time=100:00:00

module load charliecloud/0.30 # Ensure this and necessary bioinformatics tools (bcftools, htslib) are available

## ——— Static paths ———
export HOMEDIR=/dss/dsshome1/06/ru36ros
export DSSDIR=/dss/dssfs02/lwp-dss-0001/pn52ra/pn52ra-dss-0000
export REFDIR=${HOMEDIR}/REFERENCE # Ensure your reference genome is here if needed by tools
export TMPDIR=/dss/lxclscratch/06/ru36ros # Make sure this exists and has space

## Directory of per-sample GVCFs
SAMPLES_DIR=${DSSDIR}/GVCF

echo "Starting gVCF processing at $(date)"

## ——— Build the initial merge list ———
MERGE_LIST_FULL=${HOMEDIR}/mergelist_glnexus_full.txt
find "${SAMPLES_DIR}" -maxdepth 1 -type f -name '*.g.vcf' | sort > "${MERGE_LIST_FULL}"
echo "Found $(wc -l < "${MERGE_LIST_FULL}") total gVCFs -> ${MERGE_LIST_FULL}"

## --- Define problematic files (basename only) ---
# List derived from the previous error log
declare -a problematic_files=(
    "19535.g.vcf" "19565.g.vcf" "19592.g.vcf" "19655.g.vcf" "19695.g.vcf"
    "19713.g.vcf" "19736.g.vcf" "19737.g.vcf" "19741.g.vcf" "19745.g.vcf"
    "19762.g.vcf" "19789.g.vcf" "19821.g.vcf" "19841.g.vcf" "19861.g.vcf"
    "26539.g.vcf" "26579.g.vcf" "26638.g.vcf" "26822.g.vcf" "26988.g.vcf"
    "27000.g.vcf"
)

## --- Create the filtered merge list ---
MERGE_LIST_FILTERED=${HOMEDIR}/mergelist_glnexus_filtered.txt
# Start with an empty file
> "${MERGE_LIST_FILTERED}"

echo "Filtering problematic files..."
filtered_count=0
excluded_count=0
while IFS= read -r gvcf_path; do
    gvcf_basename=$(basename "${gvcf_path}")
    exclude=0
    for bad_file in "${problematic_files[@]}"; do
        if [[ "${gvcf_basename}" == "${bad_file}" ]]; then
            exclude=1
            break
        fi
    done

    if [[ $exclude -eq 0 ]]; then
        echo "${gvcf_path}" >> "${MERGE_LIST_FILTERED}"
        ((filtered_count++))
    else
        echo "Excluding: ${gvcf_basename}"
        ((excluded_count++))
    fi
done < "${MERGE_LIST_FULL}"

echo "Created filtered list with ${filtered_count} gVCFs -> ${MERGE_LIST_FILTERED}"
echo "Excluded ${excluded_count} problematic gVCFs."

# --- Check if bcftools and htslib paths are correct ---
BCFTOOLS_PATH="${HOMEDIR}/bcftools-1.19/bcftools"
HTSLIB_BGZIP_PATH="${HOMEDIR}/htslib-1.19.1/bgzip"

# --- Check if executables exist ---
if [ ! -x "${HOMEDIR}/glnexus_cli" ]; then echo "Error: glnexus_cli not found or not executable"; exit 1; fi
if [ ! -x "${BCFTOOLS_PATH}" ]; then echo "Error: bcftools not found or not executable"; exit 1; fi
if [ ! -x "${HTSLIB_BGZIP_PATH}" ]; then echo "Error: bgzip not found or not executable"; exit 1; fi
if [ ! -f "${HOMEDIR}/twistCustomExomeV2_2021_12_14.bed" ]; then echo "Error: BED file not found"; exit 1; fi
if [ ! -s "${MERGE_LIST_FILTERED}" ]; then echo "Error: Filtered merge list is empty!"; exit 1; fi # Check if list is not empty

echo "Launching GLnexus joint calling at $(date)"
# Use the filtered list
# Consider adding --threads option if glnexus supports it and you have multiple cores allocated
# Make sure the output directory (${HOMEDIR}) is writable and has enough space
${HOMEDIR}/glnexus_cli \
  -c DeepVariantWES \
  -b ${HOMEDIR}/twistCustomExomeV2_2021_12_14.bed \
  --list "${MERGE_LIST_FILTERED}" \
  > ${HOMEDIR}/joint.bcf

# Check if GLnexus succeeded before proceeding
if [ $? -ne 0 ]; then
    echo "Error: GLnexus failed at $(date)"
    exit 1 # Exit if GLnexus failed
fi
if [ ! -s "${HOMEDIR}/joint.bcf" ]; then
    echo "Error: GLnexus finished but output file ${HOMEDIR}/joint.bcf is empty or missing at $(date)"
    exit 1 # Exit if output is missing/empty
fi


echo "Converting and compressing to VCF.GZ at $(date)"
# Pipe BCF output directly to bgzip
"${BCFTOOLS_PATH}" view ${HOMEDIR}/joint.bcf \
  | "${HTSLIB_BGZIP_PATH}" -@ 4 -c \
  > ${HOMEDIR}/joint.vcf.gz

# Check if bcftools/bgzip succeeded
if [ $? -ne 0 ]; then
    echo "Error: bcftools view or bgzip failed at $(date)"
    # Optionally clean up intermediate file if desired
    # rm -f ${HOMEDIR}/joint.bcf
    exit 1
fi

echo "Successfully created ${HOMEDIR}/joint.vcf.gz"

# Optional: Clean up the intermediate BCF file
# echo "Cleaning up intermediate BCF file at $(date)"
# rm -f ${HOMEDIR}/joint.bcf

echo "Done at $(date)"
