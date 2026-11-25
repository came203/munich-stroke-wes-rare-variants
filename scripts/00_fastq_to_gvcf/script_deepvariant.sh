#!/bin/bash
#SBATCH -J test_deep
#SBATCH -o ./slurm_deep.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=teramem_inter
#SBATCH --mem=120gb
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=all
#SBATCH --mail-user=<your.email>
#SBATCH --export=NONE
#SBATCH --time=8:00:00


module load charliecloud/0.30
export PROJECTNAME=TEST
export BASE=<SAMPLE_ID>
#####MOVE FILES TO ${BASE}_1.fastq.gz
#export PATH=/cluster/rmalik/:$PATH
export HOMEDIR=/dss/dsshome1/06/ru36ros
export DSSDIR=/dss/dssfs02/lwp-dss-0001/pn52ra/pn52ra-dss-0000
export PROJECTDIR=${DSSDIR}/${BASE}
#export OLDDIR=/cluster2/rmalik2/EXIS
#export TOOLDIR=${OLDDIR}/TOOLS
export REFDIR=${HOMEDIR}/REFERENCE
export REFERENCE=${REFDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
export REFERENCE2=${REFDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
export ALDIR=${PROJECTDIR}/ALIGNMENTS
#export RESOURCES_DIR=/cluster2/rmalik2/EXIS/RESOURCES
#export GATK_DIR=${PROJECTDIR}/GATK
export TMPDIR=/dss/lxclscratch/06/ru36ros

ch-run ${HOMEDIR}/deep -w -b /dss/:/dss/ -- /opt/deepvariant/bin/run_deepvariant --model_type=WES --ref=${REFERENCE2} --reads=${ALDIR}/${BASE}_abra.bam --output_vcf=${ALDIR}/${BASE}.vcf --output_gvcf=${ALDIR}/${BASE}.g.vcf --num_shards=$(nproc) --regions=${HOMEDIR}/twistCustomExomeV2_2021_12_14.bed   --intermediate_results_dir ${TMPDIR}
