#!/bin/bash
#SBATCH -J test
#SBATCH --output=/dss/dsshome1/06/ru36ros/LOGS/R-%x.%j.out
#SBATCH --error=/dss/dsshome1/06/ru36ros/LOGS/R-%x.%j.err
#SBATCH -D ./
#SBATCH --get-user-env

#SBATCH --clusters=serial
#SBATCH --partition=serial_long
#SBATCH --mem=120G

#SBATCH --mail-type=fail
#SBATCH --mail-user=<your.email>
#SBATCH --export=NONE
#SBATCH --time=10:00:00

module load slurm_setup
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load charliecloud/0.30
#module load charliecloud/0.22
#####to run sbatch --export=ALL,var_name=${i}
#####MOVE FILES TO ${BASE}_1.fastq.gz
################ONLY CHANGE THIS

export PROJECTNAME=TEST
export BASE=${var_name}
export TMPDIR=/dss/dssfs03/pn52re/pn52re-dss-0001/TMP/${BASE}
mkdir -p /dss/dssfs03/pn52re/pn52re-dss-0001/TMP/${BASE}
#####MOVE FILES TO ${BASE}_1.fastq.gz
export HOMEDIR=/dss/dsshome1/06/ru36ros
#export DSSDIR=/dss/dssfs03/pn52re/pn52re-dss-0001
export DSSDIR=/dss/dssfs02/lwp-dss-0001/pn52ra/pn52ra-dss-0000
export PROJECTDIR=${DSSDIR}/${BASE}
#export PROJECTDIR=${DSSDIR}/${BASE}
export REFDIR=${HOMEDIR}/REFERENCE
export REFERENCE=${REFDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
export REFERENCE2=${REFDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
export ALDIR=${DSSDIR}/${BASE}/ALIGNMENTS

export BWA_CMD=/dss/dsshome1/06/ru36ros/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
export SAM_CMD=/dss/dsshome1/06/ru36ros/samtools-1.19.2/samtools
export BCF_CMD=/dss/dsshome1/06/ru36ros/bcftools-1.19/bcftools
export NGS_CMD='ch-run /dss/dsshome1/06/ru36ros/ngsb -w -b /dss/dssfs02/lwp-dss-0001/pn52ra/:/dss/dssfs02/lwp-dss-0001/pn52ra/ -- SeqPurge'
export JAVA_CMD='ch-run ./java8/ -w -b /dss/:/dss/ -- /opt/java/openjdk/bin/java'
mkdir -p ${PROJECTDIR}
mkdir -p ${ALDIR}


###trim and combine reads
${NGS_CMD} -in1 ${PROJECTDIR}/${BASE}_L001_R1_001.fastq.gz ${PROJECTDIR}/${BASE}_L002_R1_001.fastq.gz ${PROJECTDIR}/${BASE}_L003_R1_001.fastq.gz ${PROJECTDIR}/${BASE}_L004_R1_001.fastq.gz -in2 ${PROJECTDIR}/${BASE}_L001_R2_001.fastq.gz ${PROJECTDIR}/${BASE}_L002_R2_001.fastq.gz ${PROJECTDIR}/${BASE}_L003_R2_001.fastq.gz ${PROJECTDIR}/${BASE}_L004_R2_001.fastq.gz -out1 ${PROJECTDIR}/${BASE}_1_val_1.fq.gz -out2 ${PROJECTDIR}/${BASE}_2_val_2.fq.gz




###ALIGNMENT 8 min

#${BWA_CMD} mem -t 224  -M -R '@RG\tID:'${BASE}'\tPU:80C3GABXX110103\tLB:Solexa-272222\tPL:ILLUMINA\tPM:HISEQ\tSM:'${BASE} ${REFERENCE} ${PROJECTDIR}/${BASE}_1_val_1.fq.gz  ${PROJECTDIR}/${BASE}_2_val_2.fq.gz | ${HOMEDIR}/samblaster -M  | ${SAM_CMD} view -Sb - > ${ALDIR}/${BASE}.bam.uns



#${SAM_CMD} sort  -@ 224 -m 1G ${ALDIR}/${BASE}.bam.uns -o ${ALDIR}/${BASE}.bam

#${SAM_CMD} index ${ALDIR}/${BASE}.bam

#${JAVA_CMD} -Xmx28G -jar ${DSSDIR}/abra2-2.19.jar --in ${ALDIR}/${BASE}.bam --out ${ALDIR}/${BASE}_abra.bam --ref ${REFERENCE2} --threads 224 --targets ${HOMEDIR}/twistCustomExomeV2_2021_12_14.bed --tmpdir ${DSSDIR}/TMP/

#${SAM_CMD} index ${ALDIR}/${BASE}_abra.bam


##ch-run ${HOMEDIR}/deep161 -w -b /dss/:/dss/ -- /opt/deepvariant/bin/run_deepvariant --model_type=WES --ref=${REFERENCE2} --reads=${ALDIR}/${BASE}_abra.bam --output_vcf=${ALDIR}/${BASE}.vcf --output_gvcf=${ALDIR}/${BASE}.g.vcf --num_shards=56 --regions=${HOMEDIR}/twistCustomExomeV2_2021_12_14.bed
#ch-run ${HOMEDIR}/deep180 -w -b /dss/:/dss/ -- /opt/deepvariant/bin/run_deepvariant --model_type=WES --ref=${REFERENCE2} --reads=${ALDIR}/${BASE}_abra.bam --output_vcf=${ALDIR}/${BASE}.vcf --output_gvcf=${ALDIR}/${BASE}.g.vcf --num_shards=224 --regions=${HOMEDIR}/twistCustomExomeV2_2021_12_14.bed --intermediate_results_dir=/dss/dssfs03/pn52re/pn52re-dss-0001/TMP/${BASE}

#echo ${ALDIR}'/'${BASE}'.g.vcf' >> /dss/dsshome1/06/ru36ros/mergelist.txt
