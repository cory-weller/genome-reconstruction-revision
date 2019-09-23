#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-0:02:00
#SBATCH -p standard
#SBATCH --account berglandlab

individual_n=${SLURM_ARRAY_TASK_ID}

population=${1}
chromosome=${2}

filestem=${individual_n}.${chromosome}

# create subfolder
cd $population && \
mkdir -p ${filestem} && \
cd ${filestem}

# run PEAR to assemble paired end reads
../../../bin/pear-0.9.11-linux-x86_64/bin/pear -k -f ../${filestem}.F.fq.bz2 -r ../${filestem}.R.fq.bz2 -o ${filestem}


module load bwa/0.7.17
module load samtools/1.9

# generate bam files for assembled and unassembled reads
bwa mem ../../../input_data/dgrp2.reference.fasta ${filestem}.assembled.fastq | samtools view -b -o ${filestem}.assembled.bam -
bwa mem ../../../input_data/dgrp2.reference.fasta ${filestem}.unassembled.forward.fastq ${filestem}.unassembled.reverse.fastq | samtools view -b -o ${filestem}.unassembled.bam -

# merge bam files with samtools
samtools merge ${filestem}.merged.bam ${filestem}.assembled.bam ${filestem}.unassembled.bam

# Sort bam file
samtools sort ${filestem}.merged.bam > ${filestem}.sorted.bam

# Add read groupsand index .bam file
module load gatk/4.0.0.0
echo "adding read groups"
gatk AddOrReplaceReadGroups --INPUT ${filestem}.sorted.bam --OUTPUT ${filestem}.bam --RGLB library1 --RGPL illumina --RGPU platform1 --RGSM ${filestem}
samtools index ${filestem}.bam

mv ${filestem}.bam* ../ && cd ../ && rm -rf ./${filestem}/
