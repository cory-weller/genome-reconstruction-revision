#!/usr/bin/env bash


#!/usr/bin/env Rscript

library(Rsubread)

# index already built
# buildindex(basename="dmel", reference="../data_external/all_dmel.fasta")

args <- commandArgs(trailingOnly=TRUE)
file_id <- args[1]

readDir <- "/mnt/internal_1/cory/IndoorFlyCages/data_raw/cDNA_from_RNA_Seq/cDNA_master_pool/"
indexLoc <- "/mnt/internal_1/cory/IndoorFlyCages/Proj_Map_Reads/dmel"

read1 <- paste( #FILL , file_id, ".F.fq.bz2", sep="")
read2 <- paste( #FILL , file_id, ".R.fq.bz2", sep="")

align(index=indexLoc,
  type="dna",
  readfile1=read1,
  readfile2=read2,
  input_format="gzFASTQ",
  output_format="BAM", output_file=paste(file_id, ".Rsubread.bam", sep=""), nthreads=1)



  ../bin/pear-0.9.11-linux-x86_64/bin/pear

#!/usr/bin/env bash

individual_n=${1}
chromosome=${2}
filestem=${individual_n}.${chromosome}

# create subfolder
mkdir -p ${filestem} && cd ${filestem}

# run PEAR to assemble paired end reads
../../../bin/pear-0.9.11-linux-x86_64/bin/pear -k -f ../${filestem}.F.fq.bz2 -r ../${filestem}.R.fq.bz2 -o ${filestem}


# map to reference genome wit Rsubread
Rscript ${filestem} - <<EOF
args <- commandArgs(trailingOnly=TRUE)
filestem <- args[1]

library(Rsubread)

index_location <- "../../../input_data/dgrp2"

# map assembled reads
assembled <- paste(filestem, ".assembled.fastq", sep="")
align(index=index_location,
  type="dna",
  readfile1=assembled,
  input_format="FASTQ",
  output_format="BAM", output_file=paste(filestem, ".assembled.bam", sep=""), nthreads=1)

# map unassembled reads
read1unassembled <- paste(filestem, ".unassembled.forward.fastq", sep="")
read2unassembled <- paste(filestem, ".unassembled.reverse.fastq", sep="")
align(index=index_location,
  type="dna",
  readfile1=read1unassembled,
  readfile2=read2unassembled,
  input_format="FASTQ",
  output_format="BAM", output_file=paste(filestem, ".unassembled.bam", sep=""), nthreads=1)
EOF

# merge bam files with samtools
module load samtools/1.9
samtools merge ${filestem}.merged.bam ${filestem}.assembled.bam ${filestem}.unassembled.bam

# Sort bam file
samtools sort ${filestem}.merged.bam > ${filestem}.sorted.bam


# Add read groupsand index .bam file
module load gatk/4.0.0.0
echo "adding read groups"
gatk AddOrReplaceReadGroups --INPUT ${filestem}.sorted.bam --OUTPUT ${filestem}.bam --RGLB library1 --RGPL illumina --RGPU platform1 --RGSM ${filestem}
samtools index ${filestem}.bam
