#!/usr/bin/env bash

# Prepare PICARD sequence dictionary
module load picard/2.20.6

java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  R=dgrp2.reference.fasta \
  O=dgrp2.reference.dict

# Prepare SAMTOOLS index
module load samtools/1.9
samtools faidx dgrp2.reference.fasta

# Prepare BWA index
module load bwa/0.7.17
bwa index dgrp2.reference.fasta
