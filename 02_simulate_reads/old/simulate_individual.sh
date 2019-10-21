#!/usr/bin/env bash

# haplotypes_gzfile=${1}
# individual_n=${2}
# reference_fasta=${3}
# haplotypes_map_file=${4}
# readLength=100
# coverage=0.05




haplotypes_gzfile="../input_data/haplotypes.vcf.gz"
individual_n="1"
reference_fasta="../input_data/dgrp2.reference.fasta"
haplotypes_map_file="../01_forward_simulator/hybrid_swarm_32_F5_1.haps"

for

sbatch --array=1-3%3 ./simulate_reads.sh  ../input_data/haplotypes.vcf.gz 2L ../input_data/dgrp2.reference.fasta ../01_forward_simulator/hybrid_swarm_32_F5_1.haps 100 0.05
