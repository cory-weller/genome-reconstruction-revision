#!/usr/bin/env bash

vcf_file="../input_data/haplotypes.vcf.gz"
chromosome="2L"
reference_fasta="../input_data/dgrp2.reference.fasta"
# haplotypes file varies
readLength="100"
coverage="0.05"
N_individuals_to_simulate="1000"

for haplotypes_file in ../01_forward_simulator/*.haps; do
  bash run_simulated_read_job_array.sh \
    ${vcf_file} \
    ${chromosome} \
    ${reference_fasta} \
    ${haplotypes_file} \
    ${readLength} \
    ${coverage} \
    ${N_individuals_to_simulate}
done
