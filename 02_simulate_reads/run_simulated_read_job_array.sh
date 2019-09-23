#!/usr/bin/env bash
vcf_file=${1}
chromosome=${2}
reference_fasta=${3}
haplotypes_file=${4}
readLength=${5}
coverage=${6}
N_individuals_to_simulate=${7}

echo $vcf_file
echo $chromosome
echo $reference_fasta
echo $haplotypes_file
echo $readLength
echo $coverage
echo $N_individuals_to_simulate

filestem=$(basename ${haplotypes_file} | (read s; echo ${s%_1.haps}))
echo filestem is: $filestem
mkdir -p $filestem && cd $filestem

sbatch --array=1-${N_individuals_to_simulate}%1 ../simulate_reads.sh \
  ../${vcf_file} \
  ${chromosome} \
  ../${reference_fasta} \
  ../${haplotypes_file} \
  ${readLength} \
  ${coverage}
