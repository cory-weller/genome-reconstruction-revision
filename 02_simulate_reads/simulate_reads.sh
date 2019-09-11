#!/usr/bin/env bash

# arguments
# input reference FASTA
# input haplotype map
# input VCF

# haplotypes_gzfile=${1}
# individual_n=${2}
# reference_fasta=${3}
# haplotypes_map_file=${4}

haplotypes_gzfile="../input_data/haplotypes.polarized.vcf.gz"
individual_n="1"
reference_fasta="../input_data/dgrp2.reference.fasta"
haplotypes_map_file="../01_forward_simulator/hybrid_swarm_32_F5_1.haps"

# map to VCF

module load gcc/7.1.0
module load java/1.8.0
module load R/3.5.1
module load tabix/0.2.6


# BEGIN conversion of haplotype map to VCF file
# Need to be two parents wide; not single parent

Rscript - <<EOF
library(data.table)
library(foreach)

haplotypes_map_file='${haplotypes_map_file}'
ind_n <- '${individual_n}'

dat <- fread(haplotypes_map_file)
# Collapse redundant breakpoints
dat[, lineID_rleid := rleid(chromosome, haplotype, ind, lineID)]
haps.collapsed <- dat[, list(chromosome,haplotype,start=min(start), stop=max(stop), ind, sex, lineID), by=lineID_rleid]
haps.collapsed[,lineID_rleid := NULL]
dat <- haps.collapsed[!duplicated(haps.collapsed)]

fillVectorNAs <- function(x) {
        na.omit(x)[cumsum(!is.na(x))]
}

fillDataTableNAs <- function(DT, cols) {
        DT[, (cols) := lapply(.SD, function(x) fillVectorNAs(x)), .SDcols=cols]
}

# subset specific individual from table
dat <- dat[ind==ind_n]

# only include X chromosome if female
if (length(unique(dat[chromosome=="X", haplotype])) == 1) {
  # male
  chromosomes <- c("2L","2R","3L","3R")
} else {
  # female
  chromosomes <- c("2L","2R","3L","3R","X")

}


o <- foreach(chromosome_n=chromosomes, .combine="rbind") %do% {
  dat.chr_subset <- dat[chromosome==chromosome_n]
  unique_stops <- unique((sort(dat.chr_subset[,stop])))
  unique_starts <- c(0,unique_stops[1:(length(unique_stops)-1)]) +1
  dat.out <- data.table(start=unique_starts, stop=unique_stops)

  setkey(dat.out, start)
  setkey(dat.chr_subset, start)

  dat.out <- dcast(dat.out[dat.chr_subset], start+stop~haplotype, value.var="lineID")
  dat.out[, chromosome := chromosome_n]

  fillDataTableNAs(dat.out, c("1","2"))

  setcolorder(dat.out, c("chromosome","start","stop","1","2"))
  setnames(dat.out, c("1","2"), c("par1","par2"))
  return(dat.out[])
}

write.table(o, file=paste(ind_n, ".wide.haps", sep=""), quote=F, row.names=F, col.names=T, sep="\t")
EOF





# Set up associative array of lineID : column, to be used by TABIX
vals=($(zgrep -m 1 "^CHROM" ${haplotypes_gzfile} | xargs | cut -f 10-))
declare -A founderIndices

for i in $(seq 0 "${#vals[@]}"); do
    let j=$i+1
    founderIndices[[${vals[$i]}]]=$j
done

# Read through diploid paths, extracting genotypes with TABIX, appending to estimate.vcf file
# Here is where you would change the bgzipped vcf filename to whatever one has all the sites you're trying to pull out based on the path
while read chromosome start stop par1 par2; do
    col1=${founderIndices[[$par1]]}
    col2=${founderIndices[[$par2]]}
    tabix $haplotypes_gzfile ${chromosome}:${start}-${stop} | awk -v col1=$col1 -v col2=$col2 '{print $1,$2,$col1,$col2}' >> ${ind_n}.estimate.vcf
done < <(awk 'NR > 1 {print}' ${ind_n}.wide.haps)
gzip ${ind_n}.estimate.vcf




# vcf to personalized fasta
# personalized fasta to .fastq
# .fastq to mapped .bam
