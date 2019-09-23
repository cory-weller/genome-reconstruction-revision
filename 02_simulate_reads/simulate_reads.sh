#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-0:02:00
#SBATCH -p standard
#SBATCH --account berglandlab


# arguments
# input reference FASTA
# input haplotype map
# input VCF

haplotypes_gzfile=${1}
individual_n=${SLURM_ARRAY_TASK_ID}
chromosome=${2}
reference_fasta=${3}
haplotypes_map_file=${4}
readLength=${5}
coverage=${6}

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
vals=($(zgrep -m 1 "^#CHROM" ${haplotypes_gzfile} | xargs | cut -f 10-))
declare -A founderIndices

for i in $(seq 0 "${#vals[@]}"); do
    let j=$i+1
    founderIndices[[${vals[$i]}]]=$j
done

# Read through diploid paths, extracting genotypes with TABIX, appending to estimate.vcf file
# Here is where you would change the bgzipped vcf filename to whatever one has all the sites you're trying to pull out based on the path
while read chromosomeN start stop par1 par2; do
    if [[ $chromosomeN == $chromosome ]]; then
    col1=${founderIndices[[$par1]]}
    col2=${founderIndices[[$par2]]}
    tabix $haplotypes_gzfile ${chromosomeN}:${start}-${stop} | awk -v col1=$col1 -v col2=$col2 '{print $1,$2,$col1,$col2}' >> ${individual_n}.vcf
    fi
done < <(awk 'NR > 1 {print}' ${individual_n}.wide.haps)
# remove inside 'if' block to instead do all chromosomes

# vcf to personalized fasta
module load python/3.6.8

python -  <<EOF

import copy

ind_n = "${individual_n}"
chromosome = "${chromosome}"
input_fasta = "../../input_data/" + chromosome + ".fa"
input_genotypes = ind_n + ".vcf"
vcf_sites = "../../input_data/" + chromosome + ".sites"

ref_allele = "0/0"
alt_allele = "1/1"
missing_data = "./."

print(ind_n, chromosome, input_fasta)

with open(input_fasta, 'r') as infile:
  fasta = infile.read().splitlines()
  haplotype1 = list(''.join(fasta[1:]))

haplotype2 = haplotype1.copy()

# Create dictionary for this chromosome
sites = {}
with open(vcf_sites, 'r') as infile:
  for line in infile:
    chromosome, pos, ref, alt = line.rstrip().split()
    sites[int(pos)] = [ref,alt]

with open(input_genotypes, 'r') as infile:
  for line in infile:
    chrom, pos, hap1, hap2 = line.split()
    if chrom == chromosome:
      pos = int(pos)
      if hap1 == alt_allele:
        haplotype1[pos+1] = sites[pos][1]
      elif hap1 == missing_data:
        haplotype1[pos+1] = "N"
      if hap2 == alt_allele:
        haplotype2[pos+1] = sites[pos][1]
      elif hap2 == missing_data:
        haplotype2[pos+1] = "N"

with open(ind_n + "." + chromosome + ".fasta", 'w') as outfile:
  outfile.write(">" + ind_n + "_" + chromosome + "_haplotype1" "\n")
  outfile.write('\n'.join([''.join(haplotype1[i:i+50]) for i in range(0,len(haplotype1),50)]))
  outfile.write(">" + ind_n + "_" + chromosome + "_haplotype2" "\n")
  outfile.write('\n'.join([''.join(haplotype2[i:i+50]) for i in range(0,len(haplotype2),50)]))

EOF


diploidGenomeSize=$(grep "^[^>]" ${individual_n}.${chromosome}.fasta | wc -c)

function getNReads {
# getNReads $readLength $coverage $diploidGenomeSize
Rscript - <<-EOF
    readLength <- as.numeric($1)
    coverage <- as.numeric($2)
    genomeSize <- as.numeric($3)/2
    nReads = round(genomeSize * coverage / (2*readLength)) # paired end
    cat(nReads)
EOF
}

export -f getNReads

nReads=$(getNReads $readLength $coverage $diploidGenomeSize)

# From personalized fasta to fastq files with wgsim
        echo generating reads with wgsim
../../bin/wgsim-master/wgsim \
  -1 $readLength \
  -2 $readLength \
  -N $nReads \
  -e 0.001 \
  -r 0 \
  -R 0 \
  ${individual_n}.${chromosome}.fasta \
  ${individual_n}.${chromosome}.F.fq \
  ${individual_n}.${chromosome}.R.fq && \
  rm ${individual_n}.${chromosome}.fasta && \
  bzip2 ${individual_n}.${chromosome}.F.fq && \
  bzip2 ${individual_n}.${chromosome}.R.fq

exit 0
