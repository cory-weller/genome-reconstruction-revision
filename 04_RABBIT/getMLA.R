#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

population <- args[1]
population_i <- args[2]
ind <- args[3]
chromosome <- args[4]

priors.Gzfile <- paste(population, "/", chromosome, ".subset.priors.csv.gz", sep="")

library(data.table)

harp.freqs <- fread(paste("./", population, "/HARP_likelihood/", ind, ".", chromosome, ".freqs", sep=""), header=FALSE, showProgress=FALSE, na.strings="-nan")

priorsHeader <- system(paste('zgrep -m 1 "^" ', priors.Gzfile, sep=""), intern=TRUE)
priorsHeader <- unlist(strsplit(priorsHeader, ","))
lineIDs <- priorsHeader[3:(length(priorsHeader)-1)]

setnames(harp.freqs, c("chromosome","start","stop", lineIDs))

# Subset to only include known founders

harp.freqs.long <- melt(harp.freqs, measure.vars = colnames(harp.freqs)[4:length(colnames(harp.freqs))], variable="lineID", value="freq")
harp.freqs.long[, q99 := quantile(freq, 0.99, na.rm=TRUE), by=chromosome]

mla <- harp.freqs.long[, .N, by=list(chromosome, lineID, freq>=q99)][freq==TRUE][order(chromosome,-N)][,c("chromosome","lineID","N")]

print(mla)

#write.table(mla, file="", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


# Subset priors
