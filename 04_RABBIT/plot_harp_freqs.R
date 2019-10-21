#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)

dat <- fread('5.2L.freqs')

priorsHeader <- system('zgrep -m 1 "^" ../2L.subset.priors.csv.gz', intern=TRUE)
priorsHeader <- unlist(strsplit(priorsHeader, ","))
lineIDs <- priorsHeader[3:(length(priorsHeader)-1)]

setnames(dat, c("chromosome","start","stop",lineIDs))

dat.long <- melt(dat, measure.vars=lineIDs)

ggplot(dat.long, aes(x=start, y=variable, fill=value)) + geom_tile()
