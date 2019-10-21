#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(foreach)

liability <- function(x) {1/(1+exp(-1*x))}

# create dummy data set
dd <- data.table(x=0.5)
g1 <- ggplot(dd) + stat_function(fun=liability) + xlim(-5,5) + ylim(0,1) +
labs(x="Liability / Risk score", y="Probability of 'case' assignment") +
theme_few(16)

ggplot(data=data.table(x=rpois(10000, lambda=20)), aes(x=x)) + geom_histogram(binwidth=1)

# Relationship between MAF and effect size
# Linear from A to -A over the interval (0,1)
# y = A - x

# single locus effect
sle <- function(maf) {
  midpoint <- 0.1 - 0.2*maf
  rnorm(1, mean=midpoint, sd=0.005)
}

# liability from 0 to 0.1 is an increase from 50% to 52.5% chance of assignment to case group
o <- foreach(maf=seq(0.01, 1, 0.01), .combine="rbind") %do% {
  data.table("maf"=maf,
            "effect"=sapply(1:1000, function(x) sle(maf))
          )
}

ggplot(o, mapping=aes(x=factor(maf), y=effect)) + geom_boxplot()

n_loci <- rpois(1, lambda=10)

vcf <- fread('zcat ../input_data/haplotypes.polarized.vcf.gz')


vcf[, nRef := apply(.SD, 1, function(x) sum(x=="1", na.rm=TRUE)), .SDcols=colnames(vcf)[colnames(vcf) %like% "RAL"] ]
vcf[, nAlt := apply(.SD, 1, function(x) sum(x=="-1", na.rm=TRUE)), .SDcols=colnames(vcf)[colnames(vcf) %like% "RAL"] ]

vcf[, nTot := nRef + nAlt]

vcf[, refFreq := nRef/nTot]
vcf[, altFreq := nAlt/nTot]
vcf[, MAF := ifelse(refFreq < altFreq, refFreq, altFreq)]

causative <- vcf[sample(.N, size=20,replace=F)]
causative[, effect := sapply(MAF, function(x) sle(x))]

# effect is multiplied by alternate allele dosage to get individual risk score for "case" assignment
