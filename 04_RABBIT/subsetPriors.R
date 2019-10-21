#!/usr/bin/env R

library(data.table)

args <- commandArgs(trailingOnly=TRUE)

population <- args[1]
population_i <- args[2]
chromosome <- args[3]

output.fn <- paste("./", population, "/", chromosome, ".subset.priors.csv", sep="")

priors <- fread(cmd=paste("zcat ", chromosome, ".priors.csv.gz", sep=""))

founders.list <- fread(paste("../01_forward_simulator/", population, "_", population_i, ".founders", sep=""), header=FALSE)
founders.list <- founders.list[,V1]
n.founders <- length(founders.list)


subset.priors <- priors[, c(chromosome, "Ref", founders.list), with=FALSE]
subset.priors[, "Coverage" := n.founders ]
fwrite(subset.priors, file=output.fn, row.names=F, col.names=T, sep=",", quote=F)
