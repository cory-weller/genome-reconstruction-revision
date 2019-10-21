#!/usr/bin/env R

if(! require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if(! require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

if(! require(foreach)) {
  install.packages("foreach")
  library(foreach)
}


# load in each accuracy table

mainstring <- "hybrid_swarm_32"
generations <- c("F1", "F2", "F5", "F50")
population_sizes <- seq(100,1000,100)

o <- foreach(generation=generations, .combine="rbind", .errorhandling="remove") %do% {
  foreach(population_size=population_sizes, .combine="rbind") %do% {
    folder_name <- paste(mainstring, "_", generation, "_", population_size, sep="")
    dat <- readRDS(paste(folder_name, "/accuracy.RDS", sep=""))
    return(dat)
  }
}

o[, population_size := NULL]
o[, N := factor(as.numeric(N))]

fwrite(o, file="STITCH_accuracy.txt", quote=F, row.names=F, col.names=T, sep="\t")

g <- ggplot(o, mapping=aes(x=N, y=percent_correct)) +
  geom_jitter(alpha=0.2, shape=21) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  facet_grid(n_generations~.) +
  labs(x="Population Size", y="Genotype Accuracy (% sites)", title="STITCH accuracy\n(pseudoHaploid model)")

ggsave(g, file="STITCH_accuracy.png", width=10, height=20, units="cm")
