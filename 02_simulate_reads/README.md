# 02_simulate_reads

## Convert haplotype map to wide

```
bash ./convert_haps_to_wide.sh \
	-p hybrid_swarm_32_F1,hybrid_swarm_32_F2,hybrid_swarm_32_F5,hybrid_swarm_32_F50 \
	-f 1 \
	-l 1000 \
	-c 2L
```

## Generate individual genotypes
```
sbatch --array=1,2,5,50%4 tabix_genotypes.slurm hybrid_swarm_32 2L
```

# Generate individual fasta files and simulate reads
```
sbatch --array=1,2,5,50 simulate_reads.slurm hybrid_swarm_32 2L 1 1000
```

# Map reads
```
sbatch --array=1,2,5,50 map_reads.slurm hybrid_swarm_32 2L 1 1000
```
