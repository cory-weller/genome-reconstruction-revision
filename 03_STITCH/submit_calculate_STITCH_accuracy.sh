#!/usr/bin/env bash

for population in hybrid_swarm*; do   # edit glob pattern to match folders names
    sbatch calculate_STITCH_accuracy.slurm ${population}
done 
