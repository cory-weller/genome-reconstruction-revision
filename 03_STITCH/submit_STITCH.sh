#!/usr/bin/env bash

sbatch --array=1-10%3 STITCH.slurm hybrid_swarm_32_F1
sbatch --array=1-10%3 STITCH.slurm hybrid_swarm_32_F2
sbatch --array=1-10%3 STITCH.slurm hybrid_swarm_32_F5
sbatch --array=1-10%3 STITCH.slurm hybrid_swarm_32_F50

