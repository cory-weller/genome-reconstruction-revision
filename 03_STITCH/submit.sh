#!/usr/bin/env bash

parallel -j 1 sbatch --array=1-10%1 {1} run.sh ::: hybrid_swarm_32_F1 hybrid_swarm_32_F2 hybrid_swarm_32_F5 hybrid_swarm_32_F50
