#!/usr/bin/env bash

sbatch --array=331-1000%1 map_bwa.sh hybrid_swarm_32_F1 2L
sbatch --array=339-1000%1 map_bwa.sh hybrid_swarm_32_F2 2L
sbatch --array=385-1000%1 map_bwa.sh hybrid_swarm_32_F5 2L
sbatch --array=207-1000%1 map_bwa.sh hybrid_swarm_32_F50 2L
