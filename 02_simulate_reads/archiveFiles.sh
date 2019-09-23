#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-2:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

folder=${1}

cd ${folder} && \
ls | tar -czv -f ${folder}.1-1000.tar.gz --files-from -
