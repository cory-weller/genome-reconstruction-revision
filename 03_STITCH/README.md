# 03_STITCH

## Installation
Install `STITCH` using `install_STITCH.sh`

## Running STITCH on simulated populations
The `STITCH.slurm` script, submitted with `submit_STITCH.sh`, is ran as a job array for each of the simulated populations.

## Calculating STITCH accuracy
The `submit_calculate_STITCH_accuracy.sh` script submits `calculate_STITCH_accuracy.slurm` for simulated populations. Each job in turn executes `calculate_STITCH_accuracy.R`. The `R` script converts the `.vcf` file from `STITCH` into `.gds` format, and compares the estimated genotypes to the true (simulated) genotypes.
