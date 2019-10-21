# RABBIT reconstruction with founder subset

## Generating HARP likelihoods
We use `HARP` to rank potential inbred lines by likelihood of being the
ancestor of a recombinant individual. `HARP` is ran via

```
sbatch --array=1,2,5,50%1 doHarp.slurm
```

which performs likelihood calculations for F1, F2, F5, and F50 populations.
The job is ran in working memory on a ramdisk (`/dev/shm/`) due to the high
number of files written; as a result, only the final `.freq` output is written
to a physical disk.

After completing, `.freqs` output are compressed into a single `.zip` e.g. via

```
zip -T hybrid_swarm_32_F2.freq.zip *.freqs
```

## Calculating Read Counts

```
sbatch --array=1,2,5,50%1 countReads.slurm
```

Which runs `ASEReadCounter` (of `GATK`) for each `.bam` file, using a custom
`VCF` file where each 2L site is heterozygous.

## Calculating most likely ancestors
