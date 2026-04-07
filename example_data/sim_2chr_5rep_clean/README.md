# Simulated 2-Chromosome Example Dataset (Clean)

Synthetic balanced dataset for testing pipeline behavior.

- two chromosomes (`chr1`, `chr2`)
- five replicates per chromosome (`1`-`5`)
- treefiles: `trees/<chrom>/<rep>.trees`
- per-chromosome mask BEDs in `masks/`
- HapMap files in `hapmap/`
- mutation-rate maps (`chr1.mut_rate.p`, `chr2.mut_rate.p`)

Signals included:
- low-access regions in mutation-rate maps
- individuals with elevated derived mutations:
  - chr1: ind1, ind2
  - chr2: ind3, ind4
