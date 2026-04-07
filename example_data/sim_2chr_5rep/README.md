# Simulated 2-Chromosome Example Dataset

Synthetic dataset for testing the pipeline and Snakemake workflow.

- trees/chr1/1.trees ... trees/chr1/5.trees
- trees/chr2/1.trees ... trees/chr2/5.trees
- masks/chr1.example.mask.bed, masks/chr2.example.mask.bed
- hapmap/chr1.hapmap.tsv, hapmap/chr2.hapmap.tsv, hapmap/all.hapmap.tsv
- chr1.mut_rate.p, chr2.mut_rate.p, sim.fai

Designed signals:
- low-accessibility windows from zero-rate segments in mut_rate maps
- mutational-load outliers: ind1, ind2 in chr1 and ind3, ind4 in chr2
