# Outputs

← [Back to README](../README.md)

The directory layout under `out_dir` is shown in [Outputs](../README.md#outputs) in the README.

## Merged output naming

The final merged outputs (one genome-wide tree sequence per replicate) are named:

```text
<base_name>.combined.<replicate>.<suffix>
```

and are written under the configured `out_dir` in a `combined/` directory.

## Locating a tree by chromosome and position

The merge lays chromosomes end-to-end along a single coordinate axis, so a within-chromosome position (e.g. position 1234 on chromosome 8) sits at `chromosome_offset + position` in the merged sequence. To spare you that arithmetic, the merge records a `chrom_offsets` table (`[{chrom, offset, length}, ...]`) in the merged tree sequence's top-level metadata.

Use [`scripts/locate_tree.py`](../scripts/locate_tree.py) to find the local tree at a `(chromosome, position)`:

```bash
python scripts/locate_tree.py --ts <out_dir>/combined/<base>.combined.<rep>.tsz --chrom 8 --position 1234
```

It prints the genome coordinate, the covering tree's index/interval, and flags when the position falls in a masked/trimmed region (a tree with no topology, `num_edges == 0`). The same mapping is available programmatically in [`scripts/argtest_common.py`](../scripts/argtest_common.py):

```python
from argtest_common import load_ts, tree_at_chrom_position, genome_position, chrom_position_from_genome
ts = load_ts("combined/run.combined.101.tsz")
tree = tree_at_chrom_position(ts, "8", 1234)   # local tree at chr8:1234
gpos = genome_position(ts, "8", 1234)          # concatenated coordinate
chrom, pos = chrom_position_from_genome(ts, gpos)  # inverse mapping
```

**Merged files produced before this feature** lack `chrom_offsets` and will raise a `KeyError` asking you to re-merge. For those, either read the local tree straight from the per-chromosome file (native coordinates, no offset needed) — `tszip.decompress(".../step5_trimmed_samples/<chrom>/<rep>.tsz").at(position)` — or add the offset by hand. Since trimming preserves each chromosome's length, the offset table is identical across replicates, so you only need to compute it once (sum the per-chromosome `sequence_length`s in natural chromosome order).

## VCF export

Set `emit_vcf: true` to write one bgzip-less `.vcf.gz` per `(chromosome, replicate)` under `<out_dir>/vcf/<chrom>/<rep>.vcf.gz`, produced by [`scripts/export_vcf.py`](../scripts/export_vcf.py) from the filtered per-chromosome tree sequences: step 5 output by default, or step 5b output when `min_samples` is enabled. Coordinates are real per-chromosome positions, not the concatenated coordinates of the merged ARG. Notes:

- **Variable sites only** — the records are the sites carried on the ARG; the pipeline does not synthesize invariant/monomorphic positions.
- **Pruned samples are missing, not dropped.** Because `trim_samples` leaves a pruned sample *isolated* over its intervals, that sample is written as a missing genotype (`.`) at any site inside those intervals (via tskit's `isolated_as_missing`), while remaining present elsewhere. The sample/site is never globally removed.
- **REF/ALT are ancestral/derived, not reference/alternate.** tskit's VCF writer sets `REF` to each site's `ancestral_state` and `ALT` to the derived states of its mutations, so genotype `0` is ancestral and `1+` is derived. `REF` need not match a reference-genome base: it is whatever ancestral state the ARG carries, which for inferred ARGs depends on how the input was polarized (e.g. if unpolarized calls were given to the inference tool, "ancestral" may in practice be the original reference allele). Allele strings are also copied verbatim, so a tree sequence storing states as `0`/`1` (e.g. simulated data) yields `REF=0`, `ALT=1` rather than nucleotides.
- **Ploidy-aware genotypes** — a haploid individual (one sample node, e.g. the `admix` data) is coded as a single allele `0`/`1`; a diploid individual as `0|1`-style. ARG genotypes are phased (`|`).
- **One VCF per replicate** — each replicate is a distinct ARG with its own topology, sites, and per-replicate trimming, so VCFs differ across replicates. The replicate set follows the pipeline `burnin` (leading replicates already dropped); restrict further with `vcf_reps`. To get one genome-wide VCF per replicate, `bcftools concat` the per-chromosome files.
