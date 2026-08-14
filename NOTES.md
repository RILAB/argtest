# ARGtest: additional notes

Supplementary details that aren't needed for day-to-day use of the pipeline or scripts. See [README.md](README.md) for install, workflow, and the Snakemake pipeline. Run any script with `--help` for its flags.

## Sandboxed environments (read-only `~/.cache`)

Some HPC or container setups mount `~/.cache` read-only. There, prefix the Snakemake command (dry-run or real run) with cache and temp-dir redirects to `/tmp`:

```bash
XDG_CACHE_HOME=/tmp/argtest-xdg-cache TMPDIR=/tmp/argtest-tmp \
    snakemake --cores 16 --configfile config/snakemake.yaml
```

On a normal machine where `~/.cache` is writable this is not needed.

## Shared module

`scripts/argtest_common.py` contains shared tree-sequence helpers used by multiple scripts:
- TS I/O (`load_ts`, `dump_ts`)
- mutational load / stat helpers
- trimming and masking helpers

Use this module for internal script imports rather than duplicating the helpers.

## Sample ID matching (`trim_samples.py`)

`trim_samples.py` matches sample/individual IDs exactly against the tree sequence's internal individual names as produced by `argtest_common.get_individual_name()` (prefers `individual.metadata['id']` when present, otherwise a synthetic `ind<id>` name).

- Matching is **exact** and **case-sensitive**.
- `--name-substring-to-remove` (default `""`) is removed via global string replacement before matching; provide names that match the normalized individual names. Despite replacing the former `--suffix-to-strip` option, the operation is not suffix-specific: every occurrence is removed.

Per-individual analyses may pool multiple sample nodes for one individual (for example, two nodes for a diploid). Uniform ploidy is assessed only among represented individuals; individual-table rows with no sample nodes are ignored. The stricter uniform-ploidy and leaf-sample contract remains audit/evidence-gated until it has been checked on the real input corpus.
