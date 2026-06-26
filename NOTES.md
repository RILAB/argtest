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
- `--suffix-to-strip` (default `""`) is removed via simple string replacement before matching; provide names that match the post-stripping individual names.
