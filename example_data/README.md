# Example Data

This directory contains small tree-sequence fixtures used by the tests and documentation examples:

- `test_100trees.tsz` - compact single-file fixture.
- `maize.tsz` - compact maize-derived fixture for script smoke tests.
- `sim_2chr_5rep/` - two-chromosome, five-replicate simulated dataset with matching HapMap, `.fai`, masks, and mutation-rate files.

For the current end-to-end example dataset, use the repository-level `argtest-realistic-example/` workflow described in the main `README.md` and `MAKE_REALISTIC_EXAMPLE.md`.

For mutational-load outputs, use `scripts/mutload_summary.py` for the HTML diagnostic and `scripts/mutload_masks.py` for BED masks.
