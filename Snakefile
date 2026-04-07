from __future__ import annotations

import fnmatch
import re
from pathlib import Path

from snakemake.exceptions import WorkflowError


VALID_SUFFIXES = {".tree", ".trees", ".tsz"}


def natural_key(text: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def cfg_required(key: str):
    if key not in config or config[key] in (None, ""):
        raise WorkflowError(
            f"Missing required config key '{key}'. "
            "Provide a config file with --configfile (see config/snakemake.example.yaml)."
        )
    return config[key]


ROOT_DIR = Path(cfg_required("root_dir")).expanduser().resolve()
if not ROOT_DIR.exists():
    raise WorkflowError(f"root_dir does not exist: {ROOT_DIR}")
if not ROOT_DIR.is_dir():
    raise WorkflowError(f"root_dir must be a directory: {ROOT_DIR}")

TREE_PATTERN = str(config.get("tree_pattern", "*"))
OUT_DIR = Path(config.get("out_dir", "snakemake_out"))
BASE_NAME = str(config.get("base_name", ROOT_DIR.name))
ALLOW_MISSING_REPLICATES = bool(config.get("allow_missing_replicates", False))
SUFFIX_TO_STRIP = str(config.get("suffix_to_strip", "_anchorwave"))

HAPMAP = Path(cfg_required("hapmap")).expanduser().resolve()
FAI = Path(cfg_required("fai")).expanduser().resolve()
REC_FRACTION = float(cfg_required("rec_fraction"))
LOW_ACCESS_WINDOW_SIZE = float(cfg_required("low_access_window_size"))
LOW_ACCESS_CUTOFF_BP = float(cfg_required("low_access_cutoff_bp"))
MUTLOAD_CUTOFF = float(config.get("mutload_cutoff", 0.25))
MUTLOAD_FRACTION = config.get("mutload_fraction", None)
if MUTLOAD_FRACTION is not None:
    MUTLOAD_FRACTION = float(MUTLOAD_FRACTION)
    if not 0 <= MUTLOAD_FRACTION <= 1:
        raise WorkflowError("mutload_fraction must be between 0 and 1")

MUTLOAD_WINDOW_SIZE = config.get("mutload_window_size", None)
MUTLOAD_SNP_WINDOW = config.get("mutload_snp_window", None)
if (MUTLOAD_WINDOW_SIZE is None) == (MUTLOAD_SNP_WINDOW is None):
    raise WorkflowError("Set exactly one of mutload_window_size or mutload_snp_window")
if MUTLOAD_WINDOW_SIZE is not None:
    MUTLOAD_WINDOW_SIZE = float(MUTLOAD_WINDOW_SIZE)
    if MUTLOAD_WINDOW_SIZE <= 0:
        raise WorkflowError("mutload_window_size must be > 0")
    MUTLOAD_WINDOW_ARG = f"--window-size {MUTLOAD_WINDOW_SIZE}"
else:
    MUTLOAD_SNP_WINDOW = int(MUTLOAD_SNP_WINDOW)
    if MUTLOAD_SNP_WINDOW <= 0:
        raise WorkflowError("mutload_snp_window must be > 0")
    MUTLOAD_WINDOW_ARG = f"--snp-window {MUTLOAD_SNP_WINDOW}"

MUTLOAD_FRACTION_ARG = f"--fraction {MUTLOAD_FRACTION}" if MUTLOAD_FRACTION is not None else ""

MERGED_OUT_SUFFIX = config.get("merged_out_suffix", None)
if MERGED_OUT_SUFFIX is not None and MERGED_OUT_SUFFIX not in VALID_SUFFIXES:
    raise WorkflowError(f"merged_out_suffix must be one of {sorted(VALID_SUFFIXES)}")
MERGE_SUFFIX_ARG = f"--out-suffix {MERGED_OUT_SUFFIX}" if MERGED_OUT_SUFFIX else ""

STEP1_DIR = OUT_DIR / "step1_low_rec"
STEP2_DIR = OUT_DIR / "step2_low_access"
STEP3_DIR = OUT_DIR / "step3_mutload"
STEP4_MASK_DIR = OUT_DIR / "step4_masks"
STEP4_TRIM_DIR = OUT_DIR / "step4_trimmed_regions"
STEP5_DIR = OUT_DIR / "step5_trimmed_samples"
MERGED_DIR = OUT_DIR / "combined"
LOG_DIR = OUT_DIR / "logs"


def discover_tree_files():
    chrom_to_rep = {}
    for chrom_dir in sorted([p for p in ROOT_DIR.iterdir() if p.is_dir()], key=lambda p: natural_key(p.name)):
        by_rep = {}
        for path in sorted(chrom_dir.iterdir(), key=lambda p: natural_key(p.name)):
            if not path.is_file():
                continue
            if path.suffix not in VALID_SUFFIXES:
                continue
            if not fnmatch.fnmatch(path.name, TREE_PATTERN):
                continue
            rep = path.stem
            if rep in by_rep:
                raise WorkflowError(
                    f"Duplicate replicate '{rep}' in chromosome directory {chrom_dir}. "
                    "Ensure one tree file per replicate per chromosome."
                )
            by_rep[rep] = path
        if by_rep:
            chrom_to_rep[chrom_dir.name] = by_rep
    if not chrom_to_rep:
        raise WorkflowError(
            f"No chromosome subdirectories with matching tree files under {ROOT_DIR} "
            f"(pattern='{TREE_PATTERN}', suffixes={sorted(VALID_SUFFIXES)})."
        )
    return chrom_to_rep


CHROM_TO_REP = discover_tree_files()
CHROMS = sorted(CHROM_TO_REP.keys(), key=natural_key)
REPLICATE_UNION = sorted(
    {rep for reps in CHROM_TO_REP.values() for rep in reps.keys()},
    key=natural_key,
)

if not ALLOW_MISSING_REPLICATES:
    missing_messages = []
    for rep in REPLICATE_UNION:
        missing = [chrom for chrom in CHROMS if rep not in CHROM_TO_REP[chrom]]
        if missing:
            missing_messages.append(f"{rep}: missing in {', '.join(missing)}")
    if missing_messages:
        raise WorkflowError(
            "Some replicates are missing chromosome files.\n"
            + "\n".join(missing_messages)
            + "\nSet allow_missing_replicates: true to allow partial concatenation."
        )

REPLICATES = REPLICATE_UNION

TS_LOOKUP = {
    (chrom, rep): path
    for chrom, reps in CHROM_TO_REP.items()
    for rep, path in reps.items()
}

if ALLOW_MISSING_REPLICATES:
    PAIRS = sorted(TS_LOOKUP.keys(), key=lambda item: (natural_key(item[1]), natural_key(item[0])))
else:
    PAIRS = [(chrom, rep) for rep in REPLICATES for chrom in CHROMS]

PAIR_EXT = {(chrom, rep): TS_LOOKUP[(chrom, rep)].suffix.lstrip(".") for chrom, rep in PAIRS}

REPLICATE_EXT = {}
for rep in REPLICATES:
    rep_paths = [TS_LOOKUP[(chrom, rep)] for chrom in CHROMS if (chrom, rep) in TS_LOOKUP]
    if not rep_paths:
        continue
    REPLICATE_EXT[rep] = rep_paths[0].suffix.lstrip(".")

STEP1_TARGETS = [str(STEP1_DIR / f"{chrom}.low_rec.mask.bed") for chrom in CHROMS]
STEP2_TARGETS = [str(STEP2_DIR / chrom / f"{chrom}.low_access.bed") for chrom in CHROMS]
STEP5_TARGETS = [str(STEP5_DIR / chrom / f"{rep}.{PAIR_EXT[(chrom, rep)]}") for chrom, rep in PAIRS]
MERGED_TARGETS = [
    str(MERGED_DIR / f"{BASE_NAME}.combined.{rep}.{REPLICATE_EXT[rep]}")
    for rep in REPLICATES
    if rep in REPLICATE_EXT
]

if not MERGED_TARGETS:
    raise WorkflowError("No merged outputs were inferred from discovered inputs.")


def tree_inputs_for_chrom(wildcards):
    return [str(path) for path in CHROM_TO_REP[wildcards.chrom].values()]


def ts_input_for_pair(wildcards):
    key = (wildcards.chrom, wildcards.rep)
    if key not in TS_LOOKUP:
        raise WorkflowError(f"Unknown chromosome/replicate pair: {key}")
    return str(TS_LOOKUP[key])


def ts_input_for_pair_ext(wildcards):
    key = (wildcards.chrom, wildcards.rep)
    if key not in TS_LOOKUP:
        raise WorkflowError(f"Unknown chromosome/replicate pair: {key}")
    expected_ext = TS_LOOKUP[key].suffix.lstrip(".")
    if wildcards.ext != expected_ext:
        raise WorkflowError(
            f"Extension mismatch for {key}: expected .{expected_ext}, got .{wildcards.ext}"
        )
    return str(TS_LOOKUP[key])


rule all:
    input:
        MERGED_TARGETS


rule step1_low_rec_masks:
    input:
        hapmap=str(HAPMAP),
        fai=str(FAI),
    output:
        STEP1_TARGETS
    params:
        rec_fraction=REC_FRACTION,
        out_dir=str(STEP1_DIR),
    shell:
        """
        python scripts/hapmap_low_rec_mask.py \
          --hapmap {input.hapmap} \
          --fai {input.fai} \
          --rec-fraction {params.rec_fraction} \
          --out-dir {params.out_dir}
        """


rule step2_low_access:
    input:
        tree_inputs_for_chrom
    output:
        str(STEP2_DIR / "{chrom}" / "{chrom}.low_access.bed")
    params:
        ts_dir=lambda wc: str(ROOT_DIR / wc.chrom),
        window_size=LOW_ACCESS_WINDOW_SIZE,
        cutoff_bp=LOW_ACCESS_CUTOFF_BP,
        pattern=TREE_PATTERN,
    shell:
        """
        python scripts/find_low_access_regions.py \
          --ts-dir {params.ts_dir} \
          --window-size {params.window_size} \
          --cutoff-bp {params.cutoff_bp} \
          --pattern '{params.pattern}' \
          --out {output}
        """


rule step3_mutload_masks:
    input:
        ts=ts_input_for_pair
    output:
        outlier=str(STEP3_DIR / "{chrom}" / "{rep}.outliers.bed"),
        masked=str(STEP3_DIR / "{chrom}" / "{rep}.mutation_masked.bed"),
    params:
        window_arg=MUTLOAD_WINDOW_ARG,
        cutoff=MUTLOAD_CUTOFF,
        fraction_arg=MUTLOAD_FRACTION_ARG,
        suffix_to_strip=SUFFIX_TO_STRIP,
    shell:
        """
        python scripts/mutload_masks.py \
          --ts {input.ts} \
          --chrom {wildcards.chrom} \
          --outlier-bed {output.outlier} \
          --masked-bed {output.masked} \
          {params.window_arg} \
          --cutoff {params.cutoff} \
          {params.fraction_arg} \
          --suffix-to-strip {params.suffix_to_strip}
        """


rule step4_combine_remove_masks:
    input:
        low_rec=str(STEP1_DIR / "{chrom}.low_rec.mask.bed"),
        low_access=str(STEP2_DIR / "{chrom}" / "{chrom}.low_access.bed"),
        masked=str(STEP3_DIR / "{chrom}" / "{rep}.mutation_masked.bed"),
    output:
        str(STEP4_MASK_DIR / "{chrom}" / "{rep}.remove_regions.bed")
    shell:
        """
        python scripts/combine_remove_masks.py \
          --chrom {wildcards.chrom} \
          --out {output} \
          --inputs {input.low_rec} {input.low_access} {input.masked}
        """


rule step4_trim_regions_single:
    input:
        ts=ts_input_for_pair_ext,
        mask_bed=str(STEP4_MASK_DIR / "{chrom}" / "{rep}.remove_regions.bed"),
    output:
        str(STEP4_TRIM_DIR / "{chrom}" / "{rep}.{ext}")
    shell:
        """
        python scripts/trim_regions_single.py \
          --ts {input.ts} \
          --remove {input.mask_bed} \
          --out {output}
        """


rule step5_trim_samples_single:
    input:
        ts=str(STEP4_TRIM_DIR / "{chrom}" / "{rep}.{ext}"),
        outlier=str(STEP3_DIR / "{chrom}" / "{rep}.outliers.bed"),
    output:
        str(STEP5_DIR / "{chrom}" / "{rep}.{ext}")
    params:
        suffix_to_strip=SUFFIX_TO_STRIP,
        log=lambda wc: str(LOG_DIR / "step5" / wc.chrom / f"{wc.rep}.trim_samples.log"),
    shell:
        """
        python scripts/trim_samples.py \
          {input.ts} \
          --remove {input.outlier} \
          --out {output} \
          --suffix-to-strip {params.suffix_to_strip} \
          --log {params.log}
        """


rule merge_replicates:
    input:
        STEP5_TARGETS
    output:
        MERGED_TARGETS
    params:
        ts_dir=str(STEP5_DIR),
        out_dir=str(MERGED_DIR),
        base_name=BASE_NAME,
        merge_suffix_arg=MERGE_SUFFIX_ARG,
    shell:
        """
        python scripts/merge_treefiles_by_replicate.py \
          --ts-dir {params.ts_dir} \
          --layout nested \
          --base-name {params.base_name} \
          --out-dir {params.out_dir} \
          {params.merge_suffix_arg}
        """
