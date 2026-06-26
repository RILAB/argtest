from __future__ import annotations

import fnmatch
import hashlib
import re
from pathlib import Path

from snakemake.exceptions import WorkflowError


VALID_SUFFIXES = {".ts", ".trees", ".tsz"}


def natural_key(text: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def cfg_required(key: str):
    if key not in config or config[key] in (None, ""):
        raise WorkflowError(
            f"Missing required config key '{key}'. "
            "Provide a config file with --configfile (see config/snakemake.yaml)."
        )
    return config[key]


ROOT_DIR = Path(cfg_required("root_dir")).expanduser().resolve()
if not ROOT_DIR.exists():
    raise WorkflowError(f"root_dir does not exist: {ROOT_DIR}")
if not ROOT_DIR.is_dir():
    raise WorkflowError(f"root_dir must be a directory: {ROOT_DIR}")

TREE_PATTERN = str(config.get("tree_pattern", "*"))
# Optional subdirectory within each chromosome dir that holds the tree files
# (e.g. "trees" for SINGER-style layouts). Empty/unset = files live directly in
# the chromosome dir.
TREE_SUBDIR = str(config.get("tree_subdir", "") or "").strip("/")
OUT_DIR = Path(config.get("out_dir", "snakemake_out")).expanduser()

# --- Cluster resources (single-config pattern; see the "resources" config block) ---
# Per-rule mem/time/threads/partition all live in the user's --configfile. Rules
# not listed under "resources:" fall back to these defaults. Used by the SLURM
# profile (profiles/slurm). The `threads` value also affects local scheduling
# with `--cores`; the SLURM-only fields are inert for local runs.
DEFAULT_MEM_MB = 4000
DEFAULT_THREADS = 1
DEFAULT_TIME = "01:00:00"
_RESOURCES = config.get("resources", {}) or {}
SLURM_PARTITION = str(config.get("slurm_partition", "low"))
RESOURCE_RULES = {
    "step1_low_rec_masks",
    "step2_low_access",
    "step3_mutload_masks",
    "step4_combine_remove_masks",
    "step4_trim_regions_single",
    "step5_trim_samples_single",
    "merge_replicates",
    "step6_validation_plots",
    "step7_summary",
}
RESOURCE_KEYS = {"mem_mb", "time", "threads", "partition"}

if not isinstance(_RESOURCES, dict):
    raise WorkflowError("resources must be a mapping of rule names to resource settings")
for rule_name, rule_resources in _RESOURCES.items():
    if rule_name not in RESOURCE_RULES:
        raise WorkflowError(
            f"Unknown resources entry '{rule_name}'. Expected one of: "
            f"{', '.join(sorted(RESOURCE_RULES))}"
        )
    if not isinstance(rule_resources, dict):
        raise WorkflowError(f"resources.{rule_name} must be a mapping")
    unknown_keys = sorted(set(rule_resources) - RESOURCE_KEYS)
    if unknown_keys:
        raise WorkflowError(
            f"Unknown resources.{rule_name} key(s): {', '.join(unknown_keys)}. "
            f"Expected only: {', '.join(sorted(RESOURCE_KEYS))}"
        )

def _res(rule_name, key, default):
    return (_RESOURCES.get(rule_name, {}) or {}).get(key, default)

def res_mem_mb(rule_name):    return int(_res(rule_name, "mem_mb", DEFAULT_MEM_MB))
def res_time(rule_name):      return str(_res(rule_name, "time", DEFAULT_TIME))
def res_threads(rule_name):   return int(_res(rule_name, "threads", DEFAULT_THREADS))
def res_partition(rule_name): return str(_res(rule_name, "partition", SLURM_PARTITION))
# Treat an empty/whitespace base_name the same as an absent key: fall back to
# the root_dir name. Resolving the concrete name here (and passing it to
# merge_treefiles_by_replicate.py) keeps the Snakefile targets and the merge
# script's output names in lockstep — an empty string would otherwise make the
# script fall back to its own (different) default and the names would mismatch.
BASE_NAME = str(config.get("base_name") or "").strip() or ROOT_DIR.name
ALLOW_MISSING_REPLICATES = bool(config.get("allow_missing_replicates", False))
BURNIN = int(config.get("burnin", 0))
if BURNIN < 0:
    raise WorkflowError("burnin must be >= 0")
SUFFIX_TO_STRIP = str(config.get("suffix_to_strip", "_anchorwave"))

HAPMAP = Path(cfg_required("hapmap")).expanduser().resolve()
FAI = Path(cfg_required("fai")).expanduser().resolve()
REC_FRACTION = float(cfg_required("rec_fraction"))
LOW_ACCESS_WINDOW_SIZE = float(cfg_required("low_access_window_size"))
LOW_ACCESS_CUTOFF_BP = float(cfg_required("low_access_cutoff_bp"))
MUTLOAD_CUTOFF = float(config.get("mutload_cutoff", 0.5))
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

# Shared organism-wide mutation rate (per bp per generation). Serves as the
# default for both the mutload sim-based expectation (step 3) and the validation
# plots (step 6). Either step may override it with its own key below; an explicit
# null on a per-step key disables that step's scalar (validation: null skips step
# 6 entirely). Only used as a fallback — an embedded ts.metadata ratemap or a
# sibling *.mut_rate.p file always takes precedence (see resolve_mu_rate).
MUTATION_RATE = config.get("mutation_rate", None)
if MUTATION_RATE is not None:
    MUTATION_RATE = float(MUTATION_RATE)
    if MUTATION_RATE <= 0:
        raise WorkflowError("mutation_rate must be > 0")

MUTLOAD_MUTATION_RATE = config.get("mutload_mutation_rate", MUTATION_RATE)
if MUTLOAD_MUTATION_RATE is not None:
    MUTLOAD_MUTATION_RATE = float(MUTLOAD_MUTATION_RATE)
    if MUTLOAD_MUTATION_RATE <= 0:
        raise WorkflowError("mutload_mutation_rate must be > 0")
MUTLOAD_MUTATION_RATE_ARG = (
    f"--mutation-rate {MUTLOAD_MUTATION_RATE}" if MUTLOAD_MUTATION_RATE is not None else ""
)

MUTLOAD_RANDOM_SEED = int(config.get("mutload_random_seed", 1))


def mutload_seed_for(chrom: str, rep: str) -> int:
    # Deterministic per-replicate seed combining base seed, chrom, and rep.
    blob = f"{MUTLOAD_RANDOM_SEED}|{chrom}|{rep}".encode()
    h = int(hashlib.sha1(blob).hexdigest()[:8], 16)
    # msprime accepts uint32; mod down and avoid 0 so seed is always valid.
    return (h % (2**31 - 1)) + 1

VALIDATION_MUTATION_RATE = config.get("validation_mutation_rate", MUTATION_RATE)
if VALIDATION_MUTATION_RATE is not None:
    VALIDATION_MUTATION_RATE = float(VALIDATION_MUTATION_RATE)
    if VALIDATION_MUTATION_RATE <= 0:
        raise WorkflowError("validation_mutation_rate must be > 0")
VALIDATION_FIRST_CHROM_ONLY = bool(config.get("validation_first_chrom_only", True))
VALIDATION_SIM_BRANCH = bool(config.get("validation_sim_branch", False))

MERGED_OUT_SUFFIX = config.get("merged_out_suffix", None)
if MERGED_OUT_SUFFIX is not None and MERGED_OUT_SUFFIX not in VALID_SUFFIXES:
    raise WorkflowError(f"merged_out_suffix must be one of {sorted(VALID_SUFFIXES)}")
MERGE_SUFFIX_ARG = f"--out-suffix {MERGED_OUT_SUFFIX}" if MERGED_OUT_SUFFIX else ""

# Optional VCF export from the filtered per-(chrom, rep) tree sequences (step 5).
EMIT_VCF = bool(config.get("emit_vcf", False))
_VCF_REPS_CFG = config.get("vcf_reps", None)
VCF_REPS_REQUESTED = (
    [str(r) for r in _VCF_REPS_CFG] if _VCF_REPS_CFG else None
)  # None/empty -> every post-burnin replicate

# Optional user-supplied sample trimming applied in step 5, IN ADDITION to the
# step-3 mutload outliers. Both apply identically to every (chrom, rep):
#   trim_individuals — IDs removed genome-wide (string or list, joined with ",")
#   trim_remove_bed  — one or more BED files of per-individual intervals
#                      (col 4 = comma-separated sample IDs; path or list of paths)
# Sample-ID matching follows the same rules as the step-3 outlier BED, including
# suffix_to_strip. Leave unset to trim only the mutload outliers.
_TRIM_IND_CFG = config.get("trim_individuals", None)
if isinstance(_TRIM_IND_CFG, (list, tuple)):
    _TRIM_IND_CFG = ",".join(str(i) for i in _TRIM_IND_CFG)
TRIM_INDIVIDUALS = str(_TRIM_IND_CFG).strip() if _TRIM_IND_CFG is not None else ""
TRIM_INDIVIDUALS_ARG = f"--individuals {TRIM_INDIVIDUALS}" if TRIM_INDIVIDUALS else ""

_TRIM_BED_CFG = config.get("trim_remove_bed", None)
if _TRIM_BED_CFG is None:
    _TRIM_BED_CFG = []
elif not isinstance(_TRIM_BED_CFG, (list, tuple)):
    _TRIM_BED_CFG = [_TRIM_BED_CFG]
TRIM_REMOVE_BEDS = [
    str(Path(str(p)).expanduser().resolve()) for p in _TRIM_BED_CFG if str(p).strip()
]
for _bed in TRIM_REMOVE_BEDS:
    if not Path(_bed).exists():
        raise WorkflowError(f"trim_remove_bed file not found: {_bed}")
TRIM_REMOVE_ARG = " ".join(f'--remove "{b}"' for b in TRIM_REMOVE_BEDS)

STEP1_DIR = OUT_DIR / "step1_low_rec"
STEP2_DIR = OUT_DIR / "step2_low_access"
STEP3_DIR = OUT_DIR / "step3_mutload"
STEP4_MASK_DIR = OUT_DIR / "step4_masks"
STEP4_TRIM_DIR = OUT_DIR / "step4_trimmed_regions"
STEP5_DIR = OUT_DIR / "step5_trimmed_samples"
MERGED_DIR = OUT_DIR / "combined"
STEP6_DIR = OUT_DIR / "step6_validation"
VCF_DIR = OUT_DIR / "vcf"
LOG_DIR = OUT_DIR / "logs"


def discover_tree_files():
    chrom_to_rep = {}
    for chrom_dir in sorted([p for p in ROOT_DIR.iterdir() if p.is_dir()], key=lambda p: natural_key(p.name)):
        search_dir = chrom_dir / TREE_SUBDIR if TREE_SUBDIR else chrom_dir
        if not search_dir.is_dir():
            continue
        by_rep = {}
        for path in sorted(search_dir.iterdir(), key=lambda p: natural_key(p.name)):
            if not path.is_file():
                continue
            if path.suffix not in VALID_SUFFIXES:
                continue
            if not fnmatch.fnmatch(path.name, TREE_PATTERN):
                continue
            rep = path.stem
            chrom_prefix = chrom_dir.name + "."
            if rep.startswith(chrom_prefix):
                rep = rep[len(chrom_prefix):]
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

if BURNIN >= len(REPLICATE_UNION):
    raise WorkflowError(
        f"burnin ({BURNIN}) must be less than the number of replicates ({len(REPLICATE_UNION)})"
    )

REPLICATES = REPLICATE_UNION[BURNIN:]

if not ALLOW_MISSING_REPLICATES:
    missing_messages = []
    for rep in REPLICATES:
        missing = [chrom for chrom in CHROMS if rep not in CHROM_TO_REP[chrom]]
        if missing:
            missing_messages.append(f"{rep}: missing in {', '.join(missing)}")
    if missing_messages:
        raise WorkflowError(
            "Some replicates are missing chromosome files.\n"
            + "\n".join(missing_messages)
            + "\nSet allow_missing_replicates: true to allow partial concatenation."
        )

TS_LOOKUP = {
    (chrom, rep): path
    for chrom, reps in CHROM_TO_REP.items()
    for rep, path in reps.items()
}

if ALLOW_MISSING_REPLICATES:
    _rep_set = set(REPLICATES)
    PAIRS = sorted(
        [(c, r) for (c, r) in TS_LOOKUP.keys() if r in _rep_set],
        key=lambda item: (natural_key(item[1]), natural_key(item[0])),
    )
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

VALIDATION_CHROMS = [CHROMS[0]] if VALIDATION_FIRST_CHROM_ONLY else CHROMS
STEP6_TARGETS = (
    [str(STEP6_DIR / chrom / "done.txt") for chrom in VALIDATION_CHROMS]
    if VALIDATION_MUTATION_RATE is not None
    else []
)

if EMIT_VCF:
    if VCF_REPS_REQUESTED is None:
        VCF_REPS = list(REPLICATES)
    else:
        _unknown = [r for r in VCF_REPS_REQUESTED if r not in set(REPLICATES)]
        if _unknown:
            raise WorkflowError(
                f"vcf_reps lists replicates not in the post-burnin set {REPLICATES}: {_unknown}"
            )
        VCF_REPS = VCF_REPS_REQUESTED
    _vcf_rep_set = set(VCF_REPS)
    VCF_TARGETS = [
        str(VCF_DIR / chrom / f"{rep}.vcf.gz")
        for (chrom, rep) in PAIRS
        if rep in _vcf_rep_set
    ]
else:
    VCF_TARGETS = []

SUMMARY_TARGET = str(OUT_DIR / "pipeline_summary.html")


def tree_inputs_for_chrom(wildcards):
    return [str(path) for path in CHROM_TO_REP[wildcards.chrom].values()]


def step5_inputs_for_rep(wildcards):
    return [
        str(STEP5_DIR / chrom / f"{wildcards.rep}.{PAIR_EXT[(chrom, wildcards.rep)]}")
        for chrom in CHROMS
        if (chrom, wildcards.rep) in TS_LOOKUP
    ]


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
        MERGED_TARGETS + STEP6_TARGETS + VCF_TARGETS + [SUMMARY_TARGET]


rule step1_low_rec_masks:
    input:
        hapmap=str(HAPMAP),
        fai=str(FAI),
    output:
        str(STEP1_DIR / "{chrom}.low_rec.mask.bed")
    threads: res_threads("step1_low_rec_masks")
    resources:
        mem_mb=res_mem_mb("step1_low_rec_masks"),
        time=res_time("step1_low_rec_masks"),
        slurm_partition=res_partition("step1_low_rec_masks"),
    log:
        str(LOG_DIR / "step1" / "{chrom}.log")
    params:
        rec_fraction=REC_FRACTION,
        out_dir=str(STEP1_DIR),
    shell:
        """
        python scripts/hapmap_low_rec_mask.py \
          --hapmap "{input.hapmap}" \
          --fai "{input.fai}" \
          --rec-fraction {params.rec_fraction} \
          --out-dir "{params.out_dir}" \
          --chrom {wildcards.chrom} \
          --log "{log}"
        """


def first_ts_for_chrom(wildcards):
    reps = sorted(CHROM_TO_REP[wildcards.chrom].keys(), key=natural_key)
    return str(CHROM_TO_REP[wildcards.chrom][reps[0]])


rule step2_low_access:
    input:
        first_ts_for_chrom
    output:
        str(STEP2_DIR / "{chrom}" / "{chrom}.low_access.bed")
    threads: res_threads("step2_low_access")
    resources:
        mem_mb=res_mem_mb("step2_low_access"),
        time=res_time("step2_low_access"),
        slurm_partition=res_partition("step2_low_access"),
    log:
        str(LOG_DIR / "step2" / "{chrom}.log")
    params:
        window_size=LOW_ACCESS_WINDOW_SIZE,
        cutoff_bp=LOW_ACCESS_CUTOFF_BP,
    shell:
        """
        python scripts/find_low_access_regions.py \
          --ts "{input}" \
          --window-size {params.window_size} \
          --cutoff-bp {params.cutoff_bp} \
          --out "{output}" \
          --log "{log}"
        """


rule step3_mutload_masks:
    input:
        ts=ts_input_for_pair
    output:
        outlier=str(STEP3_DIR / "{chrom}" / "{rep}.outliers.bed"),
        masked=str(STEP3_DIR / "{chrom}" / "{rep}.mutation_masked.bed"),
    threads: res_threads("step3_mutload_masks")
    resources:
        mem_mb=res_mem_mb("step3_mutload_masks"),
        time=res_time("step3_mutload_masks"),
        slurm_partition=res_partition("step3_mutload_masks"),
    log:
        str(LOG_DIR / "step3" / "{chrom}" / "{rep}.log")
    params:
        window_arg=MUTLOAD_WINDOW_ARG,
        cutoff=MUTLOAD_CUTOFF,
        fraction_arg=MUTLOAD_FRACTION_ARG,
        suffix_to_strip=SUFFIX_TO_STRIP,
        mutation_rate_arg=MUTLOAD_MUTATION_RATE_ARG,
        seed=lambda wildcards: mutload_seed_for(wildcards.chrom, wildcards.rep),
    shell:
        """
        python scripts/mutload_masks.py \
          --ts "{input.ts}" \
          --chrom {wildcards.chrom} \
          --outlier-bed "{output.outlier}" \
          --masked-bed "{output.masked}" \
          {params.window_arg} \
          --cutoff {params.cutoff} \
          {params.fraction_arg} \
          --random-seed {params.seed} \
          {params.mutation_rate_arg} \
          --suffix-to-strip "{params.suffix_to_strip}" \
          --log "{log}"
        """


rule step4_combine_remove_masks:
    input:
        low_rec=str(STEP1_DIR / "{chrom}.low_rec.mask.bed"),
        low_access=str(STEP2_DIR / "{chrom}" / "{chrom}.low_access.bed"),
        masked=str(STEP3_DIR / "{chrom}" / "{rep}.mutation_masked.bed"),
    output:
        str(STEP4_MASK_DIR / "{chrom}" / "{rep}.remove_regions.bed")
    threads: res_threads("step4_combine_remove_masks")
    resources:
        mem_mb=res_mem_mb("step4_combine_remove_masks"),
        time=res_time("step4_combine_remove_masks"),
        slurm_partition=res_partition("step4_combine_remove_masks"),
    log:
        str(LOG_DIR / "step4_masks" / "{chrom}" / "{rep}.log")
    shell:
        """
        python scripts/combine_remove_masks.py \
          --chrom {wildcards.chrom} \
          --out "{output}" \
          --inputs "{input.low_rec}" "{input.low_access}" "{input.masked}" \
          --log "{log}"
        """


rule step4_trim_regions_single:
    input:
        ts=ts_input_for_pair_ext,
        mask_bed=str(STEP4_MASK_DIR / "{chrom}" / "{rep}.remove_regions.bed"),
    output:
        str(STEP4_TRIM_DIR / "{chrom}" / "{rep}.{ext}")
    threads: res_threads("step4_trim_regions_single")
    resources:
        mem_mb=res_mem_mb("step4_trim_regions_single"),
        time=res_time("step4_trim_regions_single"),
        slurm_partition=res_partition("step4_trim_regions_single"),
    wildcard_constraints:
        rep="|".join(re.escape(r) for r in REPLICATES),
        ext="|".join(re.escape(s.lstrip(".")) for s in VALID_SUFFIXES),
    log:
        str(LOG_DIR / "step4_trim" / "{chrom}" / "{rep}.{ext}.log")
    shell:
        """
        python scripts/trim_regions_single.py \
          --ts "{input.ts}" \
          --remove "{input.mask_bed}" \
          --out "{output}" \
          --log "{log}"
        """


rule step5_trim_samples_single:
    input:
        ts=str(STEP4_TRIM_DIR / "{chrom}" / "{rep}.{ext}"),
        outlier=str(STEP3_DIR / "{chrom}" / "{rep}.outliers.bed"),
        extra_beds=TRIM_REMOVE_BEDS,
    output:
        str(STEP5_DIR / "{chrom}" / "{rep}.{ext}")
    threads: res_threads("step5_trim_samples_single")
    resources:
        mem_mb=res_mem_mb("step5_trim_samples_single"),
        time=res_time("step5_trim_samples_single"),
        slurm_partition=res_partition("step5_trim_samples_single"),
    wildcard_constraints:
        rep="|".join(re.escape(r) for r in REPLICATES),
        ext="|".join(re.escape(s.lstrip(".")) for s in VALID_SUFFIXES),
    log:
        str(LOG_DIR / "step5" / "{chrom}" / "{rep}.{ext}.log")
    params:
        suffix_to_strip=SUFFIX_TO_STRIP,
        extra_remove=TRIM_REMOVE_ARG,
        individuals_arg=TRIM_INDIVIDUALS_ARG,
    shell:
        """
        python scripts/trim_samples.py \
          "{input.ts}" \
          --remove "{input.outlier}" \
          {params.extra_remove} \
          {params.individuals_arg} \
          --out "{output}" \
          --suffix-to-strip "{params.suffix_to_strip}" \
          --log "{log}"
        """


rule merge_replicates:
    input:
        step5_inputs_for_rep
    output:
        str(MERGED_DIR / f"{BASE_NAME}.combined.{{rep}}.{{ext}}")
    threads: res_threads("merge_replicates")
    resources:
        mem_mb=res_mem_mb("merge_replicates"),
        time=res_time("merge_replicates"),
        slurm_partition=res_partition("merge_replicates"),
    wildcard_constraints:
        rep="|".join(re.escape(r) for r in REPLICATES),
        ext="|".join(re.escape(s.lstrip(".")) for s in VALID_SUFFIXES),
    log:
        str(LOG_DIR / "merge" / "{rep}.{ext}.log")
    params:
        ts_dir=str(STEP5_DIR),
        out_dir=str(MERGED_DIR),
        base_name=BASE_NAME,
        merge_suffix_arg=MERGE_SUFFIX_ARG,
    shell:
        """
        python scripts/merge_treefiles_by_replicate.py \
          --ts-dir "{params.ts_dir}" \
          --layout nested \
          --base-name "{params.base_name}" \
          --out-dir "{params.out_dir}" \
          {params.merge_suffix_arg} \
          --replicate "{wildcards.rep}" \
          >> "{log}" 2>&1
        """


def vcf_input_ts(wildcards):
    return str(STEP5_DIR / wildcards.chrom / f"{wildcards.rep}.{PAIR_EXT[(wildcards.chrom, wildcards.rep)]}")


if EMIT_VCF:
    rule export_vcf:
        input:
            vcf_input_ts
        output:
            str(VCF_DIR / "{chrom}" / "{rep}.vcf.gz")
        threads: res_threads("export_vcf")
        resources:
            mem_mb=res_mem_mb("export_vcf"),
            time=res_time("export_vcf"),
            slurm_partition=res_partition("export_vcf"),
        wildcard_constraints:
            rep="|".join(re.escape(r) for r in REPLICATES),
        log:
            str(LOG_DIR / "vcf" / "{chrom}" / "{rep}.log")
        params:
            suffix_to_strip=SUFFIX_TO_STRIP,
        shell:
            """
            python scripts/export_vcf.py \
              --ts "{input}" \
              --out "{output}" \
              --chrom "{wildcards.chrom}" \
              --suffix-to-strip "{params.suffix_to_strip}" \
              --log "{log}"
            """


def step6_inputs_for_chrom(wildcards):
    return [
        str(STEP5_DIR / wildcards.chrom / f"{rep}.{PAIR_EXT[(wildcards.chrom, rep)]}")
        for rep in REPLICATES
        if (wildcards.chrom, rep) in TS_LOOKUP
    ]


def step6_original_inputs_for_chrom(wildcards):
    return [
        str(TS_LOOKUP[(wildcards.chrom, rep)])
        for rep in REPLICATES
        if (wildcards.chrom, rep) in TS_LOOKUP
    ]


if VALIDATION_MUTATION_RATE is not None:
    rule step6_validation_plots:
        input:
            cleaned=step6_inputs_for_chrom,
            original=step6_original_inputs_for_chrom,
        output:
            str(STEP6_DIR / "{chrom}" / "done.txt")
        threads: res_threads("step6_validation_plots")
        resources:
            mem_mb=res_mem_mb("step6_validation_plots"),
            time=res_time("step6_validation_plots"),
            slurm_partition=res_partition("step6_validation_plots"),
        log:
            str(LOG_DIR / "step6" / "{chrom}.log")
        params:
            cleaned_out=lambda wildcards: str(STEP6_DIR / wildcards.chrom / "cleaned"),
            original_out=lambda wildcards: str(STEP6_DIR / wildcards.chrom / "original"),
            mutation_rate=VALIDATION_MUTATION_RATE,
            sim_branch_flag="--sim-branch" if VALIDATION_SIM_BRANCH else "",
        shell:
            """
            cleaned_stage_root="$(mktemp -d /tmp/argtest-step6-cleaned.XXXXXX)"
            original_stage_root="$(mktemp -d /tmp/argtest-step6-original.XXXXXX)"
            trap 'rm -rf "$cleaned_stage_root" "$original_stage_root"' EXIT

            cleaned_stage="$cleaned_stage_root/{wildcards.chrom}"
            original_stage="$original_stage_root/{wildcards.chrom}"
            mkdir -p "$cleaned_stage" "$original_stage"

            for f in {input.cleaned:q}; do
              ln -s "$(realpath "$f")" "$cleaned_stage/$(basename "$f")"
            done

            original_files=({input.original:q})
            for f in "${{original_files[@]}}"; do
              ln -s "$(realpath "$f")" "$original_stage/$(basename "$f")"
            done

            if [ "${{#original_files[@]}}" -gt 0 ]; then
              mu_path="$(python - "${{original_files[0]}}" <<'PY'
from pathlib import Path
import sys
sys.path.insert(0, "scripts")
from argtest_common import infer_mu_path
try:
    print(infer_mu_path(Path(sys.argv[1])))
except Exception:
    print("")
PY
)"
              if [ -n "$mu_path" ] && [ -f "$mu_path" ]; then
                ln -s "$mu_path" "$original_stage_root/$(basename "$mu_path")" || true
                ln -s "$mu_path" "$original_stage/$(basename "$mu_path")" || true
              fi
            fi

            python scripts/validation_plots_from_ts.py \
              --ts-dir "$cleaned_stage" \
              --pattern "*" \
              --burnin-frac 0 \
              --mutation-rate {params.mutation_rate} \
              --out-dir "{params.cleaned_out}" \
              {params.sim_branch_flag} \
              >> "{log}" 2>&1
            python scripts/validation_plots_from_ts.py \
              --ts-dir "$original_stage" \
              --pattern "*" \
              --burnin-frac 0 \
              --mutation-rate {params.mutation_rate} \
              --out-dir "{params.original_out}" \
              {params.sim_branch_flag} \
              >> "{log}" 2>&1
            touch "{output}"
            """


rule step7_summary:
    input:
        fai=str(FAI),
        targets=MERGED_TARGETS + STEP6_TARGETS,
    output:
        SUMMARY_TARGET
    threads: res_threads("step7_summary")
    resources:
        mem_mb=res_mem_mb("step7_summary"),
        time=res_time("step7_summary"),
        slurm_partition=res_partition("step7_summary"),
    log:
        str(LOG_DIR / "step7_summary.log")
    params:
        out_dir=str(OUT_DIR),
        chroms=" ".join(CHROMS),
        replicates=" ".join(REPLICATES),
        rec_fraction=REC_FRACTION,
        low_access_window=int(LOW_ACCESS_WINDOW_SIZE),
        low_access_cutoff=int(LOW_ACCESS_CUTOFF_BP),
        mutload_cutoff=MUTLOAD_CUTOFF,
        mutation_rate=VALIDATION_MUTATION_RATE if VALIDATION_MUTATION_RATE is not None else "null",
        sim_branch="true" if VALIDATION_SIM_BRANCH else "false",
    shell:
        """
        python scripts/pipeline_summary.py \
          --out-dir "{params.out_dir}" \
          --fai "{input.fai}" \
          --chroms {params.chroms} \
          --replicates {params.replicates} \
          --out "{output}" \
          --rec-fraction {params.rec_fraction} \
          --low-access-window {params.low_access_window} \
          --low-access-cutoff {params.low_access_cutoff} \
          --mutload-cutoff {params.mutload_cutoff} \
          --mutation-rate {params.mutation_rate} \
          --sim-branch {params.sim_branch} \
          >> "{log}" 2>&1
        """
