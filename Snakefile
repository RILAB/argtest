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
OUT_DIR = Path(config.get("out_dir", "snakemake_out")).expanduser()
BASE_NAME = str(config.get("base_name", ROOT_DIR.name))
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

MUTLOAD_MUTATION_RATE = config.get("mutload_mutation_rate", None)
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

VALIDATION_MUTATION_RATE = config.get("validation_mutation_rate", None)
if VALIDATION_MUTATION_RATE is not None:
    VALIDATION_MUTATION_RATE = float(VALIDATION_MUTATION_RATE)
VALIDATION_FIRST_CHROM_ONLY = bool(config.get("validation_first_chrom_only", True))
VALIDATION_SIM_BRANCH = bool(config.get("validation_sim_branch", False))

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
STEP6_DIR = OUT_DIR / "step6_validation"
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
        MERGED_TARGETS + STEP6_TARGETS + [SUMMARY_TARGET]


rule step1_low_rec_masks:
    input:
        hapmap=str(HAPMAP),
        fai=str(FAI),
    output:
        str(STEP1_DIR / "{chrom}.low_rec.mask.bed")
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
    output:
        str(STEP5_DIR / "{chrom}" / "{rep}.{ext}")
    wildcard_constraints:
        rep="|".join(re.escape(r) for r in REPLICATES),
        ext="|".join(re.escape(s.lstrip(".")) for s in VALID_SUFFIXES),
    log:
        str(LOG_DIR / "step5" / "{chrom}" / "{rep}.{ext}.log")
    params:
        suffix_to_strip=SUFFIX_TO_STRIP,
    shell:
        """
        python scripts/trim_samples.py \
          "{input.ts}" \
          --remove "{input.outlier}" \
          --out "{output}" \
          --suffix-to-strip "{params.suffix_to_strip}" \
          --log "{log}"
        """


rule merge_replicates:
    input:
        step5_inputs_for_rep
    output:
        str(MERGED_DIR / f"{BASE_NAME}.combined.{{rep}}.{{ext}}")
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
              ln -s "$f" "$cleaned_stage/$(basename "$f")"
            done

            original_files=({input.original:q})
            for f in "${{original_files[@]}}"; do
              ln -s "$f" "$original_stage/$(basename "$f")"
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
