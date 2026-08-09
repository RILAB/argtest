from __future__ import annotations

import importlib.util
import os
import pickle
import subprocess
import sys
from pathlib import Path

import msprime
import pytest
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from argtest_common import load_ts


REPO_ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO_ROOT / "Snakefile"
SNK_AVAILABLE = importlib.util.find_spec("snakemake") is not None
SNK_SKIP_REASON = "snakemake is not installed in this environment"


def make_tree_sequence(outlier_name: str, length: int = 4000) -> tskit.TreeSequence:
    tables = tskit.TableCollection(sequence_length=length)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()

    names = ["A", "B", "C", "D"]
    node_ids = {}
    for name in names:
        ind_id = tables.individuals.add_row(metadata={"id": name})
        node_id = tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE,
            time=0,
            individual=ind_id,
            population=pop,
        )
        node_ids[name] = node_id

    anc = tables.nodes.add_row(time=1, population=pop)
    for node_id in node_ids.values():
        tables.edges.add_row(left=0, right=length, parent=anc, child=node_id)

    def add_site(position: float, sample_name: str):
        site = tables.sites.add_row(position=position, ancestral_state="0")
        tables.mutations.add_row(site=site, node=node_ids[sample_name], derived_state="1")

    # Baseline windows that survive trimming.
    for position, sample_name in [
        (100, "A"),
        (200, "B"),
        (300, "C"),
        (400, "D"),
        (3100, "A"),
        (3200, "B"),
        (3300, "C"),
        (3400, "D"),
    ]:
        add_site(position, sample_name)

    # Extra mutations in the low-accessibility windows so step 3 flags them as outliers.
    for position, sample_name in [
        (1100, "A"),
        (1200, "B"),
        (1300, "C"),
        (1400, "D"),
        (1500, outlier_name),
        (1600, outlier_name),
        (2100, "A"),
        (2200, "B"),
        (2300, "C"),
        (2400, "D"),
        (2500, outlier_name),
        (2600, outlier_name),
    ]:
        add_site(position, sample_name)

    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def write_mut_rate_map(path: Path):
    mu = msprime.RateMap(
        position=[0, 800, 3000, 4000],
        rate=[1.0, 0.0, 1.0],
    )
    with open(path, "wb") as fh:
        pickle.dump(mu, fh)


def write_hapmap(path: Path, chrom: str):
    rows = [
        f"{chrom}\t0\t1.0\t0.0",
        f"{chrom}\t500\t0.1\t0.1",
        f"{chrom}\t700\t1.2\t0.2",
        f"{chrom}\t3600\t0.05\t0.3",
    ]
    path.write_text(
        "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n" + "\n".join(rows) + "\n"
    )


def write_mask_file(path: Path, chrom: str):
    path.write_text(
        f"{chrom}\t50\t75\tmask\n"
        f"{chrom}\t3800\t3900\tmask\n"
    )


def build_dataset(root: Path) -> dict[str, Path]:
    tree_root = root / "trees"
    tree_root.mkdir(parents=True)
    out_dir = root / "out"
    out_dir.mkdir()
    hapmap = root / "hapmap.tsv"
    fai = root / "genome.fai"

    chroms = ["chr1", "chr2"]
    replicates = [1, 2]
    outlier_targets = {
        ("chr1", 1): "A",
        ("chr1", 2): "C",
        ("chr2", 1): "B",
        ("chr2", 2): "D",
    }

    fai.write_text("\n".join(f"{chrom}\t4000\t0\t0\t0" for chrom in chroms) + "\n")

    hapmap_rows = ["Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)"]
    for chrom in chroms:
        write_hapmap(root / f"{chrom}.hapmap.tsv", chrom)
        write_mask_file(root / f"{chrom}.mask.bed", chrom)
        (tree_root / chrom).mkdir(parents=True, exist_ok=True)
        write_mut_rate_map(tree_root / chrom / f"{chrom}.mut_rate.p")
        for rep in replicates:
            ts = make_tree_sequence(outlier_targets[(chrom, rep)])
            ts.dump(tree_root / chrom / f"{rep}.trees")

    for chrom in chroms:
        chrom_hapmap = root / f"{chrom}.hapmap.tsv"
        hapmap_rows.extend(chrom_hapmap.read_text().splitlines()[1:])
    hapmap.write_text("\n".join(hapmap_rows) + "\n")

    config = {
        "root_dir": str(tree_root),
        "hapmap": str(hapmap),
        "fai": str(fai),
        "tree_pattern": "*.trees",
        "rec_fraction": 0.5,
        "low_access_window_size": 1000,
        "low_access_cutoff_bp": 600,
        "mutload_window_size": 1000,
        # Sim-based expected per individual >> observed on this tiny fixture, so
        # set the cutoff wide enough that pipeline plumbing is tested without
        # mutload masking everything. Outlier-detection accuracy is covered by
        # the unit tests in tests/test_mutload_masks.py.
        "mutload_cutoff": 2.0,
        "name_substring_to_remove": "_anchorwave",
        "allow_missing_replicates": False,
        "base_name": "demo",
        "out_dir": str(out_dir),
    }
    config_path = root / "snakemake.yaml"
    config_lines = []
    for key, value in config.items():
        if isinstance(value, str):
            config_lines.append(f"{key}: {value!r}")
        else:
            config_lines.append(f"{key}: {value}")
    config_path.write_text("\n".join(config_lines) + "\n")

    return {
        "root": root,
        "tree_root": tree_root,
        "out_dir": out_dir,
        "hapmap": hapmap,
        "fai": fai,
        "config": config_path,
    }


def run_snakemake(repo_root: Path, config_path: Path, *extra_args: str):
    cache_root = config_path.parent / ".snakemake-test-cache"
    tmp_root = config_path.parent / ".snakemake-test-tmp"
    cache_root.mkdir(exist_ok=True)
    tmp_root.mkdir(exist_ok=True)
    env = os.environ.copy()
    env.update({
        "XDG_CACHE_HOME": str(cache_root),
        "TMPDIR": str(tmp_root),
    })
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(config_path),
        *extra_args,
    ]
    return subprocess.run(
        cmd,
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )


@pytest.mark.skipif(not SNK_AVAILABLE, reason=SNK_SKIP_REASON)
def test_snakemake_dry_run(tmp_path):
    dataset = build_dataset(tmp_path)
    result = run_snakemake(REPO_ROOT, dataset["config"], "-n", "-p")

    assert "rule all" in result.stdout
    assert "demo.combined.1.trees" in result.stdout
    assert "step1_low_rec_masks" in result.stdout


@pytest.mark.skipif(not SNK_AVAILABLE, reason=SNK_SKIP_REASON)
def test_snakemake_ratemap_only_validation_sim_branch_dry_run(tmp_path):
    dataset = build_dataset(tmp_path)
    with open(dataset["config"], "a") as fh:
        fh.write("validation_sim_branch: True\n")
        fh.write("validation_first_chrom_only: True\n")

    result = run_snakemake(REPO_ROOT, dataset["config"], "-n", "-p")

    assert "step6_validation_plots" in result.stdout
    assert "step6_validation/chr1/done.txt" in result.stdout
    assert "--window-size 100000" in result.stdout


@pytest.mark.skipif(not SNK_AVAILABLE, reason=SNK_SKIP_REASON)
def test_snakemake_real_run(tmp_path):
    dataset = build_dataset(tmp_path)
    result = run_snakemake(REPO_ROOT, dataset["config"], "--cores", "1", "--rerun-incomplete")

    combined_dir = dataset["out_dir"] / "combined"
    combined = sorted(combined_dir.glob("demo.combined.*.trees"))
    assert len(combined) == 2
    assert (dataset["out_dir"] / "step1_low_rec" / "chr1.low_rec.mask.bed").exists()
    assert (dataset["out_dir"] / "step2_low_access" / "chr1" / "chr1.low_access.bed").exists()
    assert (dataset["out_dir"] / "step3_mutload" / "chr1" / "1.outliers.bed").exists()
    assert (dataset["out_dir"] / "step4_masks" / "chr1" / "1.remove_regions.bed").exists()
    assert (dataset["out_dir"] / "step5_trimmed_samples" / "chr1" / "1.trees").exists()

    merged_ts = load_ts(combined[0])
    assert 0 < merged_ts.sequence_length <= 8000  # 2 chromosomes × 4000 bp each
    assert merged_ts.num_sites > 0
