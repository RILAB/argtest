from pathlib import Path
import sys

import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from argtest_common import (  # noqa: E402
    audit_individual_contract,
    get_individual_name,
    load_ts,
    name_to_nodes_map,
    sample_names,
)
from export_vcf import resolve_samples  # noqa: E402


def make_ts(ploidies=(1, 1), names=None, unassigned=0, sample_parent=False):
    tables = tskit.TableCollection(sequence_length=10)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    names = names or [f"ind{i}" for i in range(len(ploidies))]
    sample_nodes = []
    for ploidy, name in zip(ploidies, names):
        ind = tables.individuals.add_row(metadata={"id": name})
        for _ in range(ploidy):
            node_time = 1 if sample_parent and not sample_nodes else 0
            sample_nodes.append(
                tables.nodes.add_row(
                    flags=tskit.NODE_IS_SAMPLE, time=node_time, individual=ind
                )
            )
    for _ in range(unassigned):
        sample_nodes.append(tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0))
    ancestor = tables.nodes.add_row(time=2)
    for node in sample_nodes:
        if not sample_parent or node != sample_nodes[0]:
            tables.edges.add_row(0, 10, parent=ancestor, child=node)
    if sample_parent:
        child = tables.nodes.add_row(time=0)
        tables.edges.add_row(0, 10, parent=sample_nodes[0], child=child)
        tables.edges.add_row(0, 10, parent=ancestor, child=sample_nodes[0])
    tables.sort()
    return tables.tree_sequence()


def test_name_normalization_is_global_and_shared():
    ts = make_ts(names=["pre_X_mid_X_tail", "plain"])
    substring = "_X"
    assert get_individual_name(
        ts.individual(0), name_substring_to_remove=substring
    ) == "pre_mid_tail"
    assert sample_names(ts, name_substring_to_remove=substring) == [
        "pre_mid_tail",
        "plain",
    ]
    assert set(name_to_nodes_map(ts, name_substring_to_remove=substring)) == {
        "pre_mid_tail",
        "plain",
    }
    _, vcf_names, _ = resolve_samples(ts, name_substring_to_remove=substring)
    assert vcf_names == ["pre_mid_tail", "plain"]


def test_empty_name_substring_leaves_name_unchanged():
    ts = make_ts(names=["sample_X", "other"])
    assert sample_names(ts, name_substring_to_remove="") == ["sample_X", "other"]


def test_audit_ignores_unrepresented_individuals_for_ploidy():
    ts = make_ts(ploidies=(1, 1, 0))
    assert audit_individual_contract(ts) == []


def test_bundled_100trees_ignores_unrepresented_individuals():
    path = Path(__file__).resolve().parents[1] / "example_data/test_100trees.tsz"
    findings = audit_individual_contract(load_ts(path))
    assert not any("mixed ploidy" in warning for warning in findings)


def test_audit_reports_all_proposed_contract_violations():
    ts = make_ts(
        ploidies=(1, 2),
        names=("sample_drop", "sample"),
        unassigned=1,
        sample_parent=True,
    )
    findings = audit_individual_contract(
        ts, name_substring_to_remove="_drop"
    )
    assert any("no individual" in warning for warning in findings)
    assert any("duplicate normalized" in warning for warning in findings)
    assert any("mixed ploidy" in warning for warning in findings)
    assert any("edge parents" in warning for warning in findings)
