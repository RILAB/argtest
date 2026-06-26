import gzip
import sys
from pathlib import Path
from types import SimpleNamespace

import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import export_vcf


def make_ts_with_isolated_sample(length=10):
    """3 haploid individuals (1 sample node each). Sample n2 is connected only
    over [0, 5), so it is isolated over [5, 10) — mimicking a trim_samples prune.
    Two sites: pos 2 (mutation on n0, in the fully-connected region) and pos 7
    (mutation on n1, where n2 is isolated)."""
    tables = tskit.TableCollection(sequence_length=length)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    iA = tables.individuals.add_row(metadata={"id": "indA"})
    iB = tables.individuals.add_row(metadata={"id": "indB"})
    iC = tables.individuals.add_row(metadata={"id": "indC"})
    n0 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop, individual=iA)
    n1 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop, individual=iB)
    n2 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop, individual=iC)
    anc = tables.nodes.add_row(time=1, population=pop)
    tables.edges.add_row(left=0, right=length, parent=anc, child=n0)
    tables.edges.add_row(left=0, right=length, parent=anc, child=n1)
    tables.edges.add_row(left=0, right=5, parent=anc, child=n2)  # isolated over [5, 10)
    s1 = tables.sites.add_row(position=2, ancestral_state="0")
    tables.mutations.add_row(site=s1, node=n0, derived_state="1")
    s2 = tables.sites.add_row(position=7, ancestral_state="0")
    tables.mutations.add_row(site=s2, node=n1, derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def _run(ts, tmp_path, out_name="out.vcf"):
    ts_path = tmp_path / "in.trees"
    ts.dump(ts_path)
    out_path = tmp_path / out_name
    args = SimpleNamespace(
        ts=ts_path, out=out_path, chrom="chr1", suffix_to_strip="", log=None
    )
    orig = export_vcf.parse_args
    export_vcf.parse_args = lambda: args
    try:
        export_vcf.main()
    finally:
        export_vcf.parse_args = orig
    return out_path


def _read(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return fh.read()


def _parse_records(text):
    header = None
    records = []
    for line in text.splitlines():
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.split("\t")
            continue
        if line.strip():
            records.append(line.split("\t"))
    return header, records


def test_haploid_gt_and_isolated_missing(tmp_path):
    ts = make_ts_with_isolated_sample()
    out = _run(ts, tmp_path)
    header, records = _parse_records(_read(out))

    # Sample columns (from index 9) carry the individual names.
    assert header[0] == "#CHROM"
    assert header[9:] == ["indA", "indB", "indC"]
    assert header[8] == "FORMAT"

    # CHROM is the contig we passed.
    assert all(rec[0] == "chr1" for rec in records)

    # Two variable sites.
    assert len(records) == 2

    def gts(rec):
        # FORMAT is GT only -> sample fields are the bare genotype strings.
        return rec[9], rec[10], rec[11]

    # Identify the site in the n2-isolated region by its missing genotype.
    isolated_site = [rec for rec in records if "." in gts(rec)]
    connected_site = [rec for rec in records if "." not in gts(rec)]
    assert len(isolated_site) == 1
    assert len(connected_site) == 1

    # Haploid coding: single allele per sample (not "0/."), missing == ".".
    iso_a, iso_b, iso_c = gts(isolated_site[0])
    assert iso_c == "."            # pruned/isolated sample -> missing
    assert iso_b == "1"            # carrier of the mutation
    assert iso_a == "0"
    for g in (iso_a, iso_b):
        assert "/" not in g and "|" not in g  # haploid, single allele

    con_a, con_b, con_c = gts(connected_site[0])
    assert (con_a, con_b, con_c) == ("1", "0", "0")


def test_gzip_output_roundtrips(tmp_path):
    ts = make_ts_with_isolated_sample()
    out = _run(ts, tmp_path, out_name="out.vcf.gz")
    text = _read(out)  # decompresses; would raise if not valid gzip
    _, records = _parse_records(text)
    assert len(records) == 2
