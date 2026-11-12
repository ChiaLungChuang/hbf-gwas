"""Tests for the locus table.

The brief requires that every value stays `verified: false` until the user
signs it off by hand, and that a test fails the real-data-ready check while
any field is unverified. That test is `test_real_data_ready` below; it carries
the `realdata` marker so the default run (and CI's required job) skips it, and
`pytest -m realdata` runs it. It is expected to FAIL until sign-off.
"""

import dataclasses
import subprocess
import sys

import pandas as pd
import pytest

from hbfgwas import loci
from hbfgwas.loci import LOCI, Locus, Sourced

EXPECTED = {
    "BCL11A": ("rs1427407", "2p16"),
    "HBS1L-MYB": ("rs9399137", "6q23"),
    "HBB cluster": ("rs7482144", "11p15"),
}


def _all_verified(locus: Locus) -> Locus:
    """Copy of a locus with every field filled and verified (test-only)."""
    changes = {}
    for f in dataclasses.fields(Locus):
        field = getattr(locus, f.name)
        value = field.value if field.value is not None else "x"
        changes[f.name] = Sourced(
            value=value, source=field.source or "test", verified=True, verified_by="test"
        )
    return dataclasses.replace(locus, **changes)


def test_three_loci_with_the_brief_sentinels():
    assert len(LOCI) == 3
    got = {loc.gene.value: (loc.rsid.value, loc.locus.value) for loc in LOCI}
    assert got == EXPECTED


def test_every_field_is_sourced():
    for loc in LOCI:
        for f in dataclasses.fields(Locus):
            field = getattr(loc, f.name)
            assert isinstance(field, Sourced), f"{loc.gene.value}.{f.name}"
            if field.value is not None:
                assert field.source, f"{loc.gene.value}.{f.name} has a value but no source"


# The user's sign-off of 2026-11-12 (issue #3), as row IDs of docs/locus-signoff.md.
VERIFIED = {
    ("BCL11A", "chromosome"),
    ("BCL11A", "position_grch37"),
    ("BCL11A", "other_allele"),
    ("BCL11A", "effect_size"),
    ("BCL11A", "effect_allele_freq_cohort"),
    ("BCL11A", "effect_allele"),
    ("HBS1L-MYB", "chromosome"),
    ("HBS1L-MYB", "position_grch37"),
    ("HBS1L-MYB", "other_allele"),
    ("HBS1L-MYB", "effect_allele"),
    ("HBS1L-MYB", "effect_size"),
    ("HBS1L-MYB", "effect_size_scale"),
    ("HBS1L-MYB", "effect_allele_freq_cohort"),
    ("HBB cluster", "chromosome"),
    ("HBB cluster", "position_grch37"),
    ("HBB cluster", "other_allele"),
    ("HBB cluster", "effect_allele"),
    ("HBB cluster", "effect_size"),
    ("HBB cluster", "effect_size_scale"),
    ("HBB cluster", "effect_allele_freq_cohort"),
}  # S10-S12, S34, S37, S39, S14-S16, S22-S24, S26, S18-S20, S28-S30, S32
SUPERSEDED = {("BCL11A", "effect_allele_freq_1kg_afr"), ("HBS1L-MYB", "effect_allele_freq_1kg_afr")}


def test_signoff_record_is_exactly_what_was_checked():
    verified, superseded = set(), set()
    for loc in LOCI:
        for f in dataclasses.fields(Locus):
            field = getattr(loc, f.name)
            if field.verified:
                verified.add((loc.gene.value, f.name))
                assert field.verified_by, (loc.gene.value, f.name)
            if field.superseded_by:
                superseded.add((loc.gene.value, f.name))
    assert verified == VERIFIED
    assert superseded == SUPERSEDED


def test_attributions_do_not_credit_the_user_for_others_checks():
    by = {(g, f): getattr({x.gene.value: x for x in LOCI}[g], f).verified_by for g, f in VERIFIED}
    user = {k for k, v in by.items() if v.startswith("the user, by hand")}
    assert user == {("BCL11A", "effect_size"), ("BCL11A", "effect_allele_freq_cohort")}
    assert "patent" in by[("BCL11A", "effect_allele")]


def test_analysis_fields_are_ready_and_guarded():
    assert loci.analysis_ready(LOCI)
    broken = (dataclasses.replace(LOCI[0], effect_size=Sourced(0.076, "x")),) + LOCI[1:]
    assert not loci.analysis_ready(broken)


def test_a_null_value_cannot_be_verified():
    with pytest.raises(ValueError):
        Sourced(value=None, source="anything", verified=True)


def test_a_verified_value_needs_a_source():
    with pytest.raises(ValueError):
        Sourced(value=1, source=None, verified=True)


def test_chromosomes_match_the_cytobands_when_filled():
    for loc in LOCI:
        if loc.chromosome.value is not None:
            assert str(loc.locus.value).startswith(str(loc.chromosome.value))


def test_problems_list_every_unresolved_field():
    problems = loci.unverified_fields(LOCI)
    n_fields = len(dataclasses.fields(Locus))
    assert len(problems) == 3 * n_fields - len(VERIFIED) - len(SUPERSEDED) == 17
    assert not loci.real_data_ready(LOCI)


def test_ready_only_when_everything_is_verified():
    fully = tuple(_all_verified(loc) for loc in LOCI)
    assert loci.real_data_ready(fully)
    assert loci.unverified_fields(fully) == []
    one_missing = (dataclasses.replace(fully[0], effect_allele=Sourced(None, None)),) + fully[1:]
    assert not loci.real_data_ready(one_missing)
    with pytest.raises(loci.LociNotReadyError):
        loci.check_real_data_ready(one_missing)


def test_to_frame_has_values_and_flags():
    df = loci.to_frame(LOCI)
    assert isinstance(df, pd.DataFrame)
    assert list(df["rsid"]) == ["rs1427407", "rs9399137", "rs7482144"]
    assert "position_grch37_verified" in df.columns
    assert "position_grch37_source" in df.columns


def test_export_tsv(tmp_path):
    out = tmp_path / "loci.tsv"
    loci.export_tsv(out, LOCI)
    df = pd.read_csv(out, sep="\t")
    assert list(df["rsid"]) == ["rs1427407", "rs9399137", "rs7482144"]
    for col in ("gene", "chromosome", "position_grch37", "effect_allele", "all_verified"):
        assert col in df.columns
    assert not df["all_verified"].any()
    with pytest.raises(loci.LociNotReadyError):
        loci.export_tsv(tmp_path / "strict.tsv", LOCI, require_verified=True)


def test_cli_check_exits_nonzero_while_unverified():
    result = subprocess.run(
        [sys.executable, "-m", "hbfgwas.loci", "check"], capture_output=True, text=True
    )
    assert result.returncode == 1
    assert "not ready" in result.stdout.lower()


@pytest.mark.realdata
def test_real_data_ready():
    """Fails by design until the user verifies every value by hand."""
    problems = loci.unverified_fields(LOCI)
    assert not problems, "Locus table is not real-data ready:\n" + "\n".join(problems)
