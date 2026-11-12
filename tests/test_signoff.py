"""The sign-off checklist covers every field once and matches the locus table."""

import dataclasses
from pathlib import Path

from hbfgwas import signoff
from hbfgwas.loci import LOCI, Locus

DOC = Path(__file__).resolve().parents[1] / "docs" / "locus-signoff.md"


def test_every_field_of_every_locus_appears_exactly_once():
    seen = [(gene, field) for _, _, gene, field, _, _ in signoff.rows()]
    expected = [(loc.gene.value, f.name) for loc in LOCI for f in dataclasses.fields(Locus)]
    assert sorted(seen) == sorted(expected)
    assert len(seen) == 39


def test_committed_checklist_matches_the_table():
    assert DOC.read_text() == signoff.markdown()


def test_status_column_matches_the_signoff_record():
    text = DOC.read_text()
    assert text.count("| ✅ |") == 20
    assert text.count("| ↪ |") == 2
    assert text.count("| ☐ |") == 17
    for rid, _, gene, field, _, _ in signoff.rows():
        f = getattr({loc.gene.value: loc for loc in LOCI}[gene], field)
        mark = "✅" if f.verified else ("↪" if f.superseded_by else "☐")
        assert f"| {mark} | {rid} |" in text, rid


def test_the_hbb_frequency_note_names_the_right_row():
    rows = {rid: (gene, field) for rid, _, gene, field, _, _ in signoff.rows()}
    assert rows["S21"] == ("HBB cluster", "effect_allele_freq_1kg_afr")
    assert "Note on S21" in DOC.read_text()


def test_every_group_has_a_link():
    for title, links, _ in signoff.GROUPS:
        assert "](" in links, title
