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


def test_nothing_is_marked_verified():
    text = DOC.read_text()
    assert "☑" not in text and "[x]" not in text.lower()
    assert text.count("| ☐ |") == 39
    assert all(not getattr(loc, f.name).verified for loc in LOCI for f in dataclasses.fields(Locus))


def test_every_group_has_a_link():
    for title, links, _ in signoff.GROUPS:
        assert "](" in links, title
