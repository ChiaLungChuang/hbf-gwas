"""The three established HbF modifier loci and their sentinel variants.

Every value carries its source and a ``verified`` flag. Claude looked values up
from dbSNP, Ensembl GRCh37 and the primary papers (see each ``source``), but
**nothing is verified until the user checks it by hand** and flips
``verified=True`` field by field. A value that could not be found is ``None``
and can never be verified.

While any field is missing or unverified:

* ``real_data_ready()`` is False,
* ``python -m hbfgwas.loci check`` exits 1,
* ``pytest -m realdata`` fails, and
* ``scripts/04_assoc.sh`` refuses to run on this table.

The phenotype used with these loci is always simulated; 1000 Genomes has no
HbF phenotype.
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

BRIEF = "docs/BRIEF.md (specified by the user)"


@dataclass(frozen=True)
class Sourced:
    """A value, where it came from, and whether the user has verified it."""

    value: Any
    source: str | None
    verified: bool = False
    note: str | None = None

    def __post_init__(self) -> None:
        if self.verified and self.value is None:
            raise ValueError("a missing value cannot be verified")
        if self.verified and not self.source:
            raise ValueError("a verified value needs a source")


@dataclass(frozen=True)
class Locus:
    """One HbF modifier locus. Field order is the column order of exports."""

    locus: Sourced  # cytoband
    gene: Sourced
    rsid: Sourced
    chromosome: Sourced
    position_grch37: Sourced  # 1-based
    effect_allele: Sourced  # the HbF-increasing allele
    other_allele: Sourced
    effect_allele_freq_1kg_afr: Sourced  # 1000 Genomes phase 3 AFR, via Ensembl GRCh37
    effect_size: Sourced
    effect_size_scale: Sourced  # what effect_size means: units, transform, beta or R^2
    effect_cohort: Sourced  # cohort the effect size comes from
    effect_allele_freq_cohort: Sourced  # effect-allele frequency in that cohort
    citation: Sourced


def _missing() -> Sourced:
    return Sourced(value=None, source=None)


LOCI: tuple[Locus, ...] = (
    Locus(
        locus=Sourced("2p16", BRIEF),
        gene=Sourced("BCL11A", BRIEF),
        rsid=Sourced("rs1427407", BRIEF),
        chromosome=_missing(),
        position_grch37=_missing(),
        effect_allele=_missing(),
        other_allele=_missing(),
        effect_allele_freq_1kg_afr=_missing(),
        effect_size=_missing(),
        effect_size_scale=_missing(),
        effect_cohort=_missing(),
        effect_allele_freq_cohort=_missing(),
        citation=_missing(),
    ),
    Locus(
        locus=Sourced("6q23", BRIEF),
        gene=Sourced("HBS1L-MYB", BRIEF),
        rsid=Sourced("rs9399137", BRIEF),
        chromosome=_missing(),
        position_grch37=_missing(),
        effect_allele=_missing(),
        other_allele=_missing(),
        effect_allele_freq_1kg_afr=_missing(),
        effect_size=_missing(),
        effect_size_scale=_missing(),
        effect_cohort=_missing(),
        effect_allele_freq_cohort=_missing(),
        citation=_missing(),
    ),
    Locus(
        locus=Sourced("11p15", BRIEF),
        gene=Sourced("HBB cluster", BRIEF),
        rsid=Sourced("rs7482144", BRIEF),
        chromosome=_missing(),
        position_grch37=_missing(),
        effect_allele=_missing(),
        other_allele=_missing(),
        effect_allele_freq_1kg_afr=_missing(),
        effect_size=_missing(),
        effect_size_scale=_missing(),
        effect_cohort=_missing(),
        effect_allele_freq_cohort=_missing(),
        citation=_missing(),
    ),
)

FIELD_NAMES: tuple[str, ...] = tuple(f.name for f in dataclasses.fields(Locus))


class LociNotReadyError(RuntimeError):
    """Raised when the locus table is used for real data before sign-off."""


def unverified_fields(table: tuple[Locus, ...] = LOCI) -> list[str]:
    """One line per field that is missing or not yet verified by the user."""
    problems: list[str] = []
    for loc in table:
        label = loc.gene.value or loc.rsid.value or "?"
        for name in FIELD_NAMES:
            field: Sourced = getattr(loc, name)
            if field.value is None:
                problems.append(f"{label}.{name}: value missing")
            elif not field.verified:
                problems.append(f"{label}.{name}: not verified ({field.value!r})")
    return problems


def real_data_ready(table: tuple[Locus, ...] = LOCI) -> bool:
    return not unverified_fields(table)


def check_real_data_ready(table: tuple[Locus, ...] = LOCI) -> None:
    problems = unverified_fields(table)
    if problems:
        raise LociNotReadyError(
            f"Locus table not ready for real data ({len(problems)} problems):\n"
            + "\n".join(problems)
        )


def to_frame(table: tuple[Locus, ...] = LOCI) -> pd.DataFrame:
    """Wide table: each field as ``name``, ``name_source``, ``name_verified``."""
    rows = []
    for loc in table:
        row: dict[str, Any] = {}
        for name in FIELD_NAMES:
            field: Sourced = getattr(loc, name)
            row[name] = field.value
            row[f"{name}_source"] = field.source
            row[f"{name}_verified"] = field.verified
        rows.append(row)
    return pd.DataFrame(rows)


def export_tsv(
    path: str | Path, table: tuple[Locus, ...] = LOCI, require_verified: bool = False
) -> Path:
    """Write the values (no sources) plus an ``all_verified`` column, for the scripts."""
    if require_verified:
        check_real_data_ready(table)
    df = to_frame(table)
    values = df[list(FIELD_NAMES)].copy()
    values["all_verified"] = [
        all(getattr(loc, name).verified for name in FIELD_NAMES) for loc in table
    ]
    path = Path(path)
    values.to_csv(path, sep="\t", index=False, na_rep="NA")
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m hbfgwas.loci", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="exit 1 unless every value is present and verified")
    sub.add_parser("show", help="print values, sources and flags")
    exp = sub.add_parser("export", help="write a TSV of values for the shell scripts")
    exp.add_argument("--out", required=True)
    exp.add_argument("--require-verified", action="store_true")
    args = parser.parse_args(argv)

    if args.cmd == "check":
        problems = unverified_fields()
        if problems:
            print(f"Locus table NOT READY for real data: {len(problems)} problems")
            print("\n".join(f"  {p}" for p in problems))
            return 1
        print("Locus table ready: every value present and verified.")
        return 0
    if args.cmd == "show":
        with pd.option_context("display.max_columns", None, "display.width", 200):
            print(to_frame().T.to_string())
        return 0
    try:
        export_tsv(args.out, require_verified=args.require_verified)
    except LociNotReadyError as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
