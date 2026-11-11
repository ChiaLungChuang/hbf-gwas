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


# Sources, as looked up on 2026-11-10. Claude read each one; none is verified.
ENSEMBL37 = "Ensembl GRCh37 REST /variation/human/{rsid} (release 116), mappings"
ENSEMBL37_AFR = (
    "Ensembl GRCh37 REST /variation/human/{rsid}?pops=1, population 1000GENOMES:phase_3:AFR"
)
DBSNP37 = "NCBI dbSNP Variation API refsnp/{num}, GRCh37.p13 placement {acc} (SPDI 0-based {spdi})"
LETTRE2008 = (
    "Lettre G et al. (2008) DNA polymorphisms at the BCL11A, HBS1L-MYB, and beta-globin "
    "loci associate with fetal hemoglobin levels and pain crises in sickle cell disease. "
    "PNAS 105(33):11869-11874. PMID 18667698. doi:10.1073/pnas.0804799105"
)
LETTRE2008_T1 = LETTRE2008 + ". Table 1, CSSCD columns"
WONKAM2025 = (
    "Wonkam A et al. (2025) FLT1 and other candidate fetal haemoglobin modifying loci in "
    "sickle cell disease in African ancestries. Nat Commun 16(1):2092. PMID 40025045. "
    "doi:10.1038/s41467-025-57413-5"
)
BAUER2013 = (
    "Bauer DE et al. (2013) An erythroid enhancer of BCL11A subject to genetic variation "
    "determines fetal hemoglobin level. Science 342(6155):253-257. PMID 24115442. "
    "doi:10.1126/science.1242088"
)


LOCI: tuple[Locus, ...] = (
    Locus(
        locus=Sourced("2p16", BRIEF),
        gene=Sourced("BCL11A", BRIEF),
        rsid=Sourced("rs1427407", BRIEF),
        chromosome=Sourced("2", ENSEMBL37.format(rsid="rs1427407") + ": 2:60718043"),
        position_grch37=Sourced(
            60718043,
            ENSEMBL37.format(rsid="rs1427407")
            + "; "
            + DBSNP37.format(num=1427407, acc="NC_000002.11", spdi=60718042),
        ),
        effect_allele=Sourced(
            "T",
            BAUER2013 + ". Fig. 2B labels the High-HbF sequence with T at rs1427407",
            note="forward strand: T is the GRCh37 reference base (Ensembl allele string T/C/G)",
        ),
        other_allele=Sourced(
            "G",
            ENSEMBL37_AFR.format(rsid="rs1427407") + ": only T and G are observed",
            note="dbSNP also lists C; 1000 Genomes AFR carries only T and G",
        ),
        effect_allele_freq_1kg_afr=Sourced(
            0.237518910741301,
            ENSEMBL37_AFR.format(rsid="rs1427407") + ": T frequency (allele count 314)",
        ),
        effect_size=Sourced(
            0.076,
            WONKAM2025 + '. Results text: "rs1427407 therefore contributed a smaller '
            'proportion (7.6%) in HbF variability" (global meta-analysis)',
        ),
        effect_size_scale=Sourced(
            "proportion of HbF variance explained (PVE); HbF cube-root transformed",
            WONKAM2025,
            note="the paper's table does not name the allele its beta refers to",
        ),
        effect_cohort=Sourced(
            "global meta-analysis of sickle cell anaemia in African ancestries, n = 3751 "
            "(Cameroon 827, Tanzania 884, African American 2040)",
            WONKAM2025,
        ),
        effect_allele_freq_cohort=Sourced(
            0.24,
            WONKAM2025 + ': "MAFrs1427407 = 0.24" (minor allele; T is minor in AFR)',
        ),
        citation=Sourced(WONKAM2025 + "; effect allele from " + BAUER2013, "this table"),
    ),
    Locus(
        locus=Sourced("6q23", BRIEF),
        gene=Sourced("HBS1L-MYB", BRIEF),
        rsid=Sourced("rs9399137", BRIEF),
        chromosome=Sourced("6", ENSEMBL37.format(rsid="rs9399137") + ": 6:135419018"),
        position_grch37=Sourced(
            135419018,
            ENSEMBL37.format(rsid="rs9399137")
            + "; "
            + DBSNP37.format(num=9399137, acc="NC_000006.11", spdi=135419017),
        ),
        effect_allele=Sourced(
            "C",
            LETTRE2008_T1 + ': minor allele "0.06 (C)", positive effect; '
            '"Minor alleles (positive strand) are given in the parentheses"',
        ),
        other_allele=Sourced("T", ENSEMBL37.format(rsid="rs9399137") + ": allele string T/C"),
        effect_allele_freq_1kg_afr=Sourced(
            0.0416036308623298,
            ENSEMBL37_AFR.format(rsid="rs9399137") + ": C frequency (allele count 55)",
        ),
        effect_size=Sourced(
            0.571,
            LETTRE2008_T1 + ': "0.571 (0.086)"; variance explained 3.5%',
        ),
        effect_size_scale=Sourced(
            "per-allele effect of the minor allele in standard deviation units of HbF",
            LETTRE2008_T1
            + ' footnote: "Effect sizes and standard errors are given in standard deviation '
            'units for the minor allele"',
        ),
        effect_cohort=Sourced(
            "CSSCD, African American sickle cell disease, N = 1,275", LETTRE2008_T1
        ),
        effect_allele_freq_cohort=Sourced(0.06, LETTRE2008_T1 + ': "0.06 (C)"'),
        citation=Sourced(LETTRE2008_T1, "this table"),
    ),
    Locus(
        locus=Sourced("11p15", BRIEF),
        gene=Sourced("HBB cluster", BRIEF),
        rsid=Sourced("rs7482144", BRIEF),
        chromosome=Sourced("11", ENSEMBL37.format(rsid="rs7482144") + ": 11:5276169"),
        position_grch37=Sourced(
            5276169,
            ENSEMBL37.format(rsid="rs7482144")
            + "; "
            + DBSNP37.format(num=7482144, acc="NC_000011.9", spdi=5276168),
        ),
        effect_allele=Sourced(
            "A",
            LETTRE2008_T1 + ': minor allele "0.07 (A)", positive effect (positive strand)',
            note="XmnI site present; HBG2 -158 C>T on the HBG2 (reverse) strand",
        ),
        other_allele=Sourced("G", ENSEMBL37.format(rsid="rs7482144") + ": allele string G/A"),
        effect_allele_freq_1kg_afr=Sourced(
            None,
            None,
            note=(
                "NOT FOUND: Ensembl GRCh37 and GRCh38 records carry no 1000GENOMES "
                "phase_3 frequencies for rs7482144 (gnomADe:afr A = 0.1238 on GRCh37). "
                "Synonyms include rs41322148, so the 1000 Genomes call set may use another "
                "ID or lack the site: look it up by position in Phase 2"
            ),
        ),
        effect_size=Sourced(
            0.407,
            LETTRE2008_T1 + ': "0.407 (0.080)"; "In the CSSCD, rs7482144 explains 2.2% of '
            'the variation in HbF levels"',
        ),
        effect_size_scale=Sourced(
            "per-allele effect of the minor allele in standard deviation units of HbF",
            LETTRE2008_T1 + " footnote (as for rs9399137)",
        ),
        effect_cohort=Sourced(
            "CSSCD, African American sickle cell disease, N = 1,275", LETTRE2008_T1
        ),
        effect_allele_freq_cohort=Sourced(0.07, LETTRE2008_T1 + ': "0.07 (A)"'),
        citation=Sourced(LETTRE2008_T1, "this table"),
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
