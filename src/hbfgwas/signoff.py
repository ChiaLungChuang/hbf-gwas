"""Generate docs/locus-signoff.md: a checklist for the user's hand sign-off (issue #3).

One row per field of the locus table (src/hbfgwas/loci.py), with the current
value, exactly where in its source to read it, and a link to that page. Rows
are grouped by source page so each page is opened once. Nothing here marks
anything verified; the user checks each row and reports the result.

    python -m hbfgwas.signoff --out docs/locus-signoff.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from hbfgwas.loci import FIELD_NAMES, LOCI

ENS = "https://grch37.ensembl.org/Homo_sapiens/Variation"
REST = "https://grch37.rest.ensembl.org/variation/human"
DBSNP = "https://www.ncbi.nlm.nih.gov/snp"


def _ensembl_links(rsid: str) -> str:
    return (
        f"[Explore]({ENS}/Explore?v={rsid}) · [Population genetics]({ENS}/Population?v={rsid}) · "
        f"[REST JSON]({REST}/{rsid}?content-type=application/json;pops=1) · "
        f"[dbSNP]({DBSNP}/{rsid})"
    )


LETTRE = (
    "[PMC2491485](https://pmc.ncbi.nlm.nih.gov/articles/PMC2491485/) · "
    "[Europe PMC](https://europepmc.org/article/PMC/PMC2491485) · "
    "[doi:10.1073/pnas.0804799105](https://doi.org/10.1073/pnas.0804799105)"
)
WONKAM = (
    "[PMC11873275](https://pmc.ncbi.nlm.nih.gov/articles/PMC11873275/) · "
    "[Europe PMC](https://europepmc.org/article/PMC/PMC11873275) · "
    "[doi:10.1038/s41467-025-57413-5](https://doi.org/10.1038/s41467-025-57413-5)"
)
BAUER = (
    "[PMC4018826](https://pmc.ncbi.nlm.nih.gov/articles/PMC4018826/) · "
    "[Europe PMC](https://europepmc.org/article/PMC/PMC4018826) · "
    "[doi:10.1126/science.1242088](https://doi.org/10.1126/science.1242088)"
)

# Source groups, in reading order: (title, links, [(gene, field, where to read it)]).
GROUPS = [
    (
        "The brief (docs/BRIEF.md)",
        "[docs/BRIEF.md](BRIEF.md)",
        [
            (g, f, w)
            for g in ("BCL11A", "HBS1L-MYB", "HBB cluster")
            for f, w in (
                ("locus", '"What this is": "BCL11A 2p16, HBS1L-MYB 6q23, HBB cluster 11p15"'),
                ("gene", "same sentence"),
                (
                    "rsid",
                    "Phase 1 item 2: sentinel variants rs1427407, rs9399137, rs7482144 (XmnI)",
                ),
            )
        ],
    ),
    (
        "Ensembl GRCh37 and dbSNP: rs1427407",
        _ensembl_links("rs1427407"),
        [
            (
                "BCL11A",
                "chromosome",
                "Explore page: Location 'Chromosome 2:60718043'; REST `mappings[].location` = 2:60718043-60718043",
            ),
            (
                "BCL11A",
                "position_grch37",
                "same; cross-check dbSNP 'GRCh37.p13' placement chr2:60718043 (NC_000002.11)",
            ),
            (
                "BCL11A",
                "other_allele",
                "Population genetics, 1000 Genomes phase 3 AFR: only T and G observed (dbSNP also lists C)",
            ),
            (
                "BCL11A",
                "effect_allele_freq_1kg_afr",
                "Population genetics, 1000GENOMES:phase_3:AFR, allele T: 0.2375 (REST: 0.237518910741301, count 314)",
            ),
        ],
    ),
    (
        "Ensembl GRCh37 and dbSNP: rs9399137",
        _ensembl_links("rs9399137"),
        [
            (
                "HBS1L-MYB",
                "chromosome",
                "Explore page: Location 'Chromosome 6:135419018'; REST `mappings[].location`",
            ),
            (
                "HBS1L-MYB",
                "position_grch37",
                "same; cross-check dbSNP 'GRCh37.p13' placement chr6:135419018 (NC_000006.11)",
            ),
            ("HBS1L-MYB", "other_allele", "Explore page: alleles T/C"),
            (
                "HBS1L-MYB",
                "effect_allele_freq_1kg_afr",
                "Population genetics, 1000GENOMES:phase_3:AFR, allele C: 0.0416 (REST: 0.0416036308623298, count 55)",
            ),
        ],
    ),
    (
        "Ensembl GRCh37 and dbSNP: rs7482144",
        _ensembl_links("rs7482144"),
        [
            (
                "HBB cluster",
                "chromosome",
                "Explore page: Location 'Chromosome 11:5276169'; REST `mappings[].location`",
            ),
            (
                "HBB cluster",
                "position_grch37",
                "same; cross-check dbSNP 'GRCh37.p13' placement chr11:5276169 (NC_000011.9)",
            ),
            ("HBB cluster", "other_allele", "Explore page: alleles G/A"),
            (
                "HBB cluster",
                "effect_allele_freq_1kg_afr",
                "CONFIRM ABSENT: Population genetics should show no 1000 Genomes phase 3 rows (only gnomAD exomes AFR A = 0.1238, which is suspect). See the note at the top.",
            ),
        ],
    ),
    (
        "Lettre et al. 2008, PNAS 105:11869, Table 1 (CSSCD columns)",
        LETTRE,
        [
            (
                g,
                f,
                w.replace("{rs}", rs)
                .replace("{maf}", maf)
                .replace("{eff}", eff)
                .replace("{ve}", ve),
            )
            for g, rs, maf, eff, ve in (
                ("HBS1L-MYB", "rs9399137", "0.06 (C)", "0.571 (0.086)", "3.5"),
                ("HBB cluster", "rs7482144", "0.07 (A)", "0.407 (0.080)", "2.2"),
            )
            for f, w in (
                (
                    "effect_allele",
                    "Table 1, row {rs}, MAF (allele) '{maf}' with a positive effect; footnote †: 'Minor alleles (positive strand) are given in the parentheses'",
                ),
                (
                    "effect_size",
                    "Table 1, row {rs}, Effect size (SE) '{eff}'; variance explained {ve}%",
                ),
                (
                    "effect_size_scale",
                    "Table 1 footnote ‡: 'Effect sizes and standard errors are given in standard deviation units for the minor allele'",
                ),
                (
                    "effect_cohort",
                    "Table 1 column header 'CSSCD (N = 1,275)'; Methods: African American sickle cell disease",
                ),
                ("effect_allele_freq_cohort", "Table 1, row {rs}, MAF (allele) '{maf}'"),
                ("citation", "title page: authors, title, journal, volume, pages, PMID, DOI"),
            )
        ],
    ),
    (
        "Wonkam et al. 2025, Nat Commun 16:2092 (global meta-analysis)",
        WONKAM,
        [
            (
                "BCL11A",
                "effect_size",
                "Results: 'rs1427407 therefore contributed a smaller proportion (7.6%) in HbF variability'",
            ),
            (
                "BCL11A",
                "effect_size_scale",
                "Methods: HbF cube-root transformed; the 7.6% is a proportion of variance explained (PVE)",
            ),
            (
                "BCL11A",
                "effect_cohort",
                "Methods / Table 1: global meta-analysis n = 3751 (Cameroon 827, Tanzania 884, African American 2040)",
            ),
            (
                "BCL11A",
                "effect_allele_freq_cohort",
                "Results: 'MAFrs1427407 = 0.24' (the paper does not name the allele)",
            ),
            (
                "BCL11A",
                "citation",
                "title page; the field also cites Bauer et al. 2013 for the effect allele",
            ),
        ],
    ),
    (
        "Bauer et al. 2013, Science 342:253, Fig. 2B",
        BAUER,
        [
            (
                "BCL11A",
                "effect_allele",
                "Fig. 2B: the 'High-HbF' sequence carries T at rs1427407 (the 'Low-HbF' sequence carries G)",
            ),
        ],
    ),
]


def _cell(value) -> str:
    text = "null (not found)" if value is None else str(value)
    return text.replace("|", "\\|").replace("\n", " ")


def rows():
    """(row id, group index, gene, field, value, where) in document order."""
    table = {loc.gene.value: loc for loc in LOCI}
    out, k = [], 0
    for gi, (_, _, items) in enumerate(GROUPS):
        for gene, field, where in items:
            k += 1
            value = getattr(table[gene], field).value
            out.append((f"S{k:02d}", gi, gene, field, value, where))
    return out


def markdown() -> str:
    lines = [
        "# Locus table sign-off checklist (issue #3)",
        "",
        "Generated from `src/hbfgwas/loci.py` by `python -m hbfgwas.signoff`;",
        "`tests/test_signoff.py` fails if this file and the table disagree.",
        "**Nothing here is verified.** Check each row against its source and tell",
        "Claude the result by row ID, e.g. `S05 OK`, or `S12 wrong: <value> (<where>)`.",
        "Only you mark values `verified: true` (D-007).",
        "",
        "One row per field of the locus table (3 loci x 13 fields = 39 rows),",
        "grouped by source so each page is opened once.",
        "",
        "**Links checked 2026-11-11.** PMC, dbSNP and the Ensembl REST JSON links",
        "load directly. The publisher DOIs and Europe PMC load in a browser (they",
        "block scripts). The Ensembl GRCh37 *website* pages (Explore, Population",
        'genetics) returned an Ensembl server error that day ("Server Exception",',
        "a database error on their side); if they still do, the REST JSON link",
        "for the same rsID holds the same data (`mappings` for location and",
        "alleles, `populations` for 1000GENOMES:phase_3:AFR frequencies).",
        "",
        "**Note on S24 (HBB cluster, `effect_allele_freq_1kg_afr`).** The value is",
        "null because no 1000 Genomes phase 3 frequency was found. The current",
        'readiness check treats a null as "missing" and can never pass with it.',
        "If you confirm the absence, decide how that should count (e.g. allow a",
        "verified absence); Claude will not change the check without your decision.",
        "On real data, the frequency is computed from the callset itself, and HBB",
        'is "not testable in this callset" if rs7482144 is absent or fails QC.',
        "",
    ]
    by_group = {}
    for r in rows():
        by_group.setdefault(r[1], []).append(r)
    for gi, (title, links, _) in enumerate(GROUPS):
        lines += [f"## {gi + 1}. {title}", "", f"Open: {links}", ""]
        lines += [
            "| ✓ | ID | locus | field | current value | where to read it |",
            "|---|---|---|---|---|---|",
        ]
        for rid, _, gene, field, value, where in by_group[gi]:
            lines.append(f"| ☐ | {rid} | {gene} | `{field}` | {_cell(value)} | {_cell(where)} |")
        lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m hbfgwas.signoff", description=__doc__)
    parser.add_argument("--out", default=None)
    args = parser.parse_args(argv)
    text = markdown()
    if args.out:
        Path(args.out).write_text(text)
        print(f"wrote {args.out} ({len(rows())} rows)")
    else:
        sys.stdout.write(text)
    return 0


assert len(FIELD_NAMES) == 13

if __name__ == "__main__":
    sys.exit(main())
