# Locus table sign-off checklist (issue #3)

Generated from `src/hbfgwas/loci.py` by `python -m hbfgwas.signoff`;
`tests/test_signoff.py` fails if this file and the table disagree.

**Sign-off recorded 2026-11-12 from the user's message** (`loci.SIGNOFF`):
20 rows verified (✅, with who checked them), 2 superseded (↪),
17 not covered (☐). Nothing is attributed to the user that the user did
not check: checks by the user's reviewer and by the step c run are labelled
as such. Report further checks by row ID, e.g. `S05 OK` or
`S12 wrong: <value> (<where>)`.

One row per field of the locus table (3 loci x 13 fields = 39 rows),
grouped by source so each page is opened once.

**Two checks.** The full check (`python -m hbfgwas.loci check`,
`pytest -m realdata`) needs every row verified or superseded and still
fails on the 17 rows marked ☐. The analysis check
(`check --analysis`) needs only the fields the real-data scripts use
(chromosome, position, effect and other allele, effect size, effect-allele
frequency in the source cohort) and passes; `03c` and `04` require it.

**Links checked 2026-11-11.** PMC, dbSNP and the Ensembl REST JSON links
load directly. The publisher DOIs and Europe PMC load in a browser (they
block scripts). The Ensembl GRCh37 *website* pages (Explore, Population
genetics) returned an Ensembl server error that day ("Server Exception",
a database error on their side); if they still do, the REST JSON link
for the same rsID holds the same data (`mappings` for location and
alleles, `populations` for 1000GENOMES:phase_3:AFR frequencies).

**Note on S21 (HBB cluster, `effect_allele_freq_1kg_afr`).** The value is
null because no 1000 Genomes phase 3 frequency exists: the step c run found
rs7482144 absent from the callset itself. (An earlier version of this
file called this row S24 by mistake; S24 is HBS1L-MYB's effect-size scale.)
The full check counts a null as missing, so it cannot pass with it; HBB is
"not testable in this callset" either way.

## 1. The brief (docs/BRIEF.md)

Open: [docs/BRIEF.md](BRIEF.md)

| status | ID | locus | field | current value | where to read it | sign-off |
|---|---|---|---|---|---|---|
| ☐ | S01 | BCL11A | `locus` | 2p16 | "What this is": "BCL11A 2p16, HBS1L-MYB 6q23, HBB cluster 11p15" | not covered by the sign-off |
| ☐ | S02 | BCL11A | `gene` | BCL11A | same sentence | not covered by the sign-off |
| ☐ | S03 | BCL11A | `rsid` | rs1427407 | Phase 1 item 2: sentinel variants rs1427407, rs9399137, rs7482144 (XmnI) | not covered by the sign-off |
| ☐ | S04 | HBS1L-MYB | `locus` | 6q23 | "What this is": "BCL11A 2p16, HBS1L-MYB 6q23, HBB cluster 11p15" | not covered by the sign-off |
| ☐ | S05 | HBS1L-MYB | `gene` | HBS1L-MYB | same sentence | not covered by the sign-off |
| ☐ | S06 | HBS1L-MYB | `rsid` | rs9399137 | Phase 1 item 2: sentinel variants rs1427407, rs9399137, rs7482144 (XmnI) | not covered by the sign-off |
| ☐ | S07 | HBB cluster | `locus` | 11p15 | "What this is": "BCL11A 2p16, HBS1L-MYB 6q23, HBB cluster 11p15" | not covered by the sign-off |
| ☐ | S08 | HBB cluster | `gene` | HBB cluster | same sentence | not covered by the sign-off |
| ☐ | S09 | HBB cluster | `rsid` | rs7482144 | Phase 1 item 2: sentinel variants rs1427407, rs9399137, rs7482144 (XmnI) | not covered by the sign-off |

## 2. Ensembl GRCh37 and dbSNP: rs1427407

Open: [Explore](https://grch37.ensembl.org/Homo_sapiens/Variation/Explore?v=rs1427407) · [Population genetics](https://grch37.ensembl.org/Homo_sapiens/Variation/Population?v=rs1427407) · [REST JSON](https://grch37.rest.ensembl.org/variation/human/rs1427407?content-type=application/json;pops=1) · [dbSNP](https://www.ncbi.nlm.nih.gov/snp/rs1427407)

| status | ID | locus | field | current value | where to read it | sign-off |
|---|---|---|---|---|---|---|
| ✅ | S10 | BCL11A | `chromosome` | 2 | Explore page: Location 'Chromosome 2:60718043'; REST `mappings[].location` = 2:60718043-60718043 | verified by the Phase 2 step c run (runs/2026-11-11/sentinels.tsv): found in the 1000 Genomes phase 3 file at this position with these alleles |
| ✅ | S11 | BCL11A | `position_grch37` | 60718043 | same; cross-check dbSNP 'GRCh37.p13' placement chr2:60718043 (NC_000002.11) | verified by the Phase 2 step c run (runs/2026-11-11/sentinels.tsv): found in the 1000 Genomes phase 3 file at this position with these alleles |
| ✅ | S12 | BCL11A | `other_allele` | G | Population genetics, 1000 Genomes phase 3 AFR: only T and G observed (dbSNP also lists C) | verified by the Phase 2 step c run (runs/2026-11-11/sentinels.tsv): found in the 1000 Genomes phase 3 file at this position with these alleles |
| ↪ | S13 | BCL11A | `effect_allele_freq_1kg_afr` | 0.237518910741301 | Population genetics, 1000GENOMES:phase_3:AFR, allele T: 0.2375 (REST: 0.237518910741301, count 314) | superseded by the value measured in runs/2026-11-11 (T = 0.2389, 651 founders); not verified against Ensembl |

## 3. Ensembl GRCh37 and dbSNP: rs9399137

Open: [Explore](https://grch37.ensembl.org/Homo_sapiens/Variation/Explore?v=rs9399137) · [Population genetics](https://grch37.ensembl.org/Homo_sapiens/Variation/Population?v=rs9399137) · [REST JSON](https://grch37.rest.ensembl.org/variation/human/rs9399137?content-type=application/json;pops=1) · [dbSNP](https://www.ncbi.nlm.nih.gov/snp/rs9399137)

| status | ID | locus | field | current value | where to read it | sign-off |
|---|---|---|---|---|---|---|
| ✅ | S14 | HBS1L-MYB | `chromosome` | 6 | Explore page: Location 'Chromosome 6:135419018'; REST `mappings[].location` | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources (GRCh37 6:135419018 T/C); and the Phase 2 step c run (runs/2026-11-11/sentinels.tsv): found in the 1000 Genomes phase 3 file at this position with these alleles |
| ✅ | S15 | HBS1L-MYB | `position_grch37` | 135419018 | same; cross-check dbSNP 'GRCh37.p13' placement chr6:135419018 (NC_000006.11) | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources (GRCh37 6:135419018 T/C); and the Phase 2 step c run (runs/2026-11-11/sentinels.tsv): found in the 1000 Genomes phase 3 file at this position with these alleles |
| ✅ | S16 | HBS1L-MYB | `other_allele` | T | Explore page: alleles T/C | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources (GRCh37 6:135419018 T/C); and the Phase 2 step c run (runs/2026-11-11/sentinels.tsv): found in the 1000 Genomes phase 3 file at this position with these alleles |
| ↪ | S17 | HBS1L-MYB | `effect_allele_freq_1kg_afr` | 0.0416036308623298 | Population genetics, 1000GENOMES:phase_3:AFR, allele C: 0.0416 (REST: 0.0416036308623298, count 55) | superseded by the value measured in runs/2026-11-11 (C = 0.0392, 651 founders); not verified against Ensembl |

## 4. Ensembl GRCh37 and dbSNP: rs7482144

Open: [Explore](https://grch37.ensembl.org/Homo_sapiens/Variation/Explore?v=rs7482144) · [Population genetics](https://grch37.ensembl.org/Homo_sapiens/Variation/Population?v=rs7482144) · [REST JSON](https://grch37.rest.ensembl.org/variation/human/rs7482144?content-type=application/json;pops=1) · [dbSNP](https://www.ncbi.nlm.nih.gov/snp/rs7482144)

| status | ID | locus | field | current value | where to read it | sign-off |
|---|---|---|---|---|---|---|
| ✅ | S18 | HBB cluster | `chromosome` | 11 | Explore page: Location 'Chromosome 11:5276169'; REST `mappings[].location` | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources (GRCh37 11:5276169 G/A) |
| ✅ | S19 | HBB cluster | `position_grch37` | 5276169 | same; cross-check dbSNP 'GRCh37.p13' placement chr11:5276169 (NC_000011.9) | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources (GRCh37 11:5276169 G/A) |
| ✅ | S20 | HBB cluster | `other_allele` | G | Explore page: alleles G/A | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources (GRCh37 11:5276169 G/A) |
| ☐ | S21 | HBB cluster | `effect_allele_freq_1kg_afr` | null (not found) | CONFIRM ABSENT: Population genetics should show no 1000 Genomes phase 3 rows (only gnomAD exomes AFR A = 0.1238, which is suspect). See the note at the top. | not covered by the sign-off |

## 5. Lettre et al. 2008, PNAS 105:11869, Table 1 (CSSCD columns)

Open: [PMC2491485](https://pmc.ncbi.nlm.nih.gov/articles/PMC2491485/) · [Europe PMC](https://europepmc.org/article/PMC/PMC2491485) · [doi:10.1073/pnas.0804799105](https://doi.org/10.1073/pnas.0804799105)

| status | ID | locus | field | current value | where to read it | sign-off |
|---|---|---|---|---|---|---|
| ✅ | S22 | HBS1L-MYB | `effect_allele` | C | Table 1, row rs9399137, MAF (allele) '0.06 (C)' with a positive effect; footnote †: 'Minor alleles (positive strand) are given in the parentheses' | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: Lettre 2008 Table 1, minor allele '0.06 (C)' with positive effect |
| ✅ | S23 | HBS1L-MYB | `effect_size` | 0.571 | Table 1, row rs9399137, Effect size (SE) '0.571 (0.086)'; variance explained 3.5% | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: Lettre 2008 Table 1, effect 0.571 |
| ✅ | S24 | HBS1L-MYB | `effect_size_scale` | per-allele effect of the minor allele in standard deviation units of HbF | Table 1 footnote ‡: 'Effect sizes and standard errors are given in standard deviation units for the minor allele' | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: Lettre 2008, effect sizes in SD units for the minor allele |
| ☐ | S25 | HBS1L-MYB | `effect_cohort` | CSSCD, African American sickle cell disease, N = 1,275 | Table 1 column header 'CSSCD (N = 1,275)'; Methods: African American sickle cell disease | not covered by the sign-off |
| ✅ | S26 | HBS1L-MYB | `effect_allele_freq_cohort` | 0.06 | Table 1, row rs9399137, MAF (allele) '0.06 (C)' | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: Lettre 2008 Table 1, MAF 0.06 (C) |
| ☐ | S27 | HBS1L-MYB | `citation` | Lettre G et al. (2008) DNA polymorphisms at the BCL11A, HBS1L-MYB, and beta-globin loci associate with fetal hemoglobin levels and pain crises in sickle cell disease. PNAS 105(33):11869-11874. PMID 18667698. doi:10.1073/pnas.0804799105. Table 1, CSSCD columns | title page: authors, title, journal, volume, pages, PMID, DOI | not covered by the sign-off |
| ✅ | S28 | HBB cluster | `effect_allele` | A | Table 1, row rs7482144, MAF (allele) '0.07 (A)' with a positive effect; footnote †: 'Minor alleles (positive strand) are given in the parentheses' | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: Lettre 2008 Table 1, minor allele '0.07 (A)' with positive effect |
| ✅ | S29 | HBB cluster | `effect_size` | 0.407 | Table 1, row rs7482144, Effect size (SE) '0.407 (0.080)'; variance explained 2.2% | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: Lettre 2008 Table 1, effect 0.407 |
| ✅ | S30 | HBB cluster | `effect_size_scale` | per-allele effect of the minor allele in standard deviation units of HbF | Table 1 footnote ‡: 'Effect sizes and standard errors are given in standard deviation units for the minor allele' | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: Lettre 2008, effect sizes in SD units for the minor allele |
| ☐ | S31 | HBB cluster | `effect_cohort` | CSSCD, African American sickle cell disease, N = 1,275 | Table 1 column header 'CSSCD (N = 1,275)'; Methods: African American sickle cell disease | not covered by the sign-off |
| ✅ | S32 | HBB cluster | `effect_allele_freq_cohort` | 0.07 | Table 1, row rs7482144, MAF (allele) '0.07 (A)' | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: Lettre 2008 Table 1, MAF 0.07 (A) |
| ☐ | S33 | HBB cluster | `citation` | Lettre G et al. (2008) DNA polymorphisms at the BCL11A, HBS1L-MYB, and beta-globin loci associate with fetal hemoglobin levels and pain crises in sickle cell disease. PNAS 105(33):11869-11874. PMID 18667698. doi:10.1073/pnas.0804799105. Table 1, CSSCD columns | title page: authors, title, journal, volume, pages, PMID, DOI | not covered by the sign-off |

## 6. Wonkam et al. 2025, Nat Commun 16:2092 (global meta-analysis)

Open: [PMC11873275](https://pmc.ncbi.nlm.nih.gov/articles/PMC11873275/) · [Europe PMC](https://europepmc.org/article/PMC/PMC11873275) · [doi:10.1038/s41467-025-57413-5](https://doi.org/10.1038/s41467-025-57413-5)

| status | ID | locus | field | current value | where to read it | sign-off |
|---|---|---|---|---|---|---|
| ✅ | S34 | BCL11A | `effect_size` | 0.076 | Results: 'rs1427407 therefore contributed a smaller proportion (7.6%) in HbF variability' | verified by the user, by hand, in Wonkam et al. 2025, on 2026-11-11: variance explained 7.6% found |
| ☐ | S35 | BCL11A | `effect_size_scale` | proportion of HbF variance explained (PVE); HbF cube-root transformed | Methods: HbF cube-root transformed; the 7.6% is a proportion of variance explained (PVE) | not covered by the sign-off |
| ☐ | S36 | BCL11A | `effect_cohort` | global meta-analysis of sickle cell anaemia in African ancestries, n = 3751 (Cameroon 827, Tanzania 884, African American 2040) | Methods / Table 1: global meta-analysis n = 3751 (Cameroon 827, Tanzania 884, African American 2040) | not covered by the sign-off |
| ✅ | S37 | BCL11A | `effect_allele_freq_cohort` | 0.24 | Results: 'MAFrs1427407 = 0.24' (the paper does not name the allele) | verified by the user, by hand, in Wonkam et al. 2025, on 2026-11-11: MAF 0.24 found |
| ☐ | S38 | BCL11A | `citation` | Wonkam A et al. (2025) FLT1 and other candidate fetal haemoglobin modifying loci in sickle cell disease in African ancestries. Nat Commun 16(1):2092. PMID 40025045. doi:10.1038/s41467-025-57413-5; effect allele from Bauer DE et al. (2013) An erythroid enhancer of BCL11A subject to genetic variation determines fetal hemoglobin level. Science 342(6155):253-257. PMID 24115442. doi:10.1126/science.1242088 | title page; the field also cites Bauer et al. 2013 for the effect allele | not covered by the sign-off |

## 7. Bauer et al. 2013, Science 342:253, Fig. 2B

Open: [PMC4018826](https://pmc.ncbi.nlm.nih.gov/articles/PMC4018826/) · [Europe PMC](https://europepmc.org/article/PMC/PMC4018826) · [doi:10.1126/science.1242088](https://doi.org/10.1126/science.1242088)

| status | ID | locus | field | current value | where to read it | sign-off |
|---|---|---|---|---|---|---|
| ✅ | S39 | BCL11A | `effect_allele` | T | Fig. 2B: the 'High-HbF' sequence carries T at rs1427407 (the 'Low-HbF' sequence carries G) | verified by the user's reviewer (Claude, in a separate chat) on 2026-11-11, against live sources: T is the HbF-raising allele, read from the authors' patent text describing the same work, not from the paper's Fig. 2B |
