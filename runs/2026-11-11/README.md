# Run 2026-11-11: Phase 2 step c on real 1000 Genomes AFR genotypes

Steps a-c only (download, AFR subset, QC, PCA, PC plot). **No phenotype has
been simulated on these genotypes and no association has been run.** Step d
(pre-registration of the real-data criteria) has not started. Genotype files
stay in `data/1kg/` (git-ignored); this folder holds the run's documentation.

- Code: branch `phase2-step-c` (built on #27 and #28), plink2
  `PLINK v2.0.0-b.1.4 M1 (7 Nov 2026)` (`plink2_version.txt`, pinned in
  `config/analysis.yaml`).
- Download URLs: the four defaults, re-checked against the plink2 resources
  page on 2026-11-11 (unchanged). Files, sizes and sha256: `download_manifest.tsv`.
- Every plink2 `.log` and each script's stdout/stderr: `logs/`.

## Counts at each step

| step | samples | variants | source |
|---|---|---|---|
| all_phase3 (input) | 2504 | 84,805,772 | `logs/afr_step1.log` line 22 ("84805772 variants in"); the counts table shows NA because the script's log pattern missed this wording (fixed in this PR) |
| SuperPop AFR | 661 | | `01_subset_counts.tsv` |
| minus the 1000 Genomes related list (9 AFR: 6 ASW, 2 LWK, 1 ACB); autosomes; biallelic ACGT SNPs | 652 | 77,818,345 | `01_subset_counts.tsv` |
| duplicates removed (162 records at shared positions, all records of 191 duplicated IDs) | 652 | 77,817,955 | `01_subset_counts.tsv` |
| `--geno 0.02` | 652 | 77,817,955 (0 removed) | `qc_filters.tsv` |
| `--mind 0.02` | 652 (0 removed) | 77,817,955 | `qc_filters.tsv` |
| `--maf 0.01` | 652 | 14,739,729 (63,078,226 removed) | `qc_filters.tsv` |
| `--hwe 1e-6 midp keep-fewhet` | 652 | 14,698,808 (40,921 removed) | `qc_filters.tsv` |
| PCA variants: outside long-range LD regions, MAF >= 0.05 | 652 | 8,078,843 | `03_pca_counts.tsv` |
| LD-pruned (`--indep-pairwise 1000kb 1 0.1`) | 652 | 215,522 | `03_pca_counts.tsv` |
| KING `--king-cutoff 0.0884` | 652 (0 removed) | 215,522 | `03_pca_counts.tsv`, `king_removed.id` |

Sample count 652 lies in the pre-registered 640-665. 1000 Genomes phase 3 is
fully called, so `--geno` and `--mind` removed nothing, as expected; the
script printed that note. KING removed nobody, and the script said so.

## Sentinels (`sentinels.tsv`)

| locus | sentinel | found | GRCh37 position, alleles | effect allele frequency in this AFR sample | survived QC |
|---|---|---|---|---|---|
| BCL11A | rs1427407 | by rsID | 2:60718043, T/G (matches the locus table) | T = **0.2389** | yes (all four filters) |
| HBS1L-MYB | rs9399137 | by rsID | 6:135419018, T/C (matches) | C = **0.0392** | yes (all four filters) |
| HBB cluster | rs7482144 | **not found** | | | **not testable in this callset** |

Frequencies are plink2 `--freq` on the 651 founders among the 652 samples
(1,302 alleles); one sample is recorded with a parent in the data.

**rs7482144 is absent from the 1000 Genomes phase 3 callset itself, not
removed by this pipeline.** The raw `all_phase3.pvar.zst` (all 2504 samples,
84.8M variants) has no record at 11:5276169 and none under rs7482144 or its
synonym rs41322148; the nearest variant is rs533892485 at 11:5276147, which
is not used (no proxy, D-028). This is consistent with Ensembl listing no 1000
Genomes frequency for rs7482144.

## PCA (`pca_eigenval.txt`, `pc1_pc2_by_population.png`)

Eigenvalues of the 10 PCs: 4.19, 3.21, 2.67, 1.66, 1.52, 1.28, 1.27, 1.24,
1.23, 1.22. In PC1 vs PC2 (652 samples: ACB 95, ASW 55, ESN 99, GWD 113,
LWK 97, MSL 85, YRI 108), ASW and ACB spread along one line towards negative
PC1 and PC2; LWK separates from the West African populations on PC1; YRI and
ESN overlap; GWD and MSL form adjacent clusters. One ASW sample, NA20314,
lies far out at (-0.25, -0.28).

## Heterozygosity (`het_summary.tsv`; reported only, nothing removed)

Mean F 0.0017 (SD 0.0146). Beyond 3 SD: 7 GWD samples (F 0.048-0.101) and
NA20314 (ASW, F 0.274), the same sample that is extreme in the PC plot. F is
computed against allele frequencies of the pooled AFR sample, so a sample
whose ancestry differs from that pool reads as excess homozygosity.
