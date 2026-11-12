# Two questions about the Phase 2 step c run, answered (2026-11-12)

## 1. Was PCA's MAF >= 0.05 in the config before the run?

**Yes, but it was Claude's addition, never approved by the user.**

- `pca.maf_for_pca: 0.05` was added in commit `f63e8bd` (2026-11-10 18:40
  UTC), the Phase 1 pre-registration (PR #9), and recorded as decision D-012
  ("PCA settings not fixed by the brief", by Claude).
- It reached `main` with #16 (2026-11-11 16:18 UTC). The step c run used
  `main` plus #27 and #28; it started downloading at 17:36 UTC and ran PCA
  around 18:10 UTC. `runs/2026-11-11/logs/pca.0_candidates.log` shows
  `--maf 0.05`, which removed 6,053,370 of the variants left after the
  long-range-LD exclusion (`03_pca_counts.tsv`: 8,078,843 remained).
- The brief's item 8 lists the PCA steps (long-range LD exclusion, LD
  pruning, KING 0.0884, exact `--pca 10`, `--het`) and does not mention a MAF
  threshold. So this is a deviation from the brief by addition. It was
  declared before the run, but it predates the rule that Claude asks before
  such changes (D-027) and the user has not explicitly approved it. Recorded
  as D-036.

## 2. The non-founder sample: where is its parent, and how does each step treat it?

**NA20320** (ASW, female) is the only non-founder among the 652 samples. Its
`.psam` row names a mother, **NA20321** (ASW), and no father.

- NA20321 is on the 1000 Genomes second-degree related-sample list
  (`deg2_phase3.king.cutoff.out.id`), so `01_subset_afr.sh` removed it. It is
  in the downloaded `.psam` but not among the 652 samples.
- NA20320 is not on the list and was kept. Since the parent field is still
  set, plink2 counts NA20320 as a non-founder: every log reports
  "652 samples ... 651 founders". The plink2 filter documentation says why:
  "if parental IDs are provided for a sample, they are not treated as a
  founder even if neither parent is in the dataset".

How each step treated NA20320 (plink2 defaults; nothing in the config or
scripts changes them):

| step | NA20320 | basis |
|---|---|---|
| 01 AFR subset, related list | kept (its mother removed) | `01_subset_afr.sh` |
| 02 `--geno` | its genotypes counted | sample-wide missingness filter (removed nothing) |
| 02 `--mind` | kept | it has no missing calls (removed nothing) |
| 02 `--maf 0.01` | **excluded from the frequency calculation**, genotypes kept | docs: "--maf and --max-maf default to ignoring nonfounders" |
| 02 `--hwe` | **excluded from the test** | docs: "Only founders are considered by this test"; log: "(founders only)" |
| 02 sentinel frequencies | **excluded** (651 founders, 1,302 alleles) | log: "Allele frequencies (founders only)" |
| 03 `--maf 0.05` for PCA variants | excluded from the frequency calculation | as for `--maf 0.01` |
| 03 `--indep-pairwise` | **excluded from the LD (r^2) calculation** | docs (LD page): "All of the following calculations only consider founders" |
| 03 `--king-cutoff` | included, not removed (its mother is absent) | `king_removed.id` is empty |
| 03 `--pca 10` | has PC values (it is in `pca.eigenvec` and the plot) | the documentation found does not say whether `--pca` leaves nonfounders out of the decomposition; not verified |
| 03 `--het` | F reported for it (not an outlier) | `het_summary.tsv` |
| step d predictions (Python) | included: all 652 samples | `scripts/03c_prereg_predictions.sh` uses `unrelated.id` |
| step e association (planned) | included: all 652 samples | `04_assoc.sh` uses `unrelated.id` |

Effect: one of 652 samples (0.15%) is left out of the frequency-based
filters, the HWE test, the LD pruning and the reported sentinel
frequencies. It is in KING, the PCs and the association. The alternatives
would be plink2's `--make-founders` (clears the parent field when the parent
is absent) or `--nonfounders`. Either is a change to the QC; nothing has been
changed (D-037).
