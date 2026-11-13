# Running the pipeline on real 1000 Genomes genotypes

This is the Phase 2 walkthrough (see `docs/BRIEF.md`). **The phenotype is
simulated even here:** 1000 Genomes has no HbF phenotype. The real genotypes
replace the simulated ones; the HbF-like phenotype is still planted by
`hbfgwas` at the three sentinel variants. Nothing in this walkthrough has
been run yet (status: README).

Genotype files go to `data/` (git-ignored, never committed). Everything that
documents the run goes to `runs/<date>/` and **is** committed: the download
manifest, `plink2_version.txt`, every plink2 `.log`, the per-filter QC table,
PCA eigenvalues, results JSON and figures.

## 0. Before anything

1. **Sign off the locus table** (issue #3). Check every value in
   `src/hbfgwas/loci.py` against its cited source and set `verified=True`
   field by field. `python -m hbfgwas.loci check` must exit 0, and
   `pytest -m realdata` must pass. `04_assoc.sh` refuses to run otherwise.
   Effect alleles must be on the GRCh37 forward strand (1000 Genomes alleles).
2. **plink2 on PATH.** The official macOS Apple Silicon build is listed on
   <https://www.cog-genomics.org/plink/2.0/> as "macOS M1". On 2026-11-10 the
   current build was Beta 1.4 (7 Nov 2026):
   `https://s3.amazonaws.com/plink2-assets/beta1/plink2_mac_arm64_20261107.zip`.
   Builds are replaced often; take whatever the page lists on the day.
3. **Disk:** >= 20 GB free in `data/` (about 3 GB compressed download, about
   7 GB after decompression, plus intermediate filesets).
4. Python environment: `pip install -e ".[dev]"`; the scripts call `python3`
   unless `PYTHON` is set (e.g. `PYTHON=.venv/bin/python`).

Common environment for every step (run from the repository root):

```bash
export PYTHON=.venv/bin/python
export RUN_DIR=runs/$(date -u +%Y-%m-%d)
export DATA_DIR=data/1kg
```

## 1. Re-check the download URLs

The Dropbox links rotate. Compare the four defaults at the top of
`scripts/00_download_1kg.sh` with
<https://www.cog-genomics.org/plink/2.0/resources> (sections "1000 Genomes
phase 3": `all_phase3.pgen.zst`, `all_phase3_noannot.pvar.zst`,
`phase3_corrected.psam`, `deg2_phase3.king.cutoff.out.id`). Override any that
changed with `PGEN_URL`, `PVAR_URL`, `PSAM_URL`, `REL_URL`.

## 2. Download, subset, QC, PCA

```bash
scripts/00_download_1kg.sh   # manifest: $RUN_DIR/download_manifest.tsv
scripts/01_subset_afr.sh     # AFR minus related list; fails unless 640-665 samples remain
scripts/02_qc.sh             # $RUN_DIR/qc_filters.tsv, $RUN_DIR/sentinels.tsv
scripts/03_pca.sh            # $RUN_DIR/pca_eigenval.txt, 03_pca_counts.tsv, het_summary.tsv
.venv/bin/python scripts/03b_pc_plot.py   # $RUN_DIR/pc1_pc2_by_population.png
```

What to expect and check:

* `01`: sample and variant counts per step in `01_subset_counts.tsv`.
* `02`: 1000 Genomes phase 3 is **fully called**, so `--geno` and `--mind`
  remove nothing; the script prints a note saying so (a warning if they
  ever do remove something). `--maf` and `--hwe` counts are in
  `qc_filters.tsv`.
* `02` also writes `sentinels.tsv`: whether each sentinel was found (by rsID
  or, failing that, by GRCh37 position), its effect-allele frequency in this
  AFR sample, and which filter removed it, if any. rs7482144 has no 1000
  Genomes frequency in Ensembl and may be stored under rs41322148 or be
  absent.
* `03`: long-range LD regions come from `resources/long_range_ld_grch37.tsv`
  (every entry unverified; see `resources/README.md`). The KING step runs
  even though the 1000 Genomes related-sample list was already removed; if
  it removes nobody, the script says so, as 02 does for `--geno`/`--mind`.
* `03b_pc_plot.py`: PC1 vs PC2 coloured by the seven AFR populations (the
  `Population` column, not SuperPop), with the count per population.

Commit `$RUN_DIR` (logs, tables, eigenvalues). **Stop and report** (Phase 2
step c).

## 3. Pre-register the real-data criteria (Phase 2 step d)

The real-data criteria follow the **sim-v3 rule** (docs/decisions.md D-024):
gated power calibration (|empirical - predicted| <= 0.05 per locus over 2000
phenotype redraws on the real genotypes), gated single-scan recovery only for
loci with conditional predicted power >= 0.99, and every other locus's
single-scan result reported next to its predicted power. The predictions are
computed from the real genotypes and PCs before any phenotype is drawn.

Before any phenotype exists, run

```bash
scripts/03c_prereg_predictions.sh   # $RUN_DIR/prereg_predictions.json
```

It writes, from the real genotypes and PCs:
- the conditional power of each testable sentinel, and which loci the sim-v3
  rule makes REQUIRED (predicted power >= 0.99);
- the **predicted unadjusted lambda_GC** (and PC-adjusted), each as a median
  with a 5th-95th percentile over draws of the population shift (review
  item H);
- the sentinels that are not testable in this callset.

Write the real-data acceptance criteria into `config/analysis.yaml` under
`acceptance.real`. Include a `predicted_power` entry per locus, which
`05_report.py` then uses instead of computing power at report time, and the
predicted lambda_GC. Open a PR; the user reviews and merges it.
Then tag the merge commit and push the tag **before** step 4:

```bash
git tag prereg-v1 <merge commit>
git push origin prereg-v1
```

## 4. Step e: one command, run once

After the step d PR is merged, `prereg-v1` is tagged and pushed, and the step e
seed is committed to `config/step_e_seed.yaml` and pushed:

```bash
RUN_DIR=runs/$(date -u +%Y-%m-%d)-stepe PYTHON=.venv/bin/python scripts/step_e.sh
```

`step_e.sh` refuses to start, before any phenotype is simulated, unless:
- `prereg-v1` exists here and on GitHub, pointing to the same commit;
- `config/analysis.yaml` is unchanged since the tag;
- `acceptance.real` is filled;
- the locus table passes the analysis check;
- the step e seed is filled;
- the population shift is declared;
- the PCA (03) and the sensitivity PCA (03d) exist;
- `$RUN_DIR/results.json` does not exist yet (step e runs once).

It then runs, in order:
- `04_assoc.sh`: simulated phenotype; plink2 `--glm` without and with the 10 PCs; concordance subset.
- `04b_power_redraws.sh`: `acceptance.real.power_calibration.n_redraws` phenotype redraws (2000). Each redraw keeps the same planted effects and the same fixed shift, with new noise from the seed's child `power_reps`. Each is tested at the sentinels only, with and without PCs. Writes `power_redraws.json`.
- `04c_sensitivity_scan.sh`: the same phenotype without `sensitivity.exclude_samples` (NA20314), with the PCs recomputed without it by 03d. Never gated.
- `05_report.py`: writes `results.json`:
  - lambda_GC both ways;
  - type I error on variants more than 1 Mb from a planted sentinel;
  - per locus: the declared prediction, the empirical power over the redraws, and the single-scan p-value;
  - the sensitivity section;
  - the Python-vs-plink2 concordance.

  It also draws `qq_real.png` and `pc1_pc2_by_population.png`.
- `python -m hbfgwas.gate --experiment real`: writes `gate.json`.

Exit status:
- 0: every gated criterion passed;
- 3: a gated criterion failed;
- 4: the gate passed but the concordance check failed;
- 1: refused, or an error.

Commit `$RUN_DIR` under `runs/<date>-stepe/`. If the gate fails a criterion,
change nothing and report it.
