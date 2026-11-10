First, save this entire message as docs/BRIEF.md and create a CLAUDE.md that tells you to re-read it at the start of every session and phase.

You are building a new project from scratch in this empty repository: hbf-gwas.

## What this is

A small, tested GWAS QC and association pipeline, built as a portfolio piece.
Its claim is narrow and must stay narrow:

  "On real 1000 Genomes AFR genotypes, with an HbF-like phenotype simulated at
  the three established HbF modifier loci (BCL11A 2p16, HBS1L-MYB 6q23, HBB
  cluster 11p15), the pipeline recovers planted effects at the power predicted
  in advance, with calibrated type I error after ancestry adjustment."

1000 Genomes has NO HbF phenotype. Never write "rediscovers", "validates
against biology", or "replicates" anywhere in code, docs, or commit messages.
The phenotype is always simulated. Say so in the README's first paragraph.

## Constraints of this session

- You are in a cloud sandbox. Do NOT download 1000 Genomes data and do NOT run
  plink2 on real data. Write those scripts, lint them (shellcheck if
  available), and leave them unrun. Say clearly in your final report what was
  executed and what was only written.
- Never commit data. .gitignore data/ and results/ except small committed
  figures under docs/figures/.
- Do not report a test as passing unless you ran it and saw it pass.
- Do not tune simulations to hit a number I have not given you. Fix seeds,
  report what comes out.

## Phase 1: build in this session

Work test-first. One PR per numbered item, small commits.

1. Scaffold
   - Python package `hbfgwas` (src layout, pyproject.toml, Python 3.11+),
     pytest, ruff, GitHub Actions CI running lint + tests on push and PR.
   - Dependencies kept minimal: numpy, scipy, pandas, matplotlib, pyyaml.

2. hbfgwas/loci.py
   - A table of the three loci with sentinel variants: BCL11A rs1427407,
     HBS1L-MYB rs9399137, HBB cluster rs7482144 (XmnI).
   - Fields: locus, gene, rsID, chromosome, GRCh37 position, effect allele,
     literature effect size or variance explained, citation.
   - Do NOT fill positions, alleles, effect sizes, or citations from memory.
     If you cannot verify a value from an authoritative source in this
     session, set it to null and mark it `verified: false`. A test must fail
     the "real-data ready" check while any field is unverified. I will fill
     and verify these by hand.

3. hbfgwas/power.py
   - Analytic power for a single-variant additive test via the noncentral
     chi-square: inputs n, variance explained (or beta + MAF), alpha.
   - Tests against hand-checkable cases (e.g. n=650, alpha=5e-8: about 50%
     power near R2 = 4.5%).

4. hbfgwas/simulate.py
   - Genotype simulation with population structure (Balding-Nichols, several
     subpopulations, configurable Fst).
   - Phenotype simulation with three separable components: planted variant
     effects, a subpopulation-correlated non-genetic shift, and noise.
     Must also accept a real genotype matrix in place of simulated genotypes.
   - Calibration tests that fail if the simulation stops producing structure.

5. hbfgwas/assoc.py
   - Additive linear regression per variant with covariates, vectorised.
   - PCA on a standardised genotype matrix for ancestry covariates.
   - lambda_GC, empirical type I error on null variants, empirical power on
     planted variants, QQ plot.
   - This module is the tested reference implementation, not the genome-wide
     engine. Include a function that compares its betas and p-values against
     a plink2 --glm output file for the same variants, with a tolerance.

6. config/analysis.yaml (pre-registration) and a gate script
   - Declares QC thresholds, number of PCs (10), significance threshold
     (5e-8), seeds, and acceptance criteria, all before any experiment runs.
   - Acceptance criteria for the simulated experiment:
       unadjusted lambda_GC > 1.2 (structure is present)
       PC-adjusted lambda_GC in [0.97, 1.05]
       PC-adjusted type I error at alpha=0.05 in [0.04, 0.06]
   - Locus recovery criterion: a locus is REQUIRED to be recovered only if
     its pre-computed analytic power is >= 0.8. Loci below that are reported
     as "underpowered, not required" with their power stated. Do not write a
     criterion that demands all three loci.
   - CI runs the gate; a statistical regression breaks the build.
   - If the first run fails a criterion, do not edit the criterion. Stop and
     report it to me.

7. Stratification experiment (simulated genotypes)
   - Script producing: lambda_GC and type I error unadjusted vs 10 PCs, a
     two-panel QQ plot, and a results JSON. Commit the figure to
     docs/figures/ and the JSON summary alongside it.

8. Real-data scripts (WRITE ONLY, do not run)
   - scripts/00_download_1kg.sh: downloads 1000 Genomes phase 3 (GRCh37)
     PLINK2 files. URLs come from environment variables with these defaults,
     which must be re-checked against
     https://www.cog-genomics.org/plink/2.0/resources before use because the
     Dropbox links rotate:
       PGEN  https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1
       PVAR  https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1
       PSAM  https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1
       REL   https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1
     Save as all_phase3.{pgen.zst,pvar.zst,psam}. Check free disk (>= 20 GB),
     fail if any file is under 1 MB (a dead link returns an HTML page), write
     a manifest with URL, date, size, sha256, then
     `plink2 --zst-decompress all_phase3.pgen.zst all_phase3.pgen`.
   - scripts/01_subset_afr.sh: keep SuperPop AFR, remove the REL list,
     autosomes, biallelic ACGT SNPs, remove duplicates. Assert 640-665
     samples remain.
   - scripts/02_qc.sh: --geno 0.02 --mind 0.02 --maf 0.01 --hwe 1e-6 midp.
     Log samples and variants in/out per filter to a TSV. Phase 3 is fully
     called, so the missingness filters will remove nothing: the script must
     print that explicitly and the README must say it.
   - scripts/03_pca.sh: exclude long-range LD regions (GRCh37; put the region
     file in resources/ with its source cited, unverified entries flagged),
     LD-prune, KING cutoff 0.0884, exact --pca 10 (not approx), --het report.
   - scripts/04_assoc.sh: phenotype simulated on the real genotypes by
     hbfgwas, then plink2 --glm with and without the 10 PCs.
   - scripts/05_report.py: lambda_GC both ways, PC1 vs PC2 coloured by the
     seven AFR populations (not by SuperPop), locus recovery vs pre-stated
     power, and the Python-vs-plink2 concordance check.

9. README
   - First paragraph: what it claims and what it does not (simulated
     phenotype, no biological discovery).
   - A status table with two columns kept honest: "run on simulated
     genotypes" and "run on real 1000 Genomes genotypes". The second column
     stays "not yet run" until I run it.
   - How to reproduce, and a docs/running-real-data.md walkthrough.

## Out of scope

No mixed models (REGENIE, SAIGE), no imputation, no fine-mapping, no
GRCh38, no other ancestries. If you think one is needed, say why and stop.

## When Phase 1 is done

Stop and give me a report: what was built, what was run with actual output
numbers, what was written but not run, every value marked unverified, and
anything in this brief you think is wrong.

## ADDENDUM: this overrides the brief above wherever they conflict

You are NOT in a cloud sandbox. You are running locally on an Apple Silicon
Mac, in a clone of github.com/ChiaLungChuang/hbf-gwas. Ignore the section
"Constraints of this session" except for these rules, which still hold:
never commit genotype data, never report a test as passing unless you ran it,
never tune a simulation to hit a number, never edit an acceptance criterion
after seeing a result.

### GitHub is the record

- Every piece of work goes through a branch and a PR with a description of
  what was done and what the output was. Push after every working step.
- Open one GitHub issue per phase and per open decision; close them from PRs.
- Genotypes stay out of git (data/ ignored). Everything that documents a run
  IS committed: the download manifest (URL, date, size, sha256), the plink2
  version string, every plink2 .log file, the per-filter QC table, PCA
  eigenvalues, results JSON, and figures. Put these under runs/<date>/.
- Keep a CHANGELOG.md and a docs/decisions.md recording each decision, who
  made it (me or you), and why.

### Phase 1: as in the brief, items 1-7 and 9

Item 8 changes: the real-data scripts are written AND tested here on a tiny
synthetic pgen fileset generated in the test suite, so the scripts are
exercised in CI before they touch real data. Use curl -L, not wget (macOS).

Locus table: you may look up positions, alleles, and effect sizes from
dbSNP, Ensembl GRCh37, and the primary papers, citing each source. Every
value still stays `verified: false` until I sign it off by hand.

When Phase 1 is done, STOP and report. Do not start Phase 2 until I reply.

### Phase 2: real data, only after I approve Phase 1

a. Setup. Check for plink2 on PATH. If missing, tell me the exact official
   macOS ARM build to download from cog-genomics.org and wait for me to
   approve or install it. Confirm >= 20 GB free disk.
b. Before any download, re-check the four URLs against
   https://www.cog-genomics.org/plink/2.0/resources and show me any that
   differ from the defaults. Then ask me before starting the download
   (about 3 GB compressed, about 7 GB after decompression).
c. Run 00 download, 01 AFR subset, 02 QC, 03 PCA. Commit the logs and QC
   table. Produce the PC1 vs PC2 plot coloured by the seven AFR populations.
   STOP and report: sample and variant counts at each step, the real allele
   frequency of each sentinel variant in this AFR sample, and whether each
   sentinel survived QC.
d. Power and pre-registration. Using the real n and the real sentinel allele
   frequencies, compute analytic power per locus. Write the resulting
   real-data acceptance criteria into config/analysis.yaml, open a PR, and
   WAIT for me to review and merge it. After I merge, tag that commit
   `prereg-v1` and push the tag. The tag must exist on GitHub before step e
   so the history proves the criteria came first.
e. Run 04 association and 05 report. Commit outputs under runs/<date>/.
   Run the gate. If a criterion fails, do not change anything: report it.
f. Update the README status table and figures from the real run, keeping
   the simulated-genotype and real-genotype results clearly separate.

Final report: what ran, the actual numbers, what failed or surprised you,
and every value still marked unverified.
