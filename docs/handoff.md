# Handoff, 2026-11-12

Where the work stands at the end of 12 November 2026, and what waits for
you. Written by Claude at the end of your second work list, which said to
write this and stop.

## Confirmations

* **No phenotype was drawn or simulated on the real genotypes**, for any
  reason. The only real-data computations were step c (download, subset, QC,
  PCA) and step d (`03d` sensitivity PCA and `03c` predictions, which form no
  phenotype).
* **Nothing was merged.** **No tag exists** on GitHub (`git ls-remote --tags`
  is empty). **Step e was not run.**
* No force-push, no history rewrite, no visibility change. No criterion,
  seed or recorded result was changed. The repository is still private.
* The only change made on GitHub outside branches and PRs: **issue #2's
  title**, as you asked.

## What is complete from the overnight work list

| # | Item | Status |
|---|---|---|
| 1 | Sign-off PR closing #3 | **done**: #30 |
| 2 | Answer the PCA MAF and non-founder questions | **done**: #31 |
| 3 | Step d up to an open PR | **done**: #36 (committed overnight; opened today, after your second message) |
| 4 | Step e as a single command, tested in CI on synthetic data, with the refusal cases | **partly; see "Step e" below.** The refusals exist and are tested on synthetic data: no tag, unsigned table, empty seed, undeclared shift. **Missing:** the 2000 redraws on real genotypes, the NA20314 sensitivity association, and a single command. |
| 5 | Reporting on synthetic data (Manhattan, QQ both ways, PC1 vs PC2, results table, sensitivity comparison) | **skipped** |
| 6 | Timing benchmark on simulated genotypes of real size | **skipped** |
| 7 | README and docs, real results as placeholders | **partly**: #34 (claim, loci, status rows, how it was built; real results stay "not yet run") |
| 8 | Reproducibility: pin dependencies, one rebuild command, run it in CI | **partly**: `scripts/rebuild.sh`, run weekly and on demand (#32). **Dependencies are not pinned:** `pyproject.toml` has lower bounds only and there is no lock file. |
| 9 | Claims audit | **skipped** as a document. The overstatements the history scan found are fixed in #34. |
| 10 | Public-readiness checklist PR | **replaced** by the history scan (reported in chat), log redaction (#33) and the publication plan (#35) |
| 11 | Reviewer summary | **skipped**. The PR descriptions and this file cover part of it. |

**Why items were skipped:** your second message (12 November) replaced the
rest of the overnight list with seven items and said to stop after this
file. Items 5, 6, 9 and 11, the unfinished half of 4, and dependency pinning
(8) are waiting for you to say whether to do them.

## What was done from today's list

| # | Item | Result |
|---|---|---|
| 1 | Reduce CI usage | #32, and a usage report (below) |
| 2 | Log hygiene before any further real-data run | #33 |
| 3 | README | #34 |
| 4 | Issue #2 title | now "Phase 2: real 1000 Genomes AFR genotypes (steps a-c done; step d proposed, awaiting review)". Its checkboxes are unchanged (a-c are done but not ticked). |
| 5 | `docs/publication-plan.md` | #35 (plan only) |
| 6 | Open the step d PR | #36 |
| 7 | This file | #37 |

## The PRs, in the exact order to merge them

They form one stack: each one is built on the one above it. Merge strictly
in this order. Before merging a **draft**, click **Ready for review**, which
starts the full CI run (see `docs/ci.md`); GitHub will not merge a draft.
After each merge I will check that `main`'s tree is what was intended.

| Order | PR | What you approve by merging |
|---|---|---|
| 1 | #30 Locus sign-off (closes #3) | The sign-off record with who checked each row (20 verified, 2 superseded, 17 not covered), and **D-034**: 03c and 04 require the analysis check of the locus table |
| 2 | #31 PCA MAF and non-founder answers | The written answers, and the two recorded deviations: **D-036** (PCA MAF ≥ 0.05, Claude's unapproved addition) and **D-037** (NA20320's treatment). Approve or reject each. |
| 3 | #36 Step d pre-registration (not a draft) | The real-data criteria **real-v1**, the fixed population shift, the step e seed slot, the NA20314 sensitivity rule and the declared predictions. **Decide its two flags first.** Then I tag the merge commit `prereg-v1`, on your go-ahead. |
| 4 | #32 CI minutes (draft) | **D-042**: drafts run lint plus tests on Python 3.12; the full checks run before every merge; caches; a weekly clean-clone rebuild |
| 5 | #33 Log redaction (draft) | **D-043**: run records redacted when written (no hostname, repo-relative paths). **Must be merged before step e runs**, so step e's records come out redacted. |
| 6 | #34 README (draft) | The new first paragraph, two of three loci, the locus row and **the "How this was built" wording** (D-044, with the two exceptions D-022 and D-036 named) |
| 7 | #35 Publication plan (draft) | The plan (nothing executed) and **D-045**: publish later as a new repository from a scrubbed history; MIT licence |
| 8 | #37 This handoff (draft) | Documentation only |

## Every decision waiting for you

**Before merging #36 (step d):**
1. **Flag 1.** The predicted unadjusted lambda_GC is 1.1999 [1.1449,
   1.2703], so "unadjusted lambda_GC > 1.2" passes or fails about at
   random. Options:
   - keep it;
   - report it, not gated;
   - a threshold inside the predicted range;
   - a larger shift, which changes what is simulated.
2. **Flag 2.** BCL11A's conditional power is 0.9896, under the 0.99 needed
   to be REQUIRED, so no single-scan result is gated. Options: keep it, or
   change the threshold.
3. **D-034, D-036, D-037** (via #30 and #31). If you reject D-036, PCA is
   rerun without the MAF filter and the step d predictions are recomputed
   before #36 can merge.

**After #36 merges:**

4. **Go-ahead to tag** the merge commit `prereg-v1` and push the tag.
5. **The step e seed.** Pick any integer. It is committed to
   `config/step_e_seed.yaml` and pushed after the tag.
6. **The step e code PR**: the redraws, the sensitivity association and a
   single command. Should it be built, and should it land before or after
   the tag? Recommendation: before running step e. It does not change
   `analysis.yaml`, so the tag stays valid either way.

**CI:**

7. **Billing.** Jobs are still refused. The last runs, at 21:46 UTC on 12
   November, all said "The job was not started because recent account
   payments have failed or your spending limit needs to be increased". The
   $5 budget does not seem active yet; check Billing & plans. After that:
   - #30, #31 and #36 carry the old workflow (a full run, about 7 minutes
     each). Their failed checks can be re-run with one click, or you can
     rely on the local results below.
8. **Post-merge runs on main.** Keep them full (the current choice), or
   switch to lint only, saving about 7 minutes per code merge.
   Post-merge runs were 34 % of the minutes so far.

**Publication** (`docs/publication-plan.md` section 6):

9. Keep or rewrite the author email.
10. Export issues and PRs to `docs/archive/` (recommended).
11. Scrub the two pinned records (recommended) or keep them byte-identical.
12. The public repository's name.
13. Keep `docs/BRIEF.md` and `CLAUDE.md` (recommended).
14. Publish before or after step e.
15. **The MIT licence** is recorded and added in the public repository's
    first commit. Say if you also want a LICENSE file in this private
    repository.

**Remaining work:**

16. Overnight items 5, 6, 9 and 11, and dependency pinning: do, drop or
    defer.
17. Issue #2's checkboxes: tick a-c?

## Step e: the single command, and why it must not be run yet

**Do not run step e yet.** The code does not do everything the criteria
need:

* The gated power calibration compares empirical power over **2000
  phenotype redraws on the real genotypes** with the declared predictions.
  The redraws are not implemented, so the gate would stop at that
  criterion: `05_report.py` writes no `empirical_power`.
* The **NA20314 sensitivity association** (D-040) is not implemented.

Running it now would use up the single pre-registered run on an incomplete
analysis. Today it cannot start by accident: `04_assoc.sh` refuses without
the `prereg-v1` tag and with the empty seed.

**The command once the step e code PR has merged.** Run it from the
repository root, after the tag and the seed are pushed and #33 is merged:

```bash
scripts/step_e.sh
```

This script is proposed in that PR; it does not exist yet. It will run the
same steps as today's three commands, plus the redraws and the sensitivity
association:

```bash
RUN_DIR=runs/$(date -u +%Y-%m-%d)-stepe PYTHON=.venv/bin/python sh -c 'scripts/04_assoc.sh && "$PYTHON" scripts/05_report.py && "$PYTHON" -m hbfgwas.gate --experiment real --results "$RUN_DIR/results.json" --out "$RUN_DIR/gate.json"'
```

If the gate fails a criterion, nothing is changed and it is reported to you.

## CI minutes so far (details in #32)

* **This repository:** about 240 billed minutes since 10 November (47
  runs, all Linux). By job: tests 153 minutes, lint 47, gate 40. Post-merge
  runs on main took 82 minutes. A full run costs about 7 minutes.
* **The account quota** is shared with `proteomics-pipeline` (about 343
  minutes this month). By this count the total is about 580; the Billing
  page has the real figure.
* **After #32**, estimated: a draft run costs about 4 minutes, a
  documentation-only run 1, a mergeable or post-merge run about 7
  (unchanged). These are estimates, because no run of the new workflow
  has started yet.

## Test results (local runs, not CI; CI refused every job)

These were run on macOS arm64, Python 3.12, plink2 Beta 1.4 (7 Nov 2026),
with `REQUIRE_PLINK2=1`. shellcheck is not installed here, so it was not
run locally.

| Branch / PR | Result |
|---|---|
| #36 at `fdeefee` | 209 passed, 1 deselected (`realdata`) |
| #32 at `77fa4c2` | `scripts/rebuild.sh` from a clean clone: 241 passed; sim-v3 GATE PASSED; summary reproduced |
| #33 | 267 passed (full suite, before the commit); after the test-value fix in `14e3f41`, `tests/test_redact.py` 24 passed |
| #34, #35, #37 | Markdown only; wording check passed |
| #30, #31 | tested locally when opened (see each PR); their CI runs were refused for billing |
