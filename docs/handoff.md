# Handoff, 2026-11-12 (after the restack)

Where the work stands, and what waits for you. This version replaces the
first handoff (#37). Since then:
- the stack was restacked so the step d pre-registration is the **last** PR;
- the missing step e code was built;
- step d was amended as you decided.

## Confirmations

* **No phenotype was drawn or simulated on the real genotypes**, for any
  reason. The real-data computations so far:
  - step c (download, subset, QC, PCA);
  - step d (sensitivity PCA, and predictions that form no phenotype);
  - the step d amendment's rerun of those predictions, which also forms no
    phenotype (`runs/2026-11-12-stepd-amend/`).
* **Nothing was merged.** **No tag exists** on GitHub. **Step e was not
  run**, and its seed is still empty.
* No force-push, no history rewrite, no visibility change. The repository is
  still private.
* **One criterion was changed, by your decision, before any real-data
  phenotype existed:** the unadjusted lambda_GC criterion (D-048). No other
  criterion, no existing seed and no recorded result was changed. The
  amendment uses a new child (2) of the existing `real_prediction` seed
  (D-049).
* GitHub changes outside branches and PRs: issue #2's title, and the
  replaced PRs (below) are now drafts with a comment, so they cannot be
  merged by mistake. Their branches are kept.

## The restack: which PRs replace which

| Old PR (draft now, replaced) | New PR | Notes |
|---|---|---|
| none (step d's code commits, `a2bbb96` and `67a8df2`) | **#38 step e code** | The branch starts at `67a8df2`, so those hashes are kept, and adds the missing step e pieces |
| #32 CI minutes | **#39** | same content |
| #33 log redaction | **#40** | one commit, never holding the hostname; covers the new step e scripts |
| #34 README | **#41** | step d rows wait for the step d PR |
| #35 publication plan | **#42** | one sentence changed |
| #37 handoff | **#43** (this file) | |
| #36 step d | **#44**, last | merges the original `f0c9a5f`, `ed3d3ce` and `fdeefee` (hashes kept), then the amendment |

## The PRs, in the exact order to merge them

They form one stack: each one is built on the one above it. Merge strictly
in this order. A draft needs **Ready for review**, which starts the full CI
run (`docs/ci.md`), before GitHub will merge it. After each merge I will
check that `main`'s tree is what was intended.

| Order | PR | What you approve by merging |
|---|---|---|
| 1 | #30 Locus sign-off (closes #3) | The sign-off record (20 verified, each with who checked it; 2 superseded; 17 not covered), and **D-034** (03c, 04 and step e require the analysis check) |
| 2 | #31 PCA MAF and non-founder answers | The written answers; approve or reject **D-036** (PCA MAF ≥ 0.05) and **D-037** (NA20320) |
| 3 | #38 Step e code | `step_e.sh` and its refusal checks, the 2000 redraws at the sentinels, the NA20314 sensitivity scan (**D-046**), and the fixed-shift code, seed slot and sensitivity PCA from step d |
| 4 | #39 CI minutes | **D-042**: draft PRs run lint and Python 3.12; full checks before every merge; caches; weekly rebuild |
| 5 | #40 Log redaction | **D-043**: records redacted when written. Must be in before step e runs. |
| 6 | #41 README | The first paragraph, two of three loci, and the **"How this was built"** wording (D-044) |
| 7 | #42 Publication plan | The plan (nothing executed), and **D-045**: publish later as a new repository from a scrubbed history; MIT |
| 8 | #43 This handoff | Documentation only, plus **D-047** (the restack) |
| 9 | #44 Step d pre-registration (**last**) | The real-data criteria **real-v1** as amended (the unadjusted lambda_GC interval [1.1248, 1.2875], **D-048**), each criterion's pass probability (**D-049**), flag 2 unchanged (**D-050**), the fixed shift, the seed slot, the sensitivity rule, and the declared predictions. Then I tag its merge commit `prereg-v1`, on your go-ahead. |

## Every decision waiting for you

The 17 from the first handoff, with where each stands now, then the new
ones.

**Before merging the step d PR (#44):**
1. **Flag 1: decided by you (12 Nov).** "unadjusted lambda_GC > 1.2" is
   replaced by "inside the pre-stated 1st-99th percentile prediction
   interval [1.1248, 1.2875]" (D-048, in #44). Left for you: approve the
   interval, and its pass probability of 0.98 (below 0.99), by merging #44.
2. **Flag 2: decided by you, no change.** BCL11A stays reported, not gated
   (D-050).
3. **D-034, D-036, D-037** (via #30 and #31): still open. If you reject
   D-036, PCA is rerun without the MAF filter and the step d predictions are
   recomputed before #44 can merge.

**After #44 merges:**

4. **Go-ahead to tag** the merge commit `prereg-v1` and push the tag: open.
5. **The step e seed**: open. Pick any integer. It is committed to
   `config/step_e_seed.yaml` and pushed after the tag.
6. **The step e code PR: done**, as #38, which now merges before step d.
   Left for you: review and merge #38.

**CI:**

7. **Billing:** you are fixing it. Every run on 12 November was refused
   ("The job was not started because recent account payments have failed or
   your spending limit needs to be increased"). When CI runs, #30 and #31
   need a re-run.
8. **Post-merge runs on main:** keep them full (current), or lint only
   (saves about 7 minutes per code merge): open.

**Publication** (`docs/publication-plan.md` section 6):

9. Keep or rewrite the author email.
10. Export issues and PRs to `docs/archive/` (recommended).
11. Scrub the two pinned records (recommended) or keep them byte-identical.
12. The public repository's name.
13. Keep `docs/BRIEF.md` and `CLAUDE.md` (recommended).
14. Publish before or after step e.
15. **The MIT licence** is recorded for the public repository's first
    commit. Do you also want a LICENSE file in this private one?

**Remaining work:**

16. Overnight items 5 (reporting figures), 6 (timing benchmark), 9 (claims
    audit) and 11 (reviewer summary), and dependency pinning: do, drop or
    defer.
17. Issue #2's checkboxes: tick a-c?

**New since the restack:**

18. **Close the replaced PRs** #32-#37? They are drafts now, kept until you
    say.
19. **The assumptions behind the interval** (D-049): the remaining noise
    acts independently on each variant, and the observed lambda_GC behaves
    like a median over 215,522 independent variants (step c's LD-pruned
    count).
    - The check on simulated genotypes cannot test real LD. In that check
      the 1st-99th interval covered the observed value 99.3% of the time.
    - A smaller M_eff would widen the interval. Its pass probability stays
      0.98 by construction.
    - Accept these assumptions, or choose another M_eff before merging #44.
20. **"Declared predictions reproduced"** as a gated check. Step e
    recomputes each locus's conditional power from the real genotypes and
    reports it next to the declared value. Gating it would add a criterion
    whose pass probability for a correct pipeline is 1 (it is
    deterministic). Currently reported, not gated (D-046).

## Step e: the single command

After #44 is merged, `prereg-v1` is tagged and pushed, and the seed is
pushed, run from the repository root (`docs/running-real-data.md` section 4):

```bash
RUN_DIR=runs/$(date -u +%Y-%m-%d)-stepe PYTHON=.venv/bin/python scripts/step_e.sh
```

It refuses to start, before any phenotype exists, unless:
- the tag is here and on GitHub, pointing to the same commit;
- the config is unchanged since the tag;
- `acceptance.real` is filled;
- the locus table passes its analysis check;
- the seed is filled;
- the population shift is declared;
- the PCA (03) and the sensitivity PCA (03d) exist;
- `RUN_DIR/results.json` does not exist yet (step e runs once).

Then it runs 04, 04b (2000 redraws), 04c (the scan without NA20314), 05 and
the gate. Exit status:
- 0: every gated criterion passed;
- 3: a gated criterion failed (nothing is changed; it is reported);
- 4: the concordance check failed;
- 1: refused, or an error.

A correct pipeline passes all gated criteria with probability 0.98 (#44).

## Test results (local runs, not CI; CI refused every job)

These were run on macOS arm64, Python 3.12, plink2 Beta 1.4 (7 Nov 2026),
with `REQUIRE_PLINK2=1`. shellcheck is not installed here, so it was not
run.

| PR | Result |
|---|---|
| #38 | 221 passed, 1 deselected (`realdata`), including step e end to end on the synthetic fileset and its six refusals |
| #39 | `tests/test_ci_scope.py` and `tests/test_config.py`: 51 passed |
| #40 | 279 passed, 1 deselected |
| #41, #42, #43 | Markdown only; wording check passed |
| #44 | 287 passed, 1 deselected |
| #30, #31 | tested locally when opened (see each PR) |
