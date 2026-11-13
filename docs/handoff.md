# Handoff, 2026-11-12 (after your review of the step d amendment)

Where the work stands. This version replaces the one in #43. Since then:
- you decided every open item about the pre-registration (D-051 to D-056);
- one small PR (#45) adds the step e integrity check, pinned dependencies and
  a step e benchmark;
- the replaced PRs are closed.

## Confirmations

* **No phenotype was drawn or simulated on the real genotypes**, for any
  reason. The real-data computations so far:
  - step c;
  - step d's sensitivity PCA and predictions;
  - the step d amendment's rerun of those predictions.

  None of these forms a phenotype.
* **I merged nothing.** You merged #30 and #31 (22:35-22:36 UTC on 12 November). `main`'s tree is identical to #31's tip `aaa2c93`, as checked, and main's CI passed on both merges. **No tag exists** on GitHub. **Step e was not
  run**, and its seed is still empty.
* No force-push, no history rewrite, no visibility change. The repository is
  still private.
* **The one criterion change** (the unadjusted lambda_GC interval, D-048,
  accepted in D-051) was your decision, made before any real-data phenotype
  existed.

## The PRs, in the exact order to merge them

Each is built on the one above it. A draft needs **Ready for review**, which
starts the full CI run, before GitHub will merge it. CI runs again since
about 23:00 UTC on 12 November. After each merge I will check that `main`'s
tree is what was intended.

| Order | PR | What you approve by merging |
|---|---|---|
| 1 | #30 Locus sign-off (closes #3) | **merged by you**; D-034 approved (D-055) |
| 2 | #31 PCA MAF and non-founder answers | **merged by you**; D-036 (D-054) and D-037 (D-056) approved |
| 3 | #38 Step e code | `step_e.sh` and its refusals, the 2000 redraws, the NA20314 sensitivity scan (D-046) |
| 4 | #39 CI minutes | D-042. Full CI on main stays, as you decided (D-060). |
| 5 | #40 Log redaction | D-043; must be in before step e runs |
| 6 | #41 README | The first paragraph, two of three loci, "How this was built" (D-044) |
| 7 | #42 Publication plan | The plan (nothing executed) and D-045. Every publication decision is deferred until after step e (D-060). |
| 8 | #43 Handoff after the restack | Documentation, and D-047 (the restack) |
| 9 | #45 Integrity check, pinned dependencies, benchmark | D-057, D-058, D-059; your decisions D-051 to D-056 and D-060; this file |
| 10 | #44 Step d pre-registration (**last**) | real-v1 as amended and accepted, the shift, the seed slot, the sensitivity rule, the declared predictions. Then I tag its merge commit `prereg-v1`, on your go-ahead. |

## Every decision from the earlier lists, and where it stands

1. Flag 1: **accepted** ([1.1248, 1.2875], P(pass) 0.98; D-051).
2. Flag 2: **unchanged** (D-053).
3. D-034, D-036, D-037: **approved** (D-055, D-054, D-056).
4. Go-ahead to tag `prereg-v1` after #44 merges: **open**.
5. The step e seed: **open**. Pick any integer, after the tag.
6. The step e code: **done** (#38).
7. CI billing: **fixed**; #30 and #31 are merged and main's CI passed.
8. Post-merge runs on main: **keep full** (D-060).
9.-15. Publication decisions (email, export, pinned records, name,
   `BRIEF.md`/`CLAUDE.md`, timing, LICENSE): **deferred until after step e**
   (D-060).
16. Remaining work:
    - **done:** the timing benchmark and dependency pinning (#45);
    - **deferred until after step e:** the reporting figures and the claims
      audit;
    - **dropped:** the reviewer summary (D-060).
17. Issue #2 checkboxes a-c: **ticked**.
18. Replaced PRs #32-#37: **closed**, with pointers; branches kept.
19. D-049's assumptions, M_eff = 215,522: **accepted** (D-052).
20. "Declared predictions reproduced": done as an **integrity check that
    stops step e before any phenotype**, not as a criterion (D-057).

**Still open: 4 and 5.** Also each merge, which is your approval of what that
PR contains.

## Step e: the single command, and what it costs

After #44 is merged, `prereg-v1` is tagged and pushed, and the seed is
pushed:

```bash
RUN_DIR=runs/$(date -u +%Y-%m-%d)-stepe PYTHON=.venv/bin/python scripts/step_e.sh
```

It refuses before any phenotype unless all of these hold:
- the tag is here and on GitHub, pointing to the same commit;
- the config is unchanged since the tag;
- `acceptance.real` is filled;
- the locus table is signed off;
- the seed is filled;
- the shift is declared;
- the PCAs exist;
- no earlier run exists;
- **the declared predictions recompute from the genotypes** (D-057).

Then it runs 04, 04b, 04c, 05 and the gate. Exit status:
- 0: every gated criterion passed;
- 3: a gated criterion failed (nothing changed; reported);
- 4: concordance failed;
- 1: refused, or an error.

A correct pipeline passes all gated criteria with probability 0.98.

**Expected cost** (`runs/2026-11-12-benchmark/`; simulated genotypes of the
real size, a 16 GB, 8-core Mac):
- **about 8 to 10 minutes**: 7.4 measured with the sensitivity scans
  assumed, plus the null-variant sampling the benchmark skipped;
- **peak memory about 2.5 GB** (05 reading the 1.24 GB result tables);
- **about 5.3 GB written** to `data/1kg/` (four genome-wide result tables).

## Test results

These are local runs (macOS arm64, Python 3.12, plink2 Beta 1.4,
`REQUIRE_PLINK2=1`). CI runs again; check each PR's run when you click Ready
for review.

| PR | Local result |
|---|---|
| #38 | 221 passed, 1 deselected (`realdata`) |
| #39 | CI and config tests: 51 passed |
| #40 | 279 passed, 1 deselected |
| #41, #42, #43 | Markdown only; wording check passed |
| #45 | **291 passed**, 1 deselected, with `HBFGWAS_REQUIRE_LOCK=1` (installed versions equal the pins) |
| #44 | see its description (rerun after #45 is merged into it) |

CI on 12 November, after billing was fixed:
- #38: full checks, 8 billed minutes, passed;
- #39-#42: draft checks, 5 billed minutes each, passed.
