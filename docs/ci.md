# CI: what runs when, and what it costs

The repository is private, so GitHub Actions minutes are billed. Linux
runners bill each job rounded **up** to a whole minute. The rules below cut
minutes without weakening what is checked before a merge to main.
`scripts/ci_scope.py` applies them and `tests/test_ci_scope.py` tests them.

## What runs when

| Event | lint + wording check | tests | gate (simulated experiment) |
|---|---|---|---|
| Pull request into main, not a draft (opened, new push, reopened, **Ready for review**) | yes | Python 3.11 and 3.12 | yes |
| Draft pull request (how stacked PRs are opened), or a PR into another branch | yes | Python 3.12 only | no |
| Push to main (after a merge) | yes | 3.11 and 3.12 | yes |
| Manual run (Actions, CI, Run workflow) | yes | 3.11 and 3.12 | yes |
| Any of the above when every changed file is documentation that no check reads | yes | skipped | skipped |
| Weekly (Mondays 06:17 UTC) or manual: **Rebuild from a clean clone** | yes | 3.12, everything | yes |

A newer push to the same pull request or branch cancels the run it replaces.
Every job has a timeout (10-30 minutes), so a hung job cannot use hours.

"Documentation that no check reads" means `*.md` files and `LICENSE`. It
excludes `docs/locus-signoff.md` (a test compares it with its generator) and
everything under `docs/figures/` (the gate compares its summary). A test fails
if a test starts reading another `.md` file without it being added to that
list. Which files count as changed:

* for a pull request that can be merged, every file it changes against
  main, however many pushes it took;
* for a draft, only the latest push.

The draft rule is safe because a draft cannot be merged without the full run.

## Merging a stacked pull request

1. Merge the pull request below it first.
2. On the next one, click **Ready for review**. This starts the full run
   against main as it now is.
3. Merge when that run is green. GitHub will not merge a draft.

Pull requests whose branch does not contain this workflow (#30 and #31)
are judged by the old one: the full checks on every run, draft or not.

## Rebuild from a clean clone

`scripts/rebuild.sh` (locally, or Actions, "Rebuild from a clean clone", Run
workflow) does the following:

1. Clones the committed HEAD into a new temporary directory.
2. Builds a new virtualenv without the pip cache.
3. Runs lint, the wording check and every test.
4. Runs the pre-registered simulated experiment with its gate, and checks
   that the result reproduces `docs/figures/sim_stratification_summary.json`.

In Actions it downloads plink2 fresh instead of using the cache. It never
touches real genotypes.

## Caches

* **pip:** keyed on `pyproject.toml`.
* **The pinned plink2 zip:** keyed on its sha256. The zip is checked against
  the sha256 on every use, cached or not.

Both live in the Actions cache, which is free up to 10 GB per repository.
Artifacts from the gate are kept 14 days instead of 90.
