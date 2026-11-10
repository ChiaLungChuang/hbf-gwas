# Decisions

Each entry records what was decided, who decided it (**user** or **Claude**),
and why. The brief itself (`docs/BRIEF.md`) records the user's original
decisions; entries here cover interpretations of the brief and choices it
leaves open. Entries are never deleted; a reversed decision gets a new entry
that points back to the old one.

| ID | Date | Decision | By | Why |
|----|------|----------|----|-----|
| D-001 | 2026-11-10 | The first commit (`docs/BRIEF.md`, `CLAUDE.md`) went directly to `main`. All later work goes through branches and PRs. | Claude | The GitHub repository was empty, so there was no base branch for a PR to target. |
| D-002 | 2026-11-10 | Local development uses Python 3.12 (Homebrew) in `.venv`; CI tests Python 3.11 and 3.12. | Claude | The brief requires 3.11+; the Mac has 3.12 and 3.13 but no 3.11. CI covers the 3.11 floor. |
| D-003 | 2026-11-10 | CI triggers on pushes to `main`, on every pull request, and manually; not on every branch push. | Claude | Covers "on push and PR" without running every PR commit twice on a private repository's limited Actions minutes. |
| D-004 | 2026-11-10 | Claude merges its own Phase 1 PRs once CI is green, with merge commits (not squash) so the small commits stay in history. The Phase 2 step d pre-registration PR is merged by the user, as the addendum requires. | Claude | Items depend on each other; stacking nine unmerged PRs would make review harder. The addendum reserves only the pre-registration merge for the user. |
| D-005 | 2026-11-10 | `scripts/check_language.py` enforces the brief's forbidden wording on every tracked file and every commit message, in CI and in the test suite. `docs/BRIEF.md` is exempt because it quotes the banned phrases. | Claude | Makes the wording rule mechanical instead of relying on review. |
| D-006 | 2026-11-10 | Supersedes D-004: Claude does **not** merge PRs. Phase 1 PRs are stacked: each item's branch is cut from the previous item's branch and its PR targets that branch, so every PR shows only its own diff. The user merges them in order (#4 first); GitHub retargets the next PR to `main` when its base branch is merged and deleted. | Claude (forced) | The permission system denied Claude's attempt to merge #4 ("merge without review"). Merging is the user's call. |
| D-007 | 2026-11-10 | Locus effect sizes: rs9399137 and rs7482144 from Lettre et al. 2008 Table 1 (CSSCD, per-minor-allele effects in SD units: 0.571 and 0.407). rs1427407 is not in Lettre 2008, so its effect is the 7.6% of HbF variance explained in the Wonkam et al. 2025 global meta-analysis (n = 3751, African ancestries, MAF 0.24); its HbF-increasing allele (T) comes from Bauer et al. 2013 Fig. 2B. | Claude | These were the only sources found that give a standardised effect (SD units or variance explained) for each sentinel in African-ancestry cohorts. Every value stays `verified: false` until the user checks it (issue #3). |
| D-008 | 2026-11-10 | The 1000 Genomes AFR frequency of rs7482144 is left null. | Claude | Ensembl GRCh37 and GRCh38 hold no 1000 Genomes phase 3 frequencies for it, and dbSNP has none either. The brief says unverifiable values stay null. Phase 2 must look the site up by position as well as by rsID (Ensembl lists rs41322148 as a synonym). |
