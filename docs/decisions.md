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
