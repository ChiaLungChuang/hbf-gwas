# Publication plan

**A plan only. Nothing here has been executed.** No history has been
rewritten, no repository created, and no visibility changed. Each step that
creates, rewrites or pushes anything needs the user's go-ahead at the time.

## The decision (user, 2026-11-12, D-045)

* This repository stays **private** and becomes the **archive**.
* At the end of the project, a **new public repository** is published,
  built from a **scrubbed copy of this history**. The scrub removes the build
  machine's hostname and local user paths from every commit.
* Licence: **MIT**, copyright "Chia-Lung Chuang", one licence for code, docs
  and figures (section 5).

## 1. What the scrub removes, and where it is

The history scan of 2026-11-12 found no credentials, no genotype data and no
file over 1 MB. It found the following, all in run records written before
`hbfgwas.redact` (D-043) existed:

| What | Where (on `main` and the open branches) |
|---|---|
| the build machine's hostname (plink2's `Hostname:` line) | 13 plink2 logs in `runs/2026-11-11/logs/`, 3 in `runs/2026-11-12-stepd/logs/`. One older commit of `tests/test_redact.py` (7dd2bb9, on #33's branch) also used it as a test value; it was replaced by a made-up name in 14e3f41. |
| the absolute path of this checkout, under the user's home directory | the same logs (about 280 occurrences, mostly plink2 `--pfile`/`--out` arguments and `Working directory:`), and `config_path` in `runs/2026-11-10-sim-stratification/results.json`, `runs/2026-11-10-sim-v2/results.json` and `runs/2026-11-11-sim-v3/results.json` |

Not found anywhere in history: the user name or hostname in a commit message;
any other user's paths (except made-up test values such as `/Users/alice` in
`tests/test_redact.py`).

**Not part of the scrub (separate decisions, section 6):**

* the author name and email on commits;
* the `Co-Authored-By: Claude` trailers;
* `docs/BRIEF.md` and `CLAUDE.md`.

This document deliberately does not spell out the hostname or the home path.
If it did, the scrub would rewrite this file too. The expressions file is
built from the live values in step 3.

## 2. What is lost, and how the README will say so

A new repository starts with commits only. Lost:

1. **Issues and pull requests.** That means their titles, descriptions,
   comments and review threads. They hold the step-by-step reports, the
   numbers each PR reported and the user's approvals.
2. **GitHub's own timestamps.** These are when each PR was opened and merged,
   when the `prereg-v1` tag was pushed, and every CI run with its logs.
   Commit dates survive, but they are set by the committer's clock. GitHub's
   record that the pre-registration was pushed before step e is the stronger
   evidence, and it stays in the archive.
3. **Commit hashes.** Every hash changes. Hashes quoted in docs and run
   records (e.g. `f63e8bd`, `c2ca787`, `git_commit` in `results.json`) will
   point to the archive. The step 5 commit map translates them.
4. **`#N` references.** The docs and commit messages mention issues and PRs
   by number (#3, #14, #30, ...). In the new repository GitHub would link
   them to unrelated issues or PRs with the same number.
5. **Unmerged branches** (only `main` and tags are published), and the
   "Verified" badge on GitHub's web merge commits (rewritten commits lose
   their signatures).

**Recommended mitigation for 1, 2 and 4:** export the text of every issue and
PR to `docs/archive/` in the new repository. That means title, body,
comments, opened, merged and closed times, and the merge commit (old hash).
Then "#14" can still be looked up. The export goes through the same scrub.
The bodies hold only `~/...` paths; check anyway.

**README note**, added in the first new commit:

> **About this repository's history.** Published on <date> from a private
> working repository, which is kept as an archive. The history is that
> repository's `main` branch, with two kinds of text removed from every
> commit: the build machine's hostname and local user paths (see
> `docs/publication-plan.md`). Nothing else was changed except the two
> checksums in `tests/test_records.py` that pin the run records whose paths
> were removed. Commit dates and messages are unchanged. Commit hashes
> differ; old to new is in `docs/archive-commit-map.tsv`. Issues, pull
> requests and CI runs, where the step-by-step reports and the user's
> approvals were recorded, are not included. Their text is in
> `docs/archive/` [if exported], and references such as #14 refer to it. The
> archive also holds GitHub's record of when the pre-registration tag was
> pushed.

## 3. Exact steps

Everything happens in a new directory **outside** this checkout. Nothing is
pushed until section 4 passes and the user says go.

**0. Preconditions:**
- All PRs are merged, and `main` on GitHub equals the intended final tree.
  Verify the tree, as after every stacked merge.
- The step e outcome is committed (or the user decides to publish before it).
- `prereg-v1` exists on GitHub, if step e has run.

**1. A bare clone of `main` and the tags only.** A bare clone does not copy
GitHub's `refs/pull/*`. Then delete every other branch.

```bash
mkdir -p ~/hbf-gwas-publish && cd ~/hbf-gwas-publish
git clone --bare https://github.com/ChiaLungChuang/hbf-gwas.git archive.git
cp -R archive.git scrubbed.git          # archive.git stays untouched for comparison
cd scrubbed.git
git for-each-ref --format='%(refname)' refs/heads | grep -vx 'refs/heads/main' \
  | xargs -n1 git update-ref -d
git for-each-ref                         # expect: refs/heads/main and the tags only
```

**2. Install git-filter-repo.** GitHub recommends this tool for this job. It
is a download, so the user approves it first. Put it in a throwaway virtual
environment:

```bash
python3 -m venv ~/hbf-gwas-publish/venv
~/hbf-gwas-publish/venv/bin/pip install git-filter-repo
```

**3. Build the expressions file from the live values**, so they are never
written into the repository. Run these from the private checkout. Longer
paths come first, so the repository path becomes repo-relative before the
home directory becomes `~`. These rules match what `hbfgwas.redact` now does
at write time.

```bash
cd ~/projects/hbf-gwas                   # the private checkout
host=$(grep -h -m1 '^Hostname:' runs/2026-11-11/logs/afr.log | cut -d' ' -f2)
repo=$(git rev-parse --show-toplevel)    # this checkout's absolute path
home=$HOME
cat > ~/hbf-gwas-publish/expressions.txt <<EOF
${repo}/==>
${repo}==>.
${home}==>~
Hostname: ${host}==>Hostname: [redacted]
${host}==>[host]
EOF
```

The `Hostname:` rule is a literal on purpose. A regex such as
`Hostname: .*` would also rewrite source and docs that contain the word
`Hostname:` (`src/hbfgwas/redact.py`, D-043).

**4. Keep the pinned record checksums true.** `tests/test_records.py` pins the
sha256 of `runs/2026-11-10-sim-stratification/results.json` and
`runs/2026-11-10-sim-v2/results.json`. The scrub changes both files (their
`config_path`). Compute each file's checksum **after** the same replacements
and add `old==>new` lines, so the test pins the scrubbed bytes:

```bash
cd ~/projects/hbf-gwas
for f in runs/2026-11-10-sim-stratification/results.json runs/2026-11-10-sim-v2/results.json; do
  old=$(shasum -a 256 "$f" | cut -d' ' -f1)
  new=$(sed -e "s#${repo}/##g" -e "s#${repo}#.#g" -e "s#${home}#~#g" "$f" | shasum -a 256 | cut -d' ' -f1)
  echo "${old}==>${new}" >> ~/hbf-gwas-publish/expressions.txt
done
```

**Alternative (user's choice):** leave these two files byte-identical and
unscrubbed. Their only personal content is the home path in `config_path`;
there is no hostname. That keeps "records are never edited" literally true,
at the cost of the user name appearing twice. Recommendation: scrub, and say
so in the README note above.

**5. Rewrite.** Messages go through the same expressions; the scan found
nothing there, so this is a safety net. If the user picks it, the author
email is changed in the same pass (section 6).

```bash
cd ~/hbf-gwas-publish/scrubbed.git
~/hbf-gwas-publish/venv/bin/git-filter-repo \
  --replace-text ../expressions.txt --replace-message ../expressions.txt
cp filter-repo/commit-map ../archive-commit-map.tsv   # old hash -> new hash
```

**6. Check** everything in section 4. Stop here and show the user the
results.

**7. Create the public repository** (user's go-ahead). The name `hbf-gwas`
is taken by this archive. Either rename this one to `hbf-gwas-archive` first
(links in the docs to `github.com/ChiaLungChuang/hbf-gwas` then reach the
public one), or pick a new public name. Push only `main` and the tags, never
`--mirror`:

```bash
gh repo create ChiaLungChuang/<public-name> --public \
  --description "A small, tested GWAS QC and association pipeline on 1000 Genomes AFR genotypes with a simulated phenotype"
git push https://github.com/ChiaLungChuang/<public-name>.git main --tags
```

**8. First new commit in the public repository:**
- `LICENSE` (MIT, "Copyright (c) 2026 Chia-Lung Chuang");
- `license = {text = "MIT"}` in `pyproject.toml`;
- the README history note;
- `docs/archive-commit-map.tsv`;
- `docs/archive/` (if exported).

Its CI runs for free: standard GitHub runners cost nothing on public
repositories.

**9. Mark the archive.** Add a README pointer to the public repository in
this one, and optionally use GitHub's "Archive this repository", which
makes it read-only.

## 4. How the user can verify the scrub before anything is pushed

Run these in `~/hbf-gwas-publish/` after step 5. Each has an expected result.
Any other result means stop.

```bash
cd ~/hbf-gwas-publish/scrubbed.git
host=$(grep -h -m1 '^Hostname:' ~/projects/hbf-gwas/runs/2026-11-11/logs/afr.log | cut -d' ' -f2)
user=$(whoami)

# a. No hostname or user name in any file of any commit. Expect: 0
git rev-list --all | while read -r c; do git grep -I -c -E "${host}|${user}" "$c" -- ; done | wc -l

# b. ... or in any commit message. Expect: 0
git log --all --format=%B | grep -c -E "${host}|${user}"

# c. Files still containing "/Users/" hold only made-up or generic values. Expect exactly:
#    docs/publication-plan.md (this file's commands), src/hbfgwas/redact.py (docstring:
#    /Users/<name>), tests/test_redact.py (alice, bob, al)
git grep -I -l '/Users/' main --

# d. Only the expected files differ at the tip. Expect: the 27 log files under runs/*/logs/,
#    the three results.json, tests/test_records.py (two lines)
git --git-dir=../archive.git fetch -q "$PWD" main:refs/scrubbed/main
git --git-dir=../archive.git diff --stat main refs/scrubbed/main

# e. Read every changed line yourself
git --git-dir=../archive.git diff main refs/scrubbed/main | less

# f. Same number of commits, same dates and messages. Expect: no output from diff
diff <(git --git-dir=../archive.git log main --format='%ad %s') <(git log main --format='%ad %s')

# g. Refs to be published. Expect: refs/heads/main and the tags only
git for-each-ref

# h. Authors as decided in section 6
git log --all --format='%an <%ae> | %cn <%ce>' | sort | uniq -c

# i. The scrubbed tree still builds and passes, including tests/test_records.py
git clone -q "$PWD" ../check && cd ../check && REQUIRE_PLINK2=1 scripts/rebuild.sh
```

Finally, re-run the six-point history scan on the scrubbed copy:
credentials, large files, personal paths, authors, overstatements and
licence.

## 5. Licence

Chosen by the user on 2026-11-12: **MIT**, copyright "Chia-Lung Chuang",
covering code, docs and figures. It is added in the public repository's first
new commit (step 8), not into this archive or the rewritten history. The
options considered:

| Licence | In short |
|---|---|
| **MIT** (chosen) | Short and permissive: reuse with attribution |
| Apache-2.0 | Permissive, with an explicit patent grant; longer |
| BSD-3-Clause | Like MIT, and bars using the author's name to endorse derived work |
| GPL-3.0 | Copyleft: modified versions that are distributed stay GPL. plink2 (GPLv3) does not require it, because it is not redistributed |
| CC-BY-4.0 for docs and figures | A second licence alongside the code licence; not chosen |

The licence covers only this repository's own files. 1000 Genomes data are
not in it; the scripts download them under IGSR's terms.

## 6. Decisions still open for publication

1. **Author email** on the commits: `chialungchuang078@gmail.com`, already
   the public email on the GitHub profile. Either keep it, or rewrite it in
   step 5 to GitHub's noreply address with `--mailmap`. Future commits would
   then use that address too.
2. **Export issues and PRs** to `docs/archive/` (recommended), or say only
   that they exist in the archive.
3. **Scrub the two pinned records** (recommended, step 4), or keep them
   byte-identical.
4. **Public repository name** (step 7).
5. **`docs/BRIEF.md` and `CLAUDE.md`** stay as they are; the README's "How
   this was built" points to them. Removing them from history would be a
   second, larger rewrite and would hide how the work was done. Not
   recommended.
6. **When to publish:** after step e and its report, or before. Before
   means publishing without the real-data result.
