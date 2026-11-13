# Redaction of this run's records

The logs and records in this directory were redacted **when they were
written** (`hbfgwas.redact`): the machine's hostname is removed (plink2's
`Hostname:` line reads `Hostname: [redacted]`), paths inside the repository
are repo-relative, `DATA_DIR` / `RUN_DIR` outside the repository read
`[DATA_DIR]` / `[RUN_DIR]`, the system temporary directory reads `[TMP]`, the
home directory reads `~`, and any other user name in a path reads `[user]`.
Nothing else is changed: counts, frequencies, sample IDs, seeds, timings and
software versions are as produced.
