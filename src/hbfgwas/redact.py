"""Redact run records at write time: no hostname, no absolute local paths.

Every log and record that a script writes into a run directory goes through
:func:`redact_text` before it reaches disk:

* the machine's hostname is removed (plink2's ``Hostname:`` line is kept as
  ``Hostname: [redacted]``; the name itself becomes ``[host]`` anywhere else);
* paths inside the repository become repo-relative
  (``/abs/path/to/hbf-gwas/data/1kg/afr`` -> ``data/1kg/afr``; the repository
  itself -> ``.``);
* named directories outside the repository (``DATA_DIR``, ``RUN_DIR`` when
  they are elsewhere) become ``[DATA_DIR]`` / ``[RUN_DIR]``, and the system
  temporary directory becomes ``[TMP]``;
* the home directory becomes ``~``;
* any remaining ``/Users/<name>`` or ``/home/<name>`` becomes
  ``/Users/[user]`` / ``/home/[user]``, and the login name inside any other
  absolute path becomes ``[user]``.

Used by ``scripts/_common.sh`` (``python -m hbfgwas.redact``) for plink2 logs
and for each script's stdout and stderr, and by the Python writers of run
records. Records committed before this module existed are not edited.
"""

from __future__ import annotations

import argparse
import contextlib
import getpass
import io
import os
import platform
import re
import socket
import sys
import tempfile
from collections.abc import Iterable, Iterator
from pathlib import Path

NOTE_FILE = "REDACTION.md"
NOTE = """# Redaction of this run's records

The logs and records in this directory were redacted **when they were
written** (`hbfgwas.redact`): the machine's hostname is removed (plink2's
`Hostname:` line reads `Hostname: [redacted]`), paths inside the repository
are repo-relative, `DATA_DIR` / `RUN_DIR` outside the repository read
`[DATA_DIR]` / `[RUN_DIR]`, the system temporary directory reads `[TMP]`, the
home directory reads `~`, and any other user name in a path reads `[user]`.
Nothing else is changed: counts, frequencies, sample IDs, seeds, timings and
software versions are as produced.
"""

MIN_NAME = 4  # shorter host or user names are too likely to occur by chance


def repo_root() -> Path:
    """The repository root: $HBFGWAS_REPO_ROOT, else the checkout holding this package."""
    env = os.environ.get("HBFGWAS_REPO_ROOT")
    if env:
        return Path(env)
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pyproject.toml").is_file() and (parent / "src" / "hbfgwas").is_dir():
            return parent
    return Path.cwd()


def local_hostnames() -> list[str]:
    names = set()
    for n in (socket.gethostname(), platform.node()):
        if n:
            names.add(n)
            names.add(n.split(".")[0])
    return sorted((n for n in names if len(n) >= MIN_NAME), key=len, reverse=True)


def local_user() -> str | None:
    try:
        user = getpass.getuser()
    except Exception:  # noqa: BLE001 - no login name is fine
        return None
    return user if user and len(user) >= MIN_NAME else None


def _variants(path: Path | str) -> list[str]:
    """The path as given and as resolved (macOS: /var vs /private/var)."""
    raw = str(path).rstrip("/") or "/"
    out = {raw}
    with contextlib.suppress(OSError):
        out.add(os.path.realpath(raw))
    return [p for p in out if p not in ("", "/")]


def _inside(path: str, root: str) -> bool:
    return path == root or path.startswith(root + "/")


def redact_text(
    text: str,
    root: Path | str | None = None,
    aliases: dict[str, Path | str] | None = None,
    hostnames: Iterable[str] | None = None,
    home: Path | str | None = None,
    user: str | None = None,
    tmp: Iterable[Path | str] | None = None,
) -> str:
    """Return ``text`` with hostnames and absolute local paths removed."""
    root = repo_root() if root is None else root
    if aliases is None:
        aliases = {k: os.environ[k] for k in ("DATA_DIR", "RUN_DIR") if os.environ.get(k)}
    hostnames = local_hostnames() if hostnames is None else list(hostnames)
    home = Path.home() if home is None else home
    user = local_user() if user is None else user
    tmp = [tempfile.gettempdir(), "/tmp"] if tmp is None else list(tmp)

    # 1. hostname
    text = re.sub(r"(?m)^(Hostname:)[ \t]*\S.*$", r"\1 [redacted]", text)
    for h in sorted(hostnames, key=len, reverse=True):
        if len(h) >= MIN_NAME:
            text = re.sub(rf"(?<![\w.-]){re.escape(h)}(?![\w-])", "[host]", text)

    # 2. path prefixes, longest first: repository, named directories, home
    roots = _variants(root)
    prefixes: list[tuple[str, str, str]] = []  # (prefix, replacement, exact replacement)
    for r in roots:
        prefixes.append((r, "", "."))
    for name, p in (aliases or {}).items():
        for v in _variants(p):
            if any(_inside(v, r) for r in roots):
                continue  # already covered: becomes repo-relative
            prefixes.append((v, f"[{name}]/", f"[{name}]"))
    for t in tmp:
        for v in _variants(t):
            if not any(_inside(v, r) for r in roots):
                prefixes.append((v, "[TMP]/", "[TMP]"))
    for v in _variants(home):
        prefixes.append((v, "~/", "~"))
    prefixes.sort(key=lambda x: len(x[0]), reverse=True)
    for prefix, sub, exact in prefixes:
        pattern = re.compile(rf"(?<![\w.~-]){re.escape(prefix)}(/|(?![\w.-]))")
        text = pattern.sub(lambda m, s=sub, e=exact: s if m.group(1) == "/" else e, text)

    # 3. any other user directory, and the login name inside absolute paths
    text = re.sub(r"(?<![\w.-])(/Users|/home)/(?!\[user\])[^/\s\"'`]+", r"\1/[user]", text)
    if user:

        def _scrub(m: re.Match) -> str:
            return re.sub(rf"(?<![A-Za-z0-9]){re.escape(user)}(?![A-Za-z0-9])", "[user]", m[0])

        text = re.sub(r"[^\s\"'`]*/[^\s\"'`]*", _scrub, text)  # any token with a slash
    return text


def write_record(path: Path | str, text: str, **kwargs) -> Path:
    """Write a run record, redacted. ``kwargs`` go to :func:`redact_text`."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(redact_text(text, **kwargs))
    return path


def write_note(run_dir: Path | str) -> Path:
    """State in the run directory that its records are redacted (idempotent)."""
    p = Path(run_dir) / NOTE_FILE
    p.parent.mkdir(parents=True, exist_ok=True)
    if not p.exists():
        p.write_text(NOTE)
    return p


class _Tee(io.TextIOBase):
    def __init__(self, stream, buf: io.StringIO):
        self._stream, self._buf = stream, buf

    def write(self, s: str) -> int:
        self._buf.write(s)
        return self._stream.write(s)

    def flush(self) -> None:
        self._stream.flush()


@contextlib.contextmanager
def capture(run_dir: Path | str, name: str, **kwargs) -> Iterator[None]:
    """Show stdout/stderr as usual and save them, redacted, as RUN_DIR/logs/<name>.*.txt."""
    out, err = io.StringIO(), io.StringIO()
    old = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = _Tee(old[0], out), _Tee(old[1], err)
    try:
        yield
    finally:
        sys.stdout, sys.stderr = old
        logs = Path(run_dir) / "logs"
        write_record(logs / f"{name}.stdout.txt", out.getvalue(), **kwargs)
        write_record(logs / f"{name}.stderr.txt", err.getvalue(), **kwargs)
        write_note(run_dir)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("src", nargs="?", help="input file (default: stdin)")
    p.add_argument("dest", nargs="?", help="output file (default: stdout)")
    p.add_argument("--repo-root", default=None)
    p.add_argument(
        "--alias", action="append", default=[], metavar="NAME=PATH", help="e.g. DATA_DIR=/x"
    )
    p.add_argument("--note", metavar="RUN_DIR", help="write RUN_DIR/REDACTION.md and exit")
    a = p.parse_args(argv)
    if a.note:
        write_note(a.note)
        return 0
    aliases = dict(x.split("=", 1) for x in a.alias)
    text = Path(a.src).read_text(errors="replace") if a.src else sys.stdin.read()
    kwargs = {"root": a.repo_root, "aliases": aliases}
    if a.dest:
        write_record(a.dest, text, **kwargs)
    else:
        sys.stdout.write(redact_text(text, **kwargs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
