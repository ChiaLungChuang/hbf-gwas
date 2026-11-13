#!/usr/bin/env python3
"""Benchmark step e on SIMULATED genotypes of the real size: run time, memory, disk.

No real genotypes are read. plink2 ``--dummy`` writes random genotypes for 652
samples and 14,698,808 variants (the real counts after QC) into a scratch
directory. Then each expensive part of ``scripts/step_e.sh`` is run on them,
under the same code, and timed:

  integrity  export of 20,000 null variants and the declared predictions
             recomputed (03c's work, primary analysis only)
  04         plink2 --glm without covariates, and with 10 PCs; concordance export
  04b        2000 phenotype redraws tested at 2 sentinels
  04c        the two scans again without one sample: assumed to cost the same as
             04's scans (not run)
  05         reading both genome-wide results, lambda_GC, type I error, QQ figure

Each step runs in its own process. Wall time is measured here; peak memory
(maximum resident set size) comes from /usr/bin/time. The scratch directory is
deleted at the end unless --keep is given.

    python scripts/benchmark_step_e.py --work /path/to/scratch --out runs/<date>/benchmark.json
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

N_REAL, M_REAL = 652, 14_698_808
COUNTS = {"ACB": 95, "ASW": 55, "ESN": 99, "GWD": 113, "LWK": 97, "MSL": 85, "YRI": 108}
SEED = 20261114  # this benchmark's own seed; it touches no real data


def _time_cmd() -> list[str]:
    return ["/usr/bin/time", "-l"] if sys.platform == "darwin" else ["/usr/bin/time", "-v"]


def _peak_mb(stderr: str) -> float | None:
    m = re.search(r"(\d+)\s+maximum resident set size", stderr)  # macOS: bytes
    if m:
        return int(m.group(1)) / 2**20
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", stderr)  # GNU time
    return int(m.group(1)) / 2**10 if m else None


def run(name: str, cmd: list[str], work: Path, outputs: list[str]) -> dict:
    t0 = time.perf_counter()
    r = subprocess.run(_time_cmd() + cmd, cwd=work, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if r.returncode != 0:
        raise SystemExit(f"{name} failed:\n{r.stdout[-2000:]}\n{r.stderr[-4000:]}")
    size = sum((work / o).stat().st_size for o in outputs if (work / o).exists())
    row = {
        "step": name,
        "seconds": round(wall, 1),
        "peak_memory_mb": None if _peak_mb(r.stderr) is None else round(_peak_mb(r.stderr)),
        "output_bytes": size,
        "outputs": outputs,
    }
    print(
        f"{name:28s} {row['seconds']:8.1f} s  {row['peak_memory_mb'] or 0:7.0f} MB  "
        f"{size / 2**20:9.1f} MB written",
        flush=True,
    )
    return row


# --- python steps, each run in its own process --------------------------------------


def step_predictions(work: Path) -> None:
    """03c's computation on 20,000 exported null variants (the integrity check)."""
    import numpy as np

    from hbfgwas import predict
    from hbfgwas.pipeline import read_raw

    rng = np.random.default_rng(SEED)
    raw = read_raw(work / "null.raw")
    g = raw.genotypes
    g = g[:, ~np.isnan(g).any(axis=0)]
    n = g.shape[0]
    pop = np.repeat(np.array(sorted(COUNTS)), [COUNTS[k] for k in sorted(COUNTS)])[:n]
    pcs = rng.standard_normal((n, 10))
    sentinels = g[:, :2]
    betas = np.array([0.456435, 0.571])
    labels, idx = np.unique(pop, return_inverse=True)
    shift = rng.standard_normal(labels.size)[idx]
    shift = (shift - shift.mean()) / shift.std() * np.sqrt(0.2)
    for cov in (pcs, None):
        predict.conditional_power(sentinels, betas, 5e-8, covariates=cov, shift=shift)
    for cov in (None, pcs):
        predict.predicted_lambda(g, pop, 0.2, covariates=cov, n_shift_draws=200, shift=shift)
    if hasattr(predict, "predicted_null_distribution"):  # the step d amendment's method
        noise = predict.noise_variance(sentinels, betas, shift_var_fraction=float(shift.var()))
        ind = (pop[:, None] == labels[None, :]).astype(float)
        predict.predicted_null_distribution(
            g,
            sentinels @ betas + shift,
            noise,
            np.hstack([ind, pcs]),
            pcs,
            n_draws=2000,
            rng=rng,
            n_effective=215_522,
        )


def step_redraws(work: Path) -> None:
    """04b's computation: 2000 redraws tested at 2 sentinels, with and without PCs."""
    import numpy as np

    from hbfgwas import assoc, simulate

    rng = np.random.default_rng(SEED)
    pop = np.repeat(np.array(sorted(COUNTS)), [COUNTS[k] for k in sorted(COUNTS)])
    g = rng.binomial(2, [0.24, 0.04], size=(pop.size, 2)).astype(float)
    shift = {k: float(v) for k, v in zip(sorted(COUNTS), rng.standard_normal(7) * 0.4, strict=True)}
    reps = simulate.simulate_phenotype(
        g, pop, rng=rng, beta=np.array([0.456435, 0.571]), pop_shift=shift, n_reps=2000
    )
    pcs = rng.standard_normal((pop.size, 10))
    assoc.linear_assoc(g, reps.y, pcs)
    assoc.linear_assoc(g, reps.y)


def step_report(work: Path) -> None:
    """05's work on the two genome-wide results: read, lambda_GC, type I error, QQ figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    from hbfgwas import assoc
    from hbfgwas.pipeline import null_mask, read_glm

    glm = {
        "unadjusted": read_glm(work / "glm_nopc.PHENO1.glm.linear"),
        "adjusted": read_glm(work / "glm_pc.PHENO1.glm.linear"),
    }
    planted = pd.DataFrame({"chrom": ["1"], "pos": [1_000_000], "id": ["x"]})
    for df in glm.values():
        ok = df["P"].notna().to_numpy()
        assoc.lambda_gc(df.loc[ok, "P"])
        assoc.type1_error(df.loc[null_mask(df, planted, 1_000_000) & ok, "P"], 0.05)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharex=True, sharey=True)
    for ax, side in zip(axes, glm, strict=True):
        assoc.qq_plot(ax, glm[side]["P"].to_numpy(), title=side, max_points=60_000)
    fig.savefig(work / "qq.png", dpi=150)
    plt.close(fig)


STEPS = {"predictions": step_predictions, "redraws": step_redraws, "report": step_report}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--work", help="scratch directory for the simulated fileset (large)")
    ap.add_argument("--out", help="benchmark JSON (a run record)")
    ap.add_argument("--samples", type=int, default=N_REAL)
    ap.add_argument("--variants", type=int, default=M_REAL)
    ap.add_argument("--keep", action="store_true", help="keep the scratch directory")
    ap.add_argument("--step", choices=sorted(STEPS), help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    if a.step:  # child process for one python step
        STEPS[a.step](Path(a.work))
        return 0

    plink = os.environ.get("PLINK2", "plink2")
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    py = [sys.executable, str(Path(__file__).resolve()), "--work", str(work), "--step"]
    rows = []

    rows.append(
        run(
            "setup: plink2 --dummy",
            [
                plink,
                "--dummy",
                str(a.samples),
                str(a.variants),
                "0",
                "0",
                "acgt",
                "scalar-pheno",
                "--seed",
                str(SEED),
                "--make-pgen",
                "--out",
                "dummy",
            ],
            work,
            ["dummy.pgen", "dummy.pvar", "dummy.psam"],
        )
    )
    import numpy as np

    rng = np.random.default_rng(SEED)
    iids = [ln.split("\t")[0] for ln in (work / "dummy.psam").read_text().splitlines()[1:]]
    pcs = rng.standard_normal((len(iids), 10))
    with open(work / "pcs.tsv", "w") as fh:
        fh.write("#IID\t" + "\t".join(f"PC{i}" for i in range(1, 11)) + "\n")
        for iid, row in zip(iids, pcs, strict=True):
            fh.write(iid + "\t" + "\t".join(f"{x:.6f}" for x in row) + "\n")
    base = [plink, "--pfile", "dummy"]
    rows.append(
        run(
            "integrity: export 20,000 null",
            base + ["--thin-count", "20000", "--seed", str(SEED), "--export", "A", "--out", "null"],
            work,
            ["null.raw"],
        )
    )
    rows.append(run("integrity: predictions", py + ["predictions"], work, []))
    rows.append(
        run(
            "04: --glm, no covariates",
            base + ["--pheno-name", "PHENO1", "--glm", "allow-no-covars", "--out", "glm_nopc"],
            work,
            ["glm_nopc.PHENO1.glm.linear"],
        )
    )
    rows.append(
        run(
            "04: --glm, 10 PCs",
            base
            + [
                "--pheno-name",
                "PHENO1",
                "--covar",
                "pcs.tsv",
                "--covar-name",
                "PC1-PC10",
                "--covar-variance-standardize",
                "--glm",
                "hide-covar",
                "--out",
                "glm_pc",
            ],
            work,
            ["glm_pc.PHENO1.glm.linear"],
        )
    )
    rows.append(
        run(
            "04: concordance export",
            base
            + [
                "--thin-count",
                "2000",
                "--seed",
                str(SEED + 1),
                "--export",
                "A",
                "--out",
                "concord",
            ],
            work,
            ["concord.raw"],
        )
    )
    rows.append(run("04b: 2000 redraws", py + ["redraws"], work, []))
    rows.append(run("05: report", py + ["report"], work, ["qq.png"]))

    by = {r["step"]: r for r in rows}
    scans = by["04: --glm, no covariates"]["seconds"] + by["04: --glm, 10 PCs"]["seconds"]
    scan_bytes = (
        by["04: --glm, no covariates"]["output_bytes"] + by["04: --glm, 10 PCs"]["output_bytes"]
    )
    step_rows = [r for r in rows if not r["step"].startswith("setup")]
    total_s = sum(r["seconds"] for r in step_rows) + scans  # + 04c, assumed equal to 04's scans
    total_bytes = sum(r["output_bytes"] for r in step_rows) + scan_bytes
    peak = max(r["peak_memory_mb"] or 0 for r in step_rows)
    result = {
        "what": "step e benchmark on SIMULATED genotypes (plink2 --dummy); no real data",
        "seed": SEED,
        "samples": a.samples,
        "variants": a.variants,
        "machine": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "cpu_count": os.cpu_count(),
        },
        "plink2": subprocess.run(
            [plink, "--version"], capture_output=True, text=True
        ).stdout.strip(),
        "steps": rows,
        "assumed_not_run": {
            "04c: two scans without one sample": {
                "seconds": round(scans, 1),
                "output_bytes": scan_bytes,
                "basis": "the same as 04's two scans (651 instead of 652 samples)",
            }
        },
        "step_e_expected": {
            "seconds": round(total_s, 1),
            "minutes": round(total_s / 60, 1),
            "peak_memory_mb": peak,
            "disk_written_bytes": total_bytes,
            "disk_written_gb": round(total_bytes / 1e9, 2),
            "excludes": "the genotype files themselves (already on disk after step c)",
        },
        "setup_fileset_bytes": by["setup: plink2 --dummy"]["output_bytes"],
    }
    if a.out:
        from hbfgwas import redact

        redact.write_record(a.out, json.dumps(result, indent=2) + "\n")
        print(f"wrote {a.out}")
    e = result["step_e_expected"]
    print(
        f"step e expected: {e['minutes']} min, peak {e['peak_memory_mb']} MB, "
        f"{e['disk_written_gb']} GB written"
    )
    if not a.keep:
        shutil.rmtree(work)
    return 0


if __name__ == "__main__":
    sys.exit(main())
