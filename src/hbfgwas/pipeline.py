"""Helpers the real-data shell scripts call (``python -m hbfgwas.pipeline <command>``).

Commands
--------
config-get        print one value from config/analysis.yaml (dotted key)
sentinels         find the sentinel variants in a .pvar (by rsID, else by GRCh37 position)
sentinel-report   allele frequency and QC survival of each sentinel
sample-ids        reproducible random subset of variant IDs (plus required IDs; optionally
                  away from the sentinels)
predictions       Phase 2 step d: conditional power and predicted lambda_GC from the real
                  genotypes and PCs, before any phenotype
simulate-phenotype  simulate the HbF-like phenotype on genotypes exported with
                  plink2 --export A (the phenotype is always simulated)
het-summary       summary of a plink2 --het report
report            lambda_GC, type I error, locus recovery, PC plot, QQ plot and the
                  Python-vs-plink2 concordance check (used by scripts/05_report.py)
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from hbfgwas import assoc, power, simulate
from hbfgwas.config import load_config

# Review item C (D-028): no proxies; a sentinel absent from the callset or failing QC is
# reported with this label (also in config/analysis.yaml, sentinels.not_testable_label).
NOT_TESTABLE = "not testable in this callset"

# --- config ----------------------------------------------------------------


def config_get(config: dict, key: str):
    node = config
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(f"{key!r} not found in the config")
        node = node[part]
    return node


def format_scalar(value) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, (dict, list)):
        return json.dumps(value)
    return str(value)


# --- file readers ------------------------------------------------------------


def read_table(path: str | Path) -> pd.DataFrame:
    """plink2-style TSV with a '#'-prefixed header line (.psam, .eigenvec, .afreq, .id)."""
    df = pd.read_csv(path, sep=r"\s+", dtype=str, comment=None)
    df.columns = [c.lstrip("#") for c in df.columns]
    return df


def read_psam(path) -> pd.DataFrame:
    return read_table(path)


def read_eigenvec(path) -> pd.DataFrame:
    df = read_table(path)
    for c in df.columns:
        if c.startswith("PC"):
            df[c] = df[c].astype(float)
    return df


def read_ids(path) -> set[str]:
    """Variant IDs from a .snplist / .prune.in / one-ID-per-line file."""
    with open(path) as fh:
        return {line.split()[0] for line in fh if line.strip() and not line.startswith("#")}


def iter_pvar(path):
    """Yield (chrom, pos, id, ref, alt) from a .pvar."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 5)
            yield f[0], int(f[1]), f[2], f[3], f[4]


@dataclass
class RawGenotypes:
    fid: np.ndarray
    iid: np.ndarray
    ids: list[str]
    counted: list[str]
    genotypes: np.ndarray  # (n, m) float, NaN = missing


def read_raw(path) -> RawGenotypes:
    """plink2 --export A: FID IID PAT MAT SEX PHENOTYPE then <ID>_<counted allele> columns."""
    df = pd.read_csv(path, sep=r"\s+", dtype=str)
    cols = list(df.columns)
    if "IID" not in cols or "PHENOTYPE" not in cols:
        raise ValueError(f"{path}: expected IID and PHENOTYPE columns in a plink2 .raw file")
    variant_cols = cols[cols.index("PHENOTYPE") + 1 :]
    ids, counted = [], []
    for col in variant_cols:
        vid, allele = col.rsplit("_", 1)
        ids.append(vid)
        counted.append(allele)
    geno = df[variant_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    fid = df["FID"].to_numpy() if "FID" in cols else df["IID"].to_numpy()
    return RawGenotypes(fid, df["IID"].to_numpy(), ids, counted, geno)


# --- sentinels ---------------------------------------------------------------


def _has_id(v) -> bool:
    return isinstance(v, str) and v != ""


def _norm_chrom(c) -> str:
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def find_sentinels(pvar_path, loci: pd.DataFrame) -> pd.DataFrame:
    """Locate each sentinel by rsID, else by (chromosome, GRCh37 position)."""
    want_ids = {str(r) for r in loci["rsid"] if pd.notna(r)}
    want_pos = {}
    for _, r in loci.iterrows():
        if pd.notna(r["chromosome"]) and pd.notna(r["position_grch37"]):
            want_pos[(_norm_chrom(r["chromosome"]), int(float(r["position_grch37"])))] = r["rsid"]
    by_id, by_pos = {}, {}
    for chrom, pos, vid, ref, alt in iter_pvar(pvar_path):
        rec = (_norm_chrom(chrom), pos, vid, ref, alt)
        if vid in want_ids:
            by_id.setdefault(vid, rec)
        if (rec[0], pos) in want_pos:
            by_pos.setdefault((rec[0], pos), rec)
    rows = []
    for _, r in loci.iterrows():
        table_pos = (
            (_norm_chrom(r["chromosome"]), int(float(r["position_grch37"])))
            if pd.notna(r["chromosome"]) and pd.notna(r["position_grch37"])
            else None
        )
        rec, found_by = None, "not found"
        if r["rsid"] in by_id:
            rec, found_by = by_id[r["rsid"]], "rsid"
        elif table_pos is not None and table_pos in by_pos:
            rec, found_by = by_pos[table_pos], "position"
        row = {
            "gene": r["gene"],
            "rsid": r["rsid"],
            "table_chromosome": None if table_pos is None else table_pos[0],
            "table_position": None if table_pos is None else table_pos[1],
            "effect_allele": r.get("effect_allele"),
            "other_allele": r.get("other_allele"),
            "found_by": found_by,
            "id": None,
            "chrom": None,
            "pos": None,
            "ref": None,
            "alt": None,
            "position_matches_table": None,
            "alleles_match_table": None,
        }
        if rec is not None:
            row.update(id=rec[2], chrom=rec[0], pos=rec[1], ref=rec[3], alt=rec[4])
            if table_pos is not None:
                row["position_matches_table"] = (rec[0], rec[1]) == table_pos
            if pd.notna(row["effect_allele"]) and pd.notna(row["other_allele"]):
                row["alleles_match_table"] = {rec[3], rec[4]} == {
                    row["effect_allele"],
                    row["other_allele"],
                }
        rows.append(row)
    return pd.DataFrame(rows)


def sentinel_report(
    smap: pd.DataFrame, afreq: pd.DataFrame | None, steps: dict[str, set[str]]
) -> pd.DataFrame:
    out = smap.copy()
    out["alt_freq"] = np.nan
    out["obs_ct"] = np.nan
    if afreq is not None and len(afreq):
        f = afreq.set_index("ID")
        for i, r in out.iterrows():
            if r["id"] in f.index:
                out.at[i, "alt_freq"] = float(f.at[r["id"], "ALT_FREQS"])
                out.at[i, "obs_ct"] = float(f.at[r["id"], "OBS_CT"])
    eff = []
    for _, r in out.iterrows():
        if pd.isna(r["alt_freq"]) or pd.isna(r["effect_allele"]):
            eff.append(np.nan)
        elif r["effect_allele"] == r["alt"]:
            eff.append(r["alt_freq"])
        elif r["effect_allele"] == r["ref"]:
            eff.append(1.0 - r["alt_freq"])
        else:
            eff.append(np.nan)
    out["effect_allele_freq"] = eff
    for name in steps:
        out[f"passed_{name}"] = False
    removed_by = []
    for i, r in out.iterrows():
        first = None
        for name, ids in steps.items():
            ok = _has_id(r["id"]) and r["id"] in ids
            out.at[i, f"passed_{name}"] = ok
            if not ok and first is None:
                first = name
        if not _has_id(r["id"]):
            first = "not in data"
        removed_by.append(first or "")
    out["survived_qc"] = [rb == "" for rb in removed_by]
    out["removed_by"] = removed_by
    out["status"] = ["testable" if ok else NOT_TESTABLE for ok in out["survived_qc"]]
    return out


# --- variant subset ------------------------------------------------------------


def sample_ids(
    pvar_path,
    n: int,
    seed: int,
    include: set[str] | None = None,
    exclude_windows: list[tuple[str, int]] | None = None,
    window_bp: int = 0,
) -> list[str]:
    """``n`` variant IDs chosen uniformly without replacement (two passes), plus ``include``.

    Variants within ``window_bp`` of any (chrom, pos) in ``exclude_windows`` are
    never chosen (used to draw null variants away from the planted sentinels).
    """
    excl = [(_norm_chrom(c), int(p)) for c, p in (exclude_windows or [])]

    def eligible(chrom, pos):
        c = _norm_chrom(chrom)
        return not any(c == ec and abs(pos - ep) <= window_bp for ec, ep in excl)

    idx = [i for i, (c, p, _, _, _) in enumerate(iter_pvar(pvar_path)) if eligible(c, p)]
    rng = np.random.default_rng(seed)
    pick = set(np.asarray(idx)[rng.choice(len(idx), size=min(n, len(idx)), replace=False)].tolist())
    chosen = [vid for i, (_, _, vid, _, _) in enumerate(iter_pvar(pvar_path)) if i in pick]
    for vid in sorted(include or ()):
        if vid not in chosen:
            chosen.append(vid)
    return chosen


# --- step e seed and the fixed population shift ------------------------------------

# One step e seed (config/step_e_seed.yaml) is spawned, in this order, into:
STEP_E_SPAWN = ("phenotype", "power_reps", "concordance")


def read_step_e_seed(path) -> int | None:
    import yaml

    with open(path) as fh:
        value = (yaml.safe_load(fh) or {}).get("seed")
    return None if value is None else int(value)


def step_e_seeds(path) -> dict:
    seed = read_step_e_seed(path)
    if seed is None:
        raise ValueError(
            f"{path}: the step e seed is empty; it is pre-registered only after the "
            "step d pre-registration is merged and tagged prereg-v1"
        )
    children = np.random.SeedSequence(seed).spawn(len(STEP_E_SPAWN))
    return {"seed": seed, **dict(zip(STEP_E_SPAWN, children, strict=True))}


def seed_child_int(path, name: str) -> int:
    """A plain integer seed from a spawned child (for tools that take an int)."""
    return int(step_e_seeds(path)[name].generate_state(1, dtype=np.uint32)[0])


def fixed_population_shift(config: dict) -> dict | None:
    """The pre-registered shift per population, or None if not declared."""
    values = config["real_data_simulation"].get("population_shift")
    return None if values is None else {str(k): float(v) for k, v in values.items()}


# --- phenotype on real genotypes ------------------------------------------------


def effect_allele_dosages(config: dict, raw: RawGenotypes, smap: pd.DataFrame):
    """Counts of each found sentinel's HbF-increasing allele, its beta, and metadata.

    Checks that the config and the locus table agree on the effect allele and
    that it is one of the data's alleles (catches strand or allele mismatches).
    """
    effects = {e["locus"]: e for e in config["planted_effects"]}
    found = smap[smap["id"].notna()].reset_index(drop=True)
    cols, betas, loci_meta = [], [], []
    for _, r in found.iterrows():
        gene = r["gene"]
        if gene not in effects:
            raise ValueError(f"no planted effect for locus {gene!r} in the config")
        if effects[gene]["effect_allele"] != r["effect_allele"]:
            raise ValueError(f"{gene}: config effect allele differs from the locus table")
        if r["effect_allele"] not in (r["ref"], r["alt"]):
            raise ValueError(
                f"{gene}: effect allele {r['effect_allele']} is neither REF {r['ref']} nor "
                f"ALT {r['alt']} in the data (strand or allele mismatch)"
            )
        if r["id"] not in raw.ids:
            raise ValueError(f"{gene}: {r['id']} missing from the exported genotypes")
        j = raw.ids.index(r["id"])
        g = raw.genotypes[:, j]
        flipped = raw.counted[j] != r["effect_allele"]
        cols.append(2.0 - g if flipped else g)
        betas.append(float(effects[gene]["beta_sd"]))
        loci_meta.append(
            {
                "gene": gene,
                "id": r["id"],
                "counted_allele_in_export": raw.counted[j],
                "effect_allele": r["effect_allele"],
                "flipped_to_effect_allele": bool(flipped),
            }
        )
    if not cols:
        raise ValueError("no sentinel genotypes found")
    geno = np.column_stack(cols)
    if np.isnan(geno).any():
        raise ValueError("missing genotypes at a sentinel; phase 3 should be fully called")
    return geno, np.array(betas), loci_meta


def simulate_phenotype_from_raw(
    config: dict,
    raw: RawGenotypes,
    smap: pd.DataFrame,
    psam: pd.DataFrame,
    seed_file,
) -> tuple[pd.DataFrame, dict]:
    """Simulate the HbF-like phenotype on exported sentinel genotypes (step e).

    Genotypes are converted to counts of each locus's HbF-increasing allele
    (the table's effect allele) and given the per-allele effect ``beta_sd``
    from ``planted_effects``. The non-genetic shift is the exact, pre-registered
    value per population (``real_data_simulation.population_shift``); the noise
    comes from the step e seed (``seed_file``, spawned child "phenotype").
    """
    sim_cfg = config["real_data_simulation"]
    shift = fixed_population_shift(config)
    if shift is None:
        raise ValueError("real_data_simulation.population_shift is not pre-registered")
    seeds = step_e_seeds(seed_file)
    geno, betas, loci_meta = effect_allele_dosages(config, raw, smap)

    pop_col = sim_cfg["population_column"]
    pops = psam.set_index("IID").reindex(raw.iid)[pop_col]
    if pops.isna().any():
        raise ValueError("some exported samples have no population label in the .psam")
    seed = seeds["seed"]
    ph = simulate.simulate_phenotype(
        geno,
        pops.to_numpy(),
        rng=np.random.default_rng(seeds["phenotype"]),
        beta=betas,
        pop_shift=shift,
    )
    betas = list(betas)
    out = pd.DataFrame({"IID": raw.iid, "HBF_SIM": ph.y})
    if "FID" in psam.columns:
        out.insert(0, "FID", raw.fid)
    for m, j in zip(loci_meta, range(geno.shape[1]), strict=True):
        m["effect_allele_freq"] = float(geno[:, j].mean() / 2)
        m["beta_sd"] = betas[j]
        m["r2_realized"] = float(ph.r2_realized[j])
    meta = {
        "phenotype": "simulated (hbfgwas.simulate.simulate_phenotype); not a measured trait",
        "n": int(len(out)),
        "step_e_seed": seed,
        "seed_child": "phenotype",
        "shift_mode": "fixed (pre-registered population_shift)",
        "shift_var_realised": float(ph.shift_var_fraction),
        "population_column": pop_col,
        "population_shift": {str(k): float(v) for k, v in ph.pop_shift.items()},
        "variance": {
            "genetic": float(ph.genetic.var()),
            "shift": float(ph.shift.var()),
            "noise": float(ph.noise.var()),
            "total": float(ph.y.var()),
        },
        "planted": loci_meta,
    }
    return out, meta


def prereg_predictions(
    config: dict,
    raw_sentinels: RawGenotypes,
    smap: pd.DataFrame,
    raw_null: RawGenotypes,
    psam: pd.DataFrame,
    eigenvec: pd.DataFrame,
    exclude: list[str] | None = None,
) -> dict:
    """Phase 2 step d: predictions from the real genotypes and PCs, before any phenotype.

    * conditional power per testable sentinel (PC-adjusted and unadjusted),
      with the planted effects, the population shift and alpha from the config;
    * predicted lambda_GC (unadjusted and PC-adjusted) over the random null
      variants, with its spread over shift draws (review item H).
    """
    from hbfgwas import predict

    if list(raw_sentinels.iid) != list(raw_null.iid):
        raise ValueError("sentinel and null exports have different samples")
    exclude = list(exclude or [])
    if exclude:
        keep = ~np.isin(raw_sentinels.iid, exclude)
        missing = sorted(set(exclude) - set(raw_sentinels.iid))
        if missing:
            raise ValueError(f"samples to exclude are not in the export: {missing}")
        raw_sentinels = RawGenotypes(
            raw_sentinels.fid[keep],
            raw_sentinels.iid[keep],
            raw_sentinels.ids,
            raw_sentinels.counted,
            raw_sentinels.genotypes[keep],
        )
        raw_null = RawGenotypes(
            raw_null.fid[keep],
            raw_null.iid[keep],
            raw_null.ids,
            raw_null.counted,
            raw_null.genotypes[keep],
        )
    sim_cfg = config["real_data_simulation"]
    pcfg = config["prediction"]
    h = float(sim_cfg["shift_var_fraction"])
    fixed = fixed_population_shift(config)
    alpha = config["association"]["genome_wide_alpha"]
    threshold = config["acceptance"]["locus_recovery"]["required_if_predicted_power_at_least"]
    pops = psam.set_index("IID").reindex(raw_sentinels.iid)[sim_cfg["population_column"]]
    if pops.isna().any():
        raise ValueError("some exported samples have no population label in the .psam")
    n_pcs = int(config["pca"]["n_pcs"])
    pcs = (
        eigenvec.set_index("IID")
        .reindex(raw_sentinels.iid)[[f"PC{i}" for i in range(1, n_pcs + 1)]]
        .to_numpy()
    )
    if np.isnan(pcs).any():
        raise ValueError("some exported samples have no PCs")
    seq = np.random.SeedSequence(int(config["seeds"]["real_prediction"]))
    power_seed, lambda_seed = seq.spawn(2)
    svec = None
    if fixed is not None:
        missing = sorted(set(pops.astype(str)) - set(fixed))
        if missing:
            raise ValueError(f"population_shift has no value for: {missing}")
        svec = np.array([fixed[str(p)] for p in pops])

    geno, betas, meta = effect_allele_dosages(config, raw_sentinels, smap)
    power_out = {}
    for side, cov in (("pc_adjusted", pcs), ("unadjusted", None)):
        res = predict.conditional_power(
            geno,
            betas,
            alpha,
            covariates=cov,
            pop=pops.to_numpy(),
            shift_var_fraction=h,
            n_shift_draws=int(pcfg["power_shift_draws"]),
            rng=np.random.default_rng(power_seed),
            shift=svec,
        )
        for m, r in zip(meta, res, strict=True):
            power_out.setdefault(m["gene"], {})[side] = r
    for m, j in zip(meta, range(geno.shape[1]), strict=True):
        p = power_out[m["gene"]]
        p["effect_allele_freq"] = float(geno[:, j].mean() / 2)
        p["required"] = bool(p["pc_adjusted"]["power"] >= threshold)

    null = raw_null.genotypes
    null = null[:, ~np.isnan(null).any(axis=0)]
    lam = {}
    for side, cov in (("unadjusted", None), ("pc_adjusted", pcs)):
        lam[side] = predict.predicted_lambda(
            null,
            pops.to_numpy(),
            h,
            covariates=cov,
            n_shift_draws=int(pcfg["lambda_shift_draws"]),
            rng=np.random.default_rng(lambda_seed),
            shift=svec,
        )
    return {
        "computed_before_any_phenotype": True,
        "n_samples": int(len(raw_sentinels.iid)),
        "excluded_samples": exclude,
        "alpha": alpha,
        "shift_mode": "fixed" if svec is not None else "random",
        "population_shift": fixed,
        "shift_var_realised": float(svec.var()) if svec is not None else None,
        "population_counts": {str(k): int(v) for k, v in pops.value_counts().sort_index().items()},
        "shift_var_fraction": h,
        "required_if_predicted_power_at_least": threshold,
        "conditional_power": power_out,
        "not_testable": [
            g for g, i in zip(smap["gene"], smap["id"], strict=True) if not _has_id(i)
        ],
        "predicted_lambda_gc": lam,
        "n_null_variants": int(null.shape[1]),
    }


def write_pheno(df: pd.DataFrame, path) -> None:
    df = df.rename(columns={df.columns[0]: "#" + df.columns[0]})
    df.to_csv(path, sep="\t", index=False, float_format="%.10g")


# --- het --------------------------------------------------------------------------


def het_summary(het: pd.DataFrame) -> pd.DataFrame:
    f = het["F"].astype(float)
    mu, sd = float(f.mean()), float(f.std(ddof=1))
    outliers = het.loc[(f - mu).abs() > 3 * sd, "IID"].tolist()
    return pd.DataFrame(
        [
            {
                "n_samples": len(f),
                "mean_F": mu,
                "sd_F": sd,
                "min_F": float(f.min()),
                "max_F": float(f.max()),
                "n_beyond_3sd": len(outliers),
                "ids_beyond_3sd": ",".join(outliers),
                "note": "reported only; no sample removed on heterozygosity",
            }
        ]
    )


# --- report -----------------------------------------------------------------------

GLM_COLS = ["#CHROM", "POS", "ID", "A1", "TEST", "OBS_CT", "BETA", "SE", "P"]


def read_glm(path) -> pd.DataFrame:
    head = pd.read_csv(path, sep="\t", nrows=0).columns
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=[c for c in GLM_COLS if c in head],
        dtype={"ID": str, "A1": str, "#CHROM": str},
    )
    df = df.rename(columns={"#CHROM": "CHROM"})
    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"]
    for c in ("BETA", "SE", "P", "OBS_CT", "POS"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.reset_index(drop=True)


def null_mask(glm: pd.DataFrame, planted: pd.DataFrame, window: int) -> np.ndarray:
    keep = np.ones(len(glm), dtype=bool)
    chrom = glm["CHROM"].map(_norm_chrom).to_numpy()
    pos = glm["POS"].to_numpy()
    for _, r in planted.iterrows():
        keep &= ~((chrom == str(r["chrom"])) & (np.abs(pos - int(r["pos"])) <= window))
    return keep


def _python_assoc(raw: RawGenotypes, pheno: pd.DataFrame, covar: np.ndarray | None):
    y = pheno.set_index("IID").reindex(raw.iid)["HBF_SIM"].to_numpy(dtype=float)
    if np.isnan(y).any():
        raise ValueError("phenotype missing for exported samples")
    res = assoc.linear_assoc(raw.genotypes, y, covar).to_frame(ids=raw.ids)
    res["A1"] = raw.counted
    return res


# --- step e: power redraws at the sentinels, and the sensitivity scan ------------


def power_redraws(
    config: dict,
    raw: RawGenotypes,
    smap: pd.DataFrame,
    psam: pd.DataFrame,
    eigenvec: pd.DataFrame,
    seed_file,
    n_redraws: int | None = None,
) -> dict:
    """Empirical power at the planted sentinels over phenotype redraws (step e).

    Each redraw is the step e phenotype model with new noise: the same planted
    effects on the same genotypes and the same fixed population shift
    (``real_data_simulation.population_shift``); only the noise is new, drawn
    from the step e seed's child "power_reps". Every redraw is tested at the
    sentinels only, with the PCs (the gated power calibration) and without them
    (reported). The number of redraws is
    ``acceptance.real.power_calibration.n_redraws`` unless given.

    The conditional power is recomputed from the same genotypes and PCs and
    reported next to the declared value; it is not a criterion.
    """
    from hbfgwas import experiment, predict

    if n_redraws is None:
        calib = (config["acceptance"].get("real") or {}).get("power_calibration") or {}
        n_redraws = calib.get("n_redraws")
    if not n_redraws or int(n_redraws) <= 0:
        raise ValueError(
            "the number of redraws is not given and acceptance.real.power_calibration "
            "has no n_redraws"
        )
    n_redraws = int(n_redraws)
    shift = fixed_population_shift(config)
    if shift is None:
        raise ValueError("real_data_simulation.population_shift is not pre-registered")
    seeds = step_e_seeds(seed_file)
    geno, betas, meta = effect_allele_dosages(config, raw, smap)
    pop_col = config["real_data_simulation"]["population_column"]
    pops = psam.set_index("IID").reindex(raw.iid)[pop_col]
    if pops.isna().any():
        raise ValueError("some exported samples have no population label in the .psam")
    n_pcs = int(config["pca"]["n_pcs"])
    pcs = (
        eigenvec.set_index("IID")
        .reindex(raw.iid)[[f"PC{i}" for i in range(1, n_pcs + 1)]]
        .to_numpy()
    )
    if np.isnan(pcs).any():
        raise ValueError("some exported samples have no PCs")
    alpha = config["association"]["genome_wide_alpha"]

    reps = simulate.simulate_phenotype(
        geno,
        pops.to_numpy(),
        rng=np.random.default_rng(seeds["power_reps"]),
        beta=betas,
        pop_shift=shift,
        n_reps=n_redraws,
    )
    adj = assoc.linear_assoc(geno, reps.y, pcs)
    unadj = assoc.linear_assoc(geno, reps.y)
    hits = (adj.p < alpha).sum(axis=1)
    hits_u = (unadj.p < alpha).sum(axis=1)

    svec = np.array([float(shift[str(p)]) for p in pops])
    recomputed = predict.conditional_power(geno, betas, alpha, covariates=pcs, shift=svec)
    loci = {}
    for j, m in enumerate(meta):
        k, ku = int(hits[j]), int(hits_u[j])
        loci[m["gene"]] = {
            "hits_pc_adjusted": k,
            "empirical_power": k / n_redraws,
            "empirical_power_ci95": experiment._binomial_ci95(k, n_redraws),
            "hits_unadjusted": ku,
            "empirical_power_unadjusted": ku / n_redraws,
            "empirical_power_unadjusted_ci95": experiment._binomial_ci95(ku, n_redraws),
            "predicted_power_recomputed": float(recomputed[j]["power"]),
            "effect_allele_freq": float(geno[:, j].mean() / 2),
            "beta_sd": float(betas[j]),
        }
    return {
        "phenotype": "simulated (hbfgwas.simulate.simulate_phenotype); not a measured trait",
        "model": "planted effects + fixed population shift + new noise in every redraw",
        "tested_at": "the planted sentinels only",
        "n_redraws": n_redraws,
        "n_samples": int(len(raw.iid)),
        "alpha": alpha,
        "step_e_seed": seeds["seed"],
        "seed_child": "power_reps",
        "shift_var_realised": float(svec.var()),
        "loci": loci,
    }


def sensitivity_section(config: dict, data_dir: Path, planted: pd.DataFrame) -> dict:
    """Results of the sensitivity scan (04c), or a 'not run' note. Never gated."""
    data_dir = Path(data_dir)
    files = {
        "unadjusted": data_dir / "glm_sens_nopc.HBF_SIM.glm.linear",
        "adjusted": data_dir / "glm_sens_pc.HBF_SIM.glm.linear",
    }
    excluded = list((config.get("sensitivity") or {}).get("exclude_samples") or [])
    if not all(f.exists() for f in files.values()):
        return {"status": "not run", "gated": False, "excluded_samples": excluded}
    t1a = config["association"]["type1_alpha"]
    window = int(config["association"]["null_exclusion_window_bp"])
    out: dict = {
        "status": "run",
        "gated": False,
        "excluded_samples": excluded,
        "pcs": "recomputed without the excluded samples (03d_sensitivity_pca.sh)",
        "phenotype": "the step e phenotype, not redrawn, without the excluded samples",
    }
    glm = {side: read_glm(f) for side, f in files.items()}
    for side, df in glm.items():
        ok = np.isfinite(df["P"].to_numpy())
        mask = null_mask(df, planted, window) & ok
        out[side] = {
            "lambda_gc": assoc.lambda_gc(df.loc[ok, "P"]),
            "type1_error": assoc.type1_error(df.loc[mask, "P"], t1a),
            "type1_alpha": t1a,
            "n_variants": int(ok.sum()),
            "n_null_variants": int(mask.sum()),
        }
    out["n_samples_analysed"] = (
        int(np.nanmax(glm["adjusted"]["OBS_CT"])) if len(glm["adjusted"]) else 0
    )
    out["loci"] = {}
    for _, r in planted.iterrows():
        rows = {side: df[df["ID"] == r["id"]] for side, df in glm.items()}
        out["loci"][r["gene"]] = {
            f"p_{side}": float(rows[side]["P"].iloc[0]) if len(rows[side]) else float("nan")
            for side in glm
        }
    return out


def build_report(config: dict, data_dir: Path, run_dir: Path, figures: bool = True) -> dict:
    data_dir, run_dir = Path(data_dir), Path(run_dir)
    alpha = config["association"]["genome_wide_alpha"]
    t1a = config["association"]["type1_alpha"]
    window = int(config["association"]["null_exclusion_window_bp"])
    n_pcs = int(config["pca"]["n_pcs"])
    tol = {
        k: config["concordance"][k]
        for k in ("beta_rtol", "beta_atol", "log10p_rtol", "log10p_atol")
    }

    glm = {
        "unadjusted": read_glm(data_dir / "glm_nopc.HBF_SIM.glm.linear"),
        "adjusted": read_glm(data_dir / "glm_pc.HBF_SIM.glm.linear"),
    }
    smap = pd.read_csv(data_dir / "planted_map.tsv", sep="\t", dtype={"id": str, "chrom": str})
    planted = smap[smap["id"].notna()]
    meta = json.loads((run_dir / "phenotype_meta.json").read_text())

    results: dict = {
        "experiment": "real",
        "description": "1000 Genomes AFR genotypes; HbF-like phenotype simulated by hbfgwas",
    }
    for side, df in glm.items():
        ok = np.isfinite(df["P"].to_numpy())
        mask = null_mask(df, planted, window) & ok
        results[side] = {
            "lambda_gc": assoc.lambda_gc(df.loc[ok, "P"]),
            "type1_error": assoc.type1_error(df.loc[mask, "P"], t1a),
            "type1_alpha": t1a,
            "n_variants": int(ok.sum()),
            "n_null_variants": int(mask.sum()),
            "null_exclusion_window_bp": window,
            "n_covariates": 0 if side == "unadjusted" else n_pcs,
        }
    n_obs = int(np.nanmax(glm["adjusted"]["OBS_CT"])) if len(glm["adjusted"]) else 0
    results["n_samples_analysed"] = n_obs

    real_criteria = config["acceptance"].get("real") or {}
    declared_power = real_criteria.get("predicted_power", {}) if real_criteria else {}
    redraws_path = run_dir / "power_redraws.json"
    redraws = json.loads(redraws_path.read_text()) if redraws_path.exists() else None
    if real_criteria.get("power_calibration") and redraws is None:
        raise FileNotFoundError(
            f"acceptance.real pre-registers a power calibration but {redraws_path} is "
            "missing (scripts/04b_power_redraws.sh)"
        )
    effects = {e["locus"]: e for e in config["planted_effects"]}
    h = float(config["real_data_simulation"]["shift_var_fraction"])
    planted_meta = {m["gene"]: m for m in meta["planted"]}
    loci = []
    for _, r in smap.iterrows():
        gene = r["gene"]
        entry = {"locus": gene, "rsid": r["rsid"], "id": r["id"], "found_by": r["found_by"]}
        if pd.isna(r["id"]) or gene not in planted_meta:
            entry["status"] = NOT_TESTABLE  # absent from the callset or failed QC; no proxy
            loci.append(entry)
            continue
        entry["status"] = "testable"
        pm = planted_meta[gene]
        p_eff = pm["effect_allele_freq"]
        if gene in declared_power:
            pred, source = float(declared_power[gene]), "pre-registered (acceptance.real)"
        else:
            pred = float(
                power.power(
                    n=n_obs,
                    alpha=alpha,
                    beta=effects[gene]["beta_sd"],
                    maf=min(p_eff, 1 - p_eff),
                    covariate_r2=h,
                )
            )
            source = "computed at report time from the real n and frequency (NOT pre-registered)"
        rows = {side: df[df["ID"] == r["id"]] for side, df in glm.items()}
        p_adj = float(rows["adjusted"]["P"].iloc[0]) if len(rows["adjusted"]) else float("nan")
        p_raw = float(rows["unadjusted"]["P"].iloc[0]) if len(rows["unadjusted"]) else float("nan")
        entry.update(
            effect_allele_freq=p_eff,
            beta_sd=effects[gene]["beta_sd"],
            predicted_power=pred,
            predicted_power_source=source,
            p_unadjusted=p_raw,
            p_adjusted=p_adj,
            recovered=bool(np.isfinite(p_adj) and p_adj < alpha),
        )
        if redraws is not None:
            if gene not in redraws["loci"]:
                raise ValueError(f"{gene} is testable but missing from {redraws_path}")
            rd = redraws["loci"][gene]
            entry.update(
                n_redraws=redraws["n_redraws"],
                empirical_power=rd["empirical_power"],
                empirical_power_ci95=rd["empirical_power_ci95"],
                empirical_power_unadjusted=rd["empirical_power_unadjusted"],
                predicted_power_recomputed=rd["predicted_power_recomputed"],
            )
        loci.append(entry)
    results["loci_all"] = loci
    results["loci"] = [x for x in loci if "p_adjusted" in x]
    results["not_testable"] = [x["locus"] for x in loci if x["status"] == NOT_TESTABLE]

    # Python reference vs plink2 on the concordance subset.
    raw = read_raw(data_dir / "concord.raw")
    pheno = read_table(data_dir / "pheno.tsv")
    pheno["HBF_SIM"] = pheno["HBF_SIM"].astype(float)
    eig = read_eigenvec(data_dir / "pca.eigenvec").set_index("IID").reindex(raw.iid)
    pcs = eig[[f"PC{i}" for i in range(1, n_pcs + 1)]].to_numpy()
    conc = {}
    for side, covar in (("unadjusted", None), ("adjusted", pcs)):
        py = _python_assoc(raw, pheno, covar)
        res = assoc.compare_to_plink2(py, glm[side], **tol)
        conc[side] = res.to_dict()
        py.to_csv(run_dir / f"concordance_python_{side}.tsv", sep="\t", index=False)
    conc["passed"] = all(conc[s]["passed"] for s in ("unadjusted", "adjusted"))
    results["concordance"] = conc
    if redraws is not None:
        results["power_redraws"] = {
            k: redraws[k] for k in ("n_redraws", "seed_child", "step_e_seed", "model")
        }
    results["sensitivity"] = sensitivity_section(config, data_dir, planted)

    if figures:
        _figures(config, data_dir, run_dir, glm, results)
    return results


def _figures(config, data_dir, run_dir, glm, results) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharex=True, sharey=True)
    for ax, side, title, color in (
        (axes[0], "unadjusted", "No ancestry adjustment", "#c0392b"),
        (axes[1], "adjusted", f"Adjusted for {config['pca']['n_pcs']} PCs", "#2471a3"),
    ):
        assoc.qq_plot(ax, glm[side]["P"].to_numpy(), title=title, color=color, max_points=60_000)
    fig.suptitle("1000 Genomes AFR genotypes; phenotype simulated", fontsize=10)
    fig.tight_layout()
    fig.savefig(run_dir / "qq_real.png", dpi=150)
    plt.close(fig)


def pc_plot(config: dict, eigenvec_path, psam_path, out_path) -> dict:
    """PC1 vs PC2 coloured by population (the .psam Population column, not SuperPop).

    Used by scripts/03b_pc_plot.py at Phase 2 step c, right after the PCA.
    Returns the sample count per population.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    eig = read_eigenvec(eigenvec_path)
    psam = read_psam(psam_path)
    pop_col = config["real_data_simulation"]["population_column"]
    df = eig.merge(psam[["IID", pop_col]], on="IID", how="left")
    if df[pop_col].isna().any():
        raise ValueError(f"some PCA samples have no {pop_col} label in the .psam")
    counts = {}
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for i, (pop, grp) in enumerate(sorted(df.groupby(pop_col), key=lambda t: t[0])):
        ax.scatter(grp["PC1"], grp["PC2"], s=10, label=f"{pop} (n={len(grp)})", color=f"C{i}")
        counts[str(pop)] = int(len(grp))
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(fontsize=8, frameon=False)
    ax.set_title(f"1000 Genomes AFR: PC1 vs PC2 by {pop_col}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return counts


# --- CLI --------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m hbfgwas.pipeline", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("config-get")
    p.add_argument("--config", default=None)
    p.add_argument("key")

    p = sub.add_parser(
        "loci-ready", help="exit 1 unless every row of a locus TSV is analysis_ready"
    )
    p.add_argument("--loci", required=True)

    p = sub.add_parser("sentinels")
    p.add_argument("--pvar", required=True)
    p.add_argument("--loci", required=True)
    p.add_argument("--out-ids", required=True)
    p.add_argument("--out-map", required=True)

    p = sub.add_parser("sentinel-report")
    p.add_argument("--map", required=True)
    p.add_argument("--afreq", default=None)
    p.add_argument("--step", action="append", default=[], help="name=snplist (in QC order)")
    p.add_argument("--out", required=True)

    p = sub.add_parser("sample-ids")
    p.add_argument("--pvar", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--include", default=None)
    p.add_argument("--exclude-around", default=None, help="sentinel map TSV (chrom, pos)")
    p.add_argument("--window-bp", type=int, default=0)
    p.add_argument("--out", required=True)

    p = sub.add_parser("predictions")
    p.add_argument("--config", default=None)
    p.add_argument("--raw-sentinels", required=True)
    p.add_argument("--map", required=True)
    p.add_argument("--raw-null", required=True)
    p.add_argument("--psam", required=True)
    p.add_argument("--eigenvec", required=True)
    p.add_argument("--exclude", default="", help="comma-separated sample IDs to leave out")
    p.add_argument("--out", required=True)

    p = sub.add_parser("seed-child", help="print an integer seed spawned from the step e seed")
    p.add_argument("--seed-file", required=True)
    p.add_argument("--name", required=True, choices=list(STEP_E_SPAWN))

    p = sub.add_parser("simulate-phenotype")
    p.add_argument("--config", default=None)
    p.add_argument("--seed-file", required=True)
    p.add_argument("--raw", required=True)
    p.add_argument("--map", required=True)
    p.add_argument("--psam", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--meta", required=True)

    p = sub.add_parser("power-redraws", help="step e: empirical power at the sentinels")
    p.add_argument("--config", default=None)
    p.add_argument("--seed-file", required=True)
    p.add_argument("--raw", required=True)
    p.add_argument("--map", required=True)
    p.add_argument("--psam", required=True)
    p.add_argument("--eigenvec", required=True)
    p.add_argument("--n-redraws", type=int, default=None)
    p.add_argument("--out", required=True)

    p = sub.add_parser("het-summary")
    p.add_argument("--het", required=True)
    p.add_argument("--out", required=True)

    p = sub.add_parser("report")
    p.add_argument("--config", default=None)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--run-dir", required=True)
    p.add_argument("--no-figures", action="store_true")

    a = parser.parse_args(argv)
    if a.cmd == "config-get":
        print(format_scalar(config_get(load_config(a.config), a.key)))
    elif a.cmd == "loci-ready":
        loci_df = pd.read_csv(a.loci, sep="\t", dtype={"rsid": str})
        if "analysis_ready" not in loci_df.columns:
            print(f"{a.loci} has no analysis_ready column")
            return 1
        ready = loci_df["analysis_ready"].astype(str).str.lower().eq("true")
        for gene, ok in zip(loci_df["gene"], ready, strict=True):
            print(f"{gene}: {'analysis-ready' if ok else 'NOT analysis-ready'}")
        return 0 if ready.all() else 1
    elif a.cmd == "sentinels":
        loci = pd.read_csv(a.loci, sep="\t", dtype={"chromosome": str, "rsid": str})
        smap = find_sentinels(a.pvar, loci)
        smap.to_csv(a.out_map, sep="\t", index=False)
        ids = smap["id"].dropna().tolist()
        Path(a.out_ids).write_text("".join(f"{i}\n" for i in ids))
        for _, r in smap.iterrows():
            print(
                f"{r['gene']}: {r['rsid']} -> {r['found_by']}"
                + (
                    f" ({r['id']} {r['chrom']}:{r['pos']} {r['ref']}/{r['alt']})"
                    if _has_id(r["id"])
                    else ""
                )
            )
    elif a.cmd == "sentinel-report":
        smap = pd.read_csv(a.map, sep="\t", dtype={"id": str, "chrom": str})
        afreq = read_table(a.afreq) if a.afreq else None
        steps = {}
        for s in a.step:
            name, path = s.split("=", 1)
            steps[name] = read_ids(path)
        rep = sentinel_report(smap, afreq, steps)
        rep.to_csv(a.out, sep="\t", index=False)
        cols = ["gene", "rsid", "found_by", "effect_allele_freq", "removed_by", "status"]
        print(rep[cols].to_string(index=False))
    elif a.cmd == "sample-ids":
        include = read_ids(a.include) if a.include else set()
        windows = []
        if a.exclude_around:
            m = pd.read_csv(a.exclude_around, sep="\t", dtype={"chrom": str, "id": str})
            windows = [(r["chrom"], int(r["pos"])) for _, r in m.iterrows() if _has_id(r["id"])]
        ids = sample_ids(a.pvar, a.n, a.seed, include, windows, a.window_bp)
        Path(a.out).write_text("".join(f"{i}\n" for i in ids))
        print(f"{len(ids)} variant IDs -> {a.out}")
    elif a.cmd == "seed-child":
        try:
            print(seed_child_int(a.seed_file, a.name))
        except ValueError as exc:
            print(exc, file=sys.stderr)
            return 1
    elif a.cmd == "simulate-phenotype":
        config = load_config(a.config)
        smap = pd.read_csv(a.map, sep="\t", dtype={"id": str, "chrom": str})
        pheno, meta = simulate_phenotype_from_raw(
            config, read_raw(a.raw), smap, read_psam(a.psam), a.seed_file
        )
        write_pheno(pheno, a.out)
        Path(a.meta).write_text(json.dumps(meta, indent=2) + "\n")
        print(f"simulated phenotype for {meta['n']} samples -> {a.out}")
    elif a.cmd == "predictions":
        config = load_config(a.config)
        smap = pd.read_csv(a.map, sep="\t", dtype={"id": str, "chrom": str})
        res = prereg_predictions(
            config,
            read_raw(a.raw_sentinels),
            smap,
            read_raw(a.raw_null),
            read_psam(a.psam),
            read_eigenvec(a.eigenvec),
            exclude=[x for x in a.exclude.split(",") if x],
        )
        Path(a.out).write_text(json.dumps(res, indent=2, default=float) + "\n")
        lam = res["predicted_lambda_gc"]
        for side in ("unadjusted", "pc_adjusted"):
            x = lam[side]
            print(
                f"predicted lambda_GC {side}: {x['median']:.4f} "
                f"[{x['q05']:.4f}, {x['q95']:.4f}] over {x['n_variants']} null variants"
            )
        for gene, p in res["conditional_power"].items():
            print(
                f"{gene}: conditional power {p['pc_adjusted']['power']:.4f} (PC-adjusted), "
                f"{p['unadjusted']['power']:.4f} (unadjusted); "
                f"{'REQUIRED' if p['required'] else 'reported, not gated'}"
            )
        for gene in res["not_testable"]:
            print(f"{gene}: {NOT_TESTABLE}")
        print(f"wrote {a.out}")
    elif a.cmd == "power-redraws":
        config = load_config(a.config)
        smap = pd.read_csv(a.map, sep="\t", dtype={"id": str, "chrom": str})
        res = power_redraws(
            config,
            read_raw(a.raw),
            smap,
            read_psam(a.psam),
            read_eigenvec(a.eigenvec),
            a.seed_file,
            n_redraws=a.n_redraws,
        )
        Path(a.out).write_text(json.dumps(res, indent=2, default=float) + "\n")
        for gene, r in res["loci"].items():
            lo, hi = r["empirical_power_ci95"]
            print(
                f"{gene}: empirical power {r['empirical_power']:.4f} [{lo:.4f}, {hi:.4f}] over "
                f"{res['n_redraws']} redraws (PC-adjusted); recomputed prediction "
                f"{r['predicted_power_recomputed']:.4f}"
            )
        print(f"wrote {a.out}")
    elif a.cmd == "het-summary":
        s = het_summary(read_table(a.het))
        s.to_csv(a.out, sep="\t", index=False)
        print(s.to_string(index=False))
    elif a.cmd == "report":
        config = load_config(a.config)
        res = build_report(config, Path(a.data_dir), Path(a.run_dir), figures=not a.no_figures)
        out = Path(a.run_dir) / "results.json"
        out.write_text(json.dumps(res, indent=2, default=float) + "\n")
        print(json.dumps({k: res[k] for k in ("unadjusted", "adjusted")}, indent=2))
        for loc in res["loci_all"]:
            print(loc)
        print("concordance:", "PASS" if res["concordance"]["passed"] else "FAIL")
        print(f"wrote {out}")
        return 0 if res["concordance"]["passed"] else 4
    return 0


if __name__ == "__main__":
    sys.exit(main())
