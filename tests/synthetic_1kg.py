"""Build a tiny synthetic stand-in for the 1000 Genomes phase 3 plink2 fileset.

Used by tests/test_scripts.py to exercise scripts/00-05 end to end in CI. It is
NOT real data. Layout mirrors what scripts/00_download_1kg.sh expects to
download: all_phase3.pgen.zst, all_phase3.pvar.zst, all_phase3.psam
(#IID PAT MAT SEX SuperPop Population) and deg2_phase3.king.cutoff.out.id.

Contents (each special record has a known fate):
* 10 populations x 20 samples: the 7 AFR populations plus CEU, GBR, CHB.
* related-sample list: 4 AFR + 1 EUR sample -> 136 AFR samples remain.
* ~3,500 Balding-Nichols SNPs on chromosomes 1, 2, 3, 6, 11 (Fst 0.02).
* the three sentinels at their GRCh37 positions; rs7482144 is stored under
  the ID rs41322148 so the lookup must fall back to the position.
* removed by 01: 5 chrX SNPs, 5 indels, 3 multiallelic records, 2 split
  multiallelic sites (4 records), 2 records sharing the ID rsDUP.
* 10 records with ID "." -> renamed chrom:pos:ref:alt by 01 and kept.
* removed by 02: 10 monomorphic (rsMONO*) and 5 singleton (rsRARE*) records by
  --maf; 5 all-heterozygous records (rsHWE*) by --hwe.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np

from hbfgwas import simulate

AFR = ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"]
OTHER = {"CEU": "EUR", "GBR": "EUR", "CHB": "EAS"}
PER_POP = 20
REL = ["SYN_YRI_00", "SYN_LWK_00", "SYN_ASW_00", "SYN_ACB_00", "SYN_CEU_00"]
CHROM_LEN = {
    "1": 249_000_000,
    "2": 243_000_000,
    "3": 198_000_000,
    "6": 171_000_000,
    "11": 135_000_000,
}
SENTINELS = [  # chrom, pos, data ID, REF, ALT, ancestral ALT frequency
    ("2", 60718043, "rs1427407", "T", "G", 0.75),
    ("6", 135419018, "rs9399137", "T", "C", 0.12),
    ("11", 5276169, "rs41322148", "G", "A", 0.15),
]
EXPECTED_AFR_AFTER_REL = 7 * PER_POP - 4
REMOVED_BY_01 = ["rsX1", "rsINDEL1", "rsMULTI1", "rsSPLIT1a", "rsSPLIT1b", "rsDUP"]
REMOVED_BY_MAF = [f"rsMONO{i}" for i in range(10)] + [f"rsRARE{i}" for i in range(5)]
REMOVED_BY_HWE = [f"rsHWE{i}" for i in range(5)]
KEPT = ["rs1427407", "rs9399137", "rs41322148"]


def _samples():
    pops = AFR + list(OTHER)
    rows = []
    for pop in pops:
        sup = "AFR" if pop in AFR else OTHER[pop]
        for i in range(PER_POP):
            rows.append((f"SYN_{pop}_{i:02d}", pop, sup, 1 + (i % 2)))
    return rows


def _gt(g: int) -> str:
    return ("0/0", "0/1", "1/1")[g]


def build(out_dir: Path, plink2: str = "plink2", zstd: str = "zstd", seed: int = 5) -> dict:
    out_dir = Path(out_dir)
    work = out_dir / "work"
    work.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    samples = _samples()
    names = [s[0] for s in samples]
    n = len(samples)
    afr_idx = [i for i, s in enumerate(samples) if s[2] == "AFR"]
    sizes = [PER_POP] * (len(AFR) + len(OTHER))
    bases = np.array(list("ACGT"))

    records = []  # (chrom, pos, id, ref, alt, [gt strings])

    def add(chrom, pos, vid, ref, alt, gts):
        records.append((chrom, int(pos), vid, ref, alt, gts))

    taken = {(c, p) for c, p, *_ in SENTINELS}

    def free_pos(chrom):
        while True:
            p = int(rng.integers(1_000_000, CHROM_LEN[chrom] - 1_000_000))
            if (chrom, p) not in taken:
                taken.add((chrom, p))
                return p

    # regular structured SNPs
    for chrom in CHROM_LEN:
        m = 700
        geno = simulate.simulate_genotypes(sizes, fst=0.02, rng=rng, n_variants=m).genotypes
        for j in range(m):
            ref, alt = rng.choice(bases, 2, replace=False)
            add(chrom, free_pos(chrom), f"rsSYN{chrom}_{j}", ref, alt, [_gt(g) for g in geno[:, j]])
    # sentinels
    for chrom, pos, vid, ref, alt, p in SENTINELS:
        geno = simulate.simulate_genotypes(
            sizes, fst=0.02, rng=rng, ancestral_freq=np.array([p])
        ).genotypes[:, 0]
        add(chrom, pos, vid, ref, alt, [_gt(g) for g in geno])
    # missing IDs (kept, renamed by 01)
    for _ in range(10):
        g = rng.binomial(2, 0.4, size=n)
        add("1", free_pos("1"), ".", "C", "T", [_gt(x) for x in g])
    # removed by 01
    for i in range(5):
        add(
            "X",
            1_000_000 + i * 50_000,
            f"rsX{i + 1}",
            "A",
            "G",
            [("0/0", "1/1")[int(x)] for x in rng.integers(0, 2, size=n)],
        )
    for i in range(5):
        g = rng.binomial(2, 0.3, size=n)
        add("3", free_pos("3"), f"rsINDEL{i + 1}", "A", "AT", [_gt(x) for x in g])
    for i in range(3):
        gts = [("0/0", "0/1", "0/2", "1/2")[int(x)] for x in rng.integers(0, 4, size=n)]
        add("3", free_pos("3"), f"rsMULTI{i + 1}", "A", "C,G", gts)
    for i in range(2):
        pos = free_pos("1")
        for suffix, alt in (("a", "C"), ("b", "G")):
            g = rng.binomial(2, 0.3, size=n)
            add("1", pos, f"rsSPLIT{i + 1}{suffix}", "A", alt, [_gt(x) for x in g])
    for _ in range(2):
        g = rng.binomial(2, 0.3, size=n)
        add("1", free_pos("1"), "rsDUP", "G", "T", [_gt(x) for x in g])
    # removed by 02
    for i in range(10):
        add("3", free_pos("3"), f"rsMONO{i}", "C", "G", ["0/0"] * n)
    het_sample = names.index("SYN_YRI_05")
    for i in range(5):
        gts = ["0/0"] * n
        gts[het_sample] = "0/1"
        add("3", free_pos("3"), f"rsRARE{i}", "G", "A", gts)
    for i in range(5):
        add("3", free_pos("3"), f"rsHWE{i}", "T", "C", ["0/1"] * n)

    order = {c: k for k, c in enumerate(["1", "2", "3", "6", "11", "X"])}
    records.sort(key=lambda r: (order[r[0]], r[1]))

    vcf = work / "synth.vcf"
    with open(vcf, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n"
        )
        for chrom, pos, vid, ref, alt, gts in records:
            fh.write(
                f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n"
            )

    prefix = work / "all_phase3"
    subprocess.run(
        [plink2, "--vcf", str(vcf), "--make-pgen", "--out", str(prefix)],
        check=True,
        capture_output=True,
        text=True,
    )
    with open(f"{prefix}.psam", "w") as fh:
        fh.write("#IID\tPAT\tMAT\tSEX\tSuperPop\tPopulation\n")
        for iid, pop, sup, sex in samples:
            fh.write(f"{iid}\t0\t0\t{sex}\t{sup}\t{pop}\n")

    serve = out_dir / "serve"
    serve.mkdir(exist_ok=True)
    for ext in ("pgen", "pvar"):
        subprocess.run(
            [zstd, "-q", "-f", f"{prefix}.{ext}", "-o", str(serve / f"all_phase3.{ext}.zst")],
            check=True,
        )
    shutil.copy(f"{prefix}.psam", serve / "all_phase3.psam")
    (serve / "deg2_phase3.king.cutoff.out.id").write_text("#IID\n" + "".join(f"{s}\n" for s in REL))
    return {
        "serve": serve,
        "n_samples": n,
        "n_afr": len(afr_idx),
        "n_records": len(records),
        "pgen_plain": Path(f"{prefix}.pgen"),
    }
