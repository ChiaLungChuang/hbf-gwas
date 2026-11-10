#!/usr/bin/env python3
"""05: report on the real-genotype association (phenotype simulated).

lambda_GC without and with the PCs, empirical type I error on variants away
from the planted sentinels, PC1 vs PC2 coloured by the seven AFR populations
(the Population column, not SuperPop), locus recovery against the predicted
power, and the Python-vs-plink2 concordance check (hbfgwas.assoc is the
tested reference). Writes RUN_DIR/results.json, qq_real.png,
pc1_pc2_by_population.png and concordance_python_*.tsv.

Uses the same environment variables as the shell scripts (CONFIG, DATA_DIR,
RUN_DIR). The acceptance gate is a separate step:
    python -m hbfgwas.gate --experiment real --results "$RUN_DIR/results.json"
Exit code 4 if the concordance check fails.
"""

import os
import sys
import time

from hbfgwas.pipeline import main

if __name__ == "__main__":
    config = os.environ.get("CONFIG", "config/analysis.yaml")
    data_dir = os.environ.get("DATA_DIR", "data/1kg")
    run_dir = os.environ.get("RUN_DIR", "runs/" + time.strftime("%Y-%m-%d", time.gmtime()))
    os.makedirs(run_dir, exist_ok=True)
    sys.exit(main(["report", "--config", config, "--data-dir", data_dir, "--run-dir", run_dir]))
