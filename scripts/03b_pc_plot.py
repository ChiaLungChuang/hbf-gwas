#!/usr/bin/env python3
"""03b: PC1 vs PC2 coloured by the seven 1000 Genomes AFR populations.

Runs at Phase 2 step c, right after 03_pca.sh. Colours by the .psam
Population column (YRI, LWK, GWD, MSL, ESN, ASW, ACB), not by SuperPop.

In:  DATA_DIR/pca.eigenvec, DATA_DIR/afr_qc.psam
Out: RUN_DIR/pc1_pc2_by_population.png, RUN_DIR/pc_plot_counts.json

Uses the same environment variables as the shell scripts (CONFIG, DATA_DIR,
RUN_DIR).
"""

import json
import os
import sys
import time
from pathlib import Path

from hbfgwas.config import load_config
from hbfgwas.pipeline import pc_plot

if __name__ == "__main__":
    config = load_config(os.environ.get("CONFIG", "config/analysis.yaml"))
    data_dir = Path(os.environ.get("DATA_DIR", "data/1kg"))
    run_dir = Path(os.environ.get("RUN_DIR", "runs/" + time.strftime("%Y-%m-%d", time.gmtime())))
    run_dir.mkdir(parents=True, exist_ok=True)
    out = run_dir / "pc1_pc2_by_population.png"
    counts = pc_plot(config, data_dir / "pca.eigenvec", data_dir / "afr_qc.psam", out)
    (run_dir / "pc_plot_counts.json").write_text(json.dumps(counts, indent=2) + "\n")
    print(f"wrote {out}")
    for pop, n in counts.items():
        print(f"  {pop}: {n}")
    sys.exit(0)
