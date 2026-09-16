#!/usr/bin/env python3
"""Summary tables: coregulon sizes per species and figure inventory (results/tables/)."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from config import CFG, results_path  # noqa: E402

rows = []
for f in sorted((CFG["paths"]["results"] / "coregulons").glob("*.members.tsv")):
    d = pd.read_csv(f, sep="\t")
    rows.append({"species": f.name.split(".")[0], "genes": len(d), "active_seeds": int(d["is_seed"].sum()),
                 "named_genes": ", ".join(sorted(set(d["name"].dropna()) - {"-"}))})
pd.DataFrame(rows).to_csv(results_path("tables", "coregulon_summary.tsv"), sep="\t", index=False)
figs = sorted(p.stem for p in (CFG["paths"]["results"] / "figures").glob("*.pdf"))
pd.DataFrame({"figure_panel": figs}).to_csv(results_path("tables", "figure_panels.tsv"), sep="\t", index=False)
print(pd.DataFrame(rows)[["species", "genes", "active_seeds"]].to_string(index=False))
