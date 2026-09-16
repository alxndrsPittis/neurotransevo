#!/usr/bin/env python3
"""Boltz-2 binding-affinity dot plots: Fig. 4B (right) and Fig. S14.

Dot size = binding probability (affinity_probability_binary); colour =
-affinity_pred_value (higher = stronger predicted binding).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, results_path  # noqa: E402

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.6, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

LIGANDS = ["ach", "dop", "ser", "try", "pea", "his"]
LIGAND_LABELS = {"ach": "ACh", "dop": "Dopamine", "ser": "Serotonin",
                 "try": "Tryptamine", "pea": "2-PEA", "his": "Histamine"}
CMAP, VMIN, VMAX = "RdBu_r", -1.5, 1.5
SIZE_SCALE = 300

# Row order follows the Fig. 4B phylogeny (top to bottom).
FIG4B_ORDER = [
    "HRH2_MOUSE", "HRH2_HUMAN", "HRH1_HUMAN", "HRH1_MOUSE", "HRH4_HUMAN", "HRH4_MOUSE", "HRH3_HUMAN", "HRH3_MOUSE",
    "Q9VHW1_DROME",
    "ACM1_DROME", "ACM4_HUMAN", "ACM4_MOUSE", "ACM2_HUMAN", "ACM2_MOUSE", "ACM5_MOUSE", "ACM5_HUMAN",
    "ACM3_HUMAN", "ACM3_MOUSE", "ACM1_MOUSE", "ACM1_HUMAN",
    "XP_012561693.1", "t17782aep", "t21045aep", "XP_012561685.1", "XP_012561731.1",
    "t18253aep", "XP_012562020.1", "t1966aep",
    "XP_012554108.1", "XP_012561998.1", "XP_012561929.1", "t19057aep",
    "XP_012560064.1", "XP_012563056.1", "t25534aep", "XP_012562068.1", "XP_012562571.1",
    "XP_002170550.2",
    "MneR", "TriR",
]
CONTROLS = ["5HT2A_HUMAN", "5HT2B_HUMAN", "5HT2C_HUMAN", "DRD2_HUMAN", "DRD3_HUMAN", "DRD4_HUMAN"]
POCKET_CONSERVED = {"XP_012561685.1", "XP_012561731.1", "t21045aep"}
ALIAS = {"XP_012561685.1": "ACh-rec1 (XP_012561685.1)", "XP_012561731.1": "ACh-rec2 (XP_012561731.1)"}


def dotplot(prob, aff, rows, figsize, stem, title=None, highlight=()):
    fig, ax = plt.subplots(figsize=figsize)
    X, Y = np.meshgrid(range(len(LIGANDS)), range(len(rows)))
    sc = ax.scatter(X.ravel(), Y.ravel(), s=(prob.loc[rows, LIGANDS].values * SIZE_SCALE).clip(min=4).ravel(),
                    c=(-aff.loc[rows, LIGANDS].values).ravel(), cmap=CMAP, vmin=VMIN, vmax=VMAX,
                    edgecolors="k", linewidth=0.25)
    ax.set_xticks(range(len(LIGANDS)), [LIGAND_LABELS[l] for l in LIGANDS], rotation=45, ha="right")
    ax.set_yticks(range(len(rows)), [ALIAS.get(r, r) for r in rows], fontsize=7.5)
    for lbl, r in zip(ax.get_yticklabels(), rows):
        if r in highlight:
            lbl.set_fontweight("bold"); lbl.set_color("#b00")
    ax.set_xlim(-0.5, len(LIGANDS) - 0.5); ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.tick_params(length=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, shrink=0.25, pad=0.02, aspect=12)
    cbar.set_label("-affinity_pred_value\n(higher = stronger binder)", fontsize=7)
    cbar.ax.tick_params(labelsize=7, length=2); cbar.outline.set_linewidth(0.5)
    handles = [ax.scatter([], [], s=p * SIZE_SCALE, c="lightgray", edgecolor="k", linewidth=0.25, label=str(p))
               for p in (0.25, 0.5, 0.75, 1.0)]
    ax.legend(handles=handles, title="P(binder)", loc="lower left", bbox_to_anchor=(1.02, 0.0),
              fontsize=7, title_fontsize=8)
    if title:
        ax.set_title(title, fontsize=9)
    for ext in ("pdf", "svg"):
        fig.savefig(results_path("figures", f"{stem}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("figure:", results_path("figures", stem))


def main():
    df = pd.read_csv(CFG["paths"]["boltz_table"])
    prob = df.pivot(index="protein", columns="ligand", values="affinity_probability_binary")
    aff = df.pivot(index="protein", columns="ligand", values="affinity_pred_value")
    missing = set(FIG4B_ORDER + CONTROLS) - set(prob.index)
    if missing:
        sys.exit(f"receptors missing from {CFG['paths']['boltz_table']}: {sorted(missing)}")
    dotplot(prob, aff, FIG4B_ORDER, (4.8, 10.0), "Fig4B_boltz_dotplot", highlight=POCKET_CONSERVED)
    dotplot(prob, aff, CONTROLS, (4.6, 3.0), "FigS14_boltz_controls", title="Control receptors (specificity check)")
    summary = prob.loc[FIG4B_ORDER + CONTROLS, LIGANDS].round(3)
    summary.insert(0, "top_ligand", prob.loc[summary.index, LIGANDS].idxmax(axis=1))
    summary.to_csv(results_path("tables", "boltz_probabilities.tsv"), sep="\t")


if __name__ == "__main__":
    main()
