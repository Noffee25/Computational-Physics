"""Shared matplotlib style for all Project 1 figures.

Usage from any script (run with /usr/local/bin/python3 -I):

    import sys
    sys.path.insert(0, "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude")
    import plotstyle as ps
    ps.apply()
    fig, ax = plt.subplots(figsize=ps.WIDE)
    ax.plot(x, y, color=ps.COLORS["RANDU"], label="RANDU")
    ps.save(fig, "fam_name")   # writes report/figures/fam_name.png

Rules (keep every figure consistent):
  * One y-axis per panel; never twin axes. Use small multiples instead.
  * A single series uses ps.SINGLE (the blue of the earlier histograms) and
    needs no legend. Two or more series always get a legend.
  * Each generator always keeps the same color (ps.COLORS), in every figure.
  * Short, plain titles ("RANDU triplets"), axis labels with units/meaning.
  * Scatter plots with more than three generators: facet into small
    multiples rather than overlaying.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/Users/noah/Downloads/ComputationalPhysicsPHYS25000/Project1/P1_Claude"
FIGDIR = os.path.join(ROOT, "report", "figures")

SINGLE = "#3b6ea5"

# Fixed generator -> color mapping (validated categorical order, slots 1-8,
# LOGISTIC folded to neutral gray).
COLORS = {
    "RANDU":    "#2a78d6",
    "RANDM":    "#eb6834",
    "RAND":     "#1baf7a",
    "MINSTD":   "#eda100",
    "R250":     "#e87ba4",
    "RANNUM":   "#008300",
    "MT19937":  "#6250d6",
    "MIDSQ":    "#e34948",
    "LOGISTIC": "#898781",
}
ORDER = ["RANDU", "RANDM", "RAND", "MINSTD", "R250", "RANNUM", "MT19937", "MIDSQ", "LOGISTIC"]

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
EXPECT = "#52514e"   # reference curves (theory, expected value): dashed, this ink

WIDE = (10, 4.2)
SQUARE = (5.2, 5.0)
TALL = (10, 7.5)


def apply():
    plt.rcParams.update({
        "figure.dpi": 100,
        "savefig.dpi": 200,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.edgecolor": "#c3c2b7",
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "lines.linewidth": 1.6,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "figure.autolayout": False,
    })


def save(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, name + ".png")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path
