from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, ScalarFormatter


# ============================================================
# FILE PATHS
# ============================================================

BASE = Path(__file__).resolve().parent

INPUT_FILE = BASE / "figure1_downstream_peakT.csv"
OUTPUT_PNG = BASE / "Figure_1.png"
OUTPUT_PDF = BASE / "Figure_1.pdf"


# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT_FILE).sort_values("distance_m")

required_columns = {
    "distance_m",
    "Tmax_Uref5",
    "Tmax_Uref10",
    "Tmax_Uref15",
}

missing_columns = required_columns.difference(df.columns)

if missing_columns:
    raise ValueError(
        f"Missing columns: {sorted(missing_columns)}"
    )


# ============================================================
# PLOT STYLE
# ============================================================

plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 300,

    "font.family": "DejaVu Sans",
    "font.size": 9,

    "axes.labelsize": 10,
    "axes.linewidth": 0.8,

    "xtick.labelsize": 9,
    "ytick.labelsize": 9,

    "legend.fontsize": 8.5,

    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


series = [
    {
        "column": "Tmax_Uref5",
        "label": r"$U_{\mathrm{ref}} = 5$ m s$^{-1}$",
        "color": "#0072B2",
        "marker": "o",
        "linestyle": "-",
    },
    {
        "column": "Tmax_Uref10",
        "label": r"$U_{\mathrm{ref}} = 10$ m s$^{-1}$",
        "color": "#D55E00",
        "marker": "s",
        "linestyle": "--",
    },
    {
        "column": "Tmax_Uref15",
        "label": r"$U_{\mathrm{ref}} = 15$ m s$^{-1}$",
        "color": "#009E73",
        "marker": "^",
        "linestyle": "-.",
    },
]


# ============================================================
# CREATE FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(6.4, 4.0),
    constrained_layout=True,
)


# ============================================================
# PLOT
# ============================================================

for item in series:

    ax.plot(
        df["distance_m"],
        df[item["column"]],
        label=item["label"],
        color=item["color"],
        marker=item["marker"],
        linestyle=item["linestyle"],
        linewidth=1.6,
        markersize=5.5,
        markerfacecolor="white",
        markeredgecolor=item["color"],
        markeredgewidth=1.1,
        zorder=3,
    )


# ============================================================
# LABELS
# ============================================================

ax.set_xlabel("Downstream distance (m)")

ax.set_ylabel(
    r"Peak normalized contaminant concentration, $T_{\max}$ (-)"
)


# ============================================================
# AXIS FORMATTING
# ============================================================

distances = df["distance_m"].to_numpy()

ax.set_xticks(distances)

x_min = distances.min()
x_max = distances.max()
x_padding = 0.04 * (x_max - x_min)

ax.set_xlim(
    x_min - x_padding,
    x_max + x_padding,
)

ax.set_ylim(bottom=0)

formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-3, 3))

ax.yaxis.set_major_formatter(formatter)

ax.xaxis.set_minor_locator(
    AutoMinorLocator(2)
)

ax.yaxis.set_minor_locator(
    AutoMinorLocator(2)
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    length=4.5,
    width=0.8,
)

ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    length=2.5,
    width=0.6,
)

ax.grid(
    which="major",
    linestyle="--",
    linewidth=0.45,
    alpha=0.18,
    zorder=0,
)


# ============================================================
# LEGEND
# ============================================================

ax.legend(
    frameon=False,
    loc="upper right",
    handlelength=2.7,
)


# ============================================================
# SAVE AND DISPLAY
# ============================================================

fig.savefig(
    OUTPUT_PNG,
    dpi=300,
    bbox_inches="tight",
)

fig.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
)

print(f"Saved PNG: {OUTPUT_PNG}")
print(f"Saved PDF: {OUTPUT_PDF}")

plt.show()
