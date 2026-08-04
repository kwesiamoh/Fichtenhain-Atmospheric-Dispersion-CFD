from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# ============================================================
# FILE PATHS
# ============================================================

BASE = Path(__file__).resolve().parent

INPUT_FILE = BASE / "section_moments_windSpeed.csv"
OUTPUT_PNG = BASE / "Figure_2.png"
OUTPUT_PDF = BASE / "Figure_2.pdf"


# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT_FILE)

required_columns = {
    "Uref_mps",
    "distance_m",
    "sigma_y_m",
    "sigma_z_m",
}

missing_columns = required_columns.difference(df.columns)

if missing_columns:
    raise ValueError(
        f"Missing required columns: {sorted(missing_columns)}\n"
        f"Available columns: {list(df.columns)}"
    )


# ============================================================
# PLOT STYLES
# ============================================================

styles = {
    5: {
        "color": "#0072B2",
        "marker": "o",
        "linestyle": "-",
    },
    10: {
        "color": "#D55E00",
        "marker": "s",
        "linestyle": "--",
    },
    15: {
        "color": "#009E73",
        "marker": "^",
        "linestyle": "-.",
    },
}


# Moderate screen DPI prevents the display window from becoming
# excessively large. High DPI is applied only when saving.
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

    "lines.linewidth": 1.6,
    "lines.markersize": 5.5,

    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ============================================================
# CREATE FIGURE
# ============================================================

fig, axes = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(6.4, 4.8),
    sharex=True,
    constrained_layout=True,
)


# ============================================================
# PLOT EACH WIND-SPEED CASE
# ============================================================

for speed in [5, 10, 15]:

    subset = (
        df[df["Uref_mps"] == speed]
        .sort_values("distance_m")
    )

    if subset.empty:
        raise ValueError(
            f"No data found for Uref = {speed} m/s."
        )

    style = styles[speed]

    label = (
        rf"$U_{{\mathrm{{ref}}}} = {speed}$ "
        rf"m s$^{{-1}}$"
    )

    # Lateral plume spread
    axes[0].plot(
        subset["distance_m"],
        subset["sigma_y_m"],
        label=label,
        color=style["color"],
        marker=style["marker"],
        linestyle=style["linestyle"],
        markerfacecolor="white",
        markeredgecolor=style["color"],
        markeredgewidth=1.1,
        linewidth=1.6,
        markersize=5.5,
        zorder=3,
    )

    # Vertical plume spread
    axes[1].plot(
        subset["distance_m"],
        subset["sigma_z_m"],
        color=style["color"],
        marker=style["marker"],
        linestyle=style["linestyle"],
        markerfacecolor="white",
        markeredgecolor=style["color"],
        markeredgewidth=1.1,
        linewidth=1.6,
        markersize=5.5,
        zorder=3,
    )


# ============================================================
# AXIS LABELS
# ============================================================

axes[0].set_ylabel(
    r"Lateral spread, $\sigma_y$ (m)"
)

axes[1].set_ylabel(
    r"Vertical spread, $\sigma_z$ (m)"
)

axes[1].set_xlabel(
    "Downstream distance (m)"
)


# ============================================================
# PANEL LABELS
# ============================================================

axes[0].text(
    0.98,
    0.10,
    "(a)",
    transform=axes[0].transAxes,
    fontsize=10,
    fontweight="bold",
    ha="right",
    va="top",
)

axes[1].text(
    0.94,
    0.10,
    "(b)",
    transform=axes[1].transAxes,
    fontsize=10,
    fontweight="bold",
    ha="left",
    va="top",
)


# ============================================================
# AXIS FORMATTING
# ============================================================

distances = sorted(
    df["distance_m"].unique()
)

axes[1].set_xticks(distances)

x_min = min(distances)
x_max = max(distances)
x_padding = 0.04 * (x_max - x_min)

for ax in axes:

    ax.set_xlim(
        x_min - x_padding,
        x_max + x_padding,
    )

    ax.set_ylim(bottom=0)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.xaxis.set_minor_locator(
        AutoMinorLocator(2)
    )

    ax.yaxis.set_minor_locator(
        AutoMinorLocator(2)
    )

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

axes[0].legend(
    frameon=False,
    loc="upper left",
    handlelength=2.7,
    handletextpad=0.7,
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
