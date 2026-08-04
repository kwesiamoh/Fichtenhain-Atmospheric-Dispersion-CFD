import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# ============================================================
# 1. LOAD DATA
# ============================================================

csv_file = "section_moments_windSpeed.csv"
df = pd.read_csv(csv_file)

required_columns = {
    "Uref_mps",
    "distance_m",
    "sigma_y_m",
}

missing = required_columns.difference(df.columns)

if missing:
    raise ValueError(
        f"Missing required column(s): {sorted(missing)}\n"
        f"Available columns: {list(df.columns)}"
    )


# ============================================================
# 2. PUBLICATION STYLE
# ============================================================

plt.rcParams.update({
    "figure.dpi": 300,
    "font.family": "DejaVu Sans",
    "font.size": 9,

    "axes.labelsize": 11,
    "axes.linewidth": 0.8,

    "xtick.labelsize": 9,
    "ytick.labelsize": 9,

    "legend.fontsize": 9,

    "lines.linewidth": 1.8,
    "lines.markersize": 5,

    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ============================================================
# 3. CREATE FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.16, 4.2),
    dpi=300
)


# ============================================================
# 4. SERIES STYLE
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


# ============================================================
# 5. PLOT
# ============================================================

for Uref in [5, 10, 15]:

    subset = (
        df[df["Uref_mps"] == Uref]
        .sort_values("distance_m")
    )

    style = styles[Uref]

    ax.plot(
        subset["distance_m"],
        subset["sigma_y_m"],

        label=rf"$U_{{\mathrm{{ref}}}} = {Uref}$ m s$^{{-1}}$",

        color=style["color"],
        marker=style["marker"],
        linestyle=style["linestyle"],

        linewidth=1.8,
        markersize=5.5,

        markerfacecolor="white",
        markeredgecolor=style["color"],
        markeredgewidth=1.2,

        zorder=3,
    )


# ============================================================
# 6. AXIS LABELS
# ============================================================

ax.set_xlabel("Downstream distance (m)")

ax.set_ylabel(
    r"Lateral plume spread, $\sigma_y$ (m)"
)


# ============================================================
# 7. X AXIS
# ============================================================

distances = sorted(df["distance_m"].unique())

ax.set_xticks(distances)

x_min = min(distances)
x_max = max(distances)

x_pad = 0.04 * (x_max - x_min)

ax.set_xlim(
    x_min - x_pad,
    x_max + x_pad
)

ax.xaxis.set_minor_locator(
    AutoMinorLocator(2)
)


# ============================================================
# 8. Y AXIS
# ============================================================

ax.set_ylim(bottom=0)

ax.yaxis.set_minor_locator(
    AutoMinorLocator(2)
)


# ============================================================
# 9. SPINES AND TICKS
# ============================================================

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.spines["left"].set_linewidth(0.8)
ax.spines["bottom"].set_linewidth(0.8)

ax.tick_params(
    axis="both",
    which="major",
    direction="in",
    length=5,
    width=0.8,
)

ax.tick_params(
    axis="both",
    which="minor",
    direction="in",
    length=3,
    width=0.6,
)


# ============================================================
# 10. GRID
# ============================================================

ax.grid(
    which="major",
    linestyle="--",
    linewidth=0.5,
    alpha=0.20,
    zorder=0,
)


# ============================================================
# 11. LEGEND
# ============================================================

ax.legend(
    loc="upper left",
    frameon=False,
    handlelength=2.8,
    handletextpad=0.7,
)


# ============================================================
# 12. LAYOUT + DISPLAY
# ============================================================

plt.tight_layout()
plt.show()
