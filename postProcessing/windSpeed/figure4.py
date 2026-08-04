from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from matplotlib.colors import LogNorm
from matplotlib.tri import Triangulation


# ============================================================
# FILE PATHS
# ============================================================

ROOT = Path(r"C:\acfd2\Fichtenhain_Dispersion")
BASE = ROOT / "Results" / "quantitative" / "windDirection"

OUTPUT_PNG = BASE / "Figure_4.png"
OUTPUT_PDF = BASE / "Figure_4.pdf"


files = {
    r"Wind from 0° (N $\rightarrow$ S)":
        ROOT / "Fichtenhain_Plume_WD0"
        / "postProcessing" / "plumeFootprint"
        / "1800" / "releaseHeight.vtp",

    r"Wind from 90° (E $\rightarrow$ W)":
        ROOT / "Fichtenhain_Plume_WD90"
        / "postProcessing" / "plumeFootprint"
        / "1800" / "releaseHeight.vtp",

    r"Wind from 270° (W $\rightarrow$ E)":
        ROOT / "Fichtenhain_Plume_WD270"
        / "postProcessing" / "plumeFootprint"
        / "1800" / "releaseHeight.vtp",

    r"Prevailing wind from 154°":
        ROOT / "Fichtenhain_Plume_NormalWind"
        / "postProcessing" / "plumeFootprint"
        / "1800" / "releaseHeight.vtp",
}


# Stack coordinates
X_STACK = 816.4
Y_STACK = 1486.0


# ============================================================
# LOAD DATA
# ============================================================

datasets = {}
positive_values = []

for label, path in files.items():

    if not path.exists():
        raise FileNotFoundError(
            f"Missing plume-plane file:\n{path}"
        )

    mesh = pv.read(path)

    if "T" not in mesh.point_data:

        if "T" in mesh.cell_data:
            mesh = mesh.cell_data_to_point_data()

        else:
            raise KeyError(
                f"Field T was not found in:\n{path}"
            )

    surface = mesh.extract_surface().triangulate()

    points = surface.points

    # Use source-relative coordinates
    x_relative = points[:, 0] - X_STACK
    y_relative = points[:, 1] - Y_STACK

    tracer = np.asarray(
        surface.point_data["T"]
    ).ravel()

    faces = surface.faces.reshape(-1, 4)[:, 1:]

    datasets[label] = {
        "x": x_relative,
        "y": y_relative,
        "T": tracer,
        "triangles": faces,
    }

    positive = tracer[
        np.isfinite(tracer) & (tracer > 0)
    ]

    if positive.size:
        positive_values.append(positive)


if not positive_values:
    raise RuntimeError(
        "No positive T values were found in any case."
    )


# ============================================================
# COMMON COLOUR SCALE
# ============================================================

all_positive = np.concatenate(positive_values)

VMIN = 1.0e-7

data_max = float(np.max(all_positive))

# Round upper scale to the next power of ten
VMAX = 10.0 ** np.ceil(
    np.log10(data_max)
)

if VMAX <= VMIN:
    VMAX = 10.0 * VMIN

print(f"Maximum T in all four cases: {data_max:.6e}")
print(f"Colour scale: {VMIN:.1e} to {VMAX:.1e}")


# ============================================================
# COMMON SPATIAL LIMITS
# ============================================================

all_x = np.concatenate([
    data["x"] for data in datasets.values()
])

all_y = np.concatenate([
    data["y"] for data in datasets.values()
])

x_limits = (
    float(np.min(all_x)),
    float(np.max(all_x)),
)

y_limits = (
    float(np.min(all_y)),
    float(np.max(all_y)),
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
    "axes.titlesize": 9.5,
    "axes.linewidth": 0.8,

    "xtick.labelsize": 8,
    "ytick.labelsize": 8,

    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ============================================================
# CREATE FIGURE
# ============================================================

fig, axes = plt.subplots(
    2,
    2,
    figsize=(7.16, 6.2),
    sharex=True,
    sharey=True,
    constrained_layout=True,
)

axes = axes.ravel()

panel_labels = [
    "(a)",
    "(b)",
    "(c)",
    "(d)",
]

norm = LogNorm(
    vmin=VMIN,
    vmax=VMAX,
)

mappable = None


# ============================================================
# DRAW PANELS
# ============================================================

for ax, panel, (title, data) in zip(
    axes,
    panel_labels,
    datasets.items(),
):

    triangulation = Triangulation(
        data["x"],
        data["y"],
        data["triangles"],
    )

    tracer_masked = np.ma.masked_less(
        data["T"],
        VMIN,
    )

    mappable = ax.tripcolor(
        triangulation,
        tracer_masked,
        shading="gouraud",
        cmap="viridis",
        norm=norm,
        rasterized=True,
    )

    # Source position in relative coordinates
    ax.scatter(
        0.0,
        0.0,
        marker="x",
        s=38,
        linewidths=1.3,
        color="red",
        zorder=5,
    )

    ax.set_title(
        title,
        pad=6,
    )

    ax.text(
        0.025,
        0.96,
        panel,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
    )

    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_aspect("equal")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        length=4,
        width=0.8,
    )

    ax.minorticks_on()

    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        length=2.5,
        width=0.6,
    )


# ============================================================
# SHARED LABELS AND COLOUR BAR
# ============================================================

fig.supxlabel(
    "Relative east-west distance from stack (m)"
)

fig.supylabel(
    "Relative north-south distance from stack (m)"
)

colour_bar = fig.colorbar(
    mappable,
    ax=axes.tolist(),
    fraction=0.035,
    pad=0.025,
)

colour_bar.set_label(
    r"Normalized gaseous contaminant concentration, $T$ (-)"
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
