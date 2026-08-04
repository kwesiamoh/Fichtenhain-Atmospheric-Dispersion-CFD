from pathlib import Path

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import AutoMinorLocator


# ============================================================
# 1. FILE
# ============================================================

ROOT = Path(r"C:\acfd2\Fichtenhain_Dispersion")

PARTICLE_FILE = (
    ROOT
    / "Fichtenhain_DPM_Particles_Prevailing"
    / "VTK"
    / "lagrangian"
    / "kinematicCloud"
    / "kinematicCloud_9636.vtp"
)

if not PARTICLE_FILE.exists():
    raise FileNotFoundError(
        f"Particle VTP not found:\n{PARTICLE_FILE}"
    )


# ============================================================
# 2. SOURCE AND PREVAILING FLOW
# ============================================================

X0 = 816.4
Y0 = 1486.0

# OpenFOAM prevailing airflow direction
FLOW_DIR = np.array(
    [-0.432672, 0.901551],
    dtype=float
)

FLOW_DIR /= np.linalg.norm(FLOW_DIR)

# Unit vector perpendicular to the prevailing airflow
CROSS_DIR = np.array([
    -FLOW_DIR[1],
    FLOW_DIR[0]
])


# ============================================================
# 3. READ PARTICLES
# ============================================================

cloud = pv.read(PARTICLE_FILE)

xyz = np.asarray(cloud.points)

if "age" not in cloud.point_data:
    raise KeyError(
        "'age' not found.\n"
        f"Available arrays: {list(cloud.point_data.keys())}"
    )

age = np.asarray(
    cloud.point_data["age"]
).ravel()

if xyz.shape[0] != age.size:
    raise ValueError(
        "Particle coordinates and age arrays "
        "have different lengths."
    )


# ============================================================
# 4. WIND-ALIGNED COORDINATES
# ============================================================

dx = xyz[:, 0] - X0
dy = xyz[:, 1] - Y0

relative_xy = np.column_stack([
    dx,
    dy
])

# Downwind distance
s = relative_xy @ FLOW_DIR

# Signed crosswind displacement
n = relative_xy @ CROSS_DIR


valid = (
    np.isfinite(s)
    & np.isfinite(n)
    & np.isfinite(age)
    & (s >= 0)
)

s = s[valid]
n = n[valid]
age = age[valid]


print(f"Parcels          : {len(s):,}")
print(f"Downwind range   : {s.min():.2f} to {s.max():.2f} m")
print(f"Crosswind range  : {n.min():.2f} to {n.max():.2f} m")
print(f"Age range        : {age.min():.2f} to {age.max():.2f} s")


# ============================================================
# 5. PUBLICATION STYLE
# ============================================================

plt.rcParams.update({
    "figure.dpi": 300,
    "font.family": "DejaVu Sans",
    "font.size": 9,

    "axes.labelsize": 11,
    "axes.titlesize": 10,
    "axes.linewidth": 0.8,

    "xtick.labelsize": 9,
    "ytick.labelsize": 9,

    "legend.fontsize": 8,

    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def style_axis(ax):
    """Apply the same journal-style formatting to both graphs."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        length=5,
        width=0.8
    )

    ax.xaxis.set_minor_locator(
        AutoMinorLocator(2)
    )

    ax.yaxis.set_minor_locator(
        AutoMinorLocator(2)
    )

    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        length=2.5,
        width=0.6
    )

    ax.grid(
        which="major",
        linestyle="--",
        linewidth=0.45,
        alpha=0.18
    )


# ============================================================
# FIGURE 6a
# LAGRANGIAN PARCEL DISTRIBUTION
# ============================================================

age_norm = Normalize(
    vmin=age.min(),
    vmax=age.max()
)

fig, ax = plt.subplots(
    figsize=(7.16, 3.8),
    dpi=300
)


sc = ax.scatter(
    s,
    n,

    c=age,

    cmap="viridis",
    norm=age_norm,

    s=2.5,
    alpha=0.48,

    edgecolors="none",
    rasterized=True
)


# Prevailing-flow centreline
ax.axhline(
    0,
    color="#666666",
    linestyle="--",
    linewidth=0.9,
    alpha=0.75
)


# Stack source
ax.scatter(
    0,
    0,

    marker="x",
    s=45,

    linewidth=1.6,
    color="#D55E00",

    zorder=10
)


ax.set_xlabel(
    "Downwind distance (m)"
)

ax.set_ylabel(
    "Crosswind displacement (m)"
)


cbar = fig.colorbar(
    sc,
    ax=ax,
    pad=0.025,
    fraction=0.035
)

cbar.set_label(
    "Parcel age (s)"
)

cbar.ax.tick_params(
    labelsize=8
)


style_axis(ax)

fig.tight_layout()

plt.show()


# ============================================================
# FIGURE 6b
# LATERAL PARCEL SPREAD
# ============================================================

N_BINS = 25

edges = np.linspace(
    s.min(),
    s.max(),
    N_BINS + 1
)

centres = 0.5 * (
    edges[:-1] + edges[1:]
)


median_n = np.full(
    N_BINS,
    np.nan
)

q10_n = np.full(
    N_BINS,
    np.nan
)

q90_n = np.full(
    N_BINS,
    np.nan
)


for i in range(N_BINS):

    if i == N_BINS - 1:
        mask = (
            (s >= edges[i])
            & (s <= edges[i + 1])
        )
    else:
        mask = (
            (s >= edges[i])
            & (s < edges[i + 1])
        )

    # Ignore bins containing too few parcels.
    if np.count_nonzero(mask) < 10:
        continue

    values = n[mask]

    median_n[i] = np.median(
        values
    )

    q10_n[i] = np.percentile(
        values,
        10
    )

    q90_n[i] = np.percentile(
        values,
        90
    )


valid_bins = (
    np.isfinite(median_n)
    & np.isfinite(q10_n)
    & np.isfinite(q90_n)
)


fig, ax = plt.subplots(
    figsize=(7.16, 3.8),
    dpi=300
)


# Parcel envelope
ax.fill_between(
    centres[valid_bins],

    q10_n[valid_bins],
    q90_n[valid_bins],

    color="#56B4E9",
    alpha=0.22,

    linewidth=0,

    label="10–90% envelope"
)


# Median parcel trajectory
ax.plot(
    centres[valid_bins],
    median_n[valid_bins],

    color="#D55E00",

    marker="o",
    markersize=4,

    markerfacecolor="white",
    markeredgecolor="#D55E00",

    linewidth=1.6,

    label="Median"
)


ax.axhline(
    0,
    color="#666666",
    linestyle="--",
    linewidth=0.9,
    alpha=0.75
)


ax.set_xlabel(
    "Downwind distance (m)"
)

ax.set_ylabel(
    "Crosswind displacement (m)"
)


ax.legend(
    frameon=False,
    loc="upper left"
)


style_axis(ax)

fig.tight_layout()

plt.show()
