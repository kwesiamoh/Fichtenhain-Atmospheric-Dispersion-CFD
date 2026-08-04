from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt


# ============================================================
# PATHS
# ============================================================

BASE = Path(__file__).resolve().parent
INPUT_DIR = BASE / "input" / "particles"
SERIES_FILE = INPUT_DIR / "kinematicCloud.vtp.series"

OUTPUT_CSV = BASE / "tables" / "plume_rise_and_spread.csv"
OUTPUT_FIGURE = BASE / "figures" / "plume_rise_and_spread.png"


# ============================================================
# SOURCE AND FLOW DIRECTION
# ============================================================

SOURCE = np.array([816.4, 1486.0, 70.0])
INJECTION_START_TIME = 1846.9969433376684

FLOW_DIR = np.array([-0.432672, 0.901551])
FLOW_DIR /= np.linalg.norm(FLOW_DIR)

CROSS_DIR = np.array([
    -FLOW_DIR[1],
    FLOW_DIR[0],
])


# ============================================================
# WEIGHTED PERCENTILE
# ============================================================

def weighted_percentile(values, weights, percentile):
    order = np.argsort(values)

    values = values[order]
    weights = weights[order]

    cumulative = np.cumsum(weights)
    target = percentile / 100.0 * cumulative[-1]

    return float(values[np.searchsorted(cumulative, target)])


# ============================================================
# READ SERIES
# ============================================================

with SERIES_FILE.open("r", encoding="utf-8") as file:
    entries = json.load(file)["files"]


rows = []


# ============================================================
# PROCESS EACH EXPORTED TIME
# ============================================================

for index, entry in enumerate(entries, start=1):

    cloud = pv.read(INPUT_DIR / entry["name"])
    xyz = np.asarray(cloud.points, dtype=float)

    n_particle = np.asarray(
        cloud.point_data["nParticle"],
        dtype=float,
    ).ravel()

    rho = np.asarray(
        cloud.point_data["rho"],
        dtype=float,
    ).ravel()

    diameter = np.asarray(
        cloud.point_data["d"],
        dtype=float,
    ).ravel()

    mass = (
        n_particle
        * rho
        * (math.pi / 6.0)
        * diameter**3
    )

    dx = xyz[:, 0] - SOURCE[0]
    dy = xyz[:, 1] - SOURCE[1]

    downwind = (
        dx * FLOW_DIR[0]
        + dy * FLOW_DIR[1]
    )

    crosswind = (
        dx * CROSS_DIR[0]
        + dy * CROSS_DIR[1]
    )

    height_relative = xyz[:, 2] - SOURCE[2]

    # Analyse the downstream plume only.
    mask = (
        np.isfinite(downwind)
        & np.isfinite(crosswind)
        & np.isfinite(height_relative)
        & np.isfinite(mass)
        & (mass > 0.0)
        & (downwind >= 0.0)
    )

    crosswind = crosswind[mask]
    height_relative = height_relative[mask]
    weights = mass[mask]

    total_weight = np.sum(weights)

    crosswind_centroid = np.sum(
        weights * crosswind
    ) / total_weight

    plume_rise = np.sum(
        weights * height_relative
    ) / total_weight

    sigma_crosswind = np.sqrt(
        np.sum(
            weights
            * (crosswind - crosswind_centroid) ** 2
        )
        / total_weight
    )

    sigma_vertical = np.sqrt(
        np.sum(
            weights
            * (height_relative - plume_rise) ** 2
        )
        / total_weight
    )

    simulation_time = float(entry["time"])
    relative_time = simulation_time - INJECTION_START_TIME

    rows.append({
        "relative_time_s": relative_time,
        "downstream_parcel_count": int(np.count_nonzero(mask)),
        "crosswind_centroid_m": crosswind_centroid,
        "plume_centroid_height_m": SOURCE[2] + plume_rise,
        "plume_rise_above_stack_m": plume_rise,
        "sigma_crosswind_m": sigma_crosswind,
        "sigma_vertical_m": sigma_vertical,
        "crosswind_p10_m": weighted_percentile(
            crosswind, weights, 10
        ),
        "crosswind_p90_m": weighted_percentile(
            crosswind, weights, 90
        ),
        "height_p10_relative_m": weighted_percentile(
            height_relative, weights, 10
        ),
        "height_p90_relative_m": weighted_percentile(
            height_relative, weights, 90
        ),
    })

    print(
        f"[{index:02d}/{len(entries):02d}] "
        f"t={relative_time:6.2f} s | "
        f"rise={plume_rise:6.2f} m | "
        f"sigma_n={sigma_crosswind:6.2f} m | "
        f"sigma_z={sigma_vertical:6.2f} m"
    )


# ============================================================
# SAVE TABLE
# ============================================================

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_CSV, index=False)


# ============================================================
# PLOT
# ============================================================

fig, axes = plt.subplots(
    2,
    1,
    figsize=(7.16, 6.0),
    dpi=300,
    sharex=True,
)

# Vertical position
axes[0].fill_between(
    df["relative_time_s"],
    df["height_p10_relative_m"],
    df["height_p90_relative_m"],
    alpha=0.20,
    label="10–90% vertical envelope",
)

axes[0].plot(
    df["relative_time_s"],
    df["plume_rise_above_stack_m"],
    linewidth=1.8,
    label="Mass-weighted centroid",
)

axes[0].axhline(
    0.0,
    linestyle="--",
    linewidth=0.8,
)

axes[0].set_ylabel("Height relative to stack exit (m)")
axes[0].legend(frameon=False)


# Resolved spread
axes[1].plot(
    df["relative_time_s"],
    df["sigma_crosswind_m"],
    linewidth=1.8,
    label="Crosswind spread",
)

axes[1].plot(
    df["relative_time_s"],
    df["sigma_vertical_m"],
    linewidth=1.8,
    label="Vertical spread",
)

axes[1].set_xlabel("Time after release began (s)")
axes[1].set_ylabel("Resolved parcel spread (m)")
axes[1].legend(frameon=False)


for ax in axes:
    ax.grid(
        linestyle="--",
        linewidth=0.5,
        alpha=0.25,
    )

axes[1].set_xlim(0, 60)

plt.tight_layout()

plt.savefig(
    OUTPUT_FIGURE,
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# ============================================================
# FINAL SUMMARY
# ============================================================

last = df.iloc[-1]

print()
print("FINAL PLUME RISE AND SPREAD")
print("=" * 45)
print(
    f"Crosswind centroid    : "
    f"{last['crosswind_centroid_m']:.2f} m"
)
print(
    f"Centroid height       : "
    f"{last['plume_centroid_height_m']:.2f} m"
)
print(
    f"Rise above stack      : "
    f"{last['plume_rise_above_stack_m']:.2f} m"
)
print(
    f"Crosswind spread      : "
    f"{last['sigma_crosswind_m']:.2f} m"
)
print(
    f"Vertical spread       : "
    f"{last['sigma_vertical_m']:.2f} m"
)
print(
    f"Crosswind 10–90%      : "
    f"{last['crosswind_p10_m']:.2f} to "
    f"{last['crosswind_p90_m']:.2f} m"
)
print(
    f"Vertical 10–90%       : "
    f"{last['height_p10_relative_m']:.2f} to "
    f"{last['height_p90_relative_m']:.2f} m "
    f"relative to stack"
)
print()
print(f"Saved table  : {OUTPUT_CSV}")
print(f"Saved figure : {OUTPUT_FIGURE}")
