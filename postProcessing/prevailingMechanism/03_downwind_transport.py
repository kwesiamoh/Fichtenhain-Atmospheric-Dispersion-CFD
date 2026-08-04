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

OUTPUT_CSV = BASE / "tables" / "downwind_transport.csv"
OUTPUT_FIGURE = BASE / "figures" / "downwind_transport.png"


# ============================================================
# SOURCE AND WIND DIRECTION
# ============================================================

SOURCE_X = 816.4
SOURCE_Y = 1486.0
INJECTION_START_TIME = 1846.9969433376684

FLOW_DIR = np.array(
    [-0.432672, 0.901551],
    dtype=float,
)

FLOW_DIR /= np.linalg.norm(FLOW_DIR)


# ============================================================
# WEIGHTED PERCENTILE
# ============================================================

def weighted_percentile(values, weights, percentile):
    order = np.argsort(values)

    values = values[order]
    weights = weights[order]

    cumulative = np.cumsum(weights)
    target = percentile / 100.0 * cumulative[-1]

    return float(
        values[np.searchsorted(cumulative, target)]
    )


# ============================================================
# READ SERIES
# ============================================================

with SERIES_FILE.open("r", encoding="utf-8") as file:
    entries = json.load(file)["files"]


rows = []


# ============================================================
# PROCESS EACH TIME
# ============================================================

for index, entry in enumerate(entries, start=1):

    path = INPUT_DIR / entry["name"]
    simulation_time = float(entry["time"])
    relative_time = simulation_time - INJECTION_START_TIME

    cloud = pv.read(path)

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

    dx = xyz[:, 0] - SOURCE_X
    dy = xyz[:, 1] - SOURCE_Y

    downwind = (
        dx * FLOW_DIR[0]
        + dy * FLOW_DIR[1]
    )

    total_mass = np.sum(mass)

    mean_distance = np.sum(
        mass * downwind
    ) / total_mass

    upstream_mass = np.sum(
        mass[downwind < 0.0]
    )

    rows.append({
        "relative_time_s": relative_time,
        "parcel_count": cloud.n_points,
        "mean_downwind_m": mean_distance,
        "median_downwind_m": weighted_percentile(
            downwind, mass, 50
        ),
        "p90_downwind_m": weighted_percentile(
            downwind, mass, 90
        ),
        "p95_downwind_m": weighted_percentile(
            downwind, mass, 95
        ),
        "maximum_downwind_m": float(
            np.max(downwind)
        ),
        "minimum_downwind_m": float(
            np.min(downwind)
        ),
        "upstream_mass_fraction": (
            upstream_mass / total_mass
        ),
    })

    print(
        f"[{index:02d}/{len(entries):02d}] "
        f"t={relative_time:6.2f} s | "
        f"mean={mean_distance:7.2f} m | "
        f"P95={rows[-1]['p95_downwind_m']:7.2f} m | "
        f"max={rows[-1]['maximum_downwind_m']:7.2f} m"
    )


# ============================================================
# SAVE TABLE
# ============================================================

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_CSV, index=False)


# ============================================================
# PLOT
# ============================================================

plt.figure(figsize=(7.16, 4.2), dpi=300)

plt.plot(
    df["relative_time_s"],
    df["mean_downwind_m"],
    label="Mass-weighted mean",
    linewidth=1.8,
)

plt.plot(
    df["relative_time_s"],
    df["median_downwind_m"],
    label="Mass-weighted median",
    linewidth=1.8,
)

plt.plot(
    df["relative_time_s"],
    df["p95_downwind_m"],
    label="95th percentile",
    linewidth=1.8,
)

plt.xlabel("Time after release began (s)")
plt.ylabel("Downwind distance from stack (m)")

plt.xlim(0, 60)
plt.ylim(bottom=0)

plt.grid(
    linestyle="--",
    linewidth=0.5,
    alpha=0.25,
)

plt.legend(frameon=False)
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
print("FINAL DOWNWIND TRANSPORT")
print("=" * 45)
print(f"Mean distance       : {last['mean_downwind_m']:.2f} m")
print(f"Median distance     : {last['median_downwind_m']:.2f} m")
print(f"90th percentile     : {last['p90_downwind_m']:.2f} m")
print(f"95th percentile     : {last['p95_downwind_m']:.2f} m")
print(f"Maximum distance    : {last['maximum_downwind_m']:.2f} m")
print(
    f"Upstream mass       : "
    f"{100*last['upstream_mass_fraction']:.4f}%"
)
print()
print(f"Saved table  : {OUTPUT_CSV}")
print(f"Saved figure : {OUTPUT_FIGURE}")
