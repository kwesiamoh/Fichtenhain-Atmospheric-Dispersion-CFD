from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt


BASE = Path(__file__).resolve().parent
INPUT_DIR = BASE / "input" / "particles"
SERIES_FILE = INPUT_DIR / "kinematicCloud.vtp.series"

OUTPUT_CSV = BASE / "tables" / "plume_geometry_by_distance.csv"
OUTPUT_FIGURE = BASE / "figures" / "plume_geometry_by_distance.png"

SOURCE = np.array([816.4, 1486.0, 70.0])

FLOW_DIR = np.array([-0.432672, 0.901551], dtype=float)
FLOW_DIR /= np.linalg.norm(FLOW_DIR)

CROSS_DIR = np.array([-FLOW_DIR[1], FLOW_DIR[0]])

BIN_WIDTH_M = 20.0
MIN_PARCELS_PER_BIN = 100


def weighted_percentile(values, weights, percentile):
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]

    cumulative = np.cumsum(weights)
    target = percentile / 100.0 * cumulative[-1]

    return float(values[np.searchsorted(cumulative, target)])


with SERIES_FILE.open("r", encoding="utf-8") as file:
    entries = json.load(file)["files"]

final_entry = entries[-1]
final_file = INPUT_DIR / final_entry["name"]

cloud = pv.read(final_file)
xyz = np.asarray(cloud.points, dtype=float)

n_particle = np.asarray(
    cloud.point_data["nParticle"], dtype=float
).ravel()

rho = np.asarray(
    cloud.point_data["rho"], dtype=float
).ravel()

diameter = np.asarray(
    cloud.point_data["d"], dtype=float
).ravel()

active = np.asarray(
    cloud.point_data["active"], dtype=float
).ravel()

mass = (
    n_particle
    * rho
    * (math.pi / 6.0)
    * diameter**3
)

dx = xyz[:, 0] - SOURCE[0]
dy = xyz[:, 1] - SOURCE[1]

downwind = dx * FLOW_DIR[0] + dy * FLOW_DIR[1]
crosswind = dx * CROSS_DIR[0] + dy * CROSS_DIR[1]
height_relative = xyz[:, 2] - SOURCE[2]

valid = (
    np.isfinite(downwind)
    & np.isfinite(crosswind)
    & np.isfinite(height_relative)
    & np.isfinite(mass)
    & (mass > 0.0)
    & (active > 0.5)
    & (downwind >= 0.0)
)

downwind = downwind[valid]
crosswind = crosswind[valid]
height_relative = height_relative[valid]
mass = mass[valid]

maximum_distance = math.ceil(
    downwind.max() / BIN_WIDTH_M
) * BIN_WIDTH_M

edges = np.arange(
    0.0,
    maximum_distance + BIN_WIDTH_M,
    BIN_WIDTH_M,
)

rows = []

for lower, upper in zip(edges[:-1], edges[1:]):

    mask = (downwind >= lower) & (downwind < upper)

    if np.count_nonzero(mask) < MIN_PARCELS_PER_BIN:
        continue

    s = downwind[mask]
    n = crosswind[mask]
    z = height_relative[mask]
    w = mass[mask]

    total_mass = np.sum(w)

    s_mean = np.sum(w * s) / total_mass
    n_mean = np.sum(w * n) / total_mass
    z_mean = np.sum(w * z) / total_mass

    sigma_n = np.sqrt(
        np.sum(w * (n - n_mean) ** 2) / total_mass
    )

    sigma_z = np.sqrt(
        np.sum(w * (z - z_mean) ** 2) / total_mass
    )

    rows.append({
        "downwind_distance_m": s_mean,
        "bin_start_m": lower,
        "bin_end_m": upper,
        "parcel_count": int(np.count_nonzero(mask)),
        "crosswind_centroid_m": n_mean,
        "plume_rise_m": z_mean,
        "plume_centroid_height_m": SOURCE[2] + z_mean,
        "sigma_crosswind_m": sigma_n,
        "sigma_vertical_m": sigma_z,
        "height_p10_m": weighted_percentile(z, w, 10),
        "height_p90_m": weighted_percentile(z, w, 90),
    })


df = pd.DataFrame(rows)
df.to_csv(OUTPUT_CSV, index=False)


fig, axes = plt.subplots(
    2,
    1,
    figsize=(7.16, 6.0),
    dpi=300,
    sharex=True,
)

axes[0].fill_between(
    df["downwind_distance_m"],
    df["height_p10_m"],
    df["height_p90_m"],
    alpha=0.20,
    label="10–90% envelope",
)

axes[0].plot(
    df["downwind_distance_m"],
    df["plume_rise_m"],
    linewidth=1.8,
    label="Mass-weighted centroid",
)

axes[0].axhline(0.0, linestyle="--", linewidth=0.8)
axes[0].set_ylabel("Height relative to stack exit (m)")
axes[0].legend(frameon=False)

axes[1].plot(
    df["downwind_distance_m"],
    df["sigma_crosswind_m"],
    linewidth=1.8,
    label="Crosswind spread",
)

axes[1].plot(
    df["downwind_distance_m"],
    df["sigma_vertical_m"],
    linewidth=1.8,
    label="Vertical spread",
)

axes[1].set_xlabel("Downwind distance from stack (m)")
axes[1].set_ylabel("Resolved parcel spread (m)")
axes[1].legend(frameon=False)

for ax in axes:
    ax.grid(
        linestyle="--",
        linewidth=0.5,
        alpha=0.25,
    )

plt.tight_layout()
plt.savefig(OUTPUT_FIGURE, dpi=300, bbox_inches="tight")
plt.close()


print("PLUME GEOMETRY BY DOWNWIND DISTANCE")
print("=" * 50)
print(f"Final VTP file       : {final_entry['name']}")
print(f"Valid parcels        : {len(downwind):,}")
print(f"Number of bins       : {len(df)}")
print(f"Maximum analysed bin : {df['bin_end_m'].max():.0f} m")
print()
print(df[
    [
        "downwind_distance_m",
        "plume_rise_m",
        "sigma_crosswind_m",
        "sigma_vertical_m",
    ]
].to_string(index=False))
print()
print(f"Saved table  : {OUTPUT_CSV}")
print(f"Saved figure : {OUTPUT_FIGURE}")
