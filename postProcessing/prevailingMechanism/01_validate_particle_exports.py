from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import pyvista as pv


# ============================================================
# PATHS
# ============================================================

BASE = Path(__file__).resolve().parent
INPUT_DIR = BASE / "input" / "particles"
TABLE_DIR = BASE / "tables"

SERIES_FILE = INPUT_DIR / "kinematicCloud.vtp.series"
OUTPUT_CSV = TABLE_DIR / "particle_export_inventory.csv"
OUTPUT_REPORT = TABLE_DIR / "particle_export_validation.txt"

TABLE_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# CASE INFORMATION
# ============================================================

INJECTION_START_TIME = 1846.9969433376684

REQUIRED_SCALAR_FIELDS = [
    "age",
    "nParticle",
    "rho",
    "d",
    "active",
]

REQUIRED_VECTOR_FIELDS = [
    "U",
    "UTurb",
    "UCorrect",
]


# ============================================================
# HELPERS
# ============================================================

def require_finite(name: str, values: np.ndarray, filename: str) -> None:
    """Raise an error when an array contains NaN or infinity."""

    if not np.all(np.isfinite(values)):
        raise ValueError(
            f"{filename}: field '{name}' contains NaN or infinity."
        )


def parcel_mass(
    n_particle: np.ndarray,
    rho: np.ndarray,
    diameter: np.ndarray,
) -> np.ndarray:
    """
    Represented mass of each computational parcel.

    mass = nParticle × rho × particle volume
    """

    volume = (math.pi / 6.0) * diameter**3
    return n_particle * rho * volume


# ============================================================
# READ SERIES FILE
# ============================================================

if not SERIES_FILE.exists():
    raise FileNotFoundError(
        f"Series file not found:\n{SERIES_FILE}"
    )

with SERIES_FILE.open("r", encoding="utf-8") as handle:
    series_data = json.load(handle)

entries = series_data.get("files", [])

if not entries:
    raise RuntimeError(
        f"No VTP entries were found in:\n{SERIES_FILE}"
    )

filenames = [entry["name"] for entry in entries]
times = np.asarray(
    [float(entry["time"]) for entry in entries],
    dtype=float,
)

errors: list[str] = []
warnings: list[str] = []


# ============================================================
# SERIES-LEVEL CHECKS
# ============================================================

if len(filenames) != len(set(filenames)):
    errors.append("Duplicate filenames are present in the series file.")

if not np.all(np.diff(times) > 0.0):
    errors.append("Simulation times are not strictly increasing.")

listed_files = {INPUT_DIR / name for name in filenames}
actual_files = set(INPUT_DIR.glob("*.vtp"))

missing_files = sorted(listed_files - actual_files)
extra_files = sorted(actual_files - listed_files)

if missing_files:
    errors.append(
        "Files listed in the series file are missing:\n"
        + "\n".join(str(path) for path in missing_files)
    )

if extra_files:
    warnings.append(
        "VTP files exist but are not listed in the series file:\n"
        + "\n".join(str(path) for path in extra_files)
    )

first_relative_time = times[0] - INJECTION_START_TIME

if first_relative_time > 0.0:
    warnings.append(
        "The first exported cloud occurs after the injection start time. "
        "This is expected when the cloud is empty at the exact start time. "
        f"First exported relative time = {first_relative_time:.6f} s."
    )


# ============================================================
# VALIDATE EACH VTP FILE
# ============================================================

rows: list[dict[str, float | int | str]] = []

for index, entry in enumerate(entries, start=1):

    filename = entry["name"]
    simulation_time = float(entry["time"])
    relative_time = simulation_time - INJECTION_START_TIME
    path = INPUT_DIR / filename

    if not path.exists():
        continue

    try:
        cloud = pv.read(path)

        number_of_parcels = int(cloud.n_points)

        if number_of_parcels <= 0:
            raise ValueError("The VTP file contains no parcel points.")

        coordinates = np.asarray(cloud.points, dtype=float)

        if coordinates.shape != (number_of_parcels, 3):
            raise ValueError(
                "Parcel coordinates do not have the expected N × 3 shape."
            )

        require_finite("coordinates", coordinates, filename)

        available_fields = set(cloud.point_data.keys())

        required_fields = set(
            REQUIRED_SCALAR_FIELDS + REQUIRED_VECTOR_FIELDS
        )

        missing_fields = sorted(
            required_fields - available_fields
        )

        if missing_fields:
            raise KeyError(
                "Missing point-data fields: "
                + ", ".join(missing_fields)
            )

        arrays: dict[str, np.ndarray] = {}

        for field in REQUIRED_SCALAR_FIELDS:
            values = np.asarray(
                cloud.point_data[field]
            ).reshape(-1)

            if values.size != number_of_parcels:
                raise ValueError(
                    f"Field '{field}' contains {values.size} values, "
                    f"but the cloud contains {number_of_parcels} parcels."
                )

            require_finite(field, values, filename)
            arrays[field] = values.astype(float)

        for field in REQUIRED_VECTOR_FIELDS:
            values = np.asarray(
                cloud.point_data[field],
                dtype=float,
            )

            if values.shape != (number_of_parcels, 3):
                raise ValueError(
                    f"Vector field '{field}' has shape {values.shape}; "
                    f"expected ({number_of_parcels}, 3)."
                )

            require_finite(field, values, filename)
            arrays[field] = values

        if np.any(arrays["age"] < -1.0e-10):
            raise ValueError("Negative parcel ages were detected.")

        if np.any(arrays["nParticle"] < 0.0):
            raise ValueError("Negative nParticle values were detected.")

        if np.any(arrays["rho"] <= 0.0):
            raise ValueError("Non-positive parcel densities were detected.")

        if np.any(arrays["d"] <= 0.0):
            raise ValueError("Non-positive parcel diameters were detected.")

        represented_mass = parcel_mass(
            arrays["nParticle"],
            arrays["rho"],
            arrays["d"],
        )

        require_finite(
            "represented parcel mass",
            represented_mass,
            filename,
        )

        velocity_magnitude = np.linalg.norm(
            arrays["U"],
            axis=1,
        )

        rows.append({
            "sequence_index": index,
            "filename": filename,
            "simulation_time_s": simulation_time,
            "relative_time_s": relative_time,
            "parcel_count": number_of_parcels,
            "active_count": int(
                np.count_nonzero(arrays["active"] > 0.5)
            ),
            "age_min_s": float(np.min(arrays["age"])),
            "age_mean_s": float(np.mean(arrays["age"])),
            "age_max_s": float(np.max(arrays["age"])),
            "velocity_mean_mps": float(
                np.mean(velocity_magnitude)
            ),
            "velocity_max_mps": float(
                np.max(velocity_magnitude)
            ),
            "nParticle_sum": float(
                np.sum(arrays["nParticle"])
            ),
            "represented_mass_kg": float(
                np.sum(represented_mass)
            ),
            "rho_min_kgm3": float(
                np.min(arrays["rho"])
            ),
            "rho_max_kgm3": float(
                np.max(arrays["rho"])
            ),
            "diameter_min_m": float(
                np.min(arrays["d"])
            ),
            "diameter_max_m": float(
                np.max(arrays["d"])
            ),
        })

        print(
            f"[{index:02d}/{len(entries):02d}] "
            f"{filename}: "
            f"{number_of_parcels:,} parcels, "
            f"mass={np.sum(represented_mass):.6e} kg"
        )

    except Exception as exc:
        errors.append(f"{filename}: {exc}")


# ============================================================
# WRITE INVENTORY CSV
# ============================================================

if rows:
    fieldnames = list(rows[0].keys())

    with OUTPUT_CSV.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# WRITE VALIDATION REPORT
# ============================================================

report_lines = [
    "CORRECTED PARTICLE EXPORT VALIDATION",
    "=" * 60,
    f"Input directory: {INPUT_DIR}",
    f"Series entries: {len(entries)}",
    f"VTP files present: {len(actual_files)}",
    f"Successfully validated: {len(rows)}",
    "",
    f"Injection start time: {INJECTION_START_TIME:.15f}",
    f"First exported time: {times[0]:.8f}",
    f"Last exported time: {times[-1]:.8f}",
    f"First relative time: {first_relative_time:.8f}",
]

if rows:
    report_lines.extend([
        "",
        f"First parcel count: {rows[0]['parcel_count']:,}",
        f"Final parcel count: {rows[-1]['parcel_count']:,}",
        (
            "Final represented mass: "
            f"{rows[-1]['represented_mass_kg']:.8e} kg"
        ),
        (
            "Final maximum age: "
            f"{rows[-1]['age_max_s']:.4f} s"
        ),
        (
            "Final maximum parcel speed: "
            f"{rows[-1]['velocity_max_mps']:.4f} m/s"
        ),
    ])

report_lines.extend([
    "",
    f"Warnings: {len(warnings)}",
])

for warning in warnings:
    report_lines.extend([
        "-" * 60,
        warning,
    ])

report_lines.extend([
    "",
    f"Errors: {len(errors)}",
])

for error in errors:
    report_lines.extend([
        "-" * 60,
        error,
    ])

status = "PASS" if not errors else "FAIL"

report_lines.extend([
    "",
    "=" * 60,
    f"VALIDATION RESULT: {status}",
])

OUTPUT_REPORT.write_text(
    "\n".join(report_lines) + "\n",
    encoding="utf-8",
)

print()
print("\n".join(report_lines))
print()
print(f"Inventory CSV: {OUTPUT_CSV}")
print(f"Report:        {OUTPUT_REPORT}")

if errors:
    sys.exit(1)
