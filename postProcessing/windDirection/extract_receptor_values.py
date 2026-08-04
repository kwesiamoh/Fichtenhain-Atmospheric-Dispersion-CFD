from pathlib import Path
import csv

import numpy as np


# ============================================================
# PATHS
# ============================================================

ROOT = Path(r"C:\acfd2\Fichtenhain_Dispersion")
OUTPUT_DIR = ROOT / "Results" / "quantitative" / "windDirection"
OUTPUT_CSV = OUTPUT_DIR / "Figure_5_receptor_values.csv"


# ============================================================
# CASE DEFINITIONS
# ============================================================

CASES = [
    {
        "scenario": "Wind from 0°",
        "scenario_order": 1,
        "path": (
            ROOT
            / "Fichtenhain_Plume_WD0"
            / "postProcessing"
            / "receptorRing"
            / "1800"
            / "T"
        ),
    },
    {
        "scenario": "Wind from 90°",
        "scenario_order": 2,
        "path": (
            ROOT
            / "Fichtenhain_Plume_WD90"
            / "postProcessing"
            / "receptorRing"
            / "1800"
            / "T"
        ),
    },
    {
        "scenario": "Wind from 180°",
        "scenario_order": 3,
        "path": (
            ROOT
            / "Fichtenhain_Plume_WD180"
            / "postProcessing"
            / "receptorRing"
            / "1800"
            / "T"
        ),
    },
    {
        "scenario": "Prevailing wind from 154°",
        "scenario_order": 4,
        "path": (
            ROOT
            / "Fichtenhain_Plume_NormalWind"
            / "postProcessing"
            / "receptorRing"
            / "1800"
            / "T"
        ),
    },
]


RECEPTORS = [
    "N",
    "NE",
    "E",
    "SE",
    "S",
    "SW",
    "W",
    "NW",
]

EXPECTED_TIME = 1800.0


# ============================================================
# PROBE FILE READER
# ============================================================

def read_final_probe_record(
    path: Path,
    number_of_receptors: int,
) -> tuple[float, np.ndarray]:
    """
    Read the final scalar record from an OpenFOAM probe file.
    """

    if not path.exists():
        raise FileNotFoundError(
            f"Probe file not found:\n{path}"
        )

    final_record = None

    with path.open(
        mode="r",
        encoding="utf-8",
        errors="ignore",
    ) as file:

        for raw_line in file:

            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            cleaned = (
                line
                .replace("(", " ")
                .replace(")", " ")
                .replace(",", " ")
            )

            parts = cleaned.split()

            try:
                numbers = [
                    float(value)
                    for value in parts
                ]
            except ValueError:
                continue

            if len(numbers) >= number_of_receptors + 1:
                final_record = numbers[
                    :number_of_receptors + 1
                ]

    if final_record is None:
        raise RuntimeError(
            f"No valid numerical records found in:\n{path}"
        )

    time = float(final_record[0])

    values = np.asarray(
        final_record[1:],
        dtype=float,
    )

    if not np.all(np.isfinite(values)):
        raise RuntimeError(
            f"Non-finite receptor values found in:\n{path}"
        )

    return time, values


# ============================================================
# EXTRACT DATA
# ============================================================

rows = []

print("Wind-direction receptor extraction")
print("=" * 78)

for case in CASES:

    time, values = read_final_probe_record(
        path=case["path"],
        number_of_receptors=len(RECEPTORS),
    )

    if not np.isclose(
        time,
        EXPECTED_TIME,
        rtol=0.0,
        atol=1.0e-6,
    ):
        raise RuntimeError(
            f"{case['scenario']} was sampled at "
            f"{time}, not {EXPECTED_TIME}."
        )

    maximum_index = int(
        np.argmax(values)
    )

    print(
        f"{case['scenario']:28s} | "
        f"time = {time:7.1f} s | "
        f"maximum receptor = "
        f"{RECEPTORS[maximum_index]:2s} | "
        f"T = {values[maximum_index]:.6e}"
    )

    for receptor_order, (
        receptor,
        value,
    ) in enumerate(
        zip(RECEPTORS, values),
        start=1,
    ):

        rows.append({
            "scenario": case["scenario"],
            "scenario_order": case["scenario_order"],
            "sample_time_s": time,
            "receptor": receptor,
            "receptor_order": receptor_order,
            "radius_m": 200.0,
            "height_m": 70.0,
            "T": value,
        })


# ============================================================
# WRITE LONG-FORMAT CSV
# ============================================================

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

fieldnames = [
    "scenario",
    "scenario_order",
    "sample_time_s",
    "receptor",
    "receptor_order",
    "radius_m",
    "height_m",
    "T",
]

with OUTPUT_CSV.open(
    mode="w",
    newline="",
    encoding="utf-8",
) as file:

    writer = csv.DictWriter(
        file,
        fieldnames=fieldnames,
    )

    writer.writeheader()

    for row in rows:

        formatted_row = row.copy()

        formatted_row["sample_time_s"] = (
            f"{row['sample_time_s']:.6f}"
        )

        formatted_row["radius_m"] = (
            f"{row['radius_m']:.1f}"
        )

        formatted_row["height_m"] = (
            f"{row['height_m']:.1f}"
        )

        formatted_row["T"] = (
            f"{row['T']:.12e}"
        )

        writer.writerow(
            formatted_row
        )


print()
print(f"Rows written: {len(rows)}")
print(f"Saved CSV: {OUTPUT_CSV}")
