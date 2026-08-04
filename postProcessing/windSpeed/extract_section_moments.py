from pathlib import Path
import csv
import math

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path("c:\\acfd2\\Fichtenhain_Dispersion")

SOURCE_Y = 1486.0
SOURCE_Z = 70.0

DISTANCES = [50, 100, 200, 300, 400, 500]

CASES = {
    5: ROOT / "Fichtenhain_Plume_Uref5"
    / "postProcessing" / "plumeSections" / "1800",

    10: ROOT / "Fichtenhain_Plume_Transient"
              / "postProcessing" / "plumeSections" / "1400",

    15: ROOT / "Fichtenhain_Plume_Uref15"
              / "postProcessing" / "plumeSections" / "1800",
}


# ============================================================
# READ ONE CROSS-SECTION
# ============================================================

def calculate_section_moments(vtp_file: Path):
    """
    Calculate area- and concentration-weighted plume statistics
    on one VTK cross-section.

    Returns:
        Tmax
        integrated_T
        y_centroid
        z_centroid
        sigma_y
        sigma_z
    """

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(vtp_file))
    reader.Update()

    poly = reader.GetOutput()

    if poly.GetNumberOfCells() == 0:
        raise RuntimeError(f"No cells found in {vtp_file}")

    # --------------------------------------------------------
    # Triangulate surface polygons
    # --------------------------------------------------------

    triangulate = vtk.vtkTriangleFilter()
    triangulate.SetInputData(poly)
    triangulate.Update()

    tri_poly = triangulate.GetOutput()

    # --------------------------------------------------------
    # Ensure T is available as cell data
    # --------------------------------------------------------

    cell_T = tri_poly.GetCellData().GetArray("T")

    if cell_T is None:
        point_T = tri_poly.GetPointData().GetArray("T")

        if point_T is None:
            raise RuntimeError(f"No T field found in {vtp_file}")

        point_to_cell = vtk.vtkPointDataToCellData()
        point_to_cell.SetInputData(tri_poly)
        point_to_cell.PassPointDataOn()
        point_to_cell.Update()

        tri_poly = point_to_cell.GetOutput()
        cell_T = tri_poly.GetCellData().GetArray("T")

    T_values = vtk_to_numpy(cell_T).astype(float)

    # --------------------------------------------------------
    # Check physical consistency
    # --------------------------------------------------------

    if np.any(T_values < -1.0e-12):
        raise ValueError(
            f"Negative scalar values detected in {vtp_file}: "
            f"minimum T = {T_values.min():.6e}"
        )

    # Remove negligible round-off negatives only.
    T_values = np.maximum(T_values, 0.0)

    # --------------------------------------------------------
    # Calculate triangle areas and centroids
    # --------------------------------------------------------

    areas = []
    y_positions = []
    z_positions = []

    for i in range(tri_poly.GetNumberOfCells()):

        cell = tri_poly.GetCell(i)

        if cell.GetNumberOfPoints() != 3:
            raise RuntimeError(
                f"Non-triangle cell found after triangulation "
                f"in {vtp_file}"
            )

        p0 = np.array(tri_poly.GetPoint(cell.GetPointId(0)))
        p1 = np.array(tri_poly.GetPoint(cell.GetPointId(1)))
        p2 = np.array(tri_poly.GetPoint(cell.GetPointId(2)))

        area = 0.5 * np.linalg.norm(
            np.cross(p1 - p0, p2 - p0)
        )

        centroid = (p0 + p1 + p2) / 3.0

        areas.append(area)
        y_positions.append(centroid[1])
        z_positions.append(centroid[2])

    areas = np.asarray(areas)
    y = np.asarray(y_positions)
    z = np.asarray(z_positions)

    if len(T_values) != len(areas):
        raise RuntimeError(
            f"T/cell mismatch in {vtp_file}: "
            f"{len(T_values)} T values vs {len(areas)} cells"
        )

    # --------------------------------------------------------
    # Concentration-weighted moments
    #
    # Weight = T * cross-sectional element area
    # --------------------------------------------------------

    weights = T_values * areas

    total_weight = np.sum(weights)

    if total_weight <= 0.0:
        return {
            "Tmax": 0.0,
            "integrated_T": 0.0,
            "y_centroid": np.nan,
            "z_centroid": np.nan,
            "sigma_y": np.nan,
            "sigma_z": np.nan,
        }

    y_centroid = np.sum(weights * y) / total_weight
    z_centroid = np.sum(weights * z) / total_weight

    variance_y = (
        np.sum(weights * (y - y_centroid) ** 2)
        / total_weight
    )

    variance_z = (
        np.sum(weights * (z - z_centroid) ** 2)
        / total_weight
    )

    sigma_y = math.sqrt(max(variance_y, 0.0))
    sigma_z = math.sqrt(max(variance_z, 0.0))

    return {
        "Tmax": float(np.max(T_values)),
        "integrated_T": float(total_weight),
        "y_centroid": float(y_centroid),
        "z_centroid": float(z_centroid),
        "sigma_y": float(sigma_y),
        "sigma_z": float(sigma_z),
    }


# ============================================================
# PROCESS ALL CASES
# ============================================================

rows = []

for Uref, folder in CASES.items():

    print()
    print("=" * 70)
    print(f"Uref = {Uref} m/s")
    print("=" * 70)

    for distance in DISTANCES:

        vtp_file = folder / f"x{distance:03d}.vtp"

        if not vtp_file.exists():
            raise FileNotFoundError(vtp_file)

        result = calculate_section_moments(vtp_file)

        y_shift = result["y_centroid"] - SOURCE_Y
        z_shift = result["z_centroid"] - SOURCE_Z

        row = {
            "Uref_mps": Uref,
            "distance_m": distance,

            "Tmax": result["Tmax"],
            "integrated_T": result["integrated_T"],

            "y_centroid_m": result["y_centroid"],
            "z_centroid_m": result["z_centroid"],

            "y_shift_m": y_shift,
            "z_shift_m": z_shift,

            "sigma_y_m": result["sigma_y"],
            "sigma_z_m": result["sigma_z"],
        }

        rows.append(row)

        print(
            f"{distance:3d} m | "
            f"Tmax={result['Tmax']:.4e} | "
            f"yc={result['y_centroid']:.2f} m | "
            f"zc={result['z_centroid']:.2f} m | "
            f"sigma_y={result['sigma_y']:.2f} m | "
            f"sigma_z={result['sigma_z']:.2f} m"
        )


# ============================================================
# WRITE CSV
# ============================================================

output = Path("section_moments_windSpeed.csv")

fieldnames = [
    "Uref_mps",
    "distance_m",

    "Tmax",
    "integrated_T",

    "y_centroid_m",
    "z_centroid_m",

    "y_shift_m",
    "z_shift_m",

    "sigma_y_m",
    "sigma_z_m",
]


with output.open("w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
    )

    writer.writeheader()
    writer.writerows(rows)


print()
print("=" * 70)
print(f"Saved quantitative data to:")
print(output.resolve())
print("=" * 70)
