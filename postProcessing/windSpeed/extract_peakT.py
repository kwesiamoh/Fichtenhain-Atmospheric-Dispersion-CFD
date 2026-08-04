from pathlib import Path
import csv
import vtk
from vtk.util.numpy_support import vtk_to_numpy


ROOT = Path("c:\\acfd2\\Fichtenhain_Dispersion")

cases = {
    5: ROOT / "Fichtenhain_Plume_Uref5" / "postProcessing" / "plumeSections" / "1800",
    10: ROOT / "Fichtenhain_Plume_Transient" / "postProcessing" / "plumeSections" / "1400",
    15: ROOT / "Fichtenhain_Plume_Uref15" / "postProcessing" / "plumeSections" / "1800",
}

distances = [50, 100, 200, 300, 400, 500]


def read_T(vtp_file):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(vtp_file))
    reader.Update()

    data = reader.GetOutput()

    # OpenFOAM surface sampling may store T either on points or cells.
    arr = data.GetPointData().GetArray("T")
    association = "point"

    if arr is None:
        arr = data.GetCellData().GetArray("T")
        association = "cell"

    if arr is None:
        raise RuntimeError(f"No T array found in {vtp_file}")

    values = vtk_to_numpy(arr)

    return {
        "max": float(values.max()),
        "min": float(values.min()),
        "mean": float(values.mean()),
        "count": len(values),
        "association": association,
    }


results = {}

for Uref, folder in cases.items():
    results[Uref] = {}

    print(f"\n===== Uref = {Uref} m/s =====")

    for distance in distances:
        file = folder / f"x{distance:03d}.vtp"

        if not file.exists():
            raise FileNotFoundError(file)

        stats = read_T(file)
        results[Uref][distance] = stats

        print(
            f"{distance:3d} m : "
            f"max T = {stats['max']:.8e}   "
            f"mean T = {stats['mean']:.8e}   "
            f"({stats['association']} data, n={stats['count']})"
        )


output = Path("figure1_downstream_peakT.csv")

with output.open("w", newline="") as f:
    writer = csv.writer(f)

    writer.writerow([
        "distance_m",
        "Tmax_Uref5",
        "Tmax_Uref10",
        "Tmax_Uref15",
    ])

    for distance in distances:
        writer.writerow([
            distance,
            results[5][distance]["max"],
            results[10][distance]["max"],
            results[15][distance]["max"],
        ])

print(f"\nSaved: {output.resolve()}")
