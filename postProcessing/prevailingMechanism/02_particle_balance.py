from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# PATHS
# ============================================================

BASE = Path(__file__).resolve().parent

INPUT_FILE = BASE / "tables" / "particle_export_inventory.csv"
OUTPUT_CSV = BASE / "tables" / "particle_mass_balance.csv"
OUTPUT_FIGURE = BASE / "figures" / "particle_mass_balance.png"


# ============================================================
# RELEASE SETTINGS
# ============================================================

# Numerical tracer weighting, not a measured emission rate.
TOTAL_RELEASE_MASS_KG = 1.0e-6
RELEASE_DURATION_S = 60.0


# ============================================================
# READ VALIDATED INVENTORY
# ============================================================

df = pd.read_csv(INPUT_FILE).sort_values("relative_time_s")

required = {
    "relative_time_s",
    "parcel_count",
    "represented_mass_kg",
}

missing = required.difference(df.columns)

if missing:
    raise ValueError(f"Missing columns: {sorted(missing)}")


# ============================================================
# MASS BALANCE
# ============================================================

elapsed = np.clip(
    df["relative_time_s"].to_numpy(),
    0.0,
    RELEASE_DURATION_S,
)

injected_mass = (
    TOTAL_RELEASE_MASS_KG
    * elapsed
    / RELEASE_DURATION_S
)

remaining_mass = df["represented_mass_kg"].to_numpy()

escaped_mass = injected_mass - remaining_mass

if escaped_mass.min() < -1.0e-10:
    raise ValueError(
        "Remaining parcel mass exceeds injected mass."
    )

escaped_mass = np.maximum(escaped_mass, 0.0)

remaining_fraction = remaining_mass / injected_mass
escaped_fraction = escaped_mass / injected_mass


# ============================================================
# SAVE TABLE
# ============================================================

output = pd.DataFrame({
    "relative_time_s": elapsed,
    "parcel_count_remaining": df["parcel_count"],
    "injected_mass_kg": injected_mass,
    "remaining_mass_kg": remaining_mass,
    "inferred_escaped_mass_kg": escaped_mass,
    "remaining_fraction": remaining_fraction,
    "escaped_fraction": escaped_fraction,
})

output.to_csv(OUTPUT_CSV, index=False)


# ============================================================
# PLOT
# ============================================================

plt.figure(figsize=(7.16, 4.2), dpi=300)

plt.plot(
    elapsed,
    injected_mass * 1.0e6,
    label="Injected",
    linewidth=1.8,
)

plt.plot(
    elapsed,
    remaining_mass * 1.0e6,
    label="Remaining in domain",
    linewidth=1.8,
)

plt.plot(
    elapsed,
    escaped_mass * 1.0e6,
    label="Inferred escaped",
    linewidth=1.8,
)

plt.xlabel("Time after release began (s)")
plt.ylabel("Represented tracer mass (mg)")

plt.xlim(0, RELEASE_DURATION_S)
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
# FINAL WRITTEN-TIME SUMMARY
# ============================================================

last = output.iloc[-1]

print("PARTICLE MASS BALANCE")
print("=" * 45)
print(f"Written relative time : {last['relative_time_s']:.4f} s")
print(f"Parcels remaining     : {int(last['parcel_count_remaining']):,}")
print(f"Injected mass         : {last['injected_mass_kg']:.8e} kg")
print(f"Remaining mass        : {last['remaining_mass_kg']:.8e} kg")
print(f"Inferred escaped mass : {last['inferred_escaped_mass_kg']:.8e} kg")
print(f"Remaining fraction    : {100*last['remaining_fraction']:.2f}%")
print(f"Escaped fraction      : {100*last['escaped_fraction']:.2f}%")
print()
print(f"Saved table  : {OUTPUT_CSV}")
print(f"Saved figure : {OUTPUT_FIGURE}")
