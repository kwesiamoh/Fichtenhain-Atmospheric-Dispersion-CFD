# Numerical Investigation of Atmospheric Pollutant Dispersion Using OpenFOAM

This project investigates atmospheric contaminant transport around a terrain- and building-resolved industrial site using OpenFOAM.

The study was developed around the ECO PORT area in Holzwickede, Germany. Terrain and LoD2 building data were used to create the surrounding environment, while a conceptual process plant and elevated release stack were added to provide a controlled source for the dispersion study.

The main aim was to examine how changes in wind speed and wind direction affect plume transport, concentration, and spread. The work combines Eulerian passive-scalar calculations with Lagrangian particle tracking so that the plume can be examined from two different viewpoints.

The project also includes the geometry-processing, CFD setup, post-processing, and visualization workflow used to move from geographic data to the final results.

---

## Project questions

The study was built around a few practical questions:

- How does reference wind speed affect peak contaminant concentration downstream?
- How does wind speed affect lateral and vertical plume spread?
- How does changing wind direction alter the route taken by the plume through the site?
- Do virtual receptors placed around the source capture the expected downwind response?
- How far do released tracer particles travel during the simulated release?
- What fraction of the represented particle mass remains inside the computational domain?

The simulations are intended to examine transport behaviour under controlled meteorological scenarios.

---

## Study area

The computational domain represents the terrain and surrounding LoD2 building environment of the ECO PORT industrial area in Holzwickede, Germany.

A conceptual process plant was placed within this environment. The model includes the main plant structures required to create meaningful flow obstruction around the release location, including tanks, buildings, pipework and an elevated stack.

The site geometry was prepared using:

- QGIS for terrain and geospatial processing
- Blender for plant modelling and geometry preparation
- OpenFOAM for meshing and CFD
- ParaView for CFD visualization
- Python and R for quantitative post-processing

### Site model

<img width="479" height="477" alt="image" src="https://github.com/user-attachments/assets/1b4b20af-4656-4b12-ae9e-6d29140c951d" />


---

## Modelling workflow

The project was developed in stages.

### 1. Terrain and building preparation

Terrain and LoD2 building data were prepared in QGIS and transferred into the geometry workflow.

The purpose of including the surrounding built environment was to retain the major terrain and building effects that can influence the atmospheric flow approaching and leaving the plant.

### 2. Plant geometry

The conceptual plant was developed in Blender and positioned within the site model.

The Blender model also serves as the presentation environment for the final CFD visualization and animation.

### 3. Computational mesh

The OpenFOAM mesh was developed using `blockMesh` and `snappyHexMesh`.

The terrain, buildings and plant geometry were incorporated into the CFD domain before the atmospheric-flow cases were developed.

Large mesh and geometry files are not included in this repository. The relevant OpenFOAM configuration dictionaries are included as configuration snapshots.

### 4. Atmospheric boundary layer

A neutral atmospheric boundary layer was used for the airflow calculations.

The workflow included an initial developed flow field followed by transient LES calculations for the prevailing-wind case.

### 5. Eulerian plume calculations

A passive scalar, `T`, was used to represent normalized gaseous contaminant concentration.

The scalar cases were used mainly for:

- wind-speed comparison
- plume-width analysis
- wind-direction comparison
- virtual receptor sampling

`T` is dimensionless and is not a calibrated concentration in mg/m³.

### 6. Lagrangian particle calculations

A separate particle case was used to examine the movement of computational tracer parcels released from the stack.

This provided information on:

- downwind transport
- particle residence within the domain
- escape from the domain
- particle-cloud position
- particle-cloud spread

The particle analysis is used as a transport study. It is not a complete aerosol deposition model.

---

## Quality control and rejected run

One of the early transient LES runs reached its requested final timestep, but post-processing showed local velocities above 1,200 m/s.

That result was not accepted.

The problem was traced back to the inlet treatment. The simulation was returned to a clean carrier-flow checkpoint and repeated with a revised atmospheric boundary condition.

The corrected LES solution remained below approximately 30 m/s and was used as the basis for the corrected particle calculation.

The rejected result is retained only as part of the project development history and is not used in the analysis presented here.

---

## Meteorological baseline

The prevailing-wind scenario was based on historical observations from the nearby DWD Dortmund-Wickede station.

The selected prevailing direction was approximately:

**154° meteorological wind direction**

This corresponds approximately to flow from the south-southeast toward the north-northwest.

The prevailing case was used for the detailed Lagrangian transport analysis.

---

## Simulation cases

### Wind-speed comparison

Three reference wind speeds were compared while keeping the general dispersion setup consistent.

| Case | Reference wind speed |
|---|---:|
| Uref5 | 5 m/s |
| Uref10 | 10 m/s |
| Uref15 | 15 m/s |

These cases were used to compare:

- peak normalized concentration
- lateral plume spread
- vertical plume spread
- plume-centroid displacement

---

### Wind-direction comparison

Four wind directions were examined.

| Case | Wind from | Main flow direction |
|---|---:|---|
| WD0 | 0° | North to south |
| WD90 | 90° | East to west |
| WD180 | 180° | South to north |
| Prevailing | 154° | SSE to NNW |

The direction names follow the meteorological convention: the angle describes where the wind comes from.

The prevailing-wind CFD scenario was selected from hourly wind observations obtained from the Deutscher Wetterdienst (DWD) Climate Data Center. DWD station 01078 (Düsseldorf) was used because it provided a long and highly complete hourly wind record in the historical CDC wind dataset used for this study
<img width="800" height="800" alt="Figure_1" src="https://github.com/user-attachments/assets/2e5b2c13-d9ca-4e9f-9517-2eb4c9a41246" />



---

### Lagrangian particle case

The prevailing-wind particle case used an elevated release from the stack.

Main particle settings included:

| Parameter | Value |
|---|---:|
| Release location | approximately (816.4, 1486.0, 70.0) m |
| Release duration | 60 s |
| Parcel release rate | 3,000 parcels/s |
| Represented total tracer mass | 1 × 10⁻⁶ kg |
| Nominal particle diameter | 1 × 10⁻⁵ m |
| Gravity | Disabled |
| Deposition | Disabled |
| Stochastic dispersion | Disabled |

The computational parcels therefore represent kinematic tracer transport rather than a fully resolved atmospheric aerosol model.

---

# Results

## 1. Effect of wind speed on peak concentration

Peak normalized gaseous contaminant concentration decreased with downstream distance for all three wind-speed cases.

The higher-wind-speed cases also produced lower peak values at a fixed downstream location.

<img width="1954" height="1234" alt="Figure_1" src="https://github.com/user-attachments/assets/07d71d72-d74e-4597-8e30-6eaa5a5f130c" />
**[Figure: Peak normalized concentration versus downstream distance]**

### Main observation

Increasing the reference wind speed from 5 to 15 m/s reduced the peak normalized concentration along the plume pathway.

The faster flow transports the scalar downstream more quickly and changes the time available for local turbulent mixing before the plume reaches a given downstream distance.

---

## 2. Lateral and vertical plume spread

The plume became wider with downstream distance in both the lateral and vertical directions.

<img width="1954" height="1474" alt="Figure_2" src="https://github.com/user-attachments/assets/f6ecae38-ca6d-4bec-af08-1bcc4a1e7a53" />
**[Figure: Lateral and vertical plume spread]**

At 500 m downstream:

| Reference wind speed | Lateral spread, σy | Vertical spread, σz |
|---:|---:|---:|
| 5 m/s | 13.91 m | 13.95 m |
| 10 m/s | 11.27 m | 11.34 m |
| 15 m/s | 10.25 m | 10.35 m |

At the same downstream distance, the 5 m/s case produced the widest plume.

This does not mean that lower wind speed always produces greater atmospheric dispersion in every situation. In these simulations, the lower-speed plume takes longer to reach the same downstream position, giving the transported scalar more time to spread before arriving there.

---

## 3. Plume-centroid displacement

The lateral and vertical positions of the plume centroid were also tracked.

<img width="4266" height="2396" alt="Figure_3" src="https://github.com/user-attachments/assets/815f1f51-6ba0-4d3e-be4d-85c3666f9cfd" />


The centroid trajectory changed less strongly with wind speed than the concentration and plume-width metrics.

This result is treated as supporting information rather than one of the main conclusions of the study.

---

## 4. Effect of wind direction

Changing wind direction changed the route followed by the plume through the site.

<img width="1920" height="1080" alt="4in1" src="https://github.com/user-attachments/assets/aeb80a05-3a52-4922-aac5-5ebc452cfe12" />

The four cases show the plume travelling approximately:

- south for wind from 0°
- west for wind from 90°
- north for wind from 180°
- north-northwest for the prevailing wind from 154°

The release-height plots are useful for comparing plume routing, but they should not be interpreted as ground-level exposure maps.

---

## 5. Virtual receptor analysis

Eight virtual receptors were positioned on a 200 m radius around the stack:

**N, NE, E, SE, S, SW, W and NW**

All receptors were sampled at the release-height plane, `z = 70 m`.

<img width="792" height="533" alt="Rplot01" src="https://github.com/user-attachments/assets/81ae4bc5-a8d3-4c97-9b92-de5330cedcbd" />

The receptor with the highest normalized concentration was located in the expected downwind direction for each case.

| Wind scenario | Peak receptor | Maximum T |
|---|---|---:|
| Wind from 0° | S | 1.624 × 10⁻⁴ |
| Wind from 90° | W | 1.709 × 10⁻⁴ |
| Wind from 180° | N | 1.381 × 10⁻⁴ |
| Prevailing wind from 154° | NW | 6.331 × 10⁻⁶ |

The receptor analysis provides a numerical check of the plume pathways seen in the CFD fields.

These receptors are computational sampling locations. They are not physical monitoring stations.

---

## 6. Lagrangian particle transport

The corrected particle case released 176,913 computational parcels during the simulated release.

<img width="2118" height="1230" alt="particle_mass_balance" src="https://github.com/user-attachments/assets/2a749954-9bc8-4a14-84db-90fa911ad255" />

By the end of the analysed period:

- approximately 45.11% of the represented tracer mass remained in the domain
- approximately 54.89% had left the domain

This shows the transition from initial accumulation near the source to progressive transport out of the computational domain.

---

## 7. Downwind particle transport

The particle cloud was also analysed in terms of its downwind travel distance.

<img width="2118" height="1230" alt="downwind_transport" src="https://github.com/user-attachments/assets/45a744b1-06dd-48f2-adef-4bb7164a1d77" />

At the final analysed state:

- mass-weighted mean distance: approximately 199 m
- median distance: approximately 203 m
- 95th-percentile distance: approximately 385 m
- leading parcels: slightly above 400 m

The difference between the centre of the cloud and the leading particles reflects the range of velocities sampled by particles moving through the resolved flow field.

---

## 8. Particle-cloud geometry

Additional post-processing was used to examine particle-cloud height and crosswind/vertical spread.

<img width="2118" height="1769" alt="plume_rise_and_spread" src="https://github.com/user-attachments/assets/6ec72af1-2602-434e-9cca-b3ccbfedf971" />

These results are retained as supporting analysis because the main project conclusions are already captured more clearly by the concentration, plume-width, wind-direction and downwind-transport results.

---

# Blender visualization

The CFD results are also being integrated into the existing Blender site model to produce a presentation-quality animation of the tracer release.

The Blender scene contains the terrain, surrounding buildings and conceptual plant used during the geometry-development stage of the project.

The purpose of this visualization is to present the CFD results in their site context. Blender is used for presentation and animation rather than for calculating the dispersion solution.

## Planned animation

The final animation will include:

- the terrain- and building-resolved site
- the conceptual process plant
- the release stack
- the prevailing wind direction
- the corrected Lagrangian tracer sequence
- particle-age colouring
- camera movement through the site

<img width="1920" height="1080" alt="first8" src="https://github.com/user-attachments/assets/db11016a-6124-4da6-a629-4f4a3741c707" />
<img width="1920" height="1080" alt="first3" src="https://github.com/user-attachments/assets/67f84871-04fb-4533-b79c-3f3dc7f5ce6c" />


---

# What the results show

Taken together, the simulations show two different aspects of atmospheric transport.

Wind speed affects the magnitude and spatial spread of the plume. In the cases tested here, increasing the reference wind speed lowered the peak normalized concentration and reduced lateral and vertical spread at the same downstream distance.

Wind direction mainly controls the route taken by the plume through the site. The virtual receptor analysis supported the directional behaviour seen in the concentration fields.

The Lagrangian calculations add a second view of the transport process by showing how the particle population moves through and eventually leaves the computational domain.

The Eulerian and Lagrangian calculations should not be treated as identical measurements. They are complementary ways of examining the same general transport problem.

---

# Repository structure

```text
Fichtenhain-Atmospheric-Dispersion-CFD/
│
├── cases/
│   ├── baseline/
│   │   └── Fichtenhain_LES_Prevailing/
│   │
│   ├── particles/
│   │   └── Fichtenhain_DPM_Particles_Prevailing/
│   │
│   ├── windSpeed/
│   │   ├── Fichtenhain_Plume_Uref5/
│   │   ├── Fichtenhain_Plume_Uref10/
│   │   └── Fichtenhain_Plume_Uref15/
│   │
│   └── windDirection/
│       ├── Fichtenhain_Plume_WD0/
│       ├── Fichtenhain_Plume_WD90/
│       ├── Fichtenhain_Plume_WD180/
│       └── Fichtenhain_Plume_WD154_Prevailing/
│
├── postProcessing/
│   ├── prevailingMechanism/
│   ├── windSpeed/
│   └── windDirection/
│
├── results/
│   ├── figures/
│   └── tables/
│
├── docs/
│   └── images/
│
└── README.md
