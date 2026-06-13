# Package index

## From fixes to movement states

Derive a movement rate, then segment migratory vs stationary fixes.

- [`compute_step_speed()`](https://EverTSZ.github.io/movepp/reference/compute_step_speed.md)
  : Compute step speed from successive GPS fixes
- [`balm_segmentation()`](https://EverTSZ.github.io/movepp/reference/balm_segmentation.md)
  : Behaviorally-Anchored Local Moran (BALM) movement-state segmentation

## From movement states to functional places

Delineate temporary habitats and assign wintering, stopover, breeding
phases.

- [`detect_habitat_params()`](https://EverTSZ.github.io/movepp/reference/detect_habitat_params.md)
  : Derive DBSCAN habitat-delineation parameters from the data
- [`dbscan_habitats()`](https://EverTSZ.github.io/movepp/reference/dbscan_habitats.md)
  : Delineate temporary habitats by DBSCAN at a behaviourally-derived
  scale
- [`detect_dominant_axis()`](https://EverTSZ.github.io/movepp/reference/detect_dominant_axis.md)
  : Detect the Dominant Migration Axis via PCA
- [`print(`*`<movepp_dominant_axis>`*`)`](https://EverTSZ.github.io/movepp/reference/print.movepp_dominant_axis.md)
  : Print method for movepp_dominant_axis objects
- [`compute_mpi()`](https://EverTSZ.github.io/movepp/reference/compute_mpi.md)
  : Compute the Migration Phase Index (MPI)
- [`classify_phases()`](https://EverTSZ.github.io/movepp/reference/classify_phases.md)
  : Classify Habitats into Migration Phases by Peak-Anchored Density
  Cuts
- [`annotate_phases()`](https://EverTSZ.github.io/movepp/reference/annotate_phases.md)
  : Interactively annotate phases on a behaviour barcode

## Reading marks back onto places

Read time-of-day, day-of-year, and year marks back at nested timescales:
diel activity, nesting, site fidelity.

- [`detect_nests()`](https://EverTSZ.github.io/movepp/reference/detect_nests.md)
  : Detect Bird Nesting Events with Stratified Spatiotemporal DBSCAN
- [`compute_isfi()`](https://EverTSZ.github.io/movepp/reference/compute_isfi.md)
  : Compute Individual Site Fidelity Index (ISFI)
- [`compute_solar_elevation()`](https://EverTSZ.github.io/movepp/reference/compute_solar_elevation.md)
  : Compute Solar Elevation Angle at Each Tracking Point
- [`classify_circadian()`](https://EverTSZ.github.io/movepp/reference/classify_circadian.md)
  : Classify Tracking Points into Day / Twilight / Night
- [`classify_habitat_circadian()`](https://EverTSZ.github.io/movepp/reference/classify_habitat_circadian.md)
  : Type each temporary habitat by its dominant diel period

## Utilities and data

- [`make_demo_track()`](https://EverTSZ.github.io/movepp/reference/make_demo_track.md)
  : Generate a synthetic demo movement track
- [`movepp_version()`](https://EverTSZ.github.io/movepp/reference/movepp_version.md)
  : Return movepp version
- [`godwit_demo`](https://EverTSZ.github.io/movepp/reference/godwit_demo.md)
  : Black-tailed Godwit Tracking Demo Data
