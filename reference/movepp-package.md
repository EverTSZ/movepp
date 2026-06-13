# movepp: Movement Point-Pattern Analysis for Animal Tracking Data

A hierarchical spatiotemporal point-pattern framework for analyzing
animal tracking data. The package complements the companion ArcGIS Pro
toolbox by providing an open-source, scriptable, and fully reproducible
alternative.

## Details

The package reads tracking data in the **Eulerian frame**: geographic
position is the analytical primitive and time enters as a *mark* on the
resulting spatial point pattern. The workflow runs as a single pipeline:

**From fixes to movement states.** A device-agnostic movement rate
([`compute_step_speed()`](https://EverTSZ.github.io/movepp/reference/compute_step_speed.md))
feeds the Behaviorally-Anchored Local Moran
([`balm_segmentation()`](https://EverTSZ.github.io/movepp/reference/balm_segmentation.md)),
which labels each fix migratory (`HH`), stationary (`LL`), or
transitional (`HL`/`LH`) by its Moran-scatterplot quadrant, anchored to
a data-driven flight-onset threshold.

**From movement states to functional places.** The stationary (`LL`)
fixes are clustered into spatially bounded temporary habitats with
DBSCAN at a data-derived scale
([`detect_habitat_params()`](https://EverTSZ.github.io/movepp/reference/detect_habitat_params.md),
[`dbscan_habitats()`](https://EverTSZ.github.io/movepp/reference/dbscan_habitats.md));
each habitat is then assigned to a wintering, stopover, or breeding
phase via the Migration Phase Index
([`compute_mpi()`](https://EverTSZ.github.io/movepp/reference/compute_mpi.md),
[`detect_dominant_axis()`](https://EverTSZ.github.io/movepp/reference/detect_dominant_axis.md))
and density-valley phase classification
([`classify_phases()`](https://EverTSZ.github.io/movepp/reference/classify_phases.md)),
with
[`annotate_phases()`](https://EverTSZ.github.io/movepp/reference/annotate_phases.md)
as an interactive behaviour-barcode correction channel.

**Reading marks back onto places.** With habitats labelled by phase, the
time-of-day, day-of-year, and year marks are read back at nested
timescales: diel-activity type from solar elevation
([`compute_solar_elevation()`](https://EverTSZ.github.io/movepp/reference/compute_solar_elevation.md),
[`classify_circadian()`](https://EverTSZ.github.io/movepp/reference/classify_circadian.md),
[`classify_habitat_circadian()`](https://EverTSZ.github.io/movepp/reference/classify_habitat_circadian.md)),
nest detection via stratified spatiotemporal DBSCAN
([`detect_nests()`](https://EverTSZ.github.io/movepp/reference/detect_nests.md)),
and the multi-year Individual Site Fidelity Index
([`compute_isfi()`](https://EverTSZ.github.io/movepp/reference/compute_isfi.md)).

Worked articles reproducing each figure of the accompanying paper are on
the package website (<https://EverTSZ.github.io/movepp/articles/>).

## See also

Useful links:

- <https://github.com/EverTSZ/movepp>

- <https://EverTSZ.github.io/movepp>

- Report bugs at <https://github.com/EverTSZ/movepp/issues>

## Author

**Maintainer**: Hengjun Xiao <evertsz@live.com>

Authors:

- Hengjun Xiao <evertsz@live.com>

- Tomohiro Ichinose <tomohiro@sfc.keio.ac.jp> \[thesis advisor\]

- Hebo Peng

- Zhengwang Zhang

- De Chen

- Weipan Lei

- Yang Wu

- Bingrun Zhu

- Yoshiaki Miyamoto

- Donglai Li
