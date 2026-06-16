# movepp

**movepp** is an R package implementing a hierarchical spatiotemporal
point-pattern framework for analyzing animal tracking data. It
complements the companion ArcGIS Pro toolbox by providing an
open-source, scriptable, and fully reproducible alternative.

The package operationalises the methods described in:

> Xiao H., Peng H., Zhang Z., et al. (in revision). *From Movement to
> Meaning: Spatial Statistics Uncover Hidden Patterns in Avian Tracking
> Data.* Movement Ecology.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("EverTSZ/movepp")
```

## Quick start

``` r

library(movepp)
data(godwit_demo)

# 1. Device-agnostic step speed (centered) from raw GPS fixes
g <- compute_step_speed(godwit_demo,
                        time_col       = "time",
                        individual_col = "individual",
                        direction      = "centered")

# 2. Movement-state segmentation with BALM (behaviorally-anchored Local Moran)
seg <- balm_segmentation(g,
                         variable_col   = "step_speed",
                         individual_col = "individual")

# 3. Keep stationary (LL) and transient-stopover (LH) points
stationary <- seg[!is.na(seg$cluster_type) &
                  seg$cluster_type %in% c("LL", "LH"), ]

# 4. Derive the DBSCAN parameters from the data (never hard-coded)
p <- detect_habitat_params(stationary, track = g,
                           individual_col = "individual", time_col = "time")

# 5. Delineate temporary habitats with DBSCAN at the derived scale
hab <- dbscan_habitats(stationary,
                       individual_col = "individual",
                       eps = p$eps, minPts = p$minPts)

# 6. Migration Phase Index + density-valley phase classification
detect_dominant_axis(hab)              # PCA-based migration-axis detection
mpi <- compute_mpi(hab, time_col = "time")
phs <- classify_phases(
  mpi, n_phases = 3, fixed_stopover = "LH",
  phase_names = c("Wintering", "staging_stopover", "Breeding"))
table(phs$phase)

# 7. (optional) interactive manual phase annotation (Shiny gadget)
# phs <- annotate_phases(hab, time_col = "time")
```

## Workflow overview

`movepp` reads tracking data in the **Eulerian frame**: geographic
position is the analytical primitive and time enters as a *mark* on the
resulting spatial point pattern. The workflow runs as a single pipeline
— from raw fixes, to movement states, to functional places, and finally
reading the time-marks back onto those places — with each stage walked
through in a worked article on the [package
website](https://evertsz.github.io/movepp/articles/):

| Stage | Functions | What it does |
|----|----|----|
| **From fixes to movement states** | [`compute_step_speed()`](https://evertsz.github.io/movepp/reference/compute_step_speed.md), [`balm_segmentation()`](https://evertsz.github.io/movepp/reference/balm_segmentation.md) | Derive a device-agnostic movement rate, then label each fix migratory (`HH`), stationary (`LL`), or transitional (`HL`/`LH`) with BALM |
| **From movement states to functional places** | [`detect_habitat_params()`](https://evertsz.github.io/movepp/reference/detect_habitat_params.md), [`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md), [`detect_dominant_axis()`](https://evertsz.github.io/movepp/reference/detect_dominant_axis.md), [`compute_mpi()`](https://evertsz.github.io/movepp/reference/compute_mpi.md), [`classify_phases()`](https://evertsz.github.io/movepp/reference/classify_phases.md), [`annotate_phases()`](https://evertsz.github.io/movepp/reference/annotate_phases.md) | Delineate temporary habitats from the stationary (`LL`) fixes with DBSCAN, then assign each habitat a migration phase (wintering / stopover / breeding) via the Migration Phase Index |
| **Reading marks back onto places** | [`compute_solar_elevation()`](https://evertsz.github.io/movepp/reference/compute_solar_elevation.md), [`classify_circadian()`](https://evertsz.github.io/movepp/reference/classify_circadian.md), [`classify_habitat_circadian()`](https://evertsz.github.io/movepp/reference/classify_habitat_circadian.md), [`detect_nests()`](https://evertsz.github.io/movepp/reference/detect_nests.md), [`compute_isfi()`](https://evertsz.github.io/movepp/reference/compute_isfi.md) | With places labelled by phase, read the time-of-day, day-of-year, and year marks back at nested timescales: diel activity type, ST-DBSCAN nest detection, and the multi-year Individual Site Fidelity Index |

## Two kinds of stopover

BALM resolves two distinct kinds of pause. Multi-day **staging
stopovers** accumulate enough co-located fixes to form a low-speed
(`LL`) cluster and surface as the central stopover mode. Brief
in-corridor touch-downs instead appear as `LH` outliers (a low-speed fix
embedded in the high-speed migratory corridor). Setting
`fixed_stopover = "LH"` in
[`classify_phases()`](https://evertsz.github.io/movepp/reference/classify_phases.md)
hard-labels these as `transient_stopover` and removes them from the
classification, so the two stopover types stay distinct and the
transient points do not bias the fitted phase components.

## Relation to the ArcGIS Pro toolbox

`movepp` and the **Animal Movement Point-Pattern Analysis** toolbox for
ArcGIS Pro are two independent implementations of the same framework:
the toolbox runs in ArcGIS Pro’s Python environment
(`scikit-learn`/`scipy`), while `movepp` is a pure-R implementation that
needs no GIS licence. Either one alone is enough to apply the framework,
and they share the same pipeline end to end — BALM movement-state
segmentation, DBSCAN temporary habitats, density-valley MPI phase
classification, circadian classification, nest detection, and site
fidelity. The shared spatial-statistics steps are algorithmically
equivalent (matching statistically rather than bit-for-bit); the one
deliberate difference is nest detection, where the toolbox wraps
Picardi’s `nestR` (in R) while `movepp` uses its own stratified
ST-DBSCAN
([`detect_nests()`](https://evertsz.github.io/movepp/reference/detect_nests.md)).

Both implementations classify movement states with **BALM**
([`balm_segmentation()`](https://evertsz.github.io/movepp/reference/balm_segmentation.md)
in R), a movement-adapted variant of the Local Moran scatterplot. BALM
anchors the high/low reference to a data-driven behavioral threshold –
the flight-onset speed – rather than the arithmetic mean, and labels
each fix deterministically by its Moran-scatterplot quadrant rather than
by a permutation significance test; spatial support for stationary sites
is instead arbitrated by density clustering
([`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md)).
Because it uses only point locations and a movement-rate mark – never
direction or the raw temporal sequence – BALM is invariant to fix
ordering and robust to sampling-rate heterogeneity across devices,
years, and species.

## Citation

If you use `movepp` in your work, please cite both the package and the
associated paper:

``` r

citation("movepp")
```

## Licence

MIT (c) 2026 movepp authors
