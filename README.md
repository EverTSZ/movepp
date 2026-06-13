# movepp

<!-- badges: start -->
[![R-CMD-check](https://github.com/EverTSZ/movepp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EverTSZ/movepp/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**movepp** is an R package implementing a hierarchical spatiotemporal
point-pattern framework for analyzing animal tracking data. It complements
the companion ArcGIS Pro toolbox by providing an open-source, scriptable,
and fully reproducible alternative.

The package operationalises the methods described in:

> Xiao H., Peng H., Zhang Z., et al. (in revision). *From Movement to
> Meaning: Spatial Statistics Uncover Hidden Patterns in Tracking Data.*
> Movement Ecology.

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

`movepp` decomposes movement analysis into three nested scales, each
backed by a Quarto vignette (`browseVignettes("movepp")`):

| Scale | Function family | Purpose |
|-------|-----------------|---------|
| **Macro** | `compute_step_speed()`, `balm_segmentation()` | Derive a device-agnostic movement rate, then split migratory from stationary fixes with BALM (a behaviorally-anchored Local Moran) |
| **Meso** | `detect_habitat_params()`, `dbscan_habitats()`, `detect_dominant_axis()`, `compute_mpi()`, `classify_phases()`, `annotate_phases()` | Delineate temporary habitats and assign wintering / stopover / breeding phases |
| **Micro** | `detect_nests()`, `compute_isfi()`, `classify_circadian()` | Detect nesting events, quantify site fidelity, and resolve circadian rhythm |


## Two kinds of stopover

BALM resolves two distinct kinds of pause. Multi-day **staging stopovers**
accumulate enough co-located fixes to form a low-speed (`LL`) cluster and
surface as the central stopover mode. Brief in-corridor
touch-downs instead appear as `LH` outliers (a low-speed fix embedded in
the high-speed migratory corridor). Setting `fixed_stopover = "LH"` in
`classify_phases()` hard-labels these as `transient_stopover` and
removes them from the classification, so the two stopover types stay distinct and
the transient points do not bias the fitted phase components.

## Relation to the ArcGIS Pro toolbox

This R package is a sister implementation to the **Analysis Toolbox 1.2**
distributed for ArcGIS Pro. Where the toolbox classified movement states
with Anselin's Local Moran's I (the ArcGIS *Cluster and Outlier Analysis*
tool), `movepp` uses **BALM** (`balm_segmentation()`), a movement-adapted
variant of the Local Moran scatterplot. BALM anchors the high/low reference
to a data-driven behavioral threshold -- the flight-onset speed -- rather
than the arithmetic mean, and labels each fix deterministically by its
Moran-scatterplot quadrant rather than by a permutation significance test;
spatial support for stationary sites is instead arbitrated by density
clustering (`dbscan_habitats()`). Because it uses only point locations and
a movement-rate mark -- never direction or the raw temporal sequence --
BALM is invariant to fix ordering and robust to sampling-rate heterogeneity
across devices, years, and species.

## Citation

If you use `movepp` in your work, please cite both the package and the
associated paper:

``` r
citation("movepp")
```

## Licence

MIT (c) 2026 movepp authors
