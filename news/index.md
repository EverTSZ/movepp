# Changelog

## movepp (development version)

### movepp 0.1.4

- New
  [`classify_habitat_circadian()`](https://EverTSZ.github.io/movepp/reference/classify_habitat_circadian.md)
  aggregates the per-fix diel labels from
  [`classify_circadian()`](https://EverTSZ.github.io/movepp/reference/classify_circadian.md)
  to the temporary-habitat level: each habitat is typed by the diel
  period holding a majority of its fixes (else `"Mixed"`), turning
  point-level day/twilight/night into a habitat-level functional type
  (day-foraging site, night roost, twilight transit). A per-habitat
  summary is attached as the `"habitat_circadian"` attribute.

### movepp 0.1.3

- Temporary habitat delineation now uses standard DBSCAN at a single,
  data-derived spatial scale
  ([`dbscan_habitats()`](https://EverTSZ.github.io/movepp/reference/dbscan_habitats.md)),
  replacing the HDBSCAN density-hierarchy extraction
  (`hdbscan_habitats`) and its `scan_min_cluster_size` helper (both
  removed). HDBSCAN’s stability extraction has no stable optimum for the
  self-similar, scale-free spatial occupancy typical of tracking data.
- New
  [`detect_habitat_params()`](https://EverTSZ.github.io/movepp/reference/detect_habitat_params.md)
  derives both DBSCAN parameters from the data: `eps` from the upper
  edge of the local-movement complex in the per-individual step-length
  distribution, and `minPts` from the knee of the habitat-count vs
  `minPts` curve (a noise floor readable as a minimum residence time).
  [`dbscan_habitats()`](https://EverTSZ.github.io/movepp/reference/dbscan_habitats.md)
  takes `eps`/`minPts` as required arguments (no defaults).

### movepp 0.1.0

Initial release accompanying:

> Xiao H., Peng H., Zhang Z., et al. (in revision). *From Movement to
> Meaning: Spatial Statistics Uncover Hidden Patterns in Tracking Data.*
> Movement Ecology.

#### New features

- Device-agnostic step speed from raw fixes
  ([`compute_step_speed()`](https://EverTSZ.github.io/movepp/reference/compute_step_speed.md)),
  including a symmetric “centered” mode that captures both take-off and
  landing without edge NAs.
- Behaviorally-Anchored Local Moran (BALM) movement-state segmentation
  ([`balm_segmentation()`](https://EverTSZ.github.io/movepp/reference/balm_segmentation.md)):
  classifies each fix into migration (`HH`), stationary (`LL`), or
  transition (`HL`/`LH`) states from a movement-rate mark. BALM anchors
  the high/low reference to a data-driven flight-onset threshold (rather
  than the arithmetic mean) and labels fixes deterministically by
  Moran-scatterplot quadrant (no significance test; spatial support is
  arbitrated downstream by
  [`dbscan_habitats()`](https://EverTSZ.github.io/movepp/reference/dbscan_habitats.md)).
  It uses only point locations and the rate mark – never direction or
  the raw temporal sequence – so it is robust to sampling-rate
  heterogeneity. `individual_col` is optional (`NULL` pools all points).
- Temporary habitat delineation via standard DBSCAN at a
  behaviourally-derived spatial scale
  ([`dbscan_habitats()`](https://EverTSZ.github.io/movepp/reference/dbscan_habitats.md)),
  with both parameters derived from the data by
  [`detect_habitat_params()`](https://EverTSZ.github.io/movepp/reference/detect_habitat_params.md).
- Migration Phase Index
  ([`compute_mpi()`](https://EverTSZ.github.io/movepp/reference/compute_mpi.md))
  with automatic PCA-based dominant-axis detection
  ([`detect_dominant_axis()`](https://EverTSZ.github.io/movepp/reference/detect_dominant_axis.md)).
- Density-valley phase classification
  ([`classify_phases()`](https://EverTSZ.github.io/movepp/reference/classify_phases.md)):
  cuts the MPI distribution at its prominence-filtered density valleys
  and assigns each habitat by a per-cluster majority vote, with phases
  ordered by mean MPI. The `fixed_stopover` argument hard-labels brief
  in-corridor `LH` touch-downs as `transient_stopover` and excludes them
  from the cut, keeping transient stopovers distinct from multi-day
  staging stopovers.
- Interactive manual phase annotation through a linked behaviour-barcode
  / map Shiny gadget
  ([`annotate_phases()`](https://EverTSZ.github.io/movepp/reference/annotate_phases.md)).
- Stratified spatiotemporal DBSCAN for nest detection
  ([`detect_nests()`](https://EverTSZ.github.io/movepp/reference/detect_nests.md)).
- Inter-annual site fidelity
  ([`compute_isfi()`](https://EverTSZ.github.io/movepp/reference/compute_isfi.md)).
- Solar-elevation circadian classification
  ([`classify_circadian()`](https://EverTSZ.github.io/movepp/reference/classify_circadian.md),
  [`compute_solar_elevation()`](https://EverTSZ.github.io/movepp/reference/compute_solar_elevation.md)).

#### Documentation

- Seven Quarto vignettes covering the macro -\> meso -\> micro pipeline:
  data overview, BALM movement-state segmentation, DBSCAN habitats (with
  an interactive linked behaviour barcode), MPI phase classification,
  circadian classification, inter-annual site fidelity, and
  spatiotemporal nest detection.
