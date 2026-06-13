# Derive DBSCAN habitat-delineation parameters from the data

Estimates the two DBSCAN parameters (`eps` and `minPts`) used by
[`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md)
from empirical tracking data, removing the need for hard-coded
assumptions. Both parameters are derived at the population level and
have an explicit behavioural meaning. This generalized version
dynamically adapts to whether behavioural state tracking data is present
or absent.

## Usage

``` r
detect_habitat_params(
  x,
  track = NULL,
  individual_col = NULL,
  time_col = NULL,
  state_col = "cluster_type",
  stationary_states = c("LL"),
  flight_states = c("HH"),
  eps = NULL,
  eps_quantile = 0.95,
  minpts_max = 60L,
  min_residence = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  An `sf` POINT object containing the stationary subset of points to be
  clustered (the data over which `minPts` is scanned).

- track:

  An `sf` POINT object containing the complete telemetry trajectory.
  Used to extract step lengths, calculate sampling intervals, and
  analyze movement states. Required unless a custom `eps` value is
  explicitly provided without specifying `min_residence`.

- individual_col:

  Character string identifying the column for unique individual IDs in
  `x` and `track`. If `NULL`, all rows are treated as a single group.

- time_col:

  Character string specifying the timestamp column in `track`.

- state_col:

  Movement-state label column in `track` (default `"cluster_type"`, the
  BALM output). If this column is **absent** from `track`, state
  filtering is bypassed and the function falls back to Generalized Mode
  (all steps pooled) – so raw, un-segmented tracks from non-migratory
  species work without change. Set `state_col = NULL` to force the
  generalized mode even when a state column exists.

- stationary_states:

  Vector of state labels treated as local/stationary for the `eps` step
  set (default `c("LL")`). If `NULL`, state-based filtering is bypassed.

- flight_states:

  Vector of state labels treated as flight/relocation (default
  `c("HH")`), used only for the regime diagnostic (`gap_ratio`).

- eps:

  Optional numeric value. If supplied, the empirical derivation of `eps`
  (and the tracking data requirements tied to it) is skipped.

- eps_quantile:

  Numeric percentile (0 to 1) extracted from the within-site complex to
  define `eps` (default `0.95`).

- minpts_max:

  Integer establishing the upper limit for the `minPts` knee scan
  (default `60`).

- min_residence:

  Optional minimum patch residency duration expressed in **days**. If
  supplied, `minPts` bypasses the knee calculation and directly computes
  fix thresholds.

- verbose:

  Logical; if `TRUE`, outputs the derived values and diagnostic logs to
  the console (default `TRUE`).

## Value

A named list containing:

- `eps`: Calculated or user-supplied neighborhood distance threshold
  (km).

- `minPts`: Calculated or converted minimum point density threshold
  (integer).

- `regime`: Diagnostic label indicating behavior pooling (`"migratory"`,
  `"continuous"`, `"all-data-pooled"`, or `"user-supplied"`).

- `gap_ratio`: The ratio of median flight step lengths to `eps`
  (diagnostic asset).

- `n_components`: Number of GMM components selected by BIC during
  clustering.

- `sampling_interval_h`: Median telemetry sampling interval expressed in
  hours.

- `min_residence_days`: The minimum habitat residence duration implied
  by `minPts`.

## eps – the within-site movement scale

`eps` represents the radius linking telemetry fixes that belong to the
same temporary habitat patch. It can be derived in two modes:

- **State-Filtered Mode (default):** When `track` carries the
  `state_col` column (e.g. BALM `cluster_type`) and `stationary_states`
  is set, consecutive great-circle steps whose *both* endpoints are
  stationary (e.g. BALM "LL") are isolated. This cleanly strips away
  high-speed flights and between-site transits before fitting.

- **Generalized Mode (automatic fallback):** If `state_col` is absent
  from `track` (or no `stationary_states` are given), all valid positive
  step lengths across the population are pooled instead. This lets
  non-migratory species be analysed directly on the raw track, without a
  prior BALM segmentation step.

The function defaults to State-Filtered Mode but falls back to
Generalized Mode automatically (with a message) when the state column is
not present, so the same call works for both BALM-segmented and raw
tracks.

In both modes, a Gaussian Mixture Model (GMM) is fitted to the
log-transformed step lengths. The mixture is partitioned at the
**largest gap** between consecutive component means, splitting the
distribution into a lower "within-site complex" (resting, foraging,
local transits) and higher relocation modes. `eps` is set as the
`eps_quantile` (default 0.95) percentile of this within-site complex.
Pooling across the population ensures statistical stability and prevents
the overfitting (spurious components) common in per-individual fits.

Note: in Generalized Mode the within-site/relocation split relies
entirely on the GMM largest-gap heuristic. This is robust for sedentary
species, but for central-place foragers (residents that commute between
a roost and feeding sites) the commute steps can inflate `eps`; consider
a light stationarity filter, or BALM, in that case.

## minPts – the minimum-residence floor

`minPts` defines the minimum number of local fixes a cluster must
contain to be recognized as a valid habitat. By default, it is
determined via the Kneedle geometry algorithm, which captures the "knee"
of the aggregate habitat-count vs. `minPts` curve at the derived `eps`.
This curve marks the physical transition from culling accidental noise
points to eroding genuine habitat clusters.

Its operational meaning is temporal:
`minPts * sampling_interval = minimum residence duration`.
Alternatively, users can explicitly supply `min_residence` (in days) to
enforce a direct behavioral floor, which is converted to fix counts
using the median sampling interval.

## See also

[`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: BALM-segmented track -- State-Filtered Mode is used by default
p1 <- detect_habitat_params(stat_sf, track = seg,
                            individual_col = "Individual", time_col = "Time")

# Example 2: Raw track from a non-migratory species (no state column) --
# the function falls back to Generalized Mode automatically
p2 <- detect_habitat_params(raw_sf, track = raw_sf,
                            individual_col = "Individual", time_col = "Time")
} # }
```
