# Compute step speed from successive GPS fixes

For each point, computes a step-based speed using the great-circle (or
Euclidean) distance to the adjacent fix(es) in time, divided by the time
interval. Operates per individual when an individual identifier is
supplied.

## Usage

``` r
compute_step_speed(
  points,
  time_col = "time",
  individual_col = NULL,
  unit = c("km/h", "m/s", "km/day"),
  direction = c("centered", "backward", "forward"),
  max_gap_hours = Inf
)
```

## Arguments

- points:

  An `sf` POINT object.

- time_col:

  Column name (string) for timestamp; must be POSIXct or Date.

- individual_col:

  Column name (string) for individual identifier, or `NULL` to treat all
  points as one trajectory.

- unit:

  Output unit for `step_speed`: `"km/h"` (default), `"m/s"`, or
  `"km/day"`.

- direction:

  Direction of the step calculation: `"centered"` (default),
  `"backward"`, or `"forward"`. See *Direction modes*.

- max_gap_hours:

  Numeric; if the time gap to the relevant adjacent fix exceeds this
  threshold, that side is treated as missing (long gaps usually indicate
  tracker dropouts, not actual movement). Default `Inf` (no filtering).
  When `direction = "centered"` the filter is applied to each side
  independently; the centered estimate then degenerates to whichever
  side remains valid (or NA if both sides are filtered).

## Value

The input `sf` object with three new columns:

- `step_distance_km`: distance used in the calculation (km). For
  `"centered"`, the sum of distances to both adjacent fixes.

- `step_dt_hours`: time interval used in the calculation (hours). For
  `"centered"`, the sum of intervals to both adjacent fixes.

- `step_speed`: `step_distance_km / step_dt_hours` in the chosen unit.

## Why use this instead of tracker-reported speed

Tracker-reported speed varies between devices: some report Doppler
instantaneous speed, others report displacement-over-interval, others
use proprietary smoothing. Many low-cost trackers do not report speed at
all. Computing step speed from raw fix coordinates gives a reproducible,
device-agnostic variable directly comparable across datasets and
consistent with the spatial scale of the analysis.

## Direction modes

The `direction` argument controls how each fix's speed is defined:

- `"centered"` (default): the total distance to both adjacent fixes
  divided by the total time, giving a window-averaged speed around the
  focal fix. Symmetrically captures both take-off and landing events,
  and produces no edge NAs (the first and last fix degenerate to forward
  and backward respectively).

- `"backward"`: distance from the previous fix to the focal fix, divided
  by their time interval. Standard in movement ecology (e.g.
  adehabitatLT, momentuHMM) but biases detection toward landing events
  (the first fix of each individual is NA).

- `"forward"`: distance from the focal fix to the next fix, divided by
  their time interval. Biases detection toward take-off events (the last
  fix of each individual is NA).

## See also

[`balm_segmentation()`](https://EverTSZ.github.io/movepp/reference/balm_segmentation.md)
for using the computed step speed in spatial autocorrelation analysis.

## Examples

``` r
if (FALSE) { # \dontrun{
data(godwit_demo)
godwit2 <- compute_step_speed(godwit_demo,
                              time_col       = "time",
                              individual_col = "individual",
                              unit           = "km/h",
                              direction      = "centered")
summary(godwit2$step_speed)

# Feed into Local Moran's I
seg <- balm_segmentation(godwit2,
                          variable_col   = "step_speed",
                          individual_col = "individual")
} # }
```
