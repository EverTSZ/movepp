# Classify Tracking Points into Day / Twilight / Night

Partitions tracking points by solar elevation angle into three circadian
categories: **Day** (sun above horizon), **Twilight** (civil twilight,
sun between the two thresholds), and **Night** (sun well below horizon).
The default thresholds follow the standard definition of civil twilight
(0 deg and -6 deg).

## Usage

``` r
classify_circadian(
  points,
  elevation_col = "solar_elevation",
  day_threshold = 0,
  twilight_threshold = -6,
  verbose = TRUE
)
```

## Arguments

- points:

  An `sf` object with a solar elevation column (typically the output of
  [`compute_solar_elevation()`](https://evertsz.github.io/movepp/reference/compute_solar_elevation.md)).

- elevation_col:

  Column name (string) of solar elevation in degrees (default
  `"solar_elevation"`).

- day_threshold:

  Numeric; angle (deg) above which a point is classified as Day (default
  0).

- twilight_threshold:

  Numeric; angle (deg) above which a point is classified as Twilight
  (default -6, the civil-twilight bound).

- verbose:

  Logical; print summary.

## Value

The input `sf` with a new `circadian` factor column (levels: `"Day"`,
`"Twilight"`, `"Night"`).

## Examples

``` r
if (FALSE) { # \dontrun{
demo <- compute_solar_elevation(demo, time_col = "time")
demo <- classify_circadian(demo)
table(demo$circadian)
} # }
```
