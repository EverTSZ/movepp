# Compute Solar Elevation Angle at Each Tracking Point

Calculates the solar elevation angle (in degrees, with negative values
below the horizon) for each tracking point, given its coordinates and
timestamp. Used downstream by
[`classify_circadian()`](https://EverTSZ.github.io/movepp/reference/classify_circadian.md)
to partition activity into day, twilight, and night periods.

## Usage

``` r
compute_solar_elevation(points, time_col, verbose = TRUE)
```

## Arguments

- points:

  An `sf` POINT object.

- time_col:

  Column name (string) of a POSIXct/POSIXt time column. Must include
  time zone information for correct local-time computation.

- verbose:

  Logical; print progress.

## Value

The input `sf` with a new `solar_elevation` numeric column in degrees
(positive = above horizon, negative = below).

## References

Whitener, B. (2025). *Pysolar: Python Library for Solar Position
Calculations*. The R equivalent computation is performed by the
`suncalc` package (Thieurmel & Elmarhraoui).

## Examples

``` r
if (FALSE) { # \dontrun{
demo <- make_demo_track(n_individuals = 1, n_points = 100)
demo <- compute_solar_elevation(demo, time_col = "time")
hist(demo$solar_elevation)
} # }
```
