# Detect the Dominant Migration Axis via PCA

Automatically determines whether a tracking dataset is dominated by
latitudinal, longitudinal, or mixed (diagonal) movement, eliminating the
need for users to know this property *a priori*. Principal Components
Analysis (PCA) is applied to the (lon, lat) coordinates; the fraction of
total variance captured by PC1 determines whether the movement is
treated as unidirectional or mixed.

## Usage

``` r
detect_dominant_axis(points, pc1_threshold = 0.85)
```

## Arguments

- points:

  An `sf` POINT object.

- pc1_threshold:

  Numeric in (0, 1\]; fraction of variance PC1 must explain to call the
  movement unidirectional (default 0.85).

## Value

A list of class `"movepp_dominant_axis"` with elements `primary_axis`,
`secondary_axis` (or `NULL`), `pc1_variance_explained`,
`pc1_loading_lon`, `pc1_loading_lat`, and the underlying `pca`
(`prcomp`) object.

## Examples

``` r
if (FALSE) { # \dontrun{
demo <- make_demo_track()
detect_dominant_axis(demo)
} # }
```
