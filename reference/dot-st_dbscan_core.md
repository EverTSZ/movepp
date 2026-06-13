# Internal: Core ST-DBSCAN clustering algorithm

Implements the spatiotemporal DBSCAN of Birant & Kut (2007). Not
exported; called by
[`detect_nests()`](https://evertsz.github.io/movepp/reference/detect_nests.md).

## Usage

``` r
.st_dbscan_core(coords, times, eps_spatial, eps_temporal, min_pts)
```

## Arguments

- coords:

  Numeric matrix or 2-column data frame of (x, y).

- times:

  POSIXct/Date vector aligned to `coords` rows.

- eps_spatial:

  Numeric; spatial neighbourhood radius (CRS units).

- eps_temporal:

  Numeric; temporal neighbourhood radius (seconds).

- min_pts:

  Integer; minimum points per cluster.

## Value

A list with `cluster` (integer, 0 = noise, 1..k = clusters) and
`n_clusters`.
