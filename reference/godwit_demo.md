# Black-tailed Godwit Tracking Demo Data

GPS tracking data from one Black-tailed Godwit (*Limosa limosa*)
individual, covering two full migration cycles between the wintering
grounds in coastal Southeast Asia and the breeding grounds in Far East
Siberia. Used as a worked example in the package vignettes and as a
reproducibility check against the companion manuscript:

## Usage

``` r
godwit_demo
```

## Format

An `sf` POINT object in EPSG:4326 (WGS 84) with columns:

- individual:

  character; individual ID (`"Black-tailed Godwit 23_02"`).

- time:

  POSIXct; UTC timestamp of each GPS fix.

- lon:

  numeric; longitude (degrees).

- lat:

  numeric; latitude (degrees).

- speed:

  numeric; instantaneous speed (m/s as reported by the tracking device).

- geometry:

  sf POINT geometry.

## Source

Subset of tracking data from the authors\\ Black-tailed Godwit project.
The full nine-individual dataset is described in the companion
manuscript.

## Details

Xiao, H., Peng, H., Zhang, Z., et al. (in revision). *From Movement to
Meaning: Spatial Statistics Uncover Hidden Patterns in Tracking Data.*
Movement Ecology.
