# Delineate temporary habitats by DBSCAN at a behaviourally-derived scale

Clusters stationary tracking fixes into temporary habitats using
standard DBSCAN (Ester et al. 1996), applied independently to each
individual. Unlike density-hierarchy methods (HDBSCAN), both parameters
here are fixed and interpretable rather than auto-extracted:

## Usage

``` r
dbscan_habitats(x, individual_col = NULL, eps, minPts, drop_noise = TRUE)
```

## Arguments

- x:

  An sf POINT object, typically the stationary (LL/LH) subset from
  [`balm_segmentation()`](https://EverTSZ.github.io/movepp/reference/balm_segmentation.md).

- individual_col:

  Name of the column identifying individuals. Clustering is performed
  independently within each individual. If `NULL`, all rows are treated
  as a single individual.

- eps:

  Neighbourhood radius in kilometres. **Required, no default** – it must
  be derived from the data, typically via
  [`detect_habitat_params()`](https://EverTSZ.github.io/movepp/reference/detect_habitat_params.md),
  not hard-coded. Two fixes join the same habitat when reachable through
  hops no larger than `eps`.

- minPts:

  Minimum number of points within `eps` for a core point. **Required, no
  default** – typically from
  [`detect_habitat_params()`](https://EverTSZ.github.io/movepp/reference/detect_habitat_params.md).
  Acts as a uniform density / noise floor; interpretable as a minimum
  residence time given the sampling interval.

- drop_noise:

  Logical; if `TRUE` (default) noise points (DBSCAN cluster `0`) are
  dropped from the returned object.

## Value

The input sf object with an added integer column `cluster_id` giving the
per-individual habitat id (`1..k`; `0` = noise, retained only when
`drop_noise = FALSE`).

## Details

- `eps` is a fixed, biologically meaningful spatial scale – the
  species-typical maximum within-habitat movement, obtained from the
  population step-length distribution (the upper mode of local,
  non-flight steps). Occupancy at this scale is self-similar
  (scale-free), so a single explicit scale is more defensible than a
  brittle auto-selected optimum.

- `minPts` is a unified noise / minimum-duration floor – the knee of the
  cluster-count vs `minPts` curve, which for the Black-tailed Godwit
  data coincides with a minimum residence of ~3 days at a 6-hour
  sampling interval (`minPts = 12`).

Time is deliberately excluded: habitats are defined purely from spatial
occupancy, leaving temporal use to be analysed downstream.

## References

Ester, M., Kriegel, H.-P., Sander, J., & Xu, X. (1996). A density-based
algorithm for discovering clusters in large spatial databases with
noise. *Proc. KDD-96*, 226-231.

## See also

[`detect_habitat_params()`](https://EverTSZ.github.io/movepp/reference/detect_habitat_params.md)
to derive `eps` and `minPts` from data.

## Examples

``` r
if (FALSE) { # \dontrun{
stationary <- seg[seg$cluster_type %in% c("LL", "LH"), ]
## 1. derive thresholds from the data (never hard-code them)
p <- detect_habitat_params(stationary, individual_col = "Individual",
                           time_col = "Time")
## 2. run transparent per-individual DBSCAN with the derived thresholds
hab <- dbscan_habitats(stationary, individual_col = "Individual",
                       eps = p$eps, minPts = p$minPts)
tapply(hab$cluster_id, hab$Individual, function(z) length(unique(z)))
} # }
```
