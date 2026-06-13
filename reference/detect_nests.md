# Detect Bird Nesting Events with Stratified Spatiotemporal DBSCAN

Identifies nesting locations and hatching periods from high-resolution
tracking data using a three-tier stratified ST-DBSCAN scheme. The
detection rests on the assumption that a true nest yields GPS fixes
approximately bivariate-normally distributed around the actual nest
location: 95% of points within \\2\sigma\\, 68% within \\1\sigma\\, and
~38% within \\0.5\sigma\\ of the centre. A point is classified as a nest
point only if it survives the density threshold at all three sigma tiers
simultaneously, yielding markedly fewer false positives than single-tier
DBSCAN.

## Usage

``` r
detect_nests(
  points,
  time_col,
  individual_col = NULL,
  min_pts = 40L,
  search_dist = 41,
  time_window = "3 weeks",
  sigma_ratios = c(1, 0.68, 0.38),
  verbose = TRUE
)
```

## Arguments

- points:

  An `sf` POINT object with high-resolution GPS fixes during the
  breeding season.

- time_col:

  Column name (string) of a POSIXct/Date time column.

- individual_col:

  Column name (string) of the individual ID (default `NULL` = pool all
  points).

- min_pts:

  Integer; baseline minimum points per cluster at the 2\\\sigma\\ tier
  (default 40).

- search_dist:

  Numeric; baseline spatial radius at the 2\\\sigma\\ tier, in CRS
  linear units (default 41, typical for 2 x 20m GPS error).

- time_window:

  Character or `difftime`; the temporal window used by ST-DBSCAN
  (default `"3 weeks"`, suitable for a single nesting attempt).

- sigma_ratios:

  Length-3 numeric vector of (min_pts, dist) multipliers for the
  (2\\\sigma\\, 1\\\sigma\\, 0.5\\\sigma\\) tiers. Default
  `c(1.0, 0.68, 0.38)` for min_pts and is paired internally with
  `c(1.0, 0.5, 0.25)` for distance.

- verbose:

  Logical; print progress messages.

## Value

A list with two `sf` elements:

- `points`: the input `sf`, augmented with `nest_id` (`NA` if the point
  is not a nest point) and `near_dist` (distance from the point to its
  nest centre).

- `nests`: one row per detected nest, with columns `nest_id`,
  `individual`, `n_points`, `start_time`, `end_time`, `mean_time`, and
  geometry (the mean centre).

## Algorithm

For each individual (or pooled dataset if `individual_col = NULL`):

1.  Run ST-DBSCAN at three nested density tiers:

    - **2\\\sigma\\**: `min_pts` x 1.0, `search_dist` x 1.0.

    - **1\\\sigma\\**: `min_pts` x 0.68, `search_dist` x 0.5.

    - **0.5\\\sigma\\**: `min_pts` x 0.38, `search_dist` x 0.25.

2.  Retain only points classified into a non-noise cluster in **all
    three tiers** (membership gate against false positives).

3.  Assign each surviving point the cluster label of the **loosest
    (2\\\sigma\\) tier**, whose large radius collapses one nest's entire
    fix cloud into a single `nest_id`. (Labelling by the strictest tier
    would split one real nest into several sub-clusters, producing many
    spurious nest records for a single nest.)

4.  Compute the mean centre, first/last timestamp, and mean timestamp of
    each detected nest.

## Coordinate units

`search_dist` is interpreted in the linear units of the input CRS. If
`points` is in WGS84 (EPSG:4326), distances are in degrees and the
function will warn. For accurate nest detection at scales of tens of
metres, project to a metric CRS first (e.g., a local UTM zone) and
supply `search_dist` in metres.

## References

Birant, D. & Kut, A. (2007). ST-DBSCAN: An algorithm for clustering
spatial-temporal data. *Data & Knowledge Engineering* 60(1): 208-221.
[doi:10.1016/j.datak.2006.01.013](https://doi.org/10.1016/j.datak.2006.01.013)

## Examples

``` r
if (FALSE) { # \dontrun{
# Pied Avocet nesting data (project to UTM first for metric distance)
nests <- detect_nests(avocet_breeding_pts,
                      time_col = "time",
                      individual_col = "individual",
                      min_pts = 40,
                      search_dist = 41,
                      time_window = "3 weeks")
nests$nests  # mean centre + time bounds of each nest
} # }
```
