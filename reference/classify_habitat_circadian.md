# Type each temporary habitat by its dominant diel period

Aggregates the per-fix circadian labels from
[`classify_circadian()`](https://evertsz.github.io/movepp/reference/classify_circadian.md)
to the temporary-habitat (cluster) level. Each habitat is assigned the
diel period holding a majority (greater than `majority`) of its fixes,
or `"Mixed"` when no period does. The per-fix day/twilight/night label
is only an intermediate; the habitat-level functional type – a daytime
foraging site, a night roost, a twilight transit area – is the
ecological object this returns.

## Usage

``` r
classify_habitat_circadian(
  points,
  circadian_col = "circadian",
  cluster_col = "cluster_id",
  individual_col = NULL,
  majority = 0.5,
  periods = c("Day", "Twilight", "Night"),
  verbose = TRUE
)
```

## Arguments

- points:

  An `sf` object carrying a per-fix circadian label (typically the
  `circadian` column from
  [`classify_circadian()`](https://evertsz.github.io/movepp/reference/classify_circadian.md))
  and a temporary-habitat cluster id (the `cluster_id` from
  [`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md)).

- circadian_col:

  Per-fix diel-label column (default `"circadian"`).

- cluster_col:

  Temporary-habitat cluster-id column (default `"cluster_id"`).

- individual_col:

  Individual-id column, or `NULL` (default) to treat `cluster_col` as
  globally unique. When supplied, habitats are typed within each
  `(individual, cluster)` pair.

- majority:

  Numeric in \[0, 1\]; the minimum fraction the dominant period must
  reach for the habitat to take its label rather than `"Mixed"` (default
  0.5).

- periods:

  Character vector of diel levels to tabulate and rank (default
  `c("Day", "Twilight", "Night")`).

- verbose:

  Logical; print the resulting habitat-type counts.

## Value

The input `sf` with an added factor column `habitat_circadian` (levels
`periods` plus `"Mixed"`), constant within each habitat. A
one-row-per-habitat summary – cluster id, fix count, per-period
proportions, and assigned type – is attached as the attribute
`"habitat_circadian"`.

## See also

[`classify_circadian()`](https://evertsz.github.io/movepp/reference/classify_circadian.md)
for the per-fix labels this aggregates.

## Examples

``` r
if (FALSE) { # \dontrun{
pts <- compute_solar_elevation(habitat_pts, time_col = "time")
pts <- classify_circadian(pts)
pts <- classify_habitat_circadian(pts, individual_col = "individual")
attr(pts, "habitat_circadian")
} # }
```
