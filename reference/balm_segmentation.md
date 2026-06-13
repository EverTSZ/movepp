# Behaviorally-Anchored Local Moran (BALM) movement-state segmentation

Classifies each tracking fix into a movement state – migration,
stationary, or transitional – with a behaviorally-anchored variant of
the Local Moran scatterplot, using only point locations and a
movement-rate mark (no direction, no raw temporal sequence).

## Usage

``` r
balm_segmentation(
  points,
  variable_col,
  individual_col = NULL,
  reference = "auto",
  k = 8L,
  kde_adjust = 1,
  verbose = TRUE
)
```

## Arguments

- points:

  An `sf` POINT object (geographic CRS recommended).

- variable_col:

  Column name (string) of the movement-rate attribute, e.g.
  `"step_speed"` from
  [`compute_step_speed()`](https://evertsz.github.io/movepp/reference/compute_step_speed.md).

- individual_col:

  Column name (string) of the individual ID, or `NULL` (default) to
  treat all points as one group. The behavioral reference is pooled
  across all individuals; the neighbourhood is computed per group. Pool
  one species at a time.

- reference:

  Either `"auto"` (default; species-level flight onset estimated from
  all individuals pooled) or a single numeric used as a fixed reference
  for all groups.

- k:

  Integer; nearest neighbours for the spatial weights (default 8;
  inverse-distance, row-standardised).

- kde_adjust:

  Numeric; passed to the flight-onset estimator.

- verbose:

  Logical; print the pooled reference and per-group counts.

## Value

The input `sf` with `balm_deviation`, `balm_lag`, and a `cluster_type`
factor (levels `HH`, `LL`, `HL`, `LH`; NA where `variable_col` is
missing). The reference is attached as attribute `"balm_reference"`.

## Method

BALM adapts the Moran scatterplot (own value vs spatial lag) in two ways
suited to movement data:

1.  **Behavioral anchoring.** The high/low reference is not the
    arithmetic mean (standard Local Moran centering) but a data-driven
    behavioral threshold `c`: the onset of the flight mode. It is
    estimated once from the pooled (log) movement-rate of ALL
    individuals by Gaussian-mixture decomposition (component count
    chosen by BIC); `c` is the lower 5\\ yields a stable species-level
    threshold, applied to every individual, while the spatial
    neighbourhood stays per individual.

2.  **Deterministic classification.** Each fix is labelled by the signs
    of its deviation `d_i = x_i - c` and its spatial lag
    `Sum_j w_ij d_j`; no permutation significance is applied. Spatial
    support is arbitrated downstream by density clustering
    ([`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md)).

The four states follow the Moran-scatterplot quadrants:

- `HH`: fast fix among fast neighbours -\> migration.

- `LL`: slow fix among slow neighbours -\> stationary site.

- `HL`: fast fix among slow neighbours -\> local movement within a site
  (e.g. commuting or foraging burst).

- `LH`: slow fix among fast neighbours -\> brief touch-down in the
  migratory corridor (transient stopover).

Because it never uses direction or fix ordering, BALM is invariant to
the temporal sequence and robust to sampling-rate heterogeneity.

## See also

[`compute_step_speed()`](https://evertsz.github.io/movepp/reference/compute_step_speed.md)
for the input rate;
[`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md)
for downstream habitat delineation.
