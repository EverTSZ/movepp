# Compute Individual Site Fidelity Index (ISFI)

Quantifies inter-annual site fidelity per individual: the probability
that a given temporary habitat (cluster) revisited across multiple
years. Defined as \$\$\mathrm{ISFI} = \frac{\left(\sum_i
n\_{i,\mathrm{years}}\right) - n}{(c - 1) \cdot n}\$\$ where \\n\\ is
the number of unique clusters, \\c\\ the number of unique years, and
\\n\_{i,\mathrm{years}}\\ the count of distinct years in which cluster
\\i\\ was visited. ISFI ranges from 0 (each cluster visited in only one
year) to 1 (every cluster revisited in every year).

## Usage

``` r
compute_isfi(
  habitat_points,
  time_col,
  individual_col = "individual",
  cluster_col = "cluster_id",
  phase_col = NULL,
  min_years = 2L,
  ref_years = NULL,
  verbose = TRUE
)
```

## Arguments

- habitat_points:

  An `sf` object with `individual_col`, `cluster_col`, and `time_col`
  columns (typically the output of
  [`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md)
  or
  [`classify_phases()`](https://evertsz.github.io/movepp/reference/classify_phases.md)).

- time_col:

  Column name (string) of a POSIXct/Date time column. The year is
  extracted from this column.

- individual_col:

  Column name (string) of the individual ID (default `"individual"`).

- cluster_col:

  Column name (string) of the temporary habitat cluster ID (default
  `"cluster_id"`).

- phase_col:

  Optional column name for phase labels. If supplied, ISFI is computed
  separately within each (individual, phase) pair.

- min_years:

  Integer; minimum number of unique years required for a meaningful ISFI
  computation (default 2). Individuals with fewer years return `NA` ISFI
  with a warning.

- ref_years:

  Optional. A single integer, or a named numeric vector keyed by
  individual ID, giving the number of observation years to use as the
  denominator \\c\\ instead of the year span of `habitat_points`. Use
  when the cluster set is a subset of a longer record (e.g. nest-area
  fidelity within the breeding phase), so that years in which the subset
  was absent still count against fidelity. Default `NULL` (derive \\c\\
  from the data; unchanged behaviour).

- verbose:

  Logical; print progress.

## Value

A data frame with one row per individual (or per individual x phase if
`phase_col` supplied), columns:

- `individual`: individual ID.

- `phase`: phase name (only when `phase_col` supplied).

- `n_clusters`: unique clusters visited.

- `n_years`: unique years tracked.

- `sum_cluster_years`: total cluster-year occupancies (\\\sum_i
  n\_{i,\mathrm{years}}\\).

- `isfi`: site fidelity index in 0..1, `NA` if undefined.

## Use with [`detect_nests()`](https://evertsz.github.io/movepp/reference/detect_nests.md) output to derive INFI

To compute the Individual Nest Fidelity Index (INFI) introduced in the
manuscript Discussion, apply this same function to the `nests` sf
produced by
[`detect_nests()`](https://evertsz.github.io/movepp/reference/detect_nests.md)
across multi-year tracking, passing the nest table's `nest_id` as the
cluster column and the nest start time as the time column. No separate
function is needed.

## Use with phase labels for stratified analysis

Pass `phase_col` to compute ISFI separately within each phase (e.g.,
breeding, stopover, wintering). This reproduces the stratified analysis
of the original manuscript (Figure 7).

## Examples

``` r
if (FALSE) { # \dontrun{
# ISFI per individual across all habitats
fidelity <- compute_isfi(habitats, time_col = "time")

# Stratified by phase
fidelity_phase <- compute_isfi(phases, time_col = "time",
                               phase_col = "phase")
} # }
```
