# Compute the Migration Phase Index (MPI)

Calculates a per-record Migration Phase Index that combines *where* a
fix sits along the population's dominant migration axis with *when* in
the seasonal cycle it occurred. MPI is a scalar in \[0, 1\] that
approaches 1 during the breeding residency at the breeding end of the
axis and approaches 0 during the wintering residency at the opposite
end. It is fully data-driven: nothing about the species' geometry,
hemisphere, or calendar is hard-coded.

## Usage

``` r
compute_mpi(
  points,
  individual_col = NULL,
  time_col = NULL,
  q = 0.1,
  anchor_q = 0.02,
  min_span_frac = 0.3,
  winter_doy = NULL,
  breed_doy = NULL,
  pc1_threshold = 0.85,
  verbose = TRUE
)
```

## Arguments

- points:

  An `sf` POINT object (typically habitat-labelled stationary fixes).

- individual_col:

  Column identifying individuals. Spatial normalisation is per
  individual; temporal anchors are pooled across individuals. If `NULL`,
  all rows are treated as one individual.

- time_col:

  Column giving fix times (POSIXct/POSIXt/Date). Required.

- q:

  Tail fraction defining which fixes count as "at the breeding /
  wintering extreme" for pooling occupancy dates (the population
  temporal anchors and polarity test). Default 0.10. This does **not**
  set the spatial normalisation anchors – see `anchor_q`.

- anchor_q:

  Near-maximal tail fraction for the per-individual spatial
  normalisation anchors (default 0.02; `p_norm = 0` at the `anchor_q`
  quantile, `1` at the `1 - anchor_q` quantile). Smaller = anchors
  closer to the true extreme (captures intermittent breeding but more
  sensitive to outlier fixes).

- min_span_frac:

  An individual contributes to the population temporal anchors only if
  its span along the axis is at least this fraction of the population
  span (default 0.30). Resident / barely moving individuals are thus
  excluded from anchor voting but still receive an MPI.

- winter_doy, breed_doy:

  Optional day-of-year anchors (1–365). If supplied, the corresponding
  date is used directly for the temporal tent **and** to orient the
  spatial polarity (which axis extreme is the breeding end is set to the
  extreme whose occupancy dates sit nearest `breed_doy` / farthest from
  `winter_doy`), so the spatial and temporal terms stay consistent. This
  overrides the automatic temporal-concentration polarity test – use it
  whenever the printed polarity disagrees with the species' phenology.

- pc1_threshold:

  Variance fraction above which movement is treated as unidirectional
  (default 0.85).

- verbose:

  Logical; print the derived axis, anchors and polarity.

## Value

The input `sf` with an added numeric `mpi` column in \[0, 1\], plus
attributes `mpi_winter_doy`, `mpi_breed_doy` (day-of-year) and
`mpi_pc1_var`.

## Spatial term (per individual)

The migration axis is found by PCA on `(lon, lat)`. If PC1 explains at
least `pc1_threshold` of the variance the movement is treated as
unidirectional (primary axis only); otherwise PC2 is also used
(mixed-axis migration). Each individual is normalised against *its own*
axis extremes, so partial migrants and individuals with different
geographic ranges are placed on a common 0–1 scale. This is the "space =
individual" anchor.

The extremes used for normalisation are the **near-maximal** quantiles
`anchor_q` and `1 - anchor_q` (default 0.02, i.e. the 2nd and 98th
percentiles), *not* a moderate tail like the 10th/90th. The breeding
extreme is a *spatial* limit that can be badly under-represented in
*time*: an individual that skips breeding in some years spends most of
its northern time at a lower, non-breeding settlement, which would drag
a 90th-percentile anchor down and clamp the true breeding extreme to 1 –
making a non-breeding stop indistinguishable from real breeding.
Near-maximal anchors reserve `p_norm = 1` for the genuine extreme. The
separate, more moderate `q` tail is used only to collect the fixes "at
the extreme" whose dates set the population temporal anchors below.

## Temporal anchors (population level)

Which end of the axis is the *breeding* end is decided by **temporal
concentration**, not by latitude or photoperiod (both have too many
counterexamples across migration geometries). For each axis extreme, the
day-of-year of every fix near that extreme is pooled across the
population and its circular concentration `R` (mean resultant length) is
computed. Breeding sites are occupied in a tight seasonal window,
wintering sites over a longer, more diffuse period, so the extreme with
the **higher `R`** is the breeding end. A poleward cross-check is used
only to break near-ties. The breeding and wintering day-of-year anchors
(`breed_doy`, `winter_doy`) are the circular means of the pooled fix
dates at each end; they may be overridden directly. These are the "time
= population" anchors, derived from the same fix-sets that define the
spatial extremes, so space and time stay consistent.

## Temporal term (triangular tent)

The seasonal term `d_norm` is a triangular function of day-of-year,
anchored at the two population day-of-year points: it is 0 at the
wintering anchor, rises linearly to 1 at the breeding anchor, then falls
linearly back to 0. The two arcs (spring up, autumn down) fold onto the
same value, so MPI is a direction-agnostic breeding-phase index – spring
and autumn migration receive the same score by design.

## MPI

Single-axis: \\\mathrm{MPI} = \sqrt{p\_{\mathrm{norm}}\\
d\_{\mathrm{norm}}}\\. Mixed-axis: \\\mathrm{MPI} =
\sqrt\[3\]{p\_{\mathrm{norm}}\\ s\_{\mathrm{norm}}\\
d\_{\mathrm{norm}}}\\.

## See also

[`detect_dominant_axis()`](https://evertsz.github.io/movepp/reference/detect_dominant_axis.md)
for the underlying axis detection.

## Examples

``` r
if (FALSE) { # \dontrun{
hab_mpi <- compute_mpi(hab, individual_col = "Individual",
                       time_col = "Time")
hist(hab_mpi$mpi)
} # }
```
