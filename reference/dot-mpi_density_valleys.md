# Density peaks and valleys of a 1-D sample

Internal helper for
[`classify_phases()`](https://evertsz.github.io/movepp/reference/classify_phases.md).
Computes a kernel density estimate and returns its interior local maxima
(peaks) and minima (valleys), plus each valley's prominence – the
**fractional** density drop from the lower of its two nearest flanking
peaks down to the valley floor, `(lower_peak - valley) / lower_peak`.
The prominence is returned for optional ranking/flooring only; valley
selection in
[`classify_phases()`](https://evertsz.github.io/movepp/reference/classify_phases.md)
is anchored to the extreme peaks and needs no prominence threshold.

## Usage

``` r
.mpi_density_valleys(v, bw_adjust = 1)
```

## Arguments

- v:

  Numeric vector.

- bw_adjust:

  Bandwidth multiplier for
  [`stats::density()`](https://rdrr.io/r/stats/density.html).

## Value

A list with `dd` (the density), `peak_x` (peak locations, ascending),
`valley_x` (valley locations, ascending) and `valley_prom` (their
fractional prominences).
