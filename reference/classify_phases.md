# Classify Habitats into Migration Phases by Peak-Anchored Density Cuts

Assigns each temporary habitat to a migration phase (typically
*Wintering*, *Stopover*, *Breeding*) from the per-point Migration Phase
Index (MPI). The phase boundaries are **density valleys anchored to the
extreme modes** of the MPI distribution; each point is cut into a phase
deterministically and each habitat then takes the **majority** phase of
its points. No distributional model is fitted and no prominence
threshold is needed.

## Usage

``` r
classify_phases(
  habitat_points,
  mpi_col = "mpi",
  individual_col = "individual",
  cluster_col = "cluster_id",
  n_phases = 3L,
  prominence = 0,
  bw_adjust = 1,
  phase_names = NULL,
  fixed_stopover = NULL,
  cluster_type_col = "cluster_type",
  fixed_label = "transient_stopover",
  verbose = TRUE
)
```

## Arguments

- habitat_points:

  An `sf` object with at least three columns: the MPI, the individual
  ID, and the cluster ID (typically the output of
  [`compute_mpi()`](https://evertsz.github.io/movepp/reference/compute_mpi.md)
  applied to
  [`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md)
  output).

- mpi_col:

  Column name of MPI values (default `"mpi"`).

- individual_col:

  Column name of the individual ID (default `"individual"`).

- cluster_col:

  Column name of the temporary habitat cluster ID (default
  `"cluster_id"`).

- n_phases:

  Integer; the target number of biological phases (default 3 for
  wintering / stopover / breeding, 2 for species without a distinct
  stopover). The resolved number may be smaller if the data structurally
  lack the modes (see *Structural phase-count reduction*).

- prominence:

  Optional valley-prominence floor in (0, 1), measured as the fractional
  density drop from the lower flanking peak,
  `(lower_peak - valley) / lower_peak`. Used **only** to ignore
  near-flat wiggle valleys before anchoring; default `0` (off). It does
  not decide the number of phases – the peak anchoring does.

- bw_adjust:

  Numeric; bandwidth multiplier passed to
  [`stats::density()`](https://rdrr.io/r/stats/density.html) for the MPI
  KDE (default 1; larger = smoother, fewer spurious wiggle
  peaks/valleys).

- phase_names:

  Character vector of length equal to the *resolved* number of phases.
  Defaults: `c("Wintering", "Stopover", "Breeding")` for three,
  `c("Wintering", "Breeding")` for two. A mismatched vector is replaced
  by the defaults with a message.

- fixed_stopover:

  Character vector of `cluster_type` codes (e.g. `"LH"`) whose points
  are hard-assigned to `fixed_label` and excluded from the cut and the
  vote. `NULL` (default) disables this.

- cluster_type_col:

  Column name holding the Local Moran's I cluster type (default
  `"cluster_type"`); required only when `fixed_stopover` is non-NULL.

- fixed_label:

  Phase label assigned to the fixed points (default
  `"transient_stopover"`).

- verbose:

  Logical; print diagnostic information.

## Value

The input `sf` with three new columns:

- `phase_label`: integer phase rank (`1..K`) from the cut (`NA` for
  fixed and NA-MPI points).

- `phase_voted_label`: integer phase rank, the per-habitat majority
  vote.

- `phase`: factor of phase labels (with an extra `fixed_label` level
  when `fixed_stopover` is used).

An attribute `"phase_fit"` records the resolved phase count, the density
peaks and valleys, the chosen cuts, the per-phase point counts and mean
MPI, and the settings used.

## Peak-anchored cutting

MPI is constructed so that wintering piles up near 0 and breeding near

1.  The **leftmost** density peak of the pooled MPI is therefore the
    wintering mode and the **rightmost** is the breeding mode – a
    structural fact of the index, not a tuned choice. The phase
    boundaries follow directly:

    - the **lower cut** is the first density valley to the *right* of
      the leftmost (wintering) peak;

    - the **upper cut** is the last density valley to the *left* of the
      rightmost (breeding) peak.

    Everything between the two cuts is the stopover band. Because the
    cuts are anchored to the extreme peaks, any number of staging sites
    in the middle – however deep their internal valleys – are absorbed
    into the stopover phase and can never be mistaken for a phase
    boundary. This needs **no prominence threshold**: the only
    ingredient is the kernel density estimate (standard data-driven
    bandwidth).

## Assignment

Each point is labelled deterministically by the band its MPI falls in
(at or below the lower cut -\> wintering; above the upper cut -\>
breeding; in between -\> stopover), and each habitat is then assigned
its **majority** phase by voting over its points.

## Structural phase-count reduction

The resolved number of phases is `length(cuts) + 1`. It is smaller than
the requested `n_phases` only when the data structurally lack the modes
to support it: two peaks (wintering and breeding, with no staging mode
between) give a single valley and two phases; a single mode gives no
valley and one phase. This is an objective property of the MPI
distribution, not a threshold decision, so the classifier never invents
a phase the data do not contain.

## Transient vs staging stopovers (`fixed_stopover`)

Local Moran's I resolves two distinct kinds of pause. Brief in-corridor
touch-downs appear as "LH" outliers (a low-speed fix embedded in the
high-speed migratory corridor), whereas multi-day staging sites
accumulate enough co-located fixes to form an "LL" cluster and surface
as a middle-MPI mode. Setting `fixed_stopover = "LH"` hard-assigns the
LH points to `fixed_label` (default `"transient_stopover"`) and removes
them from both the cut and the per-habitat vote, so the phases are
resolved only on the genuinely stationary (LL) points.

## See also

[`compute_mpi()`](https://evertsz.github.io/movepp/reference/compute_mpi.md)
for the upstream MPI computation,
[`annotate_phases()`](https://evertsz.github.io/movepp/reference/annotate_phases.md)
for interactive manual refinement.

## Examples

``` r
if (FALSE) { # \dontrun{
demo <- make_demo_track(n_individuals = 2, n_points = 300)
seg  <- balm_segmentation(demo, "speed", "individual", verbose = FALSE)
stat <- seg[!is.na(seg$cluster_type) &
            seg$cluster_type %in% c("LL", "LH"), ]
hab  <- dbscan_habitats(stat, individual_col = "individual",
                         eps = 2, minPts = 12)
mpi  <- compute_mpi(hab, time_col = "time", verbose = FALSE)
phs  <- classify_phases(mpi, n_phases = 3, fixed_stopover = "LH")
table(phs$phase)
} # }
```
