# Interactively annotate phases on a behaviour barcode

Launches a Shiny gadget with a habitat-occupancy *behaviour barcode*
(time x habitat `rank_col`) linked to a clean coastline map (no country
borders). Box-select on EITHER panel; the selection is ringed in black
on both and can be assigned a life-history phase or a free-text custom
label. With `individual_col` set, a dropdown restricts the view to one
individual. Returns the input with an added `phase_manual` column.

## Usage

``` r
annotate_phases(
  points,
  time_col = "time",
  rank_col = "cluster_rank",
  individual_col = NULL,
  phase_col = NULL,
  basemap = c("light", "natural", "gray"),
  phases = c("Wintering", "Stopover", "Breeding"),
  viewer = c("browser", "pane", "dialog")
)
```

## Arguments

- points:

  An `sf` POINT object with at least `time_col` and `rank_col`.

- time_col:

  POSIXct/Date column name (string).

- rank_col:

  Discrete habitat-rank column (default `"cluster_rank"`).

- individual_col:

  Optional individual-id column; adds a selector. `NULL` (default) shows
  all points.

- phase_col:

  Optional column of starting phase labels to seed the annotation (e.g.
  `"phase"` from
  [`classify_phases()`](https://evertsz.github.io/movepp/reference/classify_phases.md));
  the gadget opens with these instead of all-"unassigned". `NULL`
  (default) starts blank.

- basemap:

  Coastline-map colour preset: `"light"` (default), `"natural"`, or
  `"gray"`. Coastlines, rivers, lakes, land and ocean are drawn; country
  borders never are.

- phases:

  Preset phase labels (default `c("Wintering","Stopover","Breeding")`).

- viewer:

  `"browser"` (default), `"pane"`, or `"dialog"`.

## Value

The input `sf` with an added factor column `phase_manual`, original row
order; `NULL` if cancelled.

## Details

If `phase_col` is supplied, the gadget **opens with those labels as the
starting classification** (e.g. the output of
[`classify_phases()`](https://evertsz.github.io/movepp/reference/classify_phases.md)),
so only the mis-assigned points need correcting rather than annotating
from scratch. Because the barcode is organised by habitat rank,
box-selecting a whole rank-row over a time span selects an entire
habitat – so corrections can be made per habitat (cluster), matching the
cluster-level classification output.

## Interaction

- With `individual_col`, pick a bird from the **Individual** dropdown
  (or `(all)`); switching clears the current selection.

- Box-select on the **barcode** (drag) or **map** (modebar box/lasso);
  selected points are ringed in black on BOTH panels.

- Click a preset phase button, or type a custom label and *Apply
  custom*; *Unassign* clears labels.

- Toggle **Colour by** between habitat rank and assigned label.

- *Done* returns the annotated data; *Cancel* aborts.

Requires shiny, miniUI and plotly (Suggests).

## See also

[`classify_phases()`](https://evertsz.github.io/movepp/reference/classify_phases.md)
for the automatic classification this can refine,
[`dbscan_habitats()`](https://evertsz.github.io/movepp/reference/dbscan_habitats.md)
for producing the habitat ranks.
