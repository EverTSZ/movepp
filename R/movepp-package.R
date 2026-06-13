#' movepp: Movement Point-Pattern Analysis for Animal Tracking Data
#'
#' A hierarchical spatiotemporal point-pattern framework for analyzing
#' animal tracking data. The package complements the companion ArcGIS Pro
#' toolbox by providing an open-source, scriptable, and fully reproducible
#' alternative.
#'
#' The package reads tracking data in the **Eulerian frame**: geographic
#' position is the analytical primitive and time enters as a *mark* on the
#' resulting spatial point pattern. The workflow runs as a single
#' pipeline:
#'
#' **From fixes to movement states.** A device-agnostic movement rate
#' ([compute_step_speed()]) feeds the Behaviorally-Anchored Local Moran
#' ([balm_segmentation()]), which labels each fix migratory (`HH`),
#' stationary (`LL`), or transitional (`HL`/`LH`) by its Moran-scatterplot
#' quadrant, anchored to a data-driven flight-onset threshold.
#'
#' **From movement states to functional places.** The stationary (`LL`)
#' fixes are clustered into spatially bounded temporary habitats with
#' DBSCAN at a data-derived scale ([detect_habitat_params()],
#' [dbscan_habitats()]); each habitat is then assigned to a wintering,
#' stopover, or breeding phase via the Migration Phase Index
#' ([compute_mpi()], [detect_dominant_axis()]) and density-valley phase
#' classification ([classify_phases()]), with [annotate_phases()] as an
#' interactive behaviour-barcode correction channel.
#'
#' **Reading marks back onto places.** With habitats labelled by phase,
#' the time-of-day, day-of-year, and year marks are read back at nested
#' timescales: diel-activity type from solar elevation
#' ([compute_solar_elevation()], [classify_circadian()],
#' [classify_habitat_circadian()]), nest detection via stratified
#' spatiotemporal DBSCAN ([detect_nests()]), and the multi-year Individual
#' Site Fidelity Index ([compute_isfi()]).
#'
#' Worked articles reproducing each figure of the accompanying paper are
#' on the package website
#' (\url{https://evertsz.github.io/movepp/articles/}).
#'
#' @keywords internal
"_PACKAGE"

#' Return movepp version
#'
#' Convenience function returning the installed package version.
#'
#' @return A `package_version` object.
#' @export
#' @examples
#' movepp_version()
movepp_version <- function() {
  utils::packageVersion("movepp")
}
