#' movepp: Movement Point-Pattern Analysis for Animal Tracking Data
#'
#' A hierarchical spatiotemporal point-pattern framework for analyzing
#' animal tracking data. The package complements the companion ArcGIS Pro
#' toolbox by providing an open-source, scriptable, and fully reproducible
#' alternative.
#'
#' The workflow is organised across three nested analytical scales:
#'
#' **Macro scale — migration / stationary segmentation.** A device-agnostic
#' movement rate ([compute_step_speed()]) feeds the Behaviorally-Anchored
#' Local Moran ([balm_segmentation()]), which labels each fix migratory
#' (`HH`), stationary (`LL`), or transitional (`HL`/`LH`) by its
#' Moran-scatterplot quadrant, anchored to a data-driven flight-onset
#' threshold.
#'
#' **Meso scale — temporary habitat delineation and phase classification.**
#' Stationary points are clustered into spatially bounded temporary
#' habitats with DBSCAN at a data-derived scale
#' ([detect_habitat_params()], [dbscan_habitats()]); each habitat is then
#' assigned to a wintering, stopover, or breeding phase via the Migration
#' Phase Index ([compute_mpi()]) and density-valley phase classification
#' ([classify_phases()]).
#'
#' **Micro scale — fine behavioural patterns.** Modules for nest detection
#' (stratified spatiotemporal DBSCAN, [detect_nests()]), inter-annual site
#' fidelity ([compute_isfi()]), and circadian rhythm classification from
#' solar elevation ([classify_circadian()]) operate on the habitat
#' assignments produced upstream.
#'
#' Use `vignette(package = "movepp")` to browse the worked examples that
#' reproduce each figure of the accompanying paper.
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
