#' Compute Individual Site Fidelity Index (ISFI)
#'
#' Quantifies inter-annual site fidelity per individual: the
#' probability that a given temporary habitat (cluster) revisited
#' across multiple years. Defined as
#' \deqn{\mathrm{ISFI} = \frac{\left(\sum_i n_{i,\mathrm{years}}\right) - n}{(c - 1) \cdot n}}
#' where \eqn{n} is the number of unique clusters,
#' \eqn{c} the number of unique years, and
#' \eqn{n_{i,\mathrm{years}}} the count of distinct years in which
#' cluster \eqn{i} was visited. ISFI ranges from 0 (each cluster
#' visited in only one year) to 1 (every cluster revisited in every
#' year).
#'
#' @section Use with [detect_nests()] output to derive INFI:
#' To compute the Individual Nest Fidelity Index (INFI) introduced in
#' the manuscript Discussion, apply this same function to the
#' `nests` sf produced by [detect_nests()] across multi-year tracking,
#' passing the nest table's `nest_id` as the cluster column and the
#' nest start time as the time column. No separate function is needed.
#'
#' @section Use with phase labels for stratified analysis:
#' Pass `phase_col` to compute ISFI separately within each phase
#' (e.g., breeding, stopover, wintering). This reproduces the
#' stratified analysis of the original manuscript (Figure 7).
#'
#' @param habitat_points An `sf` object with `individual_col`,
#'   `cluster_col`, and `time_col` columns (typically the output of
#'   [dbscan_habitats()] or [classify_phases()]).
#' @param individual_col Column name (string) of the individual ID
#'   (default `"individual"`).
#' @param cluster_col Column name (string) of the temporary habitat
#'   cluster ID (default `"cluster_id"`).
#' @param time_col Column name (string) of a POSIXct/Date time column.
#'   The year is extracted from this column.
#' @param phase_col Optional column name for phase labels. If supplied,
#'   ISFI is computed separately within each (individual, phase) pair.
#' @param min_years Integer; minimum number of unique years required
#'   for a meaningful ISFI computation (default 2). Individuals with
#'   fewer years return `NA` ISFI with a warning.
#' @param ref_years Optional. A single integer, or a named numeric vector
#'   keyed by individual ID, giving the number of observation years to use
#'   as the denominator \eqn{c} instead of the year span of
#'   `habitat_points`. Use when the cluster set is a subset of a longer
#'   record (e.g. nest-area fidelity within the breeding phase), so that
#'   years in which the subset was absent still count against fidelity.
#'   Default `NULL` (derive \eqn{c} from the data; unchanged behaviour).
#' @param verbose Logical; print progress.
#'
#' @return A data frame with one row per individual (or per
#'   individual x phase if `phase_col` supplied), columns:
#'   \itemize{
#'     \item `individual`: individual ID.
#'     \item `phase`: phase name (only when `phase_col` supplied).
#'     \item `n_clusters`: unique clusters visited.
#'     \item `n_years`: unique years tracked.
#'     \item `sum_cluster_years`: total cluster-year occupancies
#'       (\eqn{\sum_i n_{i,\mathrm{years}}}).
#'     \item `isfi`: site fidelity index in 0..1, `NA` if undefined.
#'   }
#'
#' @examples
#' \dontrun{
#' # ISFI per individual across all habitats
#' fidelity <- compute_isfi(habitats, time_col = "time")
#'
#' # Stratified by phase
#' fidelity_phase <- compute_isfi(phases, time_col = "time",
#'                                phase_col = "phase")
#' }
#'
#' @export
compute_isfi <- function(habitat_points,
                         time_col,
                         individual_col = "individual",
                         cluster_col    = "cluster_id",
                         phase_col      = NULL,
                         min_years      = 2L,
                         ref_years      = NULL,
                         verbose        = TRUE) {

  if (!inherits(habitat_points, "sf")) {
    # Allow data.frame too for INFI use case (sf nest tables)
    if (!is.data.frame(habitat_points)) {
      stop("`habitat_points` must be an `sf` or `data.frame`.",
           call. = FALSE)
    }
  }

  for (col in c(time_col, individual_col, cluster_col)) {
    if (!col %in% names(habitat_points)) {
      stop("Column `", col, "` not found in input.", call. = FALSE)
    }
  }
  if (!is.null(phase_col) && !phase_col %in% names(habitat_points)) {
    stop("Column `", phase_col, "` not found in input.", call. = FALSE)
  }

  tvec <- habitat_points[[time_col]]
  if (!inherits(tvec, c("POSIXct", "POSIXt", "Date"))) {
    stop("`time_col` must be a POSIXct or Date column.", call. = FALSE)
  }

  # Drop geometry if sf
  df <- if (inherits(habitat_points, "sf")) {
    sf::st_drop_geometry(habitat_points)
  } else {
    habitat_points
  }
  df$.year <- as.integer(format(tvec, "%Y"))

  # Determine grouping keys
  if (is.null(phase_col)) {
    df$.group <- df[[individual_col]]
    group_label <- "individual"
  } else {
    df$.group <- paste(df[[individual_col]], df[[phase_col]],
                       sep = "::")
    group_label <- "individual_phase"
  }

  groups <- unique(df$.group[!is.na(df$.group)])

  # Compute ISFI for each group
  rows <- lapply(groups, function(g) {
    sub <- df[df$.group == g, , drop = FALSE]
    sub <- sub[!is.na(sub[[cluster_col]]) & !is.na(sub$.year), ,
               drop = FALSE]

    n        <- length(unique(sub[[cluster_col]]))
    c_years  <- length(unique(sub$.year))

    # Optional external denominator: measure fidelity against a fixed number
    # of observation years (e.g. the individual's total breeding-phase years)
    # rather than the year span of this (possibly subset) input. Without this,
    # restricting to a single nest-bearing cluster collapses c to that
    # cluster's own years and forces ISFI to 1.
    if (!is.null(ref_years)) {
      ind_key <- as.character(sub[[individual_col]][1])
      rc <- if (!is.null(names(ref_years))) unname(ref_years[ind_key]) else ref_years[1]
      if (length(rc) == 1L && !is.na(rc)) c_years <- as.integer(rc)
    }

    if (c_years < min_years) {
      if (verbose) {
        message(sprintf(
          "[compute_isfi] '%s' has only %d year(s); ISFI undefined.",
          g, c_years))
      }
      isfi <- NA_real_
      sum_cluster_years <- NA_integer_
    } else if (n == 0L) {
      isfi <- NA_real_
      sum_cluster_years <- NA_integer_
    } else {
      # For each cluster, how many distinct years did it appear in?
      cluster_year_counts <- tapply(
        sub$.year, sub[[cluster_col]],
        function(y) length(unique(y))
      )
      sum_cluster_years <- sum(cluster_year_counts)

      denom <- (c_years - 1L) * n
      isfi <- if (denom == 0) NA_real_
              else (sum_cluster_years - n) / denom
    }

    out <- data.frame(
      individual         = sub[[individual_col]][1],
      n_clusters         = n,
      n_years            = c_years,
      sum_cluster_years  = sum_cluster_years,
      isfi               = isfi,
      stringsAsFactors   = FALSE
    )
    if (!is.null(phase_col)) {
      out$phase <- sub[[phase_col]][1]
      out <- out[, c("individual", "phase", "n_clusters", "n_years",
                     "sum_cluster_years", "isfi")]
    }
    out
  })

  result <- do.call(rbind, rows)

  if (verbose) {
    message(sprintf("[compute_isfi] Computed for %d %s(s).",
                    nrow(result), group_label))
  }

  result
}
