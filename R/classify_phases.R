#' Classify Habitats into Migration Phases by Peak-Anchored Density Cuts
#'
#' Assigns each temporary habitat to a migration phase (typically
#' *Wintering*, *Stopover*, *Breeding*) from the per-point Migration
#' Phase Index (MPI). The phase boundaries are **density valleys anchored
#' to the extreme modes** of the MPI distribution; each point is cut into
#' a phase deterministically and each habitat then takes the **majority**
#' phase of its points. No distributional model is fitted and no
#' prominence threshold is needed.
#'
#' @section Peak-anchored cutting:
#' MPI is constructed so that wintering piles up near 0 and breeding near
#' 1. The **leftmost** density peak of the pooled MPI is therefore the
#' wintering mode and the **rightmost** is the breeding mode -- a
#' structural fact of the index, not a tuned choice. The phase boundaries
#' follow directly:
#' \itemize{
#'   \item the **lower cut** is the first density valley to the *right* of
#'     the leftmost (wintering) peak;
#'   \item the **upper cut** is the last density valley to the *left* of
#'     the rightmost (breeding) peak.
#' }
#' Everything between the two cuts is the stopover band. Because the cuts
#' are anchored to the extreme peaks, any number of staging sites in the
#' middle -- however deep their internal valleys -- are absorbed into the
#' stopover phase and can never be mistaken for a phase boundary. This
#' needs **no prominence threshold**: the only ingredient is the kernel
#' density estimate (standard data-driven bandwidth).
#'
#' @section Assignment:
#' Each point is labelled deterministically by the band its MPI falls in
#' (at or below the lower cut -> wintering; above the upper cut ->
#' breeding; in between -> stopover), and each habitat is then assigned
#' its **majority** phase by voting over its points.
#'
#' @section Structural phase-count reduction:
#' The resolved number of phases is `length(cuts) + 1`. It is smaller than
#' the requested `n_phases` only when the data structurally lack the modes
#' to support it: two peaks (wintering and breeding, with no staging mode
#' between) give a single valley and two phases; a single mode gives no
#' valley and one phase. This is an objective property of the MPI
#' distribution, not a threshold decision, so the classifier never invents
#' a phase the data do not contain.
#'
#' @section Transient vs staging stopovers (`fixed_stopover`):
#' Local Moran's I resolves two distinct kinds of pause. Brief in-corridor
#' touch-downs appear as "LH" outliers (a low-speed fix embedded in the
#' high-speed migratory corridor), whereas multi-day staging sites
#' accumulate enough co-located fixes to form an "LL" cluster and surface
#' as a middle-MPI mode. Setting `fixed_stopover = "LH"` hard-assigns the
#' LH points to `fixed_label` (default `"transient_stopover"`) and removes
#' them from both the cut and the per-habitat vote, so the phases are
#' resolved only on the genuinely stationary (LL) points.
#'
#' @param habitat_points An `sf` object with at least three columns: the
#'   MPI, the individual ID, and the cluster ID (typically the output of
#'   [compute_mpi()] applied to [dbscan_habitats()] output).
#' @param mpi_col Column name of MPI values (default `"mpi"`).
#' @param individual_col Column name of the individual ID (default
#'   `"individual"`).
#' @param cluster_col Column name of the temporary habitat cluster ID
#'   (default `"cluster_id"`).
#' @param n_phases Integer; the target number of biological phases
#'   (default 3 for wintering / stopover / breeding, 2 for species without
#'   a distinct stopover). The resolved number may be smaller if the data
#'   structurally lack the modes (see *Structural phase-count reduction*).
#' @param prominence Optional valley-prominence floor in (0, 1), measured
#'   as the fractional density drop from the lower flanking peak,
#'   `(lower_peak - valley) / lower_peak`. Used **only** to ignore
#'   near-flat wiggle valleys before anchoring; default `0` (off). It does
#'   not decide the number of phases -- the peak anchoring does.
#' @param bw_adjust Numeric; bandwidth multiplier passed to
#'   [stats::density()] for the MPI KDE (default 1; larger = smoother,
#'   fewer spurious wiggle peaks/valleys).
#' @param phase_names Character vector of length equal to the *resolved*
#'   number of phases. Defaults: `c("Wintering", "Stopover", "Breeding")`
#'   for three, `c("Wintering", "Breeding")` for two. A mismatched vector
#'   is replaced by the defaults with a message.
#' @param fixed_stopover Character vector of `cluster_type` codes (e.g.
#'   `"LH"`) whose points are hard-assigned to `fixed_label` and excluded
#'   from the cut and the vote. `NULL` (default) disables this.
#' @param cluster_type_col Column name holding the Local Moran's I cluster
#'   type (default `"cluster_type"`); required only when `fixed_stopover`
#'   is non-NULL.
#' @param fixed_label Phase label assigned to the fixed points (default
#'   `"transient_stopover"`).
#' @param verbose Logical; print diagnostic information.
#'
#' @return The input `sf` with three new columns:
#'   \itemize{
#'     \item `phase_label`: integer phase rank (`1..K`) from the cut
#'       (`NA` for fixed and NA-MPI points).
#'     \item `phase_voted_label`: integer phase rank, the per-habitat
#'       majority vote.
#'     \item `phase`: factor of phase labels (with an extra `fixed_label`
#'       level when `fixed_stopover` is used).
#'   }
#'   An attribute `"phase_fit"` records the resolved phase count, the
#'   density peaks and valleys, the chosen cuts, the per-phase point
#'   counts and mean MPI, and the settings used.
#'
#' @examples
#' \dontrun{
#' demo <- make_demo_track(n_individuals = 2, n_points = 300)
#' seg  <- balm_segmentation(demo, "speed", "individual", verbose = FALSE)
#' stat <- seg[!is.na(seg$cluster_type) &
#'             seg$cluster_type %in% c("LL", "LH"), ]
#' hab  <- dbscan_habitats(stat, individual_col = "individual",
#'                          eps = 2, minPts = 12)
#' mpi  <- compute_mpi(hab, time_col = "time", verbose = FALSE)
#' phs  <- classify_phases(mpi, n_phases = 3, fixed_stopover = "LH")
#' table(phs$phase)
#' }
#'
#' @seealso [compute_mpi()] for the upstream MPI computation,
#'   [annotate_phases()] for interactive manual refinement.
#'
#' @importFrom stats density
#' @export
classify_phases <- function(habitat_points,
                            mpi_col          = "mpi",
                            individual_col   = "individual",
                            cluster_col      = "cluster_id",
                            n_phases         = 3L,
                            prominence       = 0,
                            bw_adjust        = 1,
                            phase_names      = NULL,
                            fixed_stopover   = NULL,
                            cluster_type_col = "cluster_type",
                            fixed_label      = "transient_stopover",
                            verbose          = TRUE) {
  if (!inherits(habitat_points, "sf")) {
    stop("`habitat_points` must be an `sf` object.", call. = FALSE)
  }
  for (col in c(mpi_col, individual_col, cluster_col)) {
    if (!col %in% names(habitat_points)) {
      stop("Column `", col, "` not found in `habitat_points`.",
           call. = FALSE)
    }
  }
  
  # ---- Identify fixed (hard-labelled) points, e.g. LH transients -----
  is_fixed <- rep(FALSE, nrow(habitat_points))
  if (!is.null(fixed_stopover)) {
    if (!cluster_type_col %in% names(habitat_points)) {
      stop("Column `", cluster_type_col, "` not found in ",
           "`habitat_points`; required when `fixed_stopover` is set.",
           call. = FALSE)
    }
    ct       <- as.character(habitat_points[[cluster_type_col]])
    is_fixed <- !is.na(ct) & ct %in% as.character(fixed_stopover)
    if (verbose) {
      message(sprintf(
        "[classify_phases] %d point(s) fixed as '%s' (cluster_type in {%s}); excluded from cut and voting.",
        sum(is_fixed), fixed_label,
        paste(fixed_stopover, collapse = ", ")))
    }
  }
  
  mpi_vals <- habitat_points[[mpi_col]]
  valid    <- !is.na(mpi_vals) & !is_fixed
  v        <- mpi_vals[valid]
  n_valid  <- length(v)
  if (n_valid < 10L) {
    stop("Need at least 10 non-NA, non-fixed MPI values (got ",
         n_valid, ").", call. = FALSE)
  }
  
  # ---- KDE peaks & valleys (no threshold) -------------------------------
  ex    <- .mpi_density_valleys(v, bw_adjust = bw_adjust)
  peaks <- ex$peak_x
  vals  <- ex$valley_x
  vprom <- ex$valley_prom
  # optional floor: drop near-flat wiggle valleys before anchoring (off by default)
  if (prominence > 0 && length(vals)) {
    keep <- vprom >= prominence
    vals <- vals[keep]; vprom <- vprom[keep]
  }
  
  target <- max(1L, as.integer(n_phases))
  
  # ---- Peak-anchored cuts ----------------------------------------------
  # leftmost peak = wintering mode, rightmost = breeding mode (by MPI
  # construction). lower cut = first valley right of the leftmost peak;
  # upper cut = last valley left of the rightmost peak; the middle (any
  # staging structure) is absorbed.
  cuts <- numeric(0)
  if (target >= 2L && length(peaks) >= 2L && length(vals) >= 1L) {
    lo_p <- min(peaks); hi_p <- max(peaks)
    sel  <- vals > lo_p & vals < hi_p            # valleys between the extreme peaks
    inb  <- vals[sel]; inb_prom <- vprom[sel]
    if (length(inb)) {
      if (target == 2L) {
        cuts <- inb[which.max(inb_prom)]         # single deepest divider
      } else {
        wc <- min(inb); bc <- max(inb)           # winter|stopover , stopover|breeding
        cuts <- if (wc == bc) wc else c(wc, bc)
        if (target > 3L) {                       # rare: add deepest interior valleys
          sel2 <- inb > wc & inb < bc
          if (any(sel2)) {
            extra <- inb[sel2][order(inb_prom[sel2], decreasing = TRUE)][
              seq_len(min(target - 3L, sum(sel2)))]
            cuts  <- sort(unique(c(cuts, extra)))
          }
        }
      }
    }
  }
  K <- length(cuts) + 1L
  if (verbose && K < target) {
    message(sprintf(
      "[classify_phases] data support %d phase(s) (target %d): only %d usable density valley(s) between the wintering and breeding peaks.",
      K, target, length(cuts)))
  }
  
  # ---- Hard assignment: phase = the band the MPI falls in ---------------
  point_phase_v <- if (length(cuts) == 0L) rep(1L, n_valid) else
    as.integer(cut(v, c(-Inf, cuts, Inf), labels = FALSE))
  point_phase <- rep(NA_integer_, nrow(habitat_points))
  point_phase[valid] <- point_phase_v
  habitat_points$phase_label <- point_phase
  
  # ---- Per-phase descriptive summaries (counts + mean MPI) --------------
  phase_sizes <- as.integer(tabulate(point_phase_v, nbins = K))
  phase_mpi   <- vapply(seq_len(K), function(k) {
    xk <- v[point_phase_v == k]
    if (length(xk)) mean(xk) else NA_real_
  }, numeric(1))
  
  if (verbose) {
    message(sprintf(
      "[classify_phases] %d phase(s); peaks {%s}; cut(s) {%s}; sizes {%s}; mean MPI {%s}",
      K,
      paste(sprintf("%.3f", peaks), collapse = ", "),
      if (length(cuts)) paste(sprintf("%.3f", cuts), collapse = ", ") else "none",
      paste(phase_sizes, collapse = ", "),
      paste(sprintf("%.3f", phase_mpi), collapse = ", ")))
  }
  
  # ---- Per-cluster majority voting (fixed points excluded) --------------
  df <- sf::st_drop_geometry(habitat_points)
  vote_key <- paste(df[[individual_col]], df[[cluster_col]], sep = "::")
  voted_map <- tapply(point_phase, vote_key, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) == 0L) return(NA_integer_)
    tab <- table(x_clean)
    as.integer(names(tab)[which.max(tab)])
  })
  habitat_points$phase_voted_label <- as.integer(voted_map[vote_key])
  
  # ---- Map to phase names -----------------------------------------------
  default_names <- function(k) {
    if (k == 3L)      c("Wintering", "Stopover", "Breeding")
    else if (k == 2L) c("Wintering", "Breeding")
    else if (k == 1L) "Stationary"
    else              paste0("Phase_", seq_len(k))
  }
  if (is.null(phase_names)) {
    phase_names <- default_names(K)
  } else if (length(phase_names) != K) {
    if (verbose) {
      message(sprintf(
        "[classify_phases] supplied phase_names has length %d but %d phase(s) resolved; using defaults.",
        length(phase_names), K))
    }
    phase_names <- default_names(K)
  }
  phase_chr <- as.character(phase_names[habitat_points$phase_voted_label])
  if (any(is_fixed)) phase_chr[is_fixed] <- fixed_label
  
  # ---- Level ordering (by ascending mean MPI when fixed pts present) ----
  if (any(is_fixed)) {
    lvl_means  <- tapply(mpi_vals, phase_chr, mean, na.rm = TRUE)
    all_levels <- names(sort(lvl_means))
  } else {
    all_levels <- phase_names
  }
  habitat_points$phase <- factor(phase_chr, levels = all_levels)
  
  # ---- Attach diagnostics -----------------------------------------------
  attr(habitat_points, "phase_fit") <- list(
    method         = "peak-anchored density-valley cut",
    K              = K,
    n_phases       = K,
    target_phases  = target,
    peaks          = peaks,
    valleys        = vals,
    phase_cuts     = if (length(cuts)) cuts else NA_real_,
    bw_adjust      = bw_adjust,
    prominence_floor = prominence,
    phase_sizes    = phase_sizes,
    phase_mpi_mean = phase_mpi,
    phase_names    = phase_names,
    fixed_stopover = fixed_stopover,
    fixed_label    = if (!is.null(fixed_stopover)) fixed_label else NULL,
    n_fixed        = sum(is_fixed)
  )
  habitat_points
}

#' Density peaks and valleys of a 1-D sample
#'
#' Internal helper for [classify_phases()]. Computes a kernel density
#' estimate and returns its interior local maxima (peaks) and minima
#' (valleys), plus each valley's prominence -- the **fractional** density
#' drop from the lower of its two nearest flanking peaks down to the
#' valley floor, `(lower_peak - valley) / lower_peak`. The prominence is
#' returned for optional ranking/flooring only; valley selection in
#' [classify_phases()] is anchored to the extreme peaks and needs no
#' prominence threshold.
#'
#' @param v Numeric vector.
#' @param bw_adjust Bandwidth multiplier for [stats::density()].
#' @return A list with `dd` (the density), `peak_x` (peak locations,
#'   ascending), `valley_x` (valley locations, ascending) and
#'   `valley_prom` (their fractional prominences).
#' @keywords internal
#' @importFrom stats density
.mpi_density_valleys <- function(v, bw_adjust = 1) {
  dd <- stats::density(v, adjust = bw_adjust)
  s  <- diff(sign(diff(dd$y)))
  pk <- which(s == -2) + 1L          # local maxima (peaks)
  vy <- which(s ==  2) + 1L          # local minima (valleys)
  prom <- if (!length(vy)) numeric(0) else vapply(vy, function(j) {
    lp <- pk[pk < j]; rp <- pk[pk > j]
    if (!length(lp) || !length(rp)) return(0)
    mp <- min(dd$y[max(lp)], dd$y[min(rp)])   # lower of the two nearest peaks
    if (mp <= 0) return(0)
    (mp - dd$y[j]) / mp                        # fractional drop from that peak
  }, numeric(1))
  vo <- order(dd$x[vy])
  list(dd          = dd,
       peak_x      = sort(dd$x[pk]),
       valley_x    = dd$x[vy][vo],
       valley_prom = if (length(prom)) prom[vo] else numeric(0))
}