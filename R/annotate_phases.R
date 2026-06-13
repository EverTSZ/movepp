#' Interactively annotate phases on a behaviour barcode
#'
#' Launches a Shiny gadget with a habitat-occupancy *behaviour barcode*
#' (time x habitat `rank_col`) linked to a clean coastline map (no country
#' borders). Box-select on EITHER panel; the selection is ringed in black on
#' both and can be assigned a life-history phase or a free-text custom label.
#' With `individual_col` set, a dropdown restricts the view to one individual.
#' Returns the input with an added `phase_manual` column.
#'
#' If `phase_col` is supplied, the gadget **opens with those labels as the
#' starting classification** (e.g. the output of [classify_phases()]), so
#' only the mis-assigned points need correcting rather than annotating from
#' scratch. Because the barcode is organised by habitat rank, box-selecting a
#' whole rank-row over a time span selects an entire habitat -- so corrections
#' can be made per habitat (cluster), matching the cluster-level classification output.
#'
#' @section Interaction:
#' \itemize{
#'   \item With `individual_col`, pick a bird from the **Individual** dropdown
#'     (or `(all)`); switching clears the current selection.
#'   \item Box-select on the **barcode** (drag) or **map** (modebar box/lasso);
#'     selected points are ringed in black on BOTH panels.
#'   \item Click a preset phase button, or type a custom label and *Apply
#'     custom*; *Unassign* clears labels.
#'   \item Toggle **Colour by** between habitat rank and assigned label.
#'   \item *Done* returns the annotated data; *Cancel* aborts.
#' }
#'
#' Requires \pkg{shiny}, \pkg{miniUI} and \pkg{plotly} (Suggests).
#'
#' @param points An `sf` POINT object with at least `time_col` and `rank_col`.
#' @param time_col POSIXct/Date column name (string).
#' @param rank_col Discrete habitat-rank column (default `"cluster_rank"`).
#' @param individual_col Optional individual-id column; adds a selector.
#'   `NULL` (default) shows all points.
#' @param phase_col Optional column of starting phase labels to seed the
#'   annotation (e.g. `"phase"` from [classify_phases()]); the gadget opens
#'   with these instead of all-"unassigned". `NULL` (default) starts blank.
#' @param basemap Coastline-map colour preset: `"light"` (default),
#'   `"natural"`, or `"gray"`. Coastlines, rivers, lakes, land and ocean are
#'   drawn; country borders never are.
#' @param phases Preset phase labels (default `c("Wintering","Stopover","Breeding")`).
#' @param viewer `"browser"` (default), `"pane"`, or `"dialog"`.
#'
#' @return The input `sf` with an added factor column `phase_manual`, original
#'   row order; `NULL` if cancelled.
#'
#' @seealso [classify_phases()] for the automatic classification this can
#'   refine, [dbscan_habitats()] for producing the habitat ranks.
#' @export
annotate_phases <- function(points,
                            time_col       = "time",
                            rank_col       = "cluster_rank",
                            individual_col = NULL,
                            phase_col      = NULL,
                            basemap        = c("light", "natural", "gray"),
                            phases         = c("Wintering", "Stopover", "Breeding"),
                            viewer         = c("browser", "pane", "dialog")) {
  
  for (pkg in c("shiny", "miniUI", "plotly")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for annotate_phases(). ",
           "Install it with install.packages('", pkg, "').", call. = FALSE)
  }
  if (!inherits(points, "sf")) stop("`points` must be an `sf` object.", call. = FALSE)
  if (!time_col %in% names(points)) stop("Column `", time_col, "` not found.", call. = FALSE)
  if (!rank_col %in% names(points)) stop("Column `", rank_col, "` not found.", call. = FALSE)
  has_ind <- !is.null(individual_col)
  if (has_ind && !individual_col %in% names(points))
    stop("Column `", individual_col, "` not found.", call. = FALSE)
  if (!is.null(phase_col) && !phase_col %in% names(points))
    stop("Column `", phase_col, "` not found.", call. = FALSE)
  viewer  <- match.arg(viewer)
  basemap <- match.arg(basemap)
  geo_style <- switch(basemap,
                      light   = list(land = "#f4f4f2", ocean = "#dde9f3", coast = "#9fb3c8",
                                     river = "#bcd4e6", lake = "#cfe2f3"),
                      natural = list(land = "#e8e2d0", ocean = "#cfe0ea", coast = "#8a8f96",
                                     river = "#a9c7dd", lake = "#bcd9ec"),
                      gray    = list(land = "#ededed", ocean = "#fbfbfb", coast = "#bdbdbd",
                                     river = "#d9d9d9", lake = "#e8e8e8"))
  
  co <- sf::st_coordinates(points)
  d  <- data.frame(
    .id  = seq_len(nrow(points)),
    time = points[[time_col]], lon = co[, 1], lat = co[, 2],
    rank = as.integer(points[[rank_col]]),
    .ind = if (has_ind) as.character(points[[individual_col]]) else "all")
  d <- d[order(d$time), ]
  d$tip <- paste0(format(d$time), "<br>Habitat ", d$rank,
                  if (has_ind) paste0("<br>", d$.ind) else "")
  
  ind_levels <- sort(unique(d$.ind))
  n_hab    <- max(d$rank, na.rm = TRUE)
  rank_pal <- grDevices::hcl.colors(n_hab, "Spectral")
  ring     <- list(color = "rgba(0,0,0,0)", size = 12,
                   line = list(color = "black", width = 2))
  
  preset_pal <- c("#2C7BB6", "#FDAE61", "#D7191C", "#7B3294",
                  "#1B7837", "#E66101", "#01665E")
  base_cols  <- stats::setNames(
    preset_pal[(seq_along(phases) - 1L) %% length(preset_pal) + 1L], phases)
  base_cols["unassigned"] <- "#BDBDBD"
  custom_pal <- c("#C51B7D", "#BF812D", "#35978F", "#762A83", "#4D9221")
  
  # seed starting labels (from phase_col) or start blank
  init_phase <- if (is.null(phase_col)) {
    rep("unassigned", nrow(points))
  } else {
    v <- as.character(points[[phase_col]]); v[is.na(v)] <- "unassigned"; v
  }
  seed_extra <- setdiff(unique(init_phase), names(base_cols))
  if (length(seed_extra))
    base_cols[seed_extra] <-
    custom_pal[(seq_along(seed_extra) - 1L) %% length(custom_pal) + 1L]
  
  phase_ids     <- paste0("phase_", seq_along(phases))
  phase_buttons <- mapply(
    function(id, lab, col) shiny::actionButton(
      id, lab, style = paste0("background:", col, ";color:#fff")),
    phase_ids, phases, base_cols[phases], SIMPLIFY = FALSE)
  
  indiv_control <- if (has_ind)
    shiny::selectInput("indiv", "Individual", choices = c("(all)", ind_levels),
                       selected = ind_levels[1], width = "220px") else NULL
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Behaviour barcode - annotate phases"),
    miniUI::miniContentPanel(
      shiny::fillCol(
        flex = c(NA, 1, 1),
        shiny::div(
          style = "padding:6px; line-height:2.2",
          indiv_control, phase_buttons,
          shiny::actionButton("clear", "Unassign"),
          shiny::tags$span(style = "margin-left:12px"),
          shiny::tags$input(id = "custom", type = "text",
                            placeholder = "custom label e.g. death", style = "width:160px"),
          shiny::actionButton("apply_custom", "Apply custom"),
          shiny::tags$span(style = "margin-left:12px"),
          shiny::radioButtons("colorby", NULL, inline = TRUE,
                              choices = c("Colour by label" = "label", "Colour by habitat rank" = "rank"),
                              selected = "label"),
          shiny::textOutput("status", inline = TRUE)
        ),
        plotly::plotlyOutput("barcode", height = "100%"),
        plotly::plotlyOutput("map",     height = "100%")
      )
    )
  )
  
  server <- function(input, output, session) {
    rv <- shiny::reactiveValues(phase = init_phase, colmap = base_cols)
    current_sel <- shiny::reactiveVal(integer(0))
    
    dsub <- shiny::reactive({
      if (!has_ind || identical(input$indiv, "(all)")) d
      else d[d$.ind == input$indiv, , drop = FALSE]
    })
    
    read_sel <- function(src) {
      ed <- plotly::event_data("plotly_selected", source = src)
      if (is.null(ed) || is.null(ed$customdata)) integer(0) else as.integer(ed$customdata)
    }
    shiny::observeEvent(plotly::event_data("plotly_selected", source = "bc"),
                        current_sel(read_sel("bc")), ignoreNULL = FALSE)
    shiny::observeEvent(plotly::event_data("plotly_selected", source = "map"),
                        current_sel(read_sel("map")), ignoreNULL = FALSE)
    if (has_ind)
      shiny::observeEvent(input$indiv, current_sel(integer(0)), ignoreInit = TRUE)
    
    cols_for <- function(sub) {
      if (shiny::isolate(input$colorby) == "rank") {
        rank_pal[sub$rank]
      } else {
        cols <- rv$colmap[rv$phase[sub$.id]]; cols[is.na(cols)] <- "#BDBDBD"; unname(cols)
      }
    }
    refresh_cols <- function() {
      cols <- cols_for(dsub())
      plotly::plotlyProxyInvoke(plotly::plotlyProxy("barcode", session), "restyle",
                                list(marker.color = list(cols)), 0)
      plotly::plotlyProxyInvoke(plotly::plotlyProxy("map", session), "restyle",
                                list(marker.color = list(cols)), 0)
    }
    assign_phase <- function(ph) {
      ids <- current_sel(); if (length(ids)) rv$phase[ids] <- ph
    }
    
    lapply(seq_along(phases), function(i)
      shiny::observeEvent(input[[phase_ids[i]]], { assign_phase(phases[i]); refresh_cols() }))
    shiny::observeEvent(input$clear, { assign_phase("unassigned"); refresh_cols() })
    shiny::observeEvent(input$apply_custom, {
      lab <- trimws(input$custom)
      if (nzchar(lab)) {
        if (!lab %in% names(rv$colmap)) {
          used <- length(setdiff(names(rv$colmap), names(base_cols)))
          rv$colmap[lab] <- custom_pal[(used %% length(custom_pal)) + 1L]
        }
        assign_phase(lab); refresh_cols()
      }
    })
    shiny::observeEvent(input$colorby, refresh_cols(), ignoreInit = TRUE)
    
    shiny::observeEvent(current_sel(), {
      sub <- dsub(); hl <- sub[sub$.id %in% current_sel(), , drop = FALSE]
      plotly::plotlyProxyInvoke(plotly::plotlyProxy("barcode", session), "restyle",
                                list(x = list(hl$time), y = list(hl$rank)), 1)
      plotly::plotlyProxyInvoke(plotly::plotlyProxy("map", session), "restyle",
                                list(lon = list(hl$lon), lat = list(hl$lat)), 1)
    }, ignoreNULL = FALSE)
    
    output$status <- shiny::renderText({
      tab <- sort(table(rv$phase), decreasing = TRUE)
      sprintf("  selected:%d  |  %s", length(current_sel()),
              paste(names(tab), tab, sep = ":", collapse = "  "))
    })
    
    output$barcode <- plotly::renderPlotly({
      sub <- dsub()
      plotly::plot_ly(source = "bc") |>
        plotly::add_trace(x = sub$time, y = sub$rank, customdata = sub$.id,
                          type = "scattergl", mode = "markers",
                          marker = list(color = shiny::isolate(cols_for(sub)), size = 7),
                          text = sub$tip, hoverinfo = "text", name = "points") |>
        plotly::add_trace(x = numeric(0), y = numeric(0),
                          type = "scattergl", mode = "markers",
                          marker = ring, hoverinfo = "skip", name = "selected") |>
        plotly::layout(dragmode = "select", showlegend = FALSE,
                       yaxis = list(title = "Habitat (rank, S->N)",
                                    tickvals = seq_len(max(sub$rank, na.rm = TRUE))),
                       xaxis = list(title = "time")) |>
        plotly::event_register("plotly_selected")
    })
    
    output$map <- plotly::renderPlotly({
      sub <- dsub()
      plotly::plot_ly(source = "map") |>
        plotly::add_trace(lon = sub$lon, lat = sub$lat, customdata = sub$.id,
                          type = "scattergeo", mode = "markers",
                          marker = list(color = shiny::isolate(cols_for(sub)), size = 7),
                          text = sub$tip, hoverinfo = "text", name = "points") |>
        plotly::add_trace(lon = numeric(0), lat = numeric(0),
                          type = "scattergeo", mode = "markers",
                          marker = ring, hoverinfo = "skip", name = "selected") |>
        plotly::layout(showlegend = FALSE, geo = list(
          showcountries  = FALSE,
          showcoastlines = TRUE,  coastlinecolor = geo_style$coast, coastlinewidth = 0.5,
          showland       = TRUE,  landcolor      = geo_style$land,
          showocean      = TRUE,  oceancolor     = geo_style$ocean,
          showrivers     = TRUE,  rivercolor     = geo_style$river, riverwidth = 0.5,
          showlakes      = TRUE,  lakecolor      = geo_style$lake,
          projection     = list(type = "natural earth"),
          fitbounds      = "locations")) |>
        plotly::event_register("plotly_selected")
    })
    
    shiny::observeEvent(input$done,   shiny::stopApp(rv$phase))
    shiny::observeEvent(input$cancel, shiny::stopApp(NULL))
  }
  
  vw <- switch(viewer,
               browser = shiny::browserViewer(),
               pane    = shiny::paneViewer(),
               dialog  = shiny::dialogViewer("Annotate phases"))
  res <- shiny::runGadget(ui, server, viewer = vw)
  if (is.null(res)) return(invisible(NULL))
  
  pm <- res
  used_custom <- setdiff(unique(pm), c(phases, "unassigned"))
  points$phase_manual <- factor(pm, levels = c(phases, used_custom, "unassigned"))
  points
}
