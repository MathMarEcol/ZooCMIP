
      # ---------------------------------------------------------------------------
      # PER-MODEL SPATIAL MAPS
      # ---------------------------------------------------------------------------
      #
      # One PDF per climate model, each containing a 4 (scenario: ssp126/245/370/585)
      # x 3 (variable: zooc/zmicro/zmeso) grid of % change maps - identical in
      # style to the ensemble-mean grid in 08_Climate_SpatialPlot.R, but showing
      # a SINGLE model's Percent_change layer rather than the ensemble mean.
      #
      # Where a model x variable/scenario combination does not exist (e.g. zooc-
      # only models have no zmicro/zmeso data), the panel is left blank (grey
      # "No data" placeholder), so the grid layout is always the full 4 x 3.
      #
      # Output: Figures/PerModel/<ModelName>_PercentChange.pdf  (~14 files)
      #
      # Reads from: /Volumes/T9/ClimateData/<var>/biomassChange/*.tif
      # (the same biomassChange.tif files written by 06_Climate_Review.R and
      #  used by 08_Climate_SpatialPlot.R for the ensemble mean).
      # ---------------------------------------------------------------------------

      library(terra)
      library(tidyterra)
      library(sf)
      library(ggplot2)
      library(patchwork)
      library(purrr)
      library(dplyr)
      library(tidyr)

      # ---------------------------------------------------------------------------
      # Shared constants (identical to 08_Climate_SpatialPlot.R)
      # ---------------------------------------------------------------------------

      scenario_order  <- c("ssp126", "ssp245", "ssp370", "ssp585")
      scenario_labels <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
      scenario_label_map <- setNames(scenario_labels, scenario_order)

      vars <- c("zooc", "zmicro", "zmeso")

      var_labels <- c(
        zooc   = "All Zooplankton",
        zmicro = "Microzooplankton",
        zmeso  = "Mesozooplankton"
      )

      biomasschange_dir <- function(variable) {
        file.path("", "Volumes", "T9", "ClimateData", variable, "biomassChange")
      }

      # Coastline
      world <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") |>
        sf::st_transform("EPSG:8857")

      # Equal Earth boundary polygon (densely sampled WGS84 bounding box edges,
      # reprojected) - used to mask raster corners to NA and draw the globe outline.
      ee_boundary <- {
        n <- 1000
        sf::st_sfc(
          sf::st_polygon(list(rbind(
            cbind(seq(-180,  180, length.out = n), rep( 90, n)),
            cbind(rep( 180, n), seq( 90, -90, length.out = n)),
            cbind(seq( 180, -180, length.out = n), rep(-90, n)),
            cbind(rep(-180, n), seq(-90,  90, length.out = n)),
            c(-180, 90)
          ))),
          crs = 4326
        ) |> sf::st_transform("EPSG:8857")
      }

      # ---------------------------------------------------------------------------
      # Filename parser (identical to 08_Climate_SpatialPlot.R)
      # ---------------------------------------------------------------------------

      #' Parse a biomassChange .tif filename into its metadata fields.
      #'
      #' Convention: <variable>_<model>_<scenario>_<variant>_biomassChange.tif
      parse_biomasschange_filename <- function(path) {
        bits <- strsplit(basename(path), "_")[[1]]
        n <- length(bits)
        stopifnot(
          "Unexpected biomassChange filename format" =
            n >= 4 && bits[n] == "biomassChange.tif"
        )
        tibble::tibble(
          Variable = bits[1],
          Model    = paste(bits[2:(n - 3)], collapse = "_"),
          Scenario = bits[n - 2],
          Variant  = bits[n - 1],
          Path     = path
        )
      }

      # ---------------------------------------------------------------------------
      # Build a catalogue of ALL available biomassChange files across all variables
      # ---------------------------------------------------------------------------

      all_meta <- purrr::map(vars, function(v) {
        dir   <- biomasschange_dir(v)
        files <- list.files(dir, pattern = "[.]tif$", full.names = TRUE)
        if (length(files) == 0) return(NULL)
        purrr::map(files, parse_biomasschange_filename) |> dplyr::bind_rows()
      }) |>
        dplyr::bind_rows()

      # Restrict to the four tier-1 SSPs used throughout this project
      all_meta <- all_meta |>
        dplyr::filter(Scenario %in% scenario_order)

      # Unique model names (sorted for reproducible PDF ordering)
      all_models <- sort(unique(all_meta$Model))
      message("Models found: ", paste(all_models, collapse = ", "))

      # ---------------------------------------------------------------------------
      # Map theme (identical to 08_Climate_SpatialPlot.R)
      # ---------------------------------------------------------------------------

      base_map_theme <- function(base_size = 8) {
        theme_minimal(base_size = base_size) +
          theme(
            panel.grid       = element_line(color = "gray90", linewidth = 0.2),
            panel.background = element_blank(),
            plot.background  = element_rect(fill = "white", color = NA),
            axis.text        = element_blank(),
            axis.title       = element_blank(),
            plot.title       = element_blank(),
            legend.key.height = unit(0.3, "cm"),
            legend.key.width  = unit(1,   "cm"),
            legend.title      = element_text(size = 8),
            legend.text       = element_text(size = 7),
            plot.margin       = margin(2, 2, 2, 2)
          )
      }

      # ---------------------------------------------------------------------------
      # Single-panel plotter
      # ---------------------------------------------------------------------------

      #' Plot one panel of the per-model % change grid.
      #'
      #' @param pct_rast  SpatRaster with the Percent_change layer for this
      #'   model/variable/scenario, already projected to EPSG:8857 and masked to
      #'   the Equal Earth boundary.  Pass NULL for a blank placeholder.
      #' @param show_row_label  Scenario label for the y-axis (left margin); used
      #'   for the first variable column only.
      #' @param show_col_title  Variable title for the plot title (top); used for
      #'   the first scenario row only.
      plot_model_panel <- function(pct_rast,
                                   show_row_label = NULL,
                                   show_col_title = NULL) {

        if (is.null(pct_rast)) {
          # Blank placeholder - same style as 08_Climate_SpatialPlot.R
          p <- ggplot() +
            annotate("text", x = 0, y = 0,
                     label = "No data", size = 3, color = "gray50") +
            theme_void() +
            theme(plot.margin = margin(2, 2, 2, 2))

          # Still attach row/column labels so the grid reads correctly
          if (!is.null(show_row_label)) {
            p <- p + labs(y = show_row_label) +
              theme(axis.title.y = element_text(size = 9, face = "bold", angle = 90))
          }
          if (!is.null(show_col_title)) {
            p <- p + ggtitle(show_col_title) +
              theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
          }
          return(p)
        }

        p <- ggplot() +
          geom_spatraster(data = pct_rast) +
          geom_sf(data = world,       color = "grey30", linewidth = 0.2, fill = "grey30") +
          geom_sf(data = ee_boundary, color = "gray40", linewidth = 0.3, fill = NA) +
          scale_fill_gradientn(
            colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                       "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
            limits    = c(-50, 50),
            na.value  = "transparent",
            name      = "%\nchange",
            oob       = scales::squish
          ) +
          coord_sf(
            crs    = "EPSG:8857",
            xlim   = c(-17243959, 17243959),
            ylim   = c(-8343134,   8343134),
            expand = FALSE
          ) +
          base_map_theme()

        if (!is.null(show_row_label)) {
          p <- p + labs(y = show_row_label) +
            theme(axis.title.y = element_text(size = 9, face = "bold", angle = 90))
        }
        if (!is.null(show_col_title)) {
          p <- p + ggtitle(show_col_title) +
            theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
        }

        p
      }

      # ---------------------------------------------------------------------------
      # Per-model PDF builder
      # ---------------------------------------------------------------------------

      #' Build and save one PDF for a single model.
      #'
      #' Layout: rows = Scenario (ssp126/245/370/585), columns = Variable
      #' (zooc/zmicro/zmeso).  Panels where the model has no data are blank.
      #' A single shared legend is collected at the bottom of the figure.
      #'
      #' @param model  Character string matching the Model field in all_meta.
      build_model_pdf <- function(model) {

        message("Building PDF for model: ", model)

        ee_mask <- terra::vect(ee_boundary)

        # Pre-load every available raster for this model (all vars x scenarios)
        model_meta <- all_meta |> dplyr::filter(Model == model)

        # Build the full 4 x 3 panel list (row = scenario, col = variable)
        panels <- vector("list", length(scenario_order) * length(vars))
        idx <- 1L

        for (j in seq_along(vars)) {
          variable  <- vars[j]

          for (i in seq_along(scenario_order)) {
            scenario  <- scenario_order[i]

            show_row  <- if (j == 1) scenario_label_map[[scenario]] else NULL
            show_col  <- if (i == 1) var_labels[[variable]]         else NULL

            # Look up this model/variable/scenario combination
            row <- model_meta |>
              dplyr::filter(Variable == variable, Scenario == scenario)

            if (nrow(row) == 0) {
              # No data for this combination - blank panel
              panels[[idx]] <- plot_model_panel(NULL,
                                                show_row_label = show_row,
                                                show_col_title = show_col)
            } else {
              # Load Percent_change layer, project, mask
              pct <- terra::rast(row$Path[[1]])[["Percent_change"]]
              pct <- terra::project(pct, "EPSG:8857")
              pct <- terra::mask(pct, ee_mask)

              panels[[idx]] <- plot_model_panel(pct,
                                                show_row_label = show_row,
                                                show_col_title = show_col)
            }

            idx <- idx + 1L
          }
        }

        # Arrange panels: patchwork fills column-by-column when ncol is specified,
        # but we built panels in column-major order (variable outer, scenario inner)
        # so wrap_plots with ncol = n_scn gives rows = scenario, cols = variable.
        grid <- wrap_plots(panels, ncol = length(scenario_order)) +
          plot_layout(guides = "collect") +
          plot_annotation(
            title = paste0(model, " — % Change in Zooplankton Biomass (1993\u20132014 to 2081\u20132100)"),
            theme = theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
            )
          ) &
          theme(
            legend.position   = "bottom",
            legend.key.width  = unit(3,   "cm"),
            legend.key.height = unit(0.4, "cm"),
            legend.title      = element_text(size = 9),
            legend.text       = element_text(size = 8)
          )

        # Safe filename: replace characters that are awkward in file paths
        safe_model <- gsub("[^A-Za-z0-9._-]", "_", model)
        out_path   <- file.path("Figures", "PerModel",
                                paste0(safe_model, "_PercentChange.pdf"))

        dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
        ggsave(out_path, plot = grid, width = 12, height = 10)
        message("  Saved: ", out_path)
      }

      # ---------------------------------------------------------------------------
      # Run for every model
      # ---------------------------------------------------------------------------

      purrr::walk(all_models, build_model_pdf)

      message("Done. PDFs written to Figures/PerModel/")
