# ---------------------------------------------------------------------------
# SPATIAL ENSEMBLE MAPS
# ---------------------------------------------------------------------------
#
# Two 4 (scenario: ssp126/245/370/585) x 3 (variable: zooc/zmicro/zmeso)
# patchwork grids built from the per-model/scenario biomassChange.tif files
# written by 06_Climate_Review.R's zoo_change() to
# /Volumes/T9/ClimateData/<var>/biomassChange/*.tif.
# Each file now has 4 layers: Hist_mean, Future_mean, Absolute_change,
# Percent_change (see calc_biomass_change() in 06_Climate_Review.R).
#
#   1. ENSEMBLE MEAN % CHANGE: average Hist_mean and Future_mean across
#      models at each pixel, then compute % change from those two ensemble
#      means. This avoids per-model division by near-zero baseline values
#      (see design note d in 06_Climate_Review.R).
#
#   2. MODEL AGREEMENT %: at each pixel, count how many models have the same
#      sign of Absolute_change as the ensemble mean Absolute_change, divided
#      by the number of models with valid (non-NA) data at that pixel, times
#      100. Using Absolute_change sign avoids any per-model % division.
#      Expressing agreement as a % (0-100) rather than a raw count makes
#      every panel directly comparable regardless of ensemble size, and
#      allows a single shared legend across all variables.
#
# CanESM5/CanESM5-CanOE are included (unlike the point-range plot in
# 07_Climate_Plot.R which still excludes CanESM5 at the user's request).
#
# LAYOUT: rows = Scenario (4), columns = Variable (3).
# ---------------------------------------------------------------------------
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(patchwork)

scenario_order <- c("ssp126", "ssp245", "ssp370", "ssp585")
scenario_labels <- c(
  "Low (SSP1-2.6)",
  "Medium (SSP2-4.5)",
  "High (SSP3-7.0)",
  "Very High (SSP5-8.5)"
)
scenario_label_map <- setNames(scenario_labels, scenario_order)

vars <- c("zooc", "zmicro", "zmeso")

var_labels <- c(
  zooc = "Total Zooplankton",
  zmicro = "Small Zooplankton",
  zmeso = "Large Zooplankton"
)

biomasschange_dir <- function(variable) {
  file.path("", "Volumes", "T9", "ClimateData", variable, "biomassChange")
}

# Coastline for context, matching the style already used in
# 05_Climate_Spatial.R.
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  sf::st_transform("EPSG:8857")

# Equal Earth projection boundary polygon: the outline of the full globe
# in EPSG:8857 is a smooth ellipse-like curve, NOT a rectangle. We build
# it by densely sampling the WGS84 bounding box edges and reprojecting,
# then use it to (a) mask raster corners to NA so they plot as transparent,
# and (b) draw the boundary outline on each panel.
ee_boundary <- {
  n <- 1000
  sf::st_sfc(
    sf::st_polygon(list(rbind(
      cbind(seq(-180, 180, length.out = n), rep(90, n)),
      cbind(rep(180, n), seq(90, -90, length.out = n)),
      cbind(seq(180, -180, length.out = n), rep(-90, n)),
      cbind(rep(-180, n), seq(-90, 90, length.out = n)),
      c(-180, 90)
    ))),
    crs = 4326
  ) |>
    sf::st_transform("EPSG:8857")
}

#' Parse a Data/<var>/biomassChange/*.tif filename into its metadata
#'
#' Filenames follow the fixed convention written by zoo_change() in
#' 06_Climate_Review.R:
#'   <variable>_<model>_<scenario>_<variant>_biomassChange.tif
#' i.e. a SINGLE trailing token ("biomassChange.tif"), unlike
#' parse_integrated_filename() in 06_Climate_Review.R which parses
#' *_biomass_integrated.tif filenames with TWO trailing tokens ("biomass",
#' "integrated.tif"). Re-using that function's index offsets here (n-1 for
#' Scenario) was the bug that broke every map in this section: it silently
#' shifted every field over by one, so `Scenario` actually held the Variant
#' string (e.g. "r1i1p1f1") instead of the SSP code, meaning
#' `Scenario == scenario` never matched anything and every ensemble map
#' below (mean_pct/agreement) was built from zero files.
#' Model is still recovered from the front and Scenario/Variant positionally
#' from the back, so models with hyphens (e.g. "CanESM5-CanOE") are handled
#' correctly without a regex.
parse_biomasschange_filename <- function(path) {
  bits <- strsplit(basename(path), "_")[[1]]
  n <- length(bits)
  stopifnot(
    "Unexpected biomassChange filename format" = n >= 4 &&
      bits[n] == "biomassChange.tif"
  )
  tibble::tibble(
    Variable = bits[1],
    Model = paste(bits[2:(n - 3)], collapse = "_"),
    Scenario = bits[n - 2],
    Variant = bits[n - 1],
    Path = path
  )
}

#' Build the per-pixel model-ensemble mean % change AND the per-pixel model
#' agreement count for one variable/scenario combination.
#'
#' APPROACH: uses the "Hist_mean" and "Future_mean" layers now stored in each
#' model's 4-layer biomassChange.tif (written by calc_biomass_change() in
#' 06_Climate_Review.R). Rather than averaging per-model pixel-wise %
#' changes, we:
#'   1. Average Hist_mean across models at each pixel   → ensemble_hist
#'   2. Average Future_mean across models at each pixel → ensemble_future
#'   3. mean_pct = (ensemble_future - ensemble_hist) / ensemble_hist * 100
#'
#' Dividing two ensemble means is far more robust than averaging per-model
#' ratios, because the denominator (ensemble mean historical biomass) is
#' smoothed across models and therefore much less susceptible to individual
#' model near-zero baseline artifacts.
#'
#' Model agreement is computed from each model's Absolute_change sign
#' (not per-model % sign), which avoids per-model division entirely.
#'
#' @param variable One of "zooc", "zmicro", "zmeso".
#' @param scenario One of the four tier-1 SSP codes.
#' @return A list with `mean_pct` (SpatRaster, ensemble mean % change),
#'   `agreement` (SpatRaster, integer count of models agreeing in sign with
#'   `mean_pct`), and `n_models` (number of models contributing).
build_ensemble_maps <- function(variable, scenario) {
  dir <- biomasschange_dir(variable)
  files <- list.files(dir, pattern = "[.]tif$", full.names = TRUE)

  meta <- purrr::map(files, parse_biomasschange_filename) |>
    dplyr::bind_rows() |>
    dplyr::filter(Scenario == scenario)

  if (nrow(meta) == 0) {
    return(list(
      mean_pct = NULL,
      agreement = NULL,
      n_models = 0L,
      models = character(0)
    ))
  }

  # Load Hist_mean, Future_mean and Absolute_change from each model's
  # 4-layer (Hist_mean, Future_mean, Absolute_change, Percent_change)
  # biomassChange.tif. All three are needed: Hist_mean + Future_mean for
  # the ensemble % change, Absolute_change for the agreement layer.
  load_layer <- function(path, layer_name) terra::rast(path)[[layer_name]]

  hist_stack <- terra::rast(purrr::map(meta$Path, load_layer, "Hist_mean"))
  fut_stack <- terra::rast(purrr::map(meta$Path, load_layer, "Future_mean"))
  abs_stack <- terra::rast(purrr::map(meta$Path, load_layer, "Absolute_change"))

  names(hist_stack) <- meta$Model
  names(fut_stack) <- meta$Model
  names(abs_stack) <- meta$Model

  # --- Ensemble means -------------------------------------------------------
  # na.rm = TRUE: a pixel masked to NA in one model (invalid baseline - see
  # design note d in 06_Climate_Review.R) simply drops out of that pixel's
  # mean rather than poisoning it.
  ensemble_hist <- mean(hist_stack, na.rm = TRUE)
  ensemble_fut <- mean(fut_stack, na.rm = TRUE)

  # Normalise NaN (all-NA pixel across every model) to NA for consistent
  # downstream handling (mean() of all-NA returns NaN, not NA - verified).
  ensemble_hist <- terra::ifel(is.nan(ensemble_hist), NA, ensemble_hist)
  ensemble_fut <- terra::ifel(is.nan(ensemble_fut), NA, ensemble_fut)

  # --- Ensemble % change from two means ------------------------------------
  # Where ensemble_hist is NA or <= 0 (land / all-artifact pixel), result
  # is NA. This avoids division by zero or a sign flip from a near-zero
  # denominator.
  mean_pct <- terra::ifel(
    ensemble_hist > 0,
    (ensemble_fut - ensemble_hist) / ensemble_hist * 100,
    NA
  )

  # --- Model agreement % ---------------------------------------------------
  # Ensemble absolute change (future - hist means), used as the reference
  # sign for agreement. Using absolute change avoids any per-model division.
  ensemble_abs <- ensemble_fut - ensemble_hist
  ensemble_abs <- terra::ifel(is.nan(ensemble_abs), NA, ensemble_abs)

  # Count models agreeing in sign with the ensemble absolute change.
  # sign() returns NA where input is NA, so masked pixels are automatically
  # excluded from both the numerator and denominator.
  agree <- sign(abs_stack) == sign(ensemble_abs)
  n_agree <- sum(terra::ifel(agree, 1, 0), na.rm = TRUE)
  n_valid <- sum(terra::ifel(!is.na(abs_stack), 1, 0), na.rm = TRUE)

  # Express as a percentage of the models with valid data at each pixel.
  # Where n_valid == 0 (no model has data), result is NA not 0/0.
  agreement <- terra::ifel(
    n_valid > 0,
    (n_agree / n_valid) * 100,
    NA
  )
  # Where the ensemble mean itself is NA (no model had valid data at all),
  # agreement must also be NA - not 0 - since "no data" and "0% agreement"
  # are different things and should not look the same on the map.
  agreement <- terra::ifel(is.na(mean_pct), NA, agreement)

  # --- Project and mask to Equal Earth boundary ----------------------------
  mean_pct <- terra::project(mean_pct, "EPSG:8857")
  agreement <- terra::project(agreement, "EPSG:8857")

  ee_mask <- terra::vect(ee_boundary)
  mean_pct <- terra::mask(mean_pct, ee_mask)
  agreement <- terra::mask(agreement, ee_mask)

  list(
    mean_pct = mean_pct,
    agreement = agreement,
    n_models = nrow(meta),
    models = meta$Model
  )
}

#' Common map theme/scaffolding shared by both ensemble-map figures below,
#' matching the visual style already established in 05_Climate_Spatial.R.
base_map_theme <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Helvetica"),
      panel.grid = element_line(color = "gray90", linewidth = 0.2),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_blank(),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width = unit(1, "cm"),
      legend.title = element_text(size = 12, hjust = 0.5),
      legend.title.position = "top",
      legend.text = element_text(size = 12),
      plot.margin = margin(2, 2, 2, 2)
    )
}

#' One panel of the ensemble MEAN % CHANGE grid.
#'
#' @param show_row_label Scenario label placed on the LEFT (used for the
#'   first column only, since rows = Scenario in the new orientation).
#' @param show_col_title Variable title placed on TOP (used for the first
#'   row only, since columns = Variable in the new orientation).
plot_mean_change_panel <- function(
  mean_rast,
  show_row_label = NULL,
  show_col_title = NULL
) {
  p <- ggplot() +
    geom_spatraster(data = mean_rast) +
    geom_sf(data = ee_boundary, color = "gray40", linewidth = 0.3, fill = NA) +
    geom_sf(data = world, color = "grey80", linewidth = 0.2, fill = "grey50") +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      limits = c(-50, 50),
      na.value = "transparent",
      name = "Mean zooplankton ensemble change (%)",
      oob = scales::squish
    ) +
    coord_sf(
      crs = "EPSG:8857",
      xlim = c(-17243959, 17243959),
      ylim = c(-8343134, 8343134),
      expand = FALSE
    ) +
    base_map_theme()

  if (!is.null(show_row_label)) {
    p <- p +
      labs(y = show_row_label) +
      theme(axis.title.y = element_text(size = 12, face = "plain", angle = 90))
  }
  if (!is.null(show_col_title)) {
    p <- p +
      ggtitle(show_col_title) +
      theme(plot.title = element_text(size = 12, face = "plain", hjust = 0.5))
  }

  p
}

#' One panel of the MODEL AGREEMENT % grid.
#'
#' Uses a sequential fill scale from 0 to 100 (%), fixed across every panel
#' and variable — since agreement is now expressed as a percentage of the
#' models with valid data at each pixel, the scale is always 0-100 regardless
#' of ensemble size. This allows a single shared legend for the whole figure.
#'
#' @param show_row_label Scenario label placed on the LEFT (first column
#'   only).
#' @param show_col_title Variable title placed on TOP (first row only).
plot_agreement_panel <- function(
  agreement_rast,
  show_row_label = NULL,
  show_col_title = NULL
) {
  p <- ggplot() +
    geom_spatraster(data = agreement_rast) +
    geom_sf(data = ee_boundary, color = "gray40", linewidth = 0.3, fill = NA) +
    geom_sf(data = world, color = "grey80", linewidth = 0.2, fill = "grey50") +
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      na.value = "transparent",
      name = "Model agreement (%)"
    ) +
    coord_sf(
      crs = "EPSG:8857",
      xlim = c(-17243959, 17243959),
      ylim = c(-8343134, 8343134),
      expand = FALSE
    ) +
    base_map_theme()

  if (!is.null(show_row_label)) {
    p <- p +
      labs(y = show_row_label) +
      theme(axis.title.y = element_text(size = 12, face = "plain", angle = 90))
  }
  if (!is.null(show_col_title)) {
    p <- p +
      ggtitle(show_col_title) +
      theme(
        plot.title = element_text(
          size = 12,
          face = "plain",
          hjust = 0.5,
          vjust = 2
        )
      )
  }

  p
}

# Build every variable x scenario ensemble map ONCE, since both the mean-
# change and the agreement-count figures need the exact same underlying
# per-pixel stacks - avoids re-reading every model's biomassChange.tif
# twice.
ensemble_grid <- tidyr::expand_grid(
  Variable = vars,
  Scenario = scenario_order
) %>%
  dplyr::mutate(
    result = purrr::map2(Variable, Scenario, build_ensemble_maps)
  )

n_var <- length(vars)
n_scn <- length(scenario_order)

#' Build one full variable "column": a 4-row (Scenario) stack of panels for
#' a single variable, for either the mean-change or agreement-% figure.
#'
#' @param variable One of "zooc", "zmicro", "zmeso".
#' @param kind "mean" or "agreement".
#' @param collect_own_legend If FALSE (default for both kinds now), panels
#'   are returned with their individual (identical) legends left INTACT so
#'   that a later `plot_layout(guides = "collect")` applied ONCE across the
#'   full combined grid can merge them into a SINGLE shared legend. Both
#'   figure types now share a fixed scale across all panels (mean: -50 to
#'   50%; agreement: 0 to 100%), so one legend per figure is sufficient.
#' @return A patchwork object: 4 stacked panels.
build_variable_column <- function(
  variable,
  kind = c("mean", "agreement"),
  collect_own_legend = FALSE
) {
  kind <- match.arg(kind)
  is_first_var <- variable == vars[1]

  panels <- vector("list", n_scn)

  for (j in seq_along(scenario_order)) {
    scenario <- scenario_order[j]
    row <- ensemble_grid |>
      dplyr::filter(Variable == variable, Scenario == scenario)
    res <- row$result[[1]]

    row_label <- if (is_first_var) scenario_label_map[[scenario]] else NULL
    col_title <- if (j == 1) var_labels[[variable]] else NULL

    if (is.null(res$mean_pct)) {
      panels[[j]] <- ggplot() +
        annotate(
          "text",
          x = 0,
          y = 0,
          label = "No data",
          size = 3,
          color = "gray50"
        ) +
        theme_void() +
        theme(plot.margin = margin(2, 2, 2, 2))
      next
    }

    panels[[j]] <- if (kind == "mean") {
      plot_mean_change_panel(
        res$mean_pct,
        show_row_label = row_label,
        show_col_title = col_title
      )
    } else {
      plot_agreement_panel(
        res$agreement,
        show_row_label = row_label,
        show_col_title = col_title
      )
    }
  }

  col <- wrap_plots(panels, ncol = 1)

  if (collect_own_legend) {
    col <- col +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "bottom",
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.5, "cm")
      )
  }

  col
}

# Both figure types now use a fixed scale across all panels, so build all
# columns WITHOUT per-column legend collection, then collect ONCE across
# the full combined grid below — patchwork merges identical guides into a
# single shared legend, giving 1 centred legend per figure.
agree_columns <- purrr::map(
  vars,
  build_variable_column,
  kind = "agreement",
  collect_own_legend = FALSE
)
mean_columns <- purrr::map(
  vars,
  build_variable_column,
  kind = "mean",
  collect_own_legend = FALSE
)

mean_change_grid <- wrap_plots(mean_columns, ncol = n_var) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14,
        vjust = 2
      )
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.key.width = unit(4, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

agreement_grid <- wrap_plots(agree_columns, ncol = n_var) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.key.width = unit(4, "cm"),
    legend.key.height = unit(0.6, "cm"),
    legend.title = element_text(size = 12, hjust = 0.5),
    legend.text = element_text(size = 10)
  )

ggsave(
  "Figures/EnsembleMeanPercentChange.pdf",
  plot = mean_change_grid,
  width = 10,
  height = 8
)
ggsave(
  "Figures/EnsembleMeanPercentChange.png",
  plot = mean_change_grid,
  width = 10,
  height = 8,
  dpi = 600
)
ggsave(
  "Figures/ModelAgreementCount.pdf",
  plot = agreement_grid,
  width = 10,
  height = 8
)
ggsave(
  "Figures/ModelAgreementCount.png",
  plot = agreement_grid,
  width = 10,
  height = 8,
  dpi = 600
)
