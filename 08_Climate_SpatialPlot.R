
# ---------------------------------------------------------------------------
# SPATIAL ENSEMBLE MAPS
# ---------------------------------------------------------------------------
#
# Two 4 (scenario: ssp126/245/370/585) x 3 (variable: zooc/zmicro/zmeso)
# patchwork grids, built directly from the per-model/scenario
# "Percent_change" layer already written by 06_Climate_Review.R's
# zoo_change() to Data/<var>/biomassChange/*.tif (via
# /Volumes/T9/ClimateData/<var>/biomassChange/*.tif - the CSVs in Data/ are
# only the tabular summary of these same rasters).
#
#   1. ENSEMBLE MEAN % CHANGE: for each variable/scenario cell, stack every
#      available model's Percent_change layer and take the per-pixel mean
#      (na.rm = TRUE, so a pixel masked as an invalid-baseline artifact in
#      one model - see design note (d) in 06_Climate_Review.R - simply drops
#      out of that pixel's mean rather than poisoning it).
#
#   2. MODEL AGREEMENT COUNT: for the same stack, count how many of the
#      available models agree in SIGN with that pixel's ensemble mean (i.e.
#      point the same direction - increase or decrease). A pixel where the
#      ensemble mean is itself NA (no model has valid data there) is left NA
#      in the agreement map too, rather than showing a misleading 0.
#
# CanESM5/CanESM5-CanOE are now included (unlike the point-range plot above,
# which still excludes CanESM5 - left untouched at the user's request): the
# 06_Climate_Review.R epsilon-baseline fix (EPS_BASELINE) already removed
# the regridding-artifact pixels that were producing implausible values for
# this model, so there is no longer a reason to exclude it here.
#
# LAYOUT: rows = Scenario (4), columns = Variable (3), saved as a tall,
# vertically-oriented page. This lets each COLUMN (i.e. each variable) share
# ONE legend at the BOTTOM of that column, rather than one legend per panel
# - important for the agreement-count figure, where each variable has a
# different total model count (n_models) and therefore a genuinely
# different colour-scale range; collecting a single legend across all
# columns would misrepresent the smaller-ensemble variables, but collecting
# per-column (one legend for all 4 scenario rows of e.g. zooc) is safe
# because n_models is constant WITHIN a variable across scenarios... except
# it isn't (zooc varies 11-13 models per scenario - see build_ensemble_maps()
# below), so each column instead keeps its own fixed 0..max_n_models(variable)
# scale and one collected legend per column via patchwork::plot_layout in
# the per-variable column-builder below.
# ---------------------------------------------------------------------------

library(terra)
library(tidyterra)
library(sf)
library(patchwork)

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

# Coastline for context, matching the style already used in
# 05_Climate_Spatial.R.
world <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") %>%
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
      cbind(seq(-180,  180, length.out = n), rep( 90, n)),
      cbind(rep( 180, n), seq( 90, -90, length.out = n)),
      cbind(seq( 180, -180, length.out = n), rep(-90, n)),
      cbind(rep(-180, n), seq(-90,  90, length.out = n)),
      c(-180, 90)
    ))),
    crs = 4326
  ) |> sf::st_transform("EPSG:8857")
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

#' Build the per-pixel model-ensemble mean % change AND the per-pixel model
#' agreement count for one variable/scenario combination.
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
    return(list(mean_pct = NULL, agreement = NULL, n_models = 0L, models = character(0)))
  }
  
  # Load only the Percent_change layer from each model's 2-layer
  # (Absolute_change, Percent_change) biomassChange.tif, one at a time -
  # cheap here (no parallelisation needed, unlike integrate_biomass() in
  # 06_Climate_Review.R) since these are already small, pre-computed maps.
  pct_stack <- terra::rast(
    purrr::map(meta$Path, function(p) terra::rast(p)[["Percent_change"]])
  )
  
  names(pct_stack) <- meta$Model
  
  pct_stack <- terra::project(pct_stack, "EPSG:8857")
  
  # na.rm = TRUE: a pixel masked to NA in one model (invalid baseline - see
  # design note d in 06_Climate_Review.R) simply drops out of that pixel's
  # ensemble mean rather than propagating NA to every model's agreement.
  mean_pct <- mean(pct_stack, na.rm = TRUE)
  # mean() of an all-NA pixel across every layer returns NaN, not NA
  # (verified) - normalise to NA so downstream NA-handling (mask
  # comparisons, plotting `na.value`) behaves consistently.
  mean_pct <- terra::ifel(is.nan(mean_pct), NA, mean_pct)
  
  # Model agreement: count how many models share the SAME SIGN as the
  # ensemble mean at that pixel. sign() returns NA where the input is NA,
  # so a model with no valid data at a pixel (its own baseline masked)
  # automatically does not count for OR against agreement there - it is
  # simply excluded, exactly like it is excluded from the mean above.
  agree <- sign(pct_stack) == sign(mean_pct)
  agreement <- sum(terra::ifel(agree, 1, 0), na.rm = TRUE)
  # Where the ensemble mean itself is NA (no model had valid data at all),
  # agreement must also be NA - not 0 - since "zero models agree" and "no
  # data available" are different things and should not look the same on
  # the map.
  agreement <- terra::ifel(is.na(mean_pct), NA, agreement)
  
  # Mask both outputs to the Equal Earth boundary polygon so that pixels
  # in the rectangular projected extent that fall OUTSIDE the elliptical
  # Equal Earth globe boundary are set to NA (transparent) rather than
  # showing as filled corners on the plot.
  ee_mask <- terra::vect(ee_boundary)
  mean_pct  <- terra::mask(mean_pct,  ee_mask)
  agreement <- terra::mask(agreement, ee_mask)
  
  list(
    mean_pct  = mean_pct,
    agreement = agreement,
    n_models  = nrow(meta),
    models    = meta$Model
  )
}

#' Common map theme/scaffolding shared by both ensemble-map figures below,
#' matching the visual style already established in 05_Climate_Spatial.R.
base_map_theme <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
  theme(
    panel.grid = element_line(color = "gray90", linewidth = 0.2),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_blank(),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(2, 2, 2, 2)
  )
}

#' One panel of the ensemble MEAN % CHANGE grid.
#'
#' @param show_row_label Scenario label placed on the LEFT (used for the
#'   first column only, since rows = Scenario in the new orientation).
#' @param show_col_title Variable title placed on TOP (used for the first
#'   row only, since columns = Variable in the new orientation).
plot_mean_change_panel <- function(mean_rast, show_row_label = NULL, show_col_title = NULL) {
  p <- ggplot() +
  geom_spatraster(data = mean_rast) +
  geom_sf(data = world, color = "grey30", linewidth = 0.2, fill = "grey30") +
  geom_sf(data = ee_boundary, color = "gray40", linewidth = 0.3, fill = NA) +
  scale_fill_gradientn(
    colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
    limits = c(-50, 50),
    na.value = "transparent",
    name = "Mean %\nchange",
    oob = scales::squish
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

#' One panel of the MODEL AGREEMENT COUNT grid.
#'
#' Uses a sequential (not diverging) fill scale from 0 to n_models_max, since
#' an agreement count has no meaningful "zero-centred" interpretation the
#' way a signed % change does - it is a pure count. `n_models_max` is the
#' MAXIMUM model count across all scenarios FOR THIS VARIABLE'S COLUMN (baked
#' in by build_variable_column() below), so every panel within a column
#' shares one consistent scale/legend - each variable still keeps its own
#' honestly-scaled column, since zooc has more models than zmicro/zmeso.
#'
#' @param show_row_label Scenario label placed on the LEFT (first column
#'   only).
#' @param show_col_title Variable title placed on TOP (first row only).
plot_agreement_panel <- function(agreement_rast, n_models_max, show_row_label = NULL, show_col_title = NULL) {
  p <- ggplot() +
  geom_spatraster(data = agreement_rast) +
  geom_sf(data = world, color = "grey30", linewidth = 0.2, fill = "grey30") +
  geom_sf(data = ee_boundary, color = "gray40", linewidth = 0.3, fill = NA) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, n_models_max),
    breaks = scales::breaks_pretty(n = min(5, n_models_max + 1)),
    na.value = "transparent",
    name = paste0("Models agreeing\n(of ", n_models_max, ")")
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

# Build every variable x scenario ensemble map ONCE, since both the mean-
# change and the agreement-count figures need the exact same underlying
# per-pixel stacks - avoids re-reading every model's biomassChange.tif
# twice.
ensemble_grid <- tidyr::expand_grid(Variable = vars, Scenario = scenario_order) |>
dplyr::mutate(
  result = purrr::map2(Variable, Scenario, build_ensemble_maps)
)

n_var <- length(vars)
n_scn <- length(scenario_order)

# Per-variable maximum model count across its 4 scenarios - used to give
# every agreement-count panel WITHIN a variable's column a single, honest,
# shared 0..max scale (n_models legitimately varies scenario-to-scenario for
# zooc - 11 to 13 models - since not every model/variant submitted every
# SSP; using the column's max keeps the legend consistent down the column
# while still reflecting that variable's true ensemble size).
n_models_max_by_var <- ensemble_grid |>
dplyr::mutate(n_models = purrr::map_int(result, ~ .x$n_models)) |>
dplyr::group_by(Variable) |>
dplyr::summarise(n_models_max = max(n_models), .groups = "drop") |>
dplyr::pull(n_models_max, name = Variable)

#' Build one full variable "column": a 4-row (Scenario) stack of panels for
#' a single variable, for either the mean-change or agreement-count figure.
#'
#' @param variable One of "zooc", "zmicro", "zmeso".
#' @param kind "mean" or "agreement".
#' @param collect_own_legend If TRUE (default), collect this column's
#'   panels into ONE legend at the BOTTOM of the column (scoped via `&` so
#'   it does not merge with other columns' legends) - used for the
#'   agreement-count figure, where each variable genuinely has a different
#'   colour-scale range (n_models_max_by_var differs by variable), so each
#'   column needs its own legend. If FALSE, panels are returned with their
#'   individual (identical) legends left INTACT, so that a later
#'   `plot_layout(guides = "collect")` applied ONCE across the full
#'   combined grid (all 3 columns) can merge them into a SINGLE shared
#'   legend - used for the mean-% -change figure, where every panel shares
#'   the exact same fixed scale (limits = c(-50, 50)), so one legend for
#'   the whole figure is sufficient and avoids repeating it 3 times.
#' @return A patchwork object: 4 stacked panels (+ 1 shared bottom legend
#'   if `collect_own_legend = TRUE`).
build_variable_column <- function(variable, kind = c("mean", "agreement"), collect_own_legend = TRUE) {
  kind <- match.arg(kind)
  is_first_var <- variable == vars[1]
  
  panels <- vector("list", n_scn)
  
  for (j in seq_along(scenario_order)) {
    scenario <- scenario_order[j]
    row <- ensemble_grid |>
    dplyr::filter(Variable == variable, Scenario == scenario)
    res <- row$result[[1]]
    
    row_label <- scenario_label_map[[scenario]]
    col_title <- if (j == 1) var_labels[[variable]] else NULL
    
    if (is.null(res$mean_pct)) {
      panels[[j]] <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No data", size = 3, color = "gray50") +
      theme_void() +
      theme(plot.margin = margin(2, 2, 2, 2))
      next
    }
    
    panels[[j]] <- if (kind == "mean") {
      plot_mean_change_panel(res$mean_pct, show_row_label = row_label, show_col_title = col_title)
    } else {
      plot_agreement_panel(res$agreement, n_models_max_by_var[[variable]],
        show_row_label = row_label, show_col_title = col_title)
      }
    }
    
    col <- wrap_plots(panels, ncol = 1)
    
    if (collect_own_legend) {
      col <- col +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom",
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.3, "cm"))
    }
    
    col
  }
  
  # Agreement-count: one column per variable, each with its OWN bottom-
  # collected legend (different n_models_max_by_var per variable), placed
  # side-by-side to form the full 4 (scenario) x 3 (variable) grid.
  agree_columns <- purrr::map(vars, build_variable_column, kind = "agreement", collect_own_legend = TRUE)
  
  # Mean % change: every panel shares the exact same fixed scale
  # (limits = c(-50, 50)), so build the 3 columns WITHOUT collecting a
  # legend per column, then collect ONCE across the whole combined grid
  # below - patchwork merges identical guides into a single shared legend,
  # giving 1 centred legend for the whole figure instead of 3 repeats.
  mean_columns <- purrr::map(vars, build_variable_column, kind = "mean", collect_own_legend = FALSE)
  
  mean_change_grid <- wrap_plots(mean_columns, ncol = n_var) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Ensemble Mean % Change in Zooplankton Biomass (1993-2014 to 2081-2100)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ) &
  theme(legend.position = "bottom",
  legend.key.width = unit(4, "cm"),
  legend.key.height = unit(0.4, "cm"),
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8))
  
  agreement_grid <- wrap_plots(agree_columns, ncol = n_var) +
  plot_annotation(
    title = "Number of Models Agreeing on Direction of Change",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  ) &
  theme(
    legend.title.position = "bottom",
    legend.title = element_text(size = 9, hjust = 0.5),
    legend.text  = element_text(size = 8)
  )
  
  ggsave("Figures/EnsembleMeanPercentChange.pdf", plot = mean_change_grid, width = 10, height = 8)
  ggsave("Figures/EnsembleMeanPercentChange.png", plot = mean_change_grid, width = 10, height = 8, dpi = 600)
  ggsave("Figures/ModelAgreementCount.pdf", plot = agreement_grid, width = 10, height = 8)
  ggsave("Figures/ModelAgreementCount.png", plot = agreement_grid, width = 10, height = 8, dpi = 600)
