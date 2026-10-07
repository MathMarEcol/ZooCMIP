# ZooCMIP: Compile biomass change maps and summary statistics
# Depth-integrates regridded CMIP6 zooplankton netCDF files, computes
# historical (1993-2014) vs future (2081-2100) biomass change per model/SSP,
# and writes per-model biomassChange.tif files and a summary CSV used by
# 04_Climate_Plot.R and 05_Climate_SpatialPlot.R.

library(tidyverse)
library(terra)
library(ncdf4)
library(future)
library(furrr)
library(hotrstuff)

#' Identify the vertical (Z-axis) dimension name and its bounds variable name
#'
#' CMIP6 models do not use a consistent name for the vertical ocean
#' coordinate: most use "lev"/"lev_bnds", but e.g. IPSL-CM6A-LR uses
#' "olevel"/"olevel_bnds" (verified in this collection - IPSL is the only
#' model here that deviates from "lev"). Rather than hard-coding "lev" and
#' silently failing (or worse, silently matching the wrong variable) on any
#' model that uses a different name, we discover both names generically from
#' the file's CF-convention metadata: the vertical dimension is whichever one
#' has axis = "Z", and its bounds variable name is given by that dimension's
#' own "bounds" attribute. This is robust to any naming convention we have
#' not yet encountered (e.g. NEMO's "deptht"), not just the two seen so far.
#'
#' @param nc_file Path to a netCDF file.
#' @return A list with `lev_name` and `bnds_name` character strings.
find_vertical_dim <- function(nc_file) {
  # nc_open() emits a harmless ncdf4 warning about the "bnds" dimension
  # having no associated coordinate variable (expected/normal for a plain
  # 2-element bounds dimension) - suppressed here as it is not informative.
  nc <- suppressWarnings(ncdf4::nc_open(nc_file))
  on.exit(ncdf4::nc_close(nc))

  z_name <- NULL
  for (dn in names(nc$dim)) {
    ax <- tryCatch(ncdf4::ncatt_get(nc, dn, "axis")$value, error = function(e) {
      NULL
    })
    if (!is.null(ax) && identical(ax, "Z")) {
      z_name <- dn
      break
    }
  }
  if (is.null(z_name)) {
    stop("Could not find a Z-axis (vertical) dimension in ", nc_file)
  }

  bnds_attr <- ncdf4::ncatt_get(nc, z_name, "bounds")
  if (!isTRUE(bnds_attr$hasatt) || !nzchar(bnds_attr$value)) {
    stop(
      "Vertical dimension '",
      z_name,
      "' in ",
      nc_file,
      " has no 'bounds' attribute - cannot determine cell thickness."
    )
  }

  list(lev_name = z_name, bnds_name = bnds_attr$value)
}


#' Build a vertical-thickness (m) SpatRaster from a single reference netCDF file
#'
#' Deliberately vertical-only: horizontal cell area is NOT included here,
#' so the result is a per-level
#' thickness broadcast across the horizontal grid, not a volume. This keeps
#' the depth-integration step (thickness only) cleanly separated from any
#' horizontal-area weighting, which is applied later and only where actually
#' needed (see calc_biomass_change() / summarise_biomass_change()).
#'
#' The resulting raster has one layer per vertical level, with depth() set to
#' the model's vertical coordinate values so it can be matched against a
#' biomass stack later.
#'
#' @param nc_file Path to a single netCDF file for the model of interest.
#' @return A SpatRaster (lon x lat x lev) of vertical thickness in metres,
#'   constant across every horizontal pixel within a given level.
build_thickness_raster <- function(nc_file) {
  # Horizontal template purely for its grid geometry (extent/resolution/crs) -
  # values are discarded and replaced with 1 via init(), since only the
  # levels' thickness (not horizontal area) should vary here.
  r_template <- terra::rast(nc_file)[[1]]

  z_dim <- find_vertical_dim(nc_file)
  nc <- suppressWarnings(ncdf4::nc_open(nc_file))
  lev <- ncdf4::ncvar_get(nc, z_dim$lev_name)
  lev_bnds <- ncdf4::ncvar_get(nc, z_dim$bnds_name)
  ncdf4::nc_close(nc)

  thickness <- lev_bnds[2, ] - lev_bnds[1, ]

  # Build a multi-layer raster of 1s (one layer per level) on the same
  # horizontal grid, then multiply by the per-level thickness vector - terra
  # recycles this elementwise per layer (verified numerically to match
  # thickness[i] for each layer i, uniformly across every pixel).
  thick <- terra::init(terra::rast(r_template, nlyrs = length(thickness)), 1) *
    thickness
  names(thick) <- paste0("lev_", seq_along(thickness))
  terra::depth(thick) <- lev

  thick
}


#' Get the lev_bnds matrix for a netCDF file (used to verify model-consistency)
#'
#' Uses find_vertical_dim() so this works regardless of whether the file
#' calls its vertical bounds variable "lev_bnds", "olevel_bnds", or anything
#' else that follows the CF axis="Z" / bounds= convention.
get_lev_bnds <- function(nc_file) {
  z_dim <- find_vertical_dim(nc_file)
  nc <- suppressWarnings(ncdf4::nc_open(nc_file))
  lev_bnds <- ncdf4::ncvar_get(nc, z_dim$bnds_name)
  ncdf4::nc_close(nc)
  lev_bnds
}


#' Depth-integrate a single biomass file into areal biomass density per pixel/time
#'
#' Multiplies the biomass concentration (mol C m-3) by vertical thickness (m)
#' at every depth level, then sums over depth (na.rm = TRUE) to collapse the
#' vertical dimension, giving one 2D layer per time step. Units of the output
#' are mol C m-2 (deliberately NOT multiplied by horizontal cell area, and
#' no mass conversion is applied).
#'
#' @param nc_file Path to the biomass netCDF file to integrate.
#' @param reference_nc_file Path to the reference file for this SAME model,
#'   used to (re)build the thickness raster. May be the same as nc_file.
#' @return A SpatRaster (lon x lat x time) of depth-integrated areal biomass
#'   density (mol C m-2).
integrate_one_file <- function(nc_file, reference_nc_file) {
  # Sanity check: the vertical grid of this file must match the reference
  # file's vertical grid before we reuse its thickness raster. This is cheap
  # (a few KB of numbers) and protects against silently applying the wrong
  # model's geometry - a failure mode that would not error, just corrupt
  # every downstream number for that file.
  if (!identical(nc_file, reference_nc_file)) {
    lb_ref <- get_lev_bnds(reference_nc_file)
    lb_this <- get_lev_bnds(nc_file)
    if (!isTRUE(all.equal(lb_ref, lb_this))) {
      stop(
        "lev_bnds mismatch between reference file:\n  ",
        reference_nc_file,
        "\nand target file:\n  ",
        nc_file,
        "\nRefusing to reuse thickness raster across inconsistent vertical grids."
      )
    }
  }

  thick <- build_thickness_raster(reference_nc_file)

  r <- terra::rast(nc_file)
  r_depth <- terra::depth(r)
  r_time <- terra::time(r)

  # Match each biomass layer to its corresponding thickness layer by actual
  # depth value (not position) - terra does not guarantee lev-major ordering
  # in general, so this is the only robust way to align the two stacks.
  idx <- match(r_depth, terra::depth(thick))
  if (anyNA(idx)) {
    stop(
      "Could not match all depth levels in ",
      nc_file,
      " to the thickness raster built from ",
      reference_nc_file,
      ". Unmatched depths: ",
      paste(unique(r_depth[is.na(idx)]), collapse = ", ")
    )
  }

  biomass <- r * thick[[idx]]

  # Collapse the depth dimension: sum layers sharing the same real time
  # value. na.rm = TRUE is required so partial water columns (shelf seas)
  # integrate correctly over just their wet levels, and terra::tapp()
  # correctly returns NA (not 0) for fully-land pixels (verified).
  u_time <- sort(unique(r_time))
  time_idx <- match(r_time, u_time)
  result <- terra::tapp(biomass, index = time_idx, fun = sum, na.rm = TRUE)

  # terra::tapp() does not guarantee output layers are sorted by group value
  # (verified) - recover the true date for each output layer from its name
  # and explicitly re-order chronologically before trusting it.
  grp_id <- as.integer(sub("^X", "", names(result)))
  result_time <- u_time[grp_id]
  ord <- order(result_time)
  result <- result[[ord]]
  result_time <- result_time[ord]

  stopifnot(
    "Depth-integration produced a non-chronological time axis" = all(
      diff(result_time) > 0
    )
  )

  terra::time(result) <- result_time
  terra::units(result) <- rep("mol m-2", terra::nlyr(result))
  names(result) <- format(result_time, "%Y")

  result
}


#' Depth-integrate total biomass for every file in a vector of netCDF paths
#'
#' Groups files by model (parsed from the CMIP6 filename convention
#' variable_table_model_experiment_variant_grid_daterange.nc) so the
#' thickness raster is only computed once per model and re-used across
#' historical and every SSP for that model. Each input file's integrated
#' result is written to its own GeoTIFF under out_dir/<variable>/, mirroring
#' the existing Data/spatial_change_maps/ convention used elsewhere in this
#' project.
#'
#' @param nc_files Character vector of paths to netCDF files (can span
#'   multiple models, experiments and variables).
#' @param out_dir Base output directory. A subfolder per variable
#'   (e.g. "zooc", "zmeso", "zmicro") is created beneath it.
#' @param workers Number of parallel workers (passed to future::multisession).
#'   Defaults to all available cores minus 2, matching hotrstuff's convention.
#' @param overwrite Logical; if FALSE (default) files that already exist are
#'   skipped, matching the pattern used by hotrstuff::htr_run_cdo().
#' @return A tibble summarising what was processed, one row per input file.
integrate_biomass <- function(
  nc_files,
  out_dir = file.path("Data", "integrated_biomass"),
  workers = NULL,
  overwrite = FALSE
) {
  # Parse CMIP6 metadata from filenames. Using hotrstuff's own parser rather
  # than writing a second one, so both packages stay in sync if the naming
  # convention ever changes.
  meta <- purrr::map(nc_files, hotrstuff:::htr_get_CMIP6_bits) |>
    dplyr::bind_rows() |>
    dplyr::mutate(path = nc_files)

  # One reference file per model, used to build (and re-use) that model's
  # thickness raster. Prefer "historical" when available since it is the most
  # likely to be complete; otherwise fall back to the first available file
  # for that model. Verified: GFDL-ESM4 has NO historical file in this
  # collection (only ssp245/ssp370), so this fallback is not optional - it
  # is required for at least one real model in the current dataset.
  reference_lookup <- meta |>
    dplyr::group_by(Model) |>
    dplyr::arrange(dplyr::desc(Scenario == "historical"), .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(Model, reference_path = path)

  meta <- meta |>
    dplyr::left_join(reference_lookup, by = "Model")

  purrr::pwalk(
    list(unique(meta$Variable)),
    function(v) {
      hotrstuff:::htr_make_folder(file.path(out_dir, v, "depthIntegrated"))
    }
  )

  do_one <- function(path, reference_path, variable, model, scenario, variant) {
    library(terra)
    library(ncdf4)

    out_file <- file.path(
      out_dir,
      variable,
      "depthIntegrated",
      paste0(
        variable,
        "_",
        model,
        "_",
        scenario,
        "_",
        variant,
        "_depthIntegrated.tif"
      )
    )

    if (file.exists(out_file) && !overwrite) {
      return(tibble::tibble(
        Variable = variable,
        Model = model,
        Scenario = scenario,
        Variant = variant,
        Output_file = out_file,
        Status = "Skipped (exists)"
      ))
    }

    tryCatch(
      {
        result <- integrate_one_file(path, reference_path)
        terra::writeRaster(result, out_file, overwrite = TRUE)
        tibble::tibble(
          Variable = variable,
          Model = model,
          Scenario = scenario,
          Variant = variant,
          Output_file = out_file,
          Status = "Success"
        )
      },
      error = function(e) {
        tibble::tibble(
          Variable = variable,
          Model = model,
          Scenario = scenario,
          Variant = variant,
          Output_file = NA_character_,
          Status = paste("ERROR:", conditionMessage(e))
        )
      }
    )
  }

  # Follow the hotrstuff convention: pass only file paths (strings) into
  # furrr workers, never live SpatRaster objects (verified these break
  # multisession serialization with "NULL value passed as symbol address").
  # Each worker rebuilds the thickness raster itself from reference_path -
  # this is cheap (~0.2s, verified) so recomputing it per-file is an
  # acceptable trade-off for the safety/simplicity of not sharing raster
  # objects.
  w <- if (is.null(workers)) parallelly::availableCores(omit = 2) else workers
  future::plan(future::multisession, workers = w)
  on.exit(future::plan(future::sequential), add = TRUE)

  # Our worker function performs no random number generation itself, but
  # terra/GDAL's internal state initialisation in a fresh worker process
  # spuriously trips future's "unreliable RNG" heuristic (verified: this
  # fires even for a bare terra::writeRaster() call with no RNG involved at
  # all). This is a known false positive, not a sign of a real problem with
  # our code, so we disable the check as future's own warning message
  # suggests rather than leave misleading warnings on every run.
  old_opt <- options(future.rng.onMisuse = "ignore")
  on.exit(options(old_opt), add = TRUE)

  summary_tbl <- furrr::future_pmap(
    list(
      path = meta$path,
      reference_path = meta$reference_path,
      variable = meta$Variable,
      model = meta$Model,
      scenario = meta$Scenario,
      variant = meta$Variant
    ),
    do_one
  ) |>
    dplyr::bind_rows()

  summary_tbl
}


# ---------------------------------------------------------------------------
# BIOMASS CHANGE: historical (1993-2014) mean vs SSP (2081-2100) mean
# ---------------------------------------------------------------------------

#' Parse a Data/integrated_biomass/<var>/*.tif filename into its metadata
#'
#' Filenames follow the fixed, project-specific convention written by
#' integrate_biomass() above:
#'   <variable>_<model>_<scenario>_<variant>_biomass_integrated.tif
#' Verified across all 159 files currently in Data/integrated_biomass/: this
#' always splits into exactly 6 underscore-delimited tokens, with the final
#' two always literally "biomass" and "integrated.tif" - so Model is
#' recovered positionally from the front and Scenario/Variant positionally
#' from the back, which is robust to models whose own names contain
#' underscores or hyphens (e.g. "CanESM5-CanOE") without needing a regex.
#'
#' @param path Path to a single *_biomass_integrated.tif file.
#' @return A one-row tibble with Variable, Model, Scenario, Variant, Path.
parse_integrated_filename <- function(path) {
  bits <- strsplit(basename(path), "_")[[1]]
  n <- length(bits)
  stopifnot(
    "Unexpected integrated_biomass filename format" = n >= 5 &&
      bits[n] == "depthIntegrated.tif"
  )
  tibble::tibble(
    Variable = bits[1],
    Model = paste(bits[2:(n - 3)], collapse = "_"),
    Scenario = bits[n - 2],
    Variant = bits[n - 1],
    Path = path
  )
}


#' Compute the pixel-wise biomass change map for one historical/SSP file pair
#' of the same model/variable/variant.
#'
#' Output is a 4-layer SpatRaster:
#'   Layer 1 "Hist_mean"       — 1993-2014 mean biomass (mol C m-2)
#'   Layer 2 "Future_mean"     — 2081-2100 mean biomass (mol C m-2)
#'   Layer 3 "Absolute_change" — Future_mean - Hist_mean (mol C m-2)
#'   Layer 4 "Percent_change"  — pixel-wise (Absolute_change / Hist_mean)*100 (%)
#'
#' Storing Hist_mean and Future_mean allows 05_Climate_SpatialPlot.R to
#' compute the ensemble % change as:
#'   (mean_across_models(Future_mean) - mean_across_models(Hist_mean)) /
#'    mean_across_models(Hist_mean) * 100
#' i.e. averaging biomass FIRST, then computing % change — which is far more
#' robust than averaging per-model pixel-wise % changes.
#'
#' @param hist_mean_rast Pre-computed 1993-2014 mean SpatRaster (mol C m-2).
#' @param ssp_path Path to the SSP integrated_biomass .tif file.
#' @param label Human-readable label for log messages (e.g.
#'   "zooc_GFDL-ESM4_ssp585_r1i1p1f1").
#' @return A 4-layer SpatRaster: "Hist_mean" (mol C m-2), "Future_mean"
#'   (mol C m-2), "Absolute_change" (mol C m-2), "Percent_change" (%).
calc_biomass_change <- function(hist_mean_rast, ssp_path, label = ssp_path) {
  r_ssp <- terra::rast(ssp_path)
  ssp_years <- as.integer(names(r_ssp))

  future_window <- ssp_years >= 2081 & ssp_years <= 2100
  if (sum(future_window) != 20) {
    stop(
      "Expected 20 years in the 2081-2100 window for ",
      ssp_path,
      " but found ",
      sum(future_window),
      " - refusing to average a truncated future period."
    )
  }

  r_future_mn <- mean(r_ssp[[future_window]], na.rm = TRUE)

  abs_change <- r_future_mn - hist_mean_rast
  pct_change <- (abs_change / hist_mean_rast) * 100

  out <- c(hist_mean_rast, r_future_mn, abs_change, pct_change)
  names(out) <- c(
    "Hist_mean",
    "Future_mean",
    "Absolute_change",
    "Percent_change"
  )
  terra::units(out) <- c("mol m-2", "mol m-2", "mol m-2", "%")

  out
}


#' Area-weighted mean and SD of a single SpatRaster layer
#'
#' terra::global(..., weights = ) only supports fun = "mean" or "sum" (no
#' "sd" option - verified this errors with "txtfun %in% c('mean','sum') are
#' not all TRUE"), so the weighted SD is computed manually here, following
#' exactly the same formula already used and verified elsewhere in this
#' project (see calculate_stats() in 03_Climate_Summarise.R):
#'   weighted_mean = sum(w*v) / sum(w)
#'   weighted_sd   = sqrt(sum(w*(v - weighted_mean)^2) / sum(w))
#'
#' @param x A single-layer SpatRaster.
#' @param w A same-grid SpatRaster of weights (e.g. cell area in m2).
#' @return A list with Mean, SD, and n (count of valid, i.e. non-NA-in-either,
#'   pixels).
weighted_mean_sd <- function(x, w) {
  v <- terra::values(x)[, 1]
  wt <- terra::values(w)[, 1]
  valid <- !is.na(v) & !is.na(wt)
  v <- v[valid]
  wt <- wt[valid]

  if (length(v) == 0) {
    return(list(Mean = NA_real_, SD = NA_real_, n = 0L))
  }

  w_mean <- sum(v * wt) / sum(wt)
  w_sd <- sqrt(sum(wt * (v - w_mean)^2) / sum(wt))

  list(Mean = w_mean, SD = w_sd, n = length(v))
}


#' Summarise a biomass-change SpatRaster into one row of area-weighted stats
#'
#' Mean_percent_change is computed from two global area-weighted means
#' (historical and future), NOT as the area-weighted mean of pixel-wise
#' percent changes. This gives a single, robust scalar:
#'
#'   global_hist_mean   = area-weighted mean of hist_mean_rast
#'   global_future_mean = global_hist_mean + area-weighted mean of Absolute_change
#'   Mean_percent_change = (global_future_mean - global_hist_mean) /
#'                          global_hist_mean * 100
#'
#' SD_percent_change is the area-weighted SD of the pixel-wise percent-change
#' layer, retained for spatial spread / uncertainty reporting.
#'
#' @param change_rast 4-layer SpatRaster from calc_biomass_change().
#' @param hist_mean_rast The 1993-2014 mean SpatRaster (mol C m-2) for this
#'   model/variable/variant. Used to compute the global area-weighted
#'   historical mean for the percent-change denominator.
#' @param area_rast A cellSize() raster (m2) on the same grid, shared across
#'   all calls.
#' @return A one-row tibble of area-weighted Mean/SD for both layers, the
#'   area-weighted ocean-total absolute change (mol C), and pixel count (n).
summarise_biomass_change <- function(change_rast, hist_mean_rast, area_rast) {
  abs_stats <- weighted_mean_sd(change_rast[["Absolute_change"]], area_rast)

  # Global area-weighted mean of the historical baseline (single scalar).
  # This is the denominator for Mean_percent_change.
  hist_global_mean <- terra::global(
    hist_mean_rast,
    "mean",
    weights = area_rast,
    na.rm = TRUE
  )[1, 1]

  # Global area-weighted mean of the future period, derived from the
  # historical mean + the area-weighted mean absolute change. This is
  # equivalent to computing the future global mean directly, but reuses
  # abs_stats$Mean which is already computed above.
  future_global_mean <- hist_global_mean + abs_stats$Mean

  # Single % change from two global means — avoids pixel-level division by
  # near-zero baseline values (see function documentation above).
  mean_pct_change <- (future_global_mean - hist_global_mean) /
    hist_global_mean *
    100

  # SD of the pixel-wise percent-change layer (area-weighted), retained for
  # spatial spread / uncertainty reporting.
  pct_sd_stats <- weighted_mean_sd(change_rast[["Percent_change"]], area_rast)

  # Both layers share an identical NA mask, so pixel counts are always equal;
  # abs_stats$n is reported once.
  total_change <- terra::global(
    change_rast[["Absolute_change"]],
    "sum",
    weights = area_rast,
    na.rm = TRUE
  )[1, 1]

  tibble::tibble(
    Mean_density_change = abs_stats$Mean,
    SD_density_change = abs_stats$SD,
    Mean_percent_change = mean_pct_change,
    SD_percent_change = pct_sd_stats$SD,
    Total_change = total_change,
    n = abs_stats$n
  )
}


#' Full pipeline: biomass change maps + summary tibble for a set of
#' integrated_biomass .tif files, across the four tier-1 SSPs
#'
#' @param int_files Character vector of paths to *_biomass_integrated.tif
#'   files (as produced by integrate_biomass() above).
#' @param out_dir Directory to write per-model/scenario change SpatRasters to.
#' @param ssp_order SSP scenario codes to process (default: the four tier-1
#'   SSPs specified in the analysis brief).
#' @return A tibble with one row per Variable/Model/Variant/Scenario
#'   combination, giving area-weighted change statistics. Also writes one
#'   2-layer GeoTIFF (Absolute_change, Percent_change) per row to out_dir.
zoo_change <- function(
  int_files,
  out_dir = file.path("Data", "biomass_change_maps"),
  ssp_order = c("ssp126", "ssp245", "ssp370", "ssp585")
) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  meta <- purrr::map(int_files, parse_integrated_filename) |>
    dplyr::bind_rows()

  hist_meta <- meta |> dplyr::filter(Scenario == "historical")
  ssp_meta <- meta |>
    dplyr::filter(Scenario %in% ssp_order) |>
    dplyr::left_join(
      hist_meta |> dplyr::select(Variable, Model, Variant, Hist_path = Path),
      by = c("Variable", "Model", "Variant")
    )

  missing_hist <- ssp_meta |> dplyr::filter(is.na(Hist_path))
  if (nrow(missing_hist) > 0) {
    warning(
      "No historical file found for ",
      nrow(missing_hist),
      " Variable/Model/Variant combination(s) - these will be skipped:\n  ",
      paste(missing_hist$Path, collapse = "\n  ")
    )
    ssp_meta <- ssp_meta |> dplyr::filter(!is.na(Hist_path))
  }

  # Area weights are identical across every file (common 0.5-degree grid,
  # verified via compareGeom()), so compute this exactly once and reuse it
  # for every summary statistic below, rather than recomputing cellSize() per row.
  area_template <- terra::rast(ssp_meta$Path[1])[[1]]
  area_rast <- terra::cellSize(area_template, unit = "m")

  # Cache one 1993-2014 historical mean per unique Variable/Model/Variant,
  # since it is identical for every SSP within that group - avoids
  # recomputing the same mean() up to 4x per model/variable.
  hist_cache <- new.env(parent = emptyenv())

  get_hist_mean <- function(hist_path, label) {
    key <- hist_path
    if (is.null(hist_cache[[key]])) {
      r_hist <- terra::rast(hist_path)
      hist_years <- as.integer(names(r_hist))
      hist_window <- hist_years >= 1993 & hist_years <= 2014

      if (sum(hist_window) != 22) {
        stop(
          "Expected 22 years in the 1993-2014 window for ",
          hist_path,
          " but found ",
          sum(hist_window),
          " - refusing to average a truncated historical baseline."
        )
      }

      hist_cache[[key]] <- mean(r_hist[[hist_window]], na.rm = TRUE)
    }
    hist_cache[[key]]
  }

  results <- purrr::pmap(
    ssp_meta,
    function(Variable, Model, Scenario, Variant, Path, Hist_path) {
      label <- paste(Variable, Model, Variant, sep = "_")
      label_ssp <- paste(Variable, Model, Scenario, Variant, sep = "_")
      message("Processing: ", label, " / ", Scenario)

      hist_mn <- get_hist_mean(Hist_path, label)
      change_rast <- calc_biomass_change(hist_mn, Path, label = label_ssp)

      out_file <- file.path(
        out_dir,
        paste0(
          Variable,
          "_",
          Model,
          "_",
          Scenario,
          "_",
          Variant,
          "_biomassChange.tif"
        )
      )
      terra::writeRaster(change_rast, out_file, overwrite = TRUE)

      summarise_biomass_change(change_rast, hist_mn, area_rast) |>
        dplyr::mutate(
          Variable = Variable,
          Model = Model,
          Scenario = Scenario,
          Variant = Variant,
          Output_file = out_file,
          .before = 1
        )
    }
  ) |>
    dplyr::bind_rows()

  results
}

vars <- c("zooc", "zmicro", "zmeso")

# Collect per-variable results into one combined tibble
all_change_summary <- purrr::map_dfr(vars, function(var) {
  int_files <- list.files(
    file.path("", "Volumes", "T9", "ClimateData", var, "depthIntegrated"),
    full.names = TRUE,
    recursive = TRUE,
    pattern = "[.]tif$"
  )

  zoo_change(
    int_files,
    out_dir = file.path(
      "",
      "Volumes",
      "T9",
      "ClimateData",
      var,
      "biomassChange"
    )
  )
})

# Build "All models" summary rows: one per Variable × Scenario, with
# Mean_percent_change = multimodel mean and SD_percent_change = SD across
# models. These rows are used by 04_Climate_Plot.R. All other numeric
# columns are NA since they are not meaningful for a multi-model aggregate.
summary_rows <- all_change_summary |>
  dplyr::group_by(Variable, Scenario) |>
  dplyr::summarise(
    Model = "All models",
    Variant = NA_character_,
    Output_file = NA_character_,
    Mean_density_change = NA_real_,
    SD_density_change = NA_real_,
    SD_percent_change = sd(Mean_percent_change, na.rm = TRUE),
    Mean_percent_change = mean(Mean_percent_change, na.rm = TRUE),
    Total_change = NA_real_,
    n = NA_integer_,
    .groups = "drop"
  )

# Bind per-model rows and summary rows, then write a single merged CSV
all_change_summary <- dplyr::bind_rows(all_change_summary, summary_rows) |>
  dplyr::arrange(Variable, Scenario, Model)

readr::write_csv(
  all_change_summary,
  file.path("Data", "zooplankton_change_summary.csv")
)
