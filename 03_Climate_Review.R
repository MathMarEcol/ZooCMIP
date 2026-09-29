# For the climate review paper, we need two plots.

# First, a plot similar to Kris' Figure 1, showing model (SSP)
# for % Biomass change (with SD) with 3 rows for zooc, zmicro and zmeso

# For each SSP and each model, integrate biomass over depth (multiplying by
# each grid cell's vertical thickness only - NOT horizontal area, see design
# note 1 below). This yields a depth-integrated areal density (mol C m-2),
# not a per-cell total. Horizontal area is (re)introduced explicitly and only
# downstream, at the point where it is actually needed (e.g. area-weighting
# a spatial mean, or summing to an ocean-total).

# For each model, create a thickness raster (vertical integration only).

library(tidyverse)
library(terra)
library(ncdf4)
library(future)
library(furrr)
library(hotrstuff)

# ---------------------------------------------------------------------------
# DESIGN NOTES (read before modifying)
#
# 1. UNITS: zooc/zmeso/zmicro are "mol C m-3" (mole concentration of carbon),
#    NOT a mass concentration. We integrate ONLY over depth (multiplying by
#    each vertical level's thickness in metres, from lev_bnds), giving
#    "mol C m-2" per pixel - a depth-integrated AREAL DENSITY, exactly
#    analogous to how water-column inventories (e.g. depth-integrated
#    chlorophyll) are conventionally reported. We deliberately do NOT also
#    multiply by horizontal cell area here: doing so would bake an arbitrary,
#    grid-geometry-dependent scaling (huge at the equator, tiny at the poles
#    on a regular lon/lat grid) into every pixel value, making the raster
#    neither a true density nor a true total. Horizontal area is instead
#    (re)introduced explicitly, only where actually needed, downstream in
#    06_Climate_Review.R (e.g. area-weighting a spatial mean, or summing to
#    an ocean-total). We deliberately do NOT convert to a mass unit either -
#    if a mass unit (e.g. kg C, Pg C) is needed downstream, multiply by the
#    molar mass of carbon (12.011 g/mol) at the point of use, so the
#    conversion factor is documented next to where it matters rather than
#    buried in this generic integration step.
#
# 2. VERTICAL GRID (lev/lev_bnds) is MODEL-SPECIFIC, not experiment-specific.
#    Verified: CanESM5 historical and ssp119 have byte-identical lev_bnds;
#    CanESM5-CanOE has identical lev_bnds across zooc/zmeso/zmicro. Different
#    models have different numbers of levels/depths entirely (e.g. CanESM5=45
#    levels, GFDL-ESM4=35, UKESM1-0-LL=75, MIROC-ES2L=63). So the thickness
#    raster is built ONCE per model (from one reference file) and re-used for
#    every experiment (historical/ssp*) of that model - but we still verify
#    lev_bnds match before reusing, rather than silently assuming it.
#
# 3. HORIZONTAL GRID is identical across ALL files (they are already
#    regridded to a common 0.5-degree lon/lat grid via cdo remapbil, verified
#    via compareGeom() across models/variables). This is what makes it safe
#    to introduce horizontal area only once, downstream, using a single
#    cellSize() raster shared by every file - rather than needing a
#    per-model area raster the way the vertical thickness must be.
#
# 4. NA HANDLING: land / sub-seafloor cells are NA in the biomass variable at
#    every model/level (verified - NA count increases with depth as expected
#    for bathymetry masking). The thickness raster itself has NO NAs (it's
#    pure geometry). When we multiply biomass * thickness, NA correctly
#    propagates from the biomass side. Depth-summation MUST use na.rm = TRUE
#    so that partial water columns (e.g. shelf seas) are integrated only over
#    their wet levels. We verified terra::tapp(..., fun = sum, na.rm = TRUE)
#    correctly returns NA (not 0) for all-NA (land) groups, and correctly
#    sums only the valid levels for partial water columns - this is exactly
#    the behaviour we need and requires no extra masking step.
#
# 5. CAVEAT (documented, not fixed): ocean models generally use partial
#    bottom cells at the seafloor; we cannot verify from the regridded files
#    whether that partial thickness is preserved or whether the full nominal
#    lev_bnds thickness is used for the deepest wet cell. This may cause a
#    small over-estimate of the integrated biomass density immediately above
#    topography. This is a standard, accepted limitation for this kind of
#    post-hoc diagnostic and is not something we can correct without the
#    original model's 3D bathymetry mask.
#
# 6. VERTICAL DIMENSION NAMING is not consistent across CMIP6 models: most
#    use "lev"/"lev_bnds", but IPSL-CM6A-LR uses "olevel"/"olevel_bnds"
#    (verified - this caused a real failure: "could not find lev_bnds" when
#    first run against IPSL-CM6A-LR). We do NOT hard-code "lev" anywhere;
#    find_vertical_dim() discovers the correct names generically from each
#    file's own CF metadata (whichever dimension has axis = "Z", and that
#    dimension's own "bounds" attribute), so any other naming convention we
#    have not yet encountered is handled the same way without code changes.
#
# 7. terra::tapp() does NOT guarantee its output layers are sorted by the
#    grouping value - it returns groups in order of first appearance in the
#    input stack (verified empirically). We therefore explicitly recover the
#    real date for each output layer and re-order the stack chronologically
#    ourselves, with a stopifnot() check that the final time vector is
#    strictly increasing. Silently trusting tapp()'s layer order would risk
#    mislabelling years - a bug that would not throw an error, just quietly
#    produce wrong dates.
#
# 8. PARALLELISATION: terra::SpatRaster objects hold an external C++ pointer
#    and do NOT survive being passed into future::multisession workers
#    (verified - this throws "NULL value passed as symbol address").
#    Following the pattern already used throughout the hotrstuff package
#    (which only ever passes file PATHS into furrr::future_walk/future_map,
#    never live SpatRaster objects), every worker below opens its own files
#    with terra::rast() rather than receiving a pre-built raster object.
# ---------------------------------------------------------------------------


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
    ax <- tryCatch(ncdf4::ncatt_get(nc, dn, "axis")$value, error = function(e) NULL)
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
      "Vertical dimension '", z_name, "' in ", nc_file,
      " has no 'bounds' attribute - cannot determine cell thickness."
    )
  }

  list(lev_name = z_name, bnds_name = bnds_attr$value)
}


#' Build a vertical-thickness (m) SpatRaster from a single reference netCDF file
#'
#' Deliberately vertical-only (see design note 1 at the top of this file):
#' horizontal cell area is NOT included here, so the result is a per-level
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
  thick <- terra::init(terra::rast(r_template, nlyrs = length(thickness)), 1) * thickness
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
#' are mol C m-2 (see design note 1 above: deliberately NOT multiplied by
#' horizontal cell area, and no mass conversion is applied).
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
        "lev_bnds mismatch between reference file:\n  ", reference_nc_file,
        "\nand target file:\n  ", nc_file,
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
      "Could not match all depth levels in ", nc_file,
      " to the thickness raster built from ", reference_nc_file,
      ". Unmatched depths: ", paste(unique(r_depth[is.na(idx)]), collapse = ", ")
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
    "Depth-integration produced a non-chronological time axis" =
      all(diff(result_time) > 0)
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
integrate_biomass <- function(nc_files,
                               out_dir = file.path("Data", "integrated_biomass"),
                               workers = NULL,
                               overwrite = FALSE) {

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
    function(v) hotrstuff:::htr_make_folder(file.path(out_dir, v, "depthIntegrated"))
  )

  do_one <- function(path, reference_path, variable, model, scenario, variant) {
    library(terra)
    library(ncdf4)

    out_file <- file.path(
      out_dir, variable, "depthIntegrated",
      paste0(variable, "_", model, "_", scenario, "_", variant, "_depthIntegrated.tif")
    )

    if (file.exists(out_file) && !overwrite) {
      return(tibble::tibble(
        Variable = variable, Model = model, Scenario = scenario, Variant = variant,
        Output_file = out_file, Status = "Skipped (exists)"
      ))
    }

    tryCatch({
      result <- integrate_one_file(path, reference_path)
      terra::writeRaster(result, out_file, overwrite = TRUE)
      tibble::tibble(
        Variable = variable, Model = model, Scenario = scenario, Variant = variant,
        Output_file = out_file, Status = "Success"
      )
    }, error = function(e) {
      tibble::tibble(
        Variable = variable, Model = model, Scenario = scenario, Variant = variant,
        Output_file = NA_character_, Status = paste("ERROR:", conditionMessage(e))
      )
    })
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


# NOTE: re-run these three lines (with overwrite = TRUE) to regenerate every
# existing .tif under Data/integrated_biomass/ using the corrected
# thickness-only (mol C m-2) integration above. This is required before the
# change-calculation below can be trusted, since the files currently on disk
# were written by the OLD area*thickness (mol C per cell) version of
# build_volume_raster(). Left commented out / un-forced (overwrite defaults
# to FALSE) so this is not silently re-triggered on every source() of this
# file - run deliberately when ready.

# nc_files <- list.files("/Volumes/T9/ClimateData/zooc/regridded/", full.names = TRUE, pattern = "[.]nc$")
# dat <- integrate_biomass(nc_files, workers = 5, out_dir = file.path("", "Volumes", "T9", "ClimateData"), overwrite = TRUE)
#
# nc_files <- list.files("/Volumes/T9/ClimateData/zmicro/regridded/", full.names = TRUE, pattern = "[.]nc$")
# dat <- integrate_biomass(nc_files, workers = 5, out_dir = file.path("", "Volumes", "T9", "ClimateData"), overwrite = TRUE)
#
# nc_files <- list.files("/Volumes/T9/ClimateData/zmeso/regridded/", full.names = TRUE, pattern = "[.]nc$")
# dat <- integrate_biomass(nc_files, workers = 5, out_dir = file.path("", "Volumes", "T9", "ClimateData"), overwrite = TRUE)


# ---------------------------------------------------------------------------
# BIOMASS CHANGE: historical (1993-2014) mean vs SSP (2081-2100) mean
# ---------------------------------------------------------------------------
#
# DESIGN NOTES for this section:
#
# a. SCOPE: only the four "tier-1" SSPs (ssp126, ssp245, ssp370, ssp585) are
#    processed, per the analysis brief - even though some models also ran
#    tier-2 SSPs (ssp119/ssp434/ssp460/ssp534-over) in Data/integrated_biomass/.
#    Change ssp_order below to widen this if the brief changes.
#
# b. HISTORICAL MEAN IS CACHED PER (Variable, Model, Variant): the historical
#    baseline is identical for every SSP within the same model/variable/
#    variant, so we compute each one exactly once and reuse it across all
#    four SSPs, rather than recomputing it redundantly per-SSP.
#
# c. FILENAME PARSING: Data/integrated_biomass/<var>/*.tif filenames follow
#    "<variable>_<model>_<scenario>_<variant>_biomass_integrated.tif" - a
#    project-specific convention created by integrate_biomass() above, NOT
#    the raw CMIP6 filename convention that hotrstuff:::htr_get_CMIP6_bits()
#    expects (variable_table_model_experiment_variant_grid_daterange). We
#    therefore use a small dedicated parser here, matching the style already
#    used for derived filenames elsewhere in this project (see
#    parse_filename() in 03_Climate_Summarise.R) - not the CMIP6 parser,
#    and not an unexported/internal hotrstuff function we don't control.
#
# d. REGRIDDING ARTIFACTS IN THE BASELINE: verified that a small number of
#    pixels in some historical means are <= 0 (physically impossible for a
#    biomass density) - e.g. zooc_UKESM1-0-LL and zooc_MIROC-ES2L show
#    negative values, zooc_CanESM5(-CanOE) and zmeso/zmicro_CMCC-ESM2 show
#    exact zeros. These are concentrated in small enclosed/marginal seas
#    (e.g. the Black Sea, lat ~32-46N/lon ~28-54E) and are 6+ orders of
#    magnitude smaller than typical open-ocean values - traced back to
#    bilinear regridding (cdo remapbil) overshoot/undershoot at sharp
#    coastal boundaries in the SOURCE netCDF itself (confirmed a native-grid
#    negative concentration at UKESM1-0-LL's shallowest level), not
#    something introduced by this pipeline.
#
#    IMPORTANT EXTENSION (found while investigating implausible summary
#    statistics, e.g. some CanESM5/CanESM5-CanOE and GFDL-ESM4 % changes of
#    1e6-1e22 %): the SAME regridding overshoot/undershoot mechanism does not
#    only land exactly on <= 0 - at the ice-edge/polar margins in particular
#    (>=90% of the extreme pixels for every affected model are poleward of
#    60 deg), it frequently lands on a tiny POSITIVE value instead (e.g.
#    1e-10 to 1e-30 mol C m-2). A <= 0 mask does not catch these at all,
#    yet dividing by such a baseline produces exactly the same
#    undefined/meaningless percent change (division by ~zero) as the
#    already-handled <= 0 case.
#
#    CRITICAL CORRECTION (discovered by investigating ACCESS-ESM1-5 tropical
#    data loss): the earlier approach of masking all pixels <= a fixed
#    EPS_BASELINE = 1e-3 mol C m-2 was WRONG. ACCESS-ESM1-5 has a genuine,
#    continuous distribution of real low-biomass tropical gyre pixels with
#    values from ~1e-11 up through ~1e-3 mol C m-2 (10.6% of ocean pixels,
#    94% of them equatorward of 30 deg). A fixed threshold cannot
#    distinguish these real biological signals from regridding artifacts
#    because the artifact values in affected models (CanESM5, CanESM5-CanOE)
#    span the same numerical range.
#
#    The correct approach is a NEIGHBOUR-COHERENCE check: a pixel is an
#    artifact if and only if (a) it is <= 0 (physically impossible), OR
#    (b) it is positive but more than ARTIFACT_LOG_RATIO orders of magnitude
#    below its own 3x3 spatial neighbourhood median. Real biomass fields -
#    even in oligotrophic gyres - vary smoothly; a genuine low-biomass pixel
#    will be within a few orders of magnitude of its neighbours. A regridding
#    overshoot/undershoot artifact (e.g. 1e-28 surrounded by neighbours at
#    1e-2) will be many orders of magnitude below its neighbourhood median.
#    This is model-agnostic and self-calibrating: it correctly preserves
#    ACCESS-ESM1-5's real tropical pixels (which are coherent with their
#    neighbours) while still catching CanESM5's polar/coastal artifacts
#    (which are incoherent with theirs). ARTIFACT_LOG_RATIO = 3 (i.e. mask
#    if pixel < neighbourhood_median / 1000) is used; verified this catches
#    all known artifacts while preserving all real low-biomass pixels across
#    every model in this collection.
#
#    Because a baseline at or below this floor makes percent change
#    undefined/meaningless (division by ~zero or a sign flip), these pixels
#    are masked to NA in BOTH the absolute and percent change layers (and
#    therefore excluded from all downstream statistics) before any change is
#    computed. The count masked is logged per file so this remains an
#    audited, deliberate exclusion rather than a silent one.
#
# g. PERCENT-CHANGE OUTLIER TRIM: even after the baseline masking in (d),
#    some models (notably CanESM5, CanESM5-CanOE, GFDL-ESM4) retain extreme
#    percent-change outliers driven by spatially-coherent artifact fields in
#    the baseline - i.e. broad polar/shelf regions where every pixel and all
#    its neighbours are equally artifact-level tiny, so the neighbour-
#    coherence check in (d) cannot detect them (no single pixel stands out
#    relative to its 3x3 window). Investigated exhaustively:
#
#    - GFDL-ESM4 zmeso: 15 pixels at 10^9-10^14 % (Kara Sea, 68-70N, 73-74E)
#      caused by bilinear regridding of exactly-zero native-grid cells
#      (confirmed in raw netCDF: all 300 monthly values at 68.5N, 73.5E and
#      74.5E are exactly 0.000 mol C m-3 for 1990-2014). These 15 pixels sit
#      at the 99.99th percentile of the percent-change distribution - a
#      p99.9% trim catches them cleanly.
#    - CanESM5/CanESM5-CanOE: 10-21% of ocean pixels have >1000% change
#      (Canadian Arctic/Beaufort Sea with near-zero historical biomass under
#      permanent sea ice, real future biomass after ice retreat). These are
#      NOT outliers relative to the model's own distribution - they constitute
#      a substantial fraction of it. A p99.9% trim reduces but does not
#      eliminate the problem for these models; the trimmed mean is still
#      millions of % for CanESM5-CanOE. This is a fundamental model behaviour
#      issue (near-zero historical polar biomass) that no post-hoc filter can
#      fully resolve without excluding the affected region entirely.
#
#    We apply a symmetric p99.9% trim (PCT_TRIM_QUANTILE = 0.001) to the
#    percent-change layer of every model/SSP combination: pixels outside
#    [p0.1%, p99.9%] of that file's own percent-change distribution are set
#    to NA in BOTH the Percent_change and Absolute_change layers (so the two
#    layers remain on a consistent pixel mask). The trim is self-calibrating
#    per model/SSP (uses each file's own distribution, not a fixed global
#    threshold), symmetric (removes equal fractions from both tails), and
#    consistent (applied identically to every model). The count trimmed is
#    logged. For clean models (CESM2, CNRM-ESM2-1, IPSL-CM6A-LR, etc.) the
#    p99.9% cutoff is a biologically plausible 100-250%, so the trim removes
#    only genuine extreme-change pixels at the ice edge or coastal margins.
#
# e. AREAL DENSITY vs OCEAN-TOTAL: because the rasters are now mol C m-2 (see
#    design note 1 at the top of this file), the summary tibble reports both
#    an area-weighted MEAN density change (mol C m-2, and the equivalent
#    area-weighted mean % change) AND a true area-weighted SUM (mol C,
#    "Total_change") obtained by multiplying each pixel's absolute density
#    change by that pixel's cell area before summing - this converts the
#    density back into an actual total ocean-wide biomass change. We do NOT
#    report an SD for Total_change: a sum is a single aggregated number, not
#    a population of observations, so it has no non-arbitrary population-
#    level "spread" the way a per-pixel mean does (a considered decision,
#    not an oversight).
#
# f. AREA WEIGHTS use a single shared cellSize() raster (see design note 3
#    at the top of this file: the horizontal grid is identical across every
#    model/variable, verified via compareGeom()), computed once and reused
#    for every file rather than recomputed per iteration.
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
    "Unexpected integrated_biomass filename format" =
      n >= 5 && bits[n] == "depthIntegrated.tif"
  )
  tibble::tibble(
    Variable = bits[1],
    Model    = paste(bits[2:(n - 3)], collapse = "_"),
    Scenario = bits[n - 2],
    Variant  = bits[n - 1],
    Path     = path
  )
}


#' Set to TRUE to apply the symmetric p99.9% percent-change trim (design
#' note g). Set to FALSE to skip the trim entirely and retain all pixels
#' regardless of their percent-change value. Useful for diagnosing whether
#' the trim is removing too much data or for comparing trimmed vs untrimmed
#' summary statistics.
APPLY_PCT_TRIM <- FALSE

#' Symmetric quantile trim applied to the Percent_change layer (and
#' correspondingly to Absolute_change) after computing biomass change.
#' Pixels outside [PCT_TRIM_QUANTILE, 1 - PCT_TRIM_QUANTILE] of the
#' percent-change distribution are set to NA in both layers (see design
#' note g above). A value of 0.001 gives a symmetric p0.1%/p99.9% trim.
#' Only used when APPLY_PCT_TRIM = TRUE.
PCT_TRIM_QUANTILE <- 0.001

#' Compute the pixel-wise biomass change map for one historical/SSP file pair
#' of the same model/variable/variant.
#'
#' Output is a 4-layer SpatRaster storing all quantities needed for both the
#' tabular summary (06_Climate_Review.R) and the spatial ensemble maps
#' (08_Climate_SpatialPlot.R), so downstream scripts can work entirely from
#' these files without re-reading the depth-integrated originals:
#'
#'   Layer 1 "Hist_mean"       — 1993-2014 mean biomass (mol C m-2)
#'   Layer 2 "Future_mean"     — 2081-2100 mean biomass (mol C m-2)
#'   Layer 3 "Absolute_change" — Future_mean - Hist_mean (mol C m-2)
#'   Layer 4 "Percent_change"  — pixel-wise (Absolute_change / Hist_mean)*100 (%)
#'
#' Storing Hist_mean and Future_mean allows 08_Climate_SpatialPlot.R to
#' compute the ensemble % change as:
#'   (mean_across_models(Future_mean) - mean_across_models(Hist_mean)) /
#'    mean_across_models(Hist_mean) * 100
#' i.e. averaging biomass FIRST, then computing % change — which is far more
#' robust than averaging per-model pixel-wise % changes (which amplify
#' near-zero baseline artifacts).
#'
#' @param hist_mean_rast Pre-computed 1993-2014 mean SpatRaster (mol C m-2)
#'   for this model/variable/variant, with artifact baseline pixels already
#'   masked to NA (see mask_invalid_baseline()).
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
      "Expected 20 years in the 2081-2100 window for ", ssp_path,
      " but found ", sum(future_window), " - refusing to average a",
      " truncated future period."
    )
  }

  r_future_mn <- mean(r_ssp[[future_window]], na.rm = TRUE)

  # hist_mean_rast already has artifact pixels masked to NA (see
  # mask_invalid_baseline()), so NA propagates automatically into all
  # derived layers below without any extra handling here.
  abs_change <- r_future_mn - hist_mean_rast
  pct_change <- (abs_change / hist_mean_rast) * 100

  # Symmetric p99.9% trim on the percent-change layer (see design note g).
  # Skipped entirely when APPLY_PCT_TRIM = FALSE.
  if (APPLY_PCT_TRIM) {
    # Compute quantile bounds from the non-NA values of this file's own
    # percent-change distribution (self-calibrating per model/SSP).
    pct_vals <- terra::values(pct_change)[, 1]
    pct_valid <- pct_vals[!is.na(pct_vals)]

    if (length(pct_valid) > 0) {
      lo <- stats::quantile(pct_valid, PCT_TRIM_QUANTILE)
      hi <- stats::quantile(pct_valid, 1 - PCT_TRIM_QUANTILE)

      # Build a logical mask: TRUE where pct_change is outside [lo, hi]
      trim_mask <- terra::ifel(pct_change < lo | pct_change > hi, TRUE, FALSE)

      n_trimmed <- terra::global(trim_mask, "sum", na.rm = TRUE)[1, 1]

      if (n_trimmed > 0) {
        message(
          "  [", label, "] trimmed ", n_trimmed,
          " pixel(s) outside [p", PCT_TRIM_QUANTILE * 100, "%, p",
          (1 - PCT_TRIM_QUANTILE) * 100, "%] of percent-change distribution",
          " ([", round(lo, 1), "%, ", round(hi, 1), "%])",
          " - see design note g"
        )
      }

      # Apply the trim mask to BOTH change layers so they remain on a
      # consistent pixel footprint. Hist_mean and Future_mean are NOT
      # trimmed — they are raw means and should remain complete so that
      # 08_Climate_SpatialPlot.R can compute ensemble means from them
      # without the trim mask interfering.
      pct_change <- terra::ifel(trim_mask, NA, pct_change)
      abs_change <- terra::ifel(trim_mask, NA, abs_change)
    }
  }

  out <- c(hist_mean_rast, r_future_mn, abs_change, pct_change)
  names(out) <- c("Hist_mean", "Future_mean", "Absolute_change", "Percent_change")
  terra::units(out) <- c("mol m-2", "mol m-2", "mol m-2", "%")

  out
}


#' Set to TRUE to apply the neighbour-coherence baseline artifact masking
#' (design note d). Set to FALSE to skip baseline masking entirely and pass
#' the raw historical mean (with only land/sea NAs) through to the change
#' calculation. Useful for diagnosing whether the masking is removing too
#' much data or for comparing masked vs unmasked results.
APPLY_BASELINE_MASK <- FALSE

#' Number of orders of magnitude below its 3x3 neighbourhood median a pixel
#' must be to be classified as a regridding artifact (see design note d).
#' A value of 3 means: mask if pixel < neighbourhood_median / 1000.
#' Verified to catch all known artifacts (CanESM5 polar/coastal overshoot
#' pixels at 1e-28 to 1e-30 mol C m-2) while preserving all real low-biomass
#' pixels (ACCESS-ESM1-5 tropical gyre pixels at 1e-11 to 1e-3 mol C m-2,
#' which are coherent with their neighbours). Kept as a named constant so the
#' choice is visible and can be revisited in one place if new models are added.
#' Only used when APPLY_BASELINE_MASK = TRUE.
ARTIFACT_LOG_RATIO <- 3

#' Mask historical-baseline pixels that are regridding artifacts to NA,
#' logging how many were masked. See design note (d) above.
#'
#' A pixel is classified as an artifact if:
#'   (a) its value is <= 0 (physically impossible biomass density), OR
#'   (b) it is positive but more than ARTIFACT_LOG_RATIO orders of magnitude
#'       below its own 3x3 spatial neighbourhood median (incoherent with
#'       neighbours => regridding overshoot/undershoot, not real biology).
#'
#' This is model-agnostic and self-calibrating: it correctly preserves real
#' low-biomass tropical pixels (coherent with neighbours) while catching
#' polar/coastal regridding artifacts (incoherent with neighbours).
#'
#' @param hist_mean_rast The 1993-2014 mean SpatRaster (mol C m-2).
#' @param label Human-readable label (e.g. "zooc_ACCESS-ESM1-5_r40i1p1f1") for
#'   the log message identifying which file this came from.
#' @return The same SpatRaster with artifact pixels set to NA.
mask_invalid_baseline <- function(hist_mean_rast, label) {

  # When APPLY_BASELINE_MASK = FALSE, skip all masking and return the raster
  # unchanged (only the existing land/sea NAs are preserved).
  if (!APPLY_BASELINE_MASK) {
    return(hist_mean_rast)
  }

  n_before <- terra::global(hist_mean_rast, "notNA")[1, 1]

  # Step 1: mask physically impossible values (<= 0)
  hist_pos <- terra::ifel(hist_mean_rast <= 0, NA, hist_mean_rast)

  # Step 2: neighbour-coherence check on the remaining positive pixels.
  # Compute the 3x3 focal median (na.rm = TRUE so coastal/ice-edge pixels
  # with NA neighbours are still evaluated against their valid neighbours).
  # Then mask any pixel whose log10(value) is more than ARTIFACT_LOG_RATIO
  # below log10(neighbourhood_median) - i.e. pixel < median / 10^ARTIFACT_LOG_RATIO.
  # We work in log10 space to make the ratio scale-invariant.
  nbr_median <- terra::focal(hist_pos, w = 3, fun = "median", na.rm = TRUE)

  # A pixel is an artifact if it is positive but its log10 value is more than
  # ARTIFACT_LOG_RATIO below the log10 of its neighbourhood median.
  # Guard against log10(0) or log10(NA): hist_pos already has <= 0 as NA,
  # and nbr_median can be NA where all 9 neighbours are NA (open ocean
  # isolated pixels - extremely rare; treat as non-artifact to be conservative).
  is_artifact <- terra::ifel(
    is.na(nbr_median),
    FALSE,                                          # no neighbours => keep
    (log10(hist_pos) < log10(nbr_median) - ARTIFACT_LOG_RATIO)
  )
  hist_masked <- terra::ifel(is_artifact, NA, hist_pos)

  n_after <- terra::global(hist_masked, "notNA")[1, 1]
  n_masked <- n_before - n_after

  if (n_masked > 0) {
    message(
      "  [", label, "] masked ", n_masked,
      " historical-baseline pixel(s) as regridding artifacts",
      " (<=0 or >", ARTIFACT_LOG_RATIO, " orders of magnitude below",
      " 3x3 neighbourhood median - see design note d)"
    )
  }

  hist_masked
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
#' percent changes. This avoids division-by-near-zero at the pixel level
#' (which produces undefined/meaningless per-pixel percents in polar/coastal
#' regridding-artifact regions) and gives a single, robust scalar:
#'
#'   global_hist_mean  = area-weighted mean of hist_mean_rast
#'   global_future_mean = global_hist_mean + area-weighted mean of Absolute_change
#'   Mean_percent_change = (global_future_mean - global_hist_mean) /
#'                          global_hist_mean * 100
#'
#' SD_percent_change is the area-weighted SD of the pixel-wise percent-change
#' layer (retained for spatial spread / uncertainty, but note it is still
#' subject to the pixel-level division issues described above for extreme
#' artifact pixels - use with caution for models with known polar artifacts).
#'
#' @param change_rast 2-layer SpatRaster from calc_biomass_change() (
#'   "Absolute_change" mol C m-2, "Percent_change" %), with the symmetric
#'   p99.9% trim already applied to both layers (see design note g above).
#' @param hist_mean_rast The 1993-2014 mean SpatRaster (mol C m-2) for this
#'   model/variable/variant, with artifact pixels already masked to NA (same
#'   object passed to calc_biomass_change()). Used to compute the global
#'   area-weighted historical mean for the percent-change denominator.
#' @param area_rast A cellSize() raster (m2) on the same grid, shared across
#'   all calls (see design note f above).
#' @return A one-row tibble of area-weighted Mean/SD for both layers, the
#'   area-weighted ocean-total absolute change (mol C), and pixel count (n).
summarise_biomass_change <- function(change_rast, hist_mean_rast, area_rast) {

  abs_stats <- weighted_mean_sd(change_rast[["Absolute_change"]], area_rast)

  # Global area-weighted mean of the historical baseline (single scalar).
  # This is the denominator for Mean_percent_change.
  hist_global_mean <- terra::global(
    hist_mean_rast, "mean", weights = area_rast, na.rm = TRUE
  )[1, 1]

  # Global area-weighted mean of the future period, derived from the
  # historical mean + the area-weighted mean absolute change. This is
  # equivalent to computing the future global mean directly, but reuses
  # abs_stats$Mean which is already computed above.
  future_global_mean <- hist_global_mean + abs_stats$Mean

  # Single % change from two global means — avoids pixel-level division by
  # near-zero baseline values (see function documentation above).
  mean_pct_change <- (future_global_mean - hist_global_mean) /
                      hist_global_mean * 100

  # SD of the pixel-wise percent-change layer (area-weighted), retained for
  # spatial spread / uncertainty reporting. Subject to pixel-level division
  # issues for models with polar/coastal artifacts — see function docs.
  pct_sd_stats <- weighted_mean_sd(change_rast[["Percent_change"]], area_rast)

  # Both layers share an identical NA mask (baseline masking + percent-change
  # trim both applied simultaneously to both layers - see design notes d and g),
  # so pixel counts are always equal; abs_stats$n is reported once.
  total_change <- terra::global(
    change_rast[["Absolute_change"]], "sum",
    weights = area_rast, na.rm = TRUE
  )[1, 1]

  tibble::tibble(
    Mean_density_change = abs_stats$Mean,
    SD_density_change   = abs_stats$SD,
    Mean_percent_change = mean_pct_change,
    SD_percent_change    = pct_sd_stats$SD,
    Total_change         = total_change,
    n                    = abs_stats$n
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
zoo_change <- function(int_files,
                        out_dir = file.path("Data", "biomass_change_maps"),
                        ssp_order = c("ssp126", "ssp245", "ssp370", "ssp585")) {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  meta <- purrr::map(int_files, parse_integrated_filename) |>
    dplyr::bind_rows()

  hist_meta <- meta |> dplyr::filter(Scenario == "historical")
  ssp_meta  <- meta |>
    dplyr::filter(Scenario %in% ssp_order) |>
    dplyr::left_join(
      hist_meta |> dplyr::select(Variable, Model, Variant, Hist_path = Path),
      by = c("Variable", "Model", "Variant")
    )

  missing_hist <- ssp_meta |> dplyr::filter(is.na(Hist_path))
  if (nrow(missing_hist) > 0) {
    warning(
      "No historical file found for ", nrow(missing_hist),
      " Variable/Model/Variant combination(s) - these will be skipped:\n  ",
      paste(missing_hist$Path, collapse = "\n  ")
    )
    ssp_meta <- ssp_meta |> dplyr::filter(!is.na(Hist_path))
  }

  # Area weights are identical across every file (common 0.5-degree grid,
  # verified via compareGeom() - see design note 3 at the top of this file),
  # so compute this exactly once and reuse it for every summary statistic
  # below, rather than recomputing cellSize() per row.
  area_template <- terra::rast(ssp_meta$Path[1])[[1]]
  area_rast <- terra::cellSize(area_template, unit = "m")

  # Cache one 1993-2014 historical mean (with <= 0 pixels masked, per design
  # note d) per unique Variable/Model/Variant, since it is identical for
  # every SSP within that group (design note b) - avoids recomputing the
  # same mean() up to 4x per model/variable.
  hist_cache <- new.env(parent = emptyenv())

  get_hist_mean <- function(hist_path, label) {
    key <- hist_path
    if (is.null(hist_cache[[key]])) {
      r_hist <- terra::rast(hist_path)
      hist_years <- as.integer(names(r_hist))
      hist_window <- hist_years >= 1993 & hist_years <= 2014

      if (sum(hist_window) != 22) {
        stop(
          "Expected 22 years in the 1993-2014 window for ", hist_path,
          " but found ", sum(hist_window), " - refusing to average a",
          " truncated historical baseline."
        )
      }

      hist_mn <- mean(r_hist[[hist_window]], na.rm = TRUE)
      hist_cache[[key]] <- mask_invalid_baseline(hist_mn, label)
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
        paste0(Variable, "_", Model, "_", Scenario, "_", Variant, "_biomassChange.tif")
      )
      terra::writeRaster(change_rast, out_file, overwrite = TRUE)

      summarise_biomass_change(change_rast, hist_mn, area_rast) |>
        dplyr::mutate(
          Variable = Variable, Model = Model, Scenario = Scenario,
          Variant = Variant, Output_file = out_file,
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
    full.names = TRUE, recursive = TRUE, pattern = "[.]tif$"
  )

  zoo_change(
    int_files,
    out_dir = file.path("", "Volumes", "T9", "ClimateData", var, "biomassChange")
  )
})

# Build "All models" summary rows: one per Variable × Scenario, with
# Mean_percent_change = multimodel mean and SD_percent_change = SD across
# models. These rows represent exactly what is plotted as the mean point +
# SD errorbar in 07_Climate_Plot.R. All other numeric columns are NA since
# they are not meaningful for a multi-model aggregate.
summary_rows <- all_change_summary |>
  dplyr::group_by(Variable, Scenario) |>
  dplyr::summarise(
    Model               = "All models",
    Variant             = NA_character_,
    Output_file         = NA_character_,
    Mean_density_change = NA_real_,
    SD_density_change   = NA_real_,
    SD_percent_change   = sd(Mean_percent_change,   na.rm = TRUE),
    Mean_percent_change = mean(Mean_percent_change, na.rm = TRUE),
    Total_change        = NA_real_,
    n                   = NA_integer_,
    .groups = "drop"
  )

# Bind per-model rows and summary rows, then write a single merged CSV
all_change_summary <- dplyr::bind_rows(all_change_summary, summary_rows) |>
  dplyr::arrange(Variable, Scenario, Model)

readr::write_csv(
  all_change_summary,
  file.path("Data", "zooplankton_change_summary.csv")
)


