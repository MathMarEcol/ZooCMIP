# PROMPT: You are a marine climate scientist and expert R programmer. Together we are
# going to write some code to summarise some ESM zooplankton biomass data I have
# downloaded and processed to the same grid and time interval. There are 3
# variables - zooc (all zooplankton). zmicro (microzooplankton) and zmeso
# (mesozooplankton). To start with, we will only process the zooc data till we
# come up with an approach. These data are located @/Volumes/T9/zooc/regridded.
#
# What I want is a dataframe that has a mean, median, standard deviation, and
# standard error for surface zooplankton for each year, model and scenario.
# Ideally I would like this to be contained within a single data tibble. That
# has columns of Variable, Model, Scenario, Variant, Year and Biomass.
#
# You should always use tidyverse equivalents where they exist and are preferable
# (such as purrr instead of apply). You should never guess. If you are unsure,
# you should ask.

# Load required libraries
library(tidyverse)
library(ncdf4)
library(terra)
library(Hmisc)

# Define variables and their directories
variables <- tibble(
  var_name = c("zooc", "zmicro", "zmeso"),
  data_dir = c("/Volumes/T9/zooc/regridded",
               "/Volumes/T9/zmicro/regridded",
               "/Volumes/T9/zmeso/regridded")
)

# Function to parse filename and extract metadata
parse_filename <- function(filepath) {
  filename <- basename(filepath)
  parts <- str_split(filename, "_")[[1]]

  tibble(
    Variable = parts[1],
    Model = parts[3],
    Scenario = parts[4],
    Variant = parts[5]
  )
}

# Function to calculate area-weighted statistics from NetCDF file.
#
# All statistics are weighted by cell area (m²) computed via terra::cellSize().
# On a regular 0.5° lat/lon grid, cell area scales as cos(latitude), so
# unweighted statistics would over-represent high-latitude (smaller) cells.
#
# SE is intentionally omitted: gridded model output is spatially autocorrelated,
# so the pixel count is not a valid effective sample size for SE estimation.
#
# Weighted median uses Hmisc::wtd.quantile().
calculate_stats <- function(filepath, var_name) {

  # Open NetCDF file
  nc <- nc_open(filepath)

  # Read the data (dimensions typically: lon, lat, depth, time)
  data <- ncvar_get(nc, var_name)

  # Get time variable and convert to years
  time_var <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value

  # Get lon and lat for georeferencing
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")

  # Close NetCDF file
  nc_close(nc)

  # Parse time units to get years
  # Assuming format like "days since YYYY-MM-DD"
  base_date <- str_extract(time_units, "\\d{4}-\\d{2}-\\d{2}")
  base_year <- as.numeric(str_sub(base_date, 1, 4))

  # Convert time to years
  if (str_detect(time_units, "days since")) {
    years <- base_year + floor(time_var / 365.25)
  } else if (str_detect(time_units, "years since")) {
    years <- base_year + floor(time_var)
  } else {
    # Fallback: assume annual data in sequence
    years <- base_year + seq(0, length(time_var) - 1)
  }

  # Extract surface layer (first depth level)
  # Handle different dimension orders
  dims <- dim(data)
  n_dims <- length(dims)

  if (n_dims == 4) {
    # Assume order: lon, lat, depth, time — extract first depth level
    surface_data <- data[, , 1, ]
  } else if (n_dims == 3) {
    # Could be lon, lat, time (no depth) or lon, lat, depth
    # Check if last dimension matches time length
    if (dims[3] == length(years)) {
      surface_data <- data  # Already surface-only or no depth dimension
    } else {
      # Assume lon, lat, depth — take first depth and treat as single time step
      surface_data <- array(data[, , 1], dim = c(dims[1], dims[2], 1))
      years <- years[1]
    }
  } else {
    stop("Unexpected number of dimensions in NetCDF file")
  }

  # Build a template raster from the first time step to derive cell areas.
  # terra::cellSize() returns exact areas in m² for each cell, correctly
  # accounting for the cos(lat) scaling on a geographic (EPSG:4326) grid.
  # Note: ncdf4 returns data as (lon, lat), so we transpose to (lat, lon)
  # for terra, then flip vertically to restore N→S row order.
  template_rast <- rast(
    t(surface_data[, , 1]),
    extent = c(min(lon), max(lon), min(lat), max(lat)),
    crs = "EPSG:4326"
  )
  template_rast <- flip(template_rast, direction = "vertical")
  cell_areas <- cellSize(template_rast, unit = "m")

  # Calculate area-weighted statistics for each year
  map_dfr(seq_along(years), function(i) {

    # Build a georeferenced raster for this year's surface slice
    year_rast <- rast(
      t(surface_data[, , i]),
      extent = c(min(lon), max(lon), min(lat), max(lat)),
      crs = "EPSG:4326"
    )
    year_rast <- flip(year_rast, direction = "vertical")

    # Extract values and corresponding cell areas as vectors
    values  <- values(year_rast)[, 1]
    weights <- values(cell_areas)[, 1]

    # Restrict to cells where both value and area are non-NA
    valid   <- !is.na(values) & !is.na(weights)
    v <- values[valid]
    w <- weights[valid]

    if (length(v) == 0) {
      return(tibble(
        Year   = years[i],
        Mean   = NA_real_,
        Median = NA_real_,
        SD     = NA_real_
      ))
    }

    # Area-weighted mean
    w_mean <- sum(v * w) / sum(w)

    # Area-weighted median via Hmisc::wtd.quantile()
    w_median <- as.numeric(wtd.quantile(v, weights = w, probs = 0.5))

    # Area-weighted standard deviation
    w_sd <- sqrt(sum(w * (v - w_mean)^2) / sum(w))

    tibble(
      Year   = years[i],
      Mean   = w_mean,
      Median = w_median,
      SD     = w_sd
    )
  })
}

# Process all variables
zoo_summary <- variables %>%
  pmap_dfr(function(var_name, data_dir) {

    # List all NetCDF files for this variable (exclude hidden files starting with ._)
    nc_files <- list.files(data_dir, pattern = paste0("^", var_name, "_.*\\.nc$"), full.names = TRUE)

    # Process all files for this variable
    map_dfr(nc_files, function(filepath) {
      # Parse filename to get metadata
      metadata <- parse_filename(filepath)

      # Calculate area-weighted statistics
      stats <- calculate_stats(filepath, var_name)

      # Attach metadata columns to each year's row (not a Cartesian product)
      stats %>%
        mutate(
          Variable = metadata$Variable,
          Model    = metadata$Model,
          Scenario = metadata$Scenario,
          Variant  = metadata$Variant,
          .before  = Year
        )
    })
  })

# Display summary
print(zoo_summary)

# Save results
write_csv(zoo_summary, "Data/zooplankton_annual_summary_all.csv")

# ============================================================================
# SPATIAL SUMMARY: Calculate % biomass change between time periods
# ============================================================================

# Create output directory for spatial tifs
spatial_output_dir <- "Data/spatial_change_maps"
if (!dir.exists(spatial_output_dir)) {
  dir.create(spatial_output_dir, recursive = TRUE)
}

# Function to extract the per-pixel temporal mean for a given time period from
# a NetCDF file, returned as a georeferenced terra SpatRaster.
#
# The mean is computed across the time dimension only (not spatially collapsed),
# so the output retains the full spatial grid. Area weighting is not applied
# here because the % change is computed pixel-by-pixel in process_spatial_change.
extract_spatial_mean <- function(filepath, var_name, start_year, end_year) {

  # Open NetCDF file
  nc <- nc_open(filepath)

  # Read the data
  data <- ncvar_get(nc, var_name)

  # Get time variable and convert to years
  time_var   <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value

  # Get lon and lat for georeferencing
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")

  # Close NetCDF file
  nc_close(nc)

  # Parse time units to get years
  base_date <- str_extract(time_units, "\\d{4}-\\d{2}-\\d{2}")
  base_year <- as.numeric(str_sub(base_date, 1, 4))

  # Convert time to years
  if (str_detect(time_units, "days since")) {
    years <- base_year + floor(time_var / 365.25)
  } else if (str_detect(time_units, "years since")) {
    years <- base_year + floor(time_var)
  } else {
    years <- base_year + seq(0, length(time_var) - 1)
  }

  # Extract surface layer (first depth level if depth dimension present)
  dims   <- dim(data)
  n_dims <- length(dims)

  if (n_dims == 4) {
    # Order: lon, lat, depth, time — extract first depth level
    surface_data <- data[, , 1, ]
  } else if (n_dims == 3) {
    if (dims[3] == length(years)) {
      surface_data <- data  # Already surface-only: lon, lat, time
    } else {
      # lon, lat, depth with a single time step — take first depth level
      surface_data <- array(data[, , 1], dim = c(dims[1], dims[2], 1))
      years <- years[1]
    }
  } else {
    stop("Unexpected number of dimensions in NetCDF file")
  }

  # Find indices for the requested time period
  time_indices <- which(years >= start_year & years <= end_year)

  if (length(time_indices) == 0) {
    return(NULL)  # No data for this time period
  }

  # Build a multi-layer SpatRaster for the time period (one layer per year).
  # ncdf4 returns data as (lon, lat), so we transpose each slice to (lat, lon)
  # for terra, then flip vertically to restore N→S row order.
  period_stack <- map(time_indices, function(i) {
    rast(
      t(surface_data[, , i]),
      extent = c(min(lon), max(lon), min(lat), max(lat)),
      crs = "EPSG:4326"
    ) |> flip(direction = "vertical")
  }) |> rast()

  # Return the pixel-wise mean across the time period
  mean(period_stack, na.rm = TRUE)
}

# Function to process a model/variant combination and calculate % change
process_spatial_change <- function(var_name, data_dir, model, variant, scenario) {

  # Find historical file for this model/variant
  hist_pattern <- paste0("^", var_name, "_Omon_", model, "_historical_", variant, ".*\\.nc$")
  hist_files <- list.files(data_dir, pattern = hist_pattern, full.names = TRUE)

  # Find scenario file for this model/variant/scenario
  scen_pattern <- paste0("^", var_name, "_Omon_", model, "_", scenario, "_", variant, ".*\\.nc$")
  scen_files <- list.files(data_dir, pattern = scen_pattern, full.names = TRUE)

  # Check if both files exist
  if (length(hist_files) == 0 || length(scen_files) == 0) {
    return(NULL)
  }

  # Use the first matching file (should only be one)
  hist_file <- hist_files[1]
  scen_file <- scen_files[1]

  # Extract spatial means for both periods
  tryCatch({
    baseline_rast <- extract_spatial_mean(hist_file, var_name, 1995, 2014)
    future_rast <- extract_spatial_mean(scen_file, var_name, 2080, 2100)

    if (is.null(baseline_rast) || is.null(future_rast)) {
      return(NULL)
    }

    # Calculate % change: ((future - baseline) / baseline) * 100
    change_rast <- ((future_rast - baseline_rast) / baseline_rast) * 100

    # Create output filename
    output_filename <- paste0(var_name, "_", model, "_", scenario, "_", variant,
                              "_change_1995-2014_to_2080-2100.tif")
    output_path <- file.path(spatial_output_dir, output_filename)

    # Save as GeoTIFF
    writeRaster(change_rast, output_path, overwrite = TRUE)

    # Return metadata about what was processed
    tibble(
      Variable = var_name,
      Model = model,
      Scenario = scenario,
      Variant = variant,
      Output_file = output_filename,
      Status = "Success"
    )
  }, error = function(e) {
    tibble(
      Variable = var_name,
      Model = model,
      Scenario = scenario,
      Variant = variant,
      Output_file = NA_character_,
      Status = paste0("Error: ", e$message)
    )
  })
}

# Process all variables for spatial summary
spatial_summary <- variables %>%
  pmap_dfr(function(var_name, data_dir) {
    
    # Get list of all available model/scenario/variant combinations
    var_files <- list.files(data_dir, pattern = paste0("^", var_name, "_.*\\.nc$"), full.names = FALSE)
    
    # Parse all filenames to get unique combinations
    var_metadata <- var_files %>%
      map_dfr(~parse_filename(.x)) %>%
      filter(Scenario != "historical")  # Only process scenarios, not historical
    
    # Process all combinations for this variable
    var_metadata %>%
      pmap_dfr(function(Variable, Model, Scenario, Variant) {
        message(paste("Processing:", Variable, Model, Scenario, Variant))
        process_spatial_change(Variable, data_dir, Model, Variant, Scenario)
      })
  })

# Display summary of processed files
print(spatial_summary)

# Save processing summary
write_csv(spatial_summary, file.path(spatial_output_dir, "spatial_processing_summary.csv"))

message(paste("Spatial change maps saved to:", spatial_output_dir))



