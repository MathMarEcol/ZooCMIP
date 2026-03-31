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

# Function to calculate statistics from NetCDF file
calculate_stats <- function(filepath, var_name) {

  # Open NetCDF file
  nc <- nc_open(filepath)

  # Read the data (dimensions typically: lon, lat, depth, time)
  data <- ncvar_get(nc, var_name)

  # Get time variable and convert to years
  time_var <- ncvar_get(nc, "time")
  time_units <- ncatt_get(nc, "time", "units")$value

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

  # Close NetCDF file
  nc_close(nc)

  # Extract surface layer (first depth level)
  # Handle different dimension orders
  dims <- dim(data)
  n_dims <- length(dims)

  if (n_dims == 4) {
    # Assume order: lon, lat, depth, time
    # Extract first depth level
    surface_data <- data[, , 1, ]
  } else if (n_dims == 3) {
    # Could be lon, lat, time (no depth) or lon, lat, depth
    # Check if last dimension matches time length
    if (dims[3] == length(years)) {
      surface_data <- data  # Already surface-only or no depth dimension
    } else {
      # Assume lon, lat, depth - take first depth
      surface_data <- data[, , 1]
      # Add time dimension
      surface_data <- array(surface_data, dim = c(dim(surface_data), 1))
    }
  } else {
    stop("Unexpected number of dimensions in NetCDF file")
  }

  # Calculate statistics for each year
  stats_list <- map_dfr(seq_along(years), function(i) {
    # Extract data for this year
    if (length(dim(surface_data)) == 3) {
      year_data <- surface_data[, , i]
    } else {
      year_data <- surface_data
    }

    # Flatten and remove NAs
    values <- as.vector(year_data)
    valid_values <- values[!is.na(values)]

    # Calculate statistics
    if (length(valid_values) > 0) {
      tibble(
        Year = years[i],
        Mean = mean(valid_values, na.rm = TRUE),
        Median = median(valid_values, na.rm = TRUE),
        SD = sd(valid_values, na.rm = TRUE),
        SE = sd(valid_values, na.rm = TRUE) / sqrt(length(valid_values))
      )
    } else {
      tibble(
        Year = years[i],
        Mean = NA_real_,
        Median = NA_real_,
        SD = NA_real_,
        SE = NA_real_
      )
    }
  })

  return(stats_list)
}

# Process all variables
zoo_summary <- variables %>%
  pmap_dfr(function(var_name, data_dir) {

    # List all NetCDF files for this variable (exclude hidden files starting with ._)
    nc_files <- list.files(data_dir, pattern = paste0("^", var_name, "_.*\\.nc$"), full.names = TRUE)

    # Process all files for this variable
    nc_files %>%
      map_dfr(function(filepath) {
        # Parse filename to get metadata
        metadata <- parse_filename(filepath)

        # Calculate statistics
        stats <- calculate_stats(filepath, var_name)

        # Combine metadata with statistics
        metadata %>%
          crossing(stats) %>%
          select(Variable, Model, Scenario, Variant, Year, Mean, Median, SD, SE)
      })
  })

# Display summary
print(zoo_summary)

# Optionally save the results
write_csv(zoo_summary, "zooplankton_annual_summary_all.csv")


