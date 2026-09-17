# ZooCMIP: Prepare Files
# Parse CMIP6 zooplankton model filenames and check coverage
# Filename structure: VARIABLE_TIMESCALE_MODEL_EXPERIMENT_VARIANT_GRID_START-ENDTIME.nc

library(tidyverse)

# List all .nc files in the raw folder
files <- list.files("/Volumes/T9/ClimateData/zmeso/regridded/", pattern = "\\.nc$", full.names = FALSE)

# Parse the filenames into components
df <- tibble(filename = files) %>%
  # Remove .nc extension and split by underscore
  mutate(
    basename = str_remove(filename, "\\.nc$"),
    parts = str_split(basename, "_")
  ) %>%
  # Extract components
  mutate(
    variable = map_chr(parts, 1),
    timescale = map_chr(parts, 2),
    model = map_chr(parts, 3),
    experiment = map_chr(parts, 4),
    variant = map_chr(parts, 5),
    grid = map_chr(parts, 6),
    timerange = map_chr(parts, 7)
  ) %>%
  # Split start and end time
  separate(timerange, into = c("start_time", "end_time"), sep = "-") %>%
  select(-basename, -parts)

# Check GRID and VARIANT consistency within each model
consistency_check <- df %>%
  group_by(model) %>%
  summarise(
    n_variants = n_distinct(variant),
    variants = paste(unique(variant), collapse = ", "),
    n_grids = n_distinct(grid),
    grids = paste(unique(grid), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_variants > 1 | n_grids > 1)

if (nrow(consistency_check) > 0) {
  warning("Inconsistent VARIANT or GRID found within models:")
  print(consistency_check)
}

# Summary by model, experiment, variant, and time range
summary_df <- df %>%
  group_by(model, experiment, variable, variant, grid) %>%
  summarise(
    start_time = min(start_time),
    end_time = max(end_time),
    n_files = n(),
    .groups = "drop"
  ) %>%
  arrange(model, variable, experiment, variant)

# View results
print(summary_df)

# ---- Check for temporal gaps in coverage ----
# Helper function to add one month to YYYYMM string
add_one_month <- function(yyyymm) {
  year <- as.integer(substr(yyyymm, 1, 4))
  month <- as.integer(substr(yyyymm, 5, 6))
  month <- month + 1
  if (month > 12) {
    month <- 1
    year <- year + 1
  }
  sprintf("%04d%02d", year, month)
}

# Check for gaps within each variable/model/experiment/variant combination
time_ranges <- df %>%
  group_by(variable, model, experiment, variant) %>%
  arrange(start_time, .by_group = TRUE) %>%
  mutate(
    row_num = row_number(),
    prev_end = lag(end_time, default = NA_character_),
    expected_start = map_chr(prev_end, ~ if (is.na(.x)) NA_character_ else add_one_month(.x)),
    has_gap = row_num > 1 & start_time != expected_start
  ) %>%
  ungroup()

# Extract gaps between files (only checking gaps between consecutive files, not first file)
gaps <- time_ranges %>%
  filter(has_gap) %>%
  select(variable, model, experiment, variant, expected_start, actual_start = start_time)

if (nrow(gaps) > 0) {
  warning("Temporal gaps detected - missing files between these periods:")
  print(gaps)
} else {
  message("No temporal gaps detected - all variable/model/experiment combinations have continuous coverage.")
}

# ---- Check for late starts ----
# Historical experiments should start by 1950, other experiments by 2015
late_starts <- df %>%
  group_by(variable, model, experiment, variant) %>%
  summarise(first_start = min(start_time), .groups = "drop") %>%
  mutate(
    required_start = if_else(experiment == "historical", "19500101", "20150101"),
    starts_late = first_start > required_start
  ) %>%
  filter(starts_late) %>%
  select(variable, model, experiment, variant, required_start, actual_start = first_start)

if (nrow(late_starts) > 0) {
  warning("Late starts detected - data begins after expected start date:")
  print(late_starts)
} else {
  message("No late starts detected - all series begin on time.")
}
