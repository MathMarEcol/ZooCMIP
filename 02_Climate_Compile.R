library(hotrstuff)
library(tidyverse)

base_dir <- file.path("", "Volumes", "T9")

var <- "zmeso"

model_list <- list.files(file.path(base_dir, var, "raw"), pattern = ".nc", full.names = TRUE) %>%
  purrr::map(hotrstuff:::htr_get_CMIP6_bits) %>%
  bind_rows() %>%
  pull(Model) %>%
  unique()

model_list

# TODO hotrstuff doesn't error if no files or directory exists
# TODO Does hotrstuff warn if dates are missing

htr_merge_files(
  indir = file.path(base_dir, var, "raw"), # input directory
  outdir = file.path(base_dir, var, "merged"), # output directory
  year_start = 1950, # earliest year across all the scenarios considered (e.g., historical, ssp126, ssp245, ssp585)
  year_end = 2100, # latest year across all the scenarios considered
  ncores = 3,
)

# Warning x17 with zmicro
# cdi  warning (set_coordinates_varids): Coordinates variable lev can't be assigned!


#
# htr_change_freq(
#   freq = "yearly",
#   indir = file.path(base_dir, var, "merged"), # input directory
#   outdir = file.path(base_dir, var, "annual"),
#   ncores = 3
# )
#
#
# htr_regrid_esm(
#   indir = file.path(base_dir, var, "annual"),
#   outdir = file.path(base_dir, var, "regridded"),
#   cell_res = 0.5,
#   layer = "annual",
#   ncores = 3
# )






# htr_create_ensemble(
#   indir = file.path(base_dir, var, "regridded"), # input directory
#   outdir = file.path(base_dir, var, "ensemble"), # output directory
#   model_list = model_list, # list of models for ensemble
#   variable = var, # variable name
#   freq = "Omon", # original frequency of data
#   scenario = "ssp245", # emission scenario
#   mean = TRUE # if false, takes the median
# )


# htr_create_ensemble(
#   indir = file.path(base_dir, var, "regridded"), # input directory
#   outdir = file.path(base_dir, var, "ensemble"), # output directory
#   model_list = model_list, # list of models for ensemble
#   variable = var, # variable name
#   freq = "Omon", # original frequency of data
#   scenario = "ssp370", # emission scenario
#   mean = TRUE # if false, takes the median
# )
#
#
# htr_create_ensemble(
#   indir = file.path(base_dir, var, "regridded"), # input directory
#   outdir = file.path(base_dir, var, "ensemble"), # output directory
#   model_list = model_list, # list of models for ensemble
#   variable = var, # variable name
#   freq = "Omon", # original frequency of data
#   scenario = "ssp585", # emission scenario
#   mean = TRUE # if false, takes the median
# )
#


