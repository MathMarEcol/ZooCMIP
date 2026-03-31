# Load required libraries
library(tidyverse)
library(terra)
library(tidyterra)
library(patchwork)
library(sf)

# Define input directory
spatial_dir <- "Data/spatial_change_maps"

# Read the processing summary to get list of available files
spatial_files <- read_csv(file.path(spatial_dir, "spatial_processing_summary.csv")) %>%
  filter(Status == "Success", !is.na(Output_file))

# Define scenario order and labels
scenario_order <- c(
  # "ssp119",
  "ssp126",
  "ssp245",
  "ssp370",
  # "ssp434",
  # "ssp460",
  "ssp534-over",
  "ssp585"
)
scenario_labels <- c(
  # "SSP1-1.9",
  "SSP1-2.6",
  "SSP2-4.5",
  "SSP3-7.0",
  # "SSP4-3.4",
  # "SSP4-6.0",
  "SSP5-3.4-over",
  "SSP5-8.5"
)


# Create a lookup for scenario labels
scenario_label_map <- setNames(scenario_labels, scenario_order)

# Load world map for overlay (simplified for faster plotting)
world <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")

# Function to create a single spatial plot
plot_spatial_change <- function(filepath, model_name, scenario,
                                show_xlab = FALSE, show_ylab = FALSE) {

  # Load raster
  r <- rast(file.path(spatial_dir, filepath))

  # Create scenario label
  scenario_label <- scenario_label_map[scenario]
  if (is.na(scenario_label)) scenario_label <- scenario

  # Create plot
  p <- ggplot() +
    geom_spatraster(data = r) +
    geom_sf(data = world, color = "gray30", linewidth = 0.2, fill = NA) +
    scale_fill_gradientn(
      colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                 "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
      limits = c(-50, 50),
      na.value = "transparent",
      name = "Change in\nbiomass (%)",
      oob = scales::squish
    ) +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_line(color = "gray90", linewidth = 0.2),
      panel.background = element_rect(fill = "gray95", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_blank(),
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.margin = margin(2, 2, 2, 2)
    )

  # Add axis labels if specified
  if (show_xlab) {
    p <- p + labs(x = scenario_label) +
      theme(axis.title.x = element_text(size = 9, face = "bold", hjust = 0.5))
  }

  if (show_ylab) {
    p <- p + labs(y = model_name) +
      theme(axis.title.y = element_text(size = 9, face = "bold", angle = 90))
  }

  return(p)
}

# Function to create a grid of plots for one variable
create_variable_grid <- function(var_name) {

  # Filter files for this variable
  var_files <- spatial_files %>%
    filter(Variable == var_name) %>%
    # Filter to only include scenarios in our order
    filter(Scenario %in% scenario_order) %>%
    # Order by model and scenario
    arrange(Model, match(Scenario, scenario_order))

  # Get unique models and scenarios that have data
  available_models <- unique(var_files$Model) %>% sort()
  available_scenarios <- scenario_order[scenario_order %in% unique(var_files$Scenario)]

  # Create a complete grid of all combinations
  grid_data <- expand_grid(
    Model = available_models,
    Scenario = available_scenarios
  ) %>%
    left_join(var_files, by = c("Model", "Scenario"))

  # Calculate grid dimensions
  n_models <- length(available_models)
  n_scenarios <- length(available_scenarios)

  # Create all plots
  plot_list <- list()

  for (i in seq_len(nrow(grid_data))) {
    row <- grid_data[i, ]

    # Determine position in grid
    model_idx <- which(available_models == row$Model)
    scenario_idx <- which(available_scenarios == row$Scenario)

    # Determine if labels should be shown
    show_xlab <- (model_idx == n_models)  # Bottom row
    show_ylab <- (scenario_idx == 1)      # Left column

    if (!is.na(row$Output_file)) {
      # Create plot
      p <- plot_spatial_change(
        row$Output_file,
        row$Model,
        row$Scenario,
        show_xlab = show_xlab,
        show_ylab = show_ylab
      )
    } else {
      # Create empty placeholder
      p <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No data", size = 3, color = "gray50") +
        theme_void() +
        theme(plot.margin = margin(2, 2, 2, 2))

      if (show_xlab) {
        scenario_label <- scenario_label_map[row$Scenario]
        p <- p + labs(x = scenario_label) +
          theme(axis.title.x = element_text(size = 9, face = "bold"))
      }

      if (show_ylab) {
        p <- p + labs(y = row$Model) +
          theme(axis.title.y = element_text(size = 9, face = "bold", angle = 90))
      }
    }

    plot_list[[i]] <- p
  }

  # Create variable label
  var_label <- case_when(
    var_name == "zooc" ~ "All Zooplankton",
    var_name == "zmicro" ~ "Microzooplankton",
    var_name == "zmeso" ~ "Mesozooplankton",
    TRUE ~ var_name
  )

  # Combine plots using patchwork
  combined <- wrap_plots(plot_list, ncol = n_scenarios, nrow = n_models) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste(var_label, "- Biomass Change (1995-2014 to 2080-2100)"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    ) &
    theme(legend.position = "right")

  return(list(
    plot = combined,
    n_models = n_models,
    n_scenarios = n_scenarios
  ))
}

# ============================================
# CREATE PLOTS FOR EACH VARIABLE
# ============================================

# Get list of variables that have data
available_variables <- unique(spatial_files$Variable)
message(paste("Available variables:", paste(available_variables, collapse = ", ")))

# Variables to process
variables <- c("zooc", "zmicro", "zmeso")

# Create and save plots for each variable that has data
for (var in variables) {
  
  # Check if this variable has data
  if (!var %in% available_variables) {
    message(paste("Skipping", var, "- no spatial data files found"))
    next
  }
  
  message(paste("Creating spatial plot for", var))

  result <- create_variable_grid(var)

  # Display plot
  print(result$plot)

  # Calculate appropriate dimensions
  plot_width <- 3 * result$n_scenarios + 1  # Extra space for legend
  plot_height <- 2.5 * result$n_models + 1  # Extra space for title

  # Save plot
  output_file <- file.path("Figures", paste0("zooplankton_spatial_change_", var, ".png"))
  ggsave(output_file, result$plot,
         width = plot_width, height = plot_height, dpi = 300)

  message(paste("Saved:", output_file))
}
