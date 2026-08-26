# Load required libraries
library(tidyverse)
library(patchwork)

# Read the summary data
zoo_data <- read_csv("Data/zooplankton_annual_summary_all.csv") %>%
  filter(Year >= 1950, Year <= 2100) %>%
  filter(Scenario %in% c("historical",
                         "ssp119",
                         "ssp126",
                         "ssp245",
                         "ssp370",
                         "ssp434",
                         "ssp460",
                         "ssp534-over",
                         "ssp585"))

# Calculate baseline (mean biomass for 1995-2014 historical period)
baseline <- zoo_data %>%
  filter(Scenario == "historical", Year >= 1995, Year <= 2014) %>%
  group_by(Variable, Model) %>%
  summarise(baseline_mean = mean(Mean, na.rm = TRUE), .groups = "drop")

# Join baseline and calculate percentage change
zoo_change <- zoo_data %>%
  left_join(baseline, by = c("Variable", "Model")) %>%
  mutate(
    pct_change = ((Mean - baseline_mean) / baseline_mean) * 100,
    scenario_type = if_else(Scenario == "historical", "Historical", "Future")
  )

# Convert Scenario to factor with all unique scenarios in sensible order
# This ensures all scenarios appear in the legend even if not in every subplot
all_scenarios <- unique(zoo_change$Scenario)
scenario_order <- c("historical", "ssp119", "ssp126", "ssp245", "ssp370",
                   "ssp434", "ssp460", "ssp534-over", "ssp585")
# Only include scenarios that exist in the data
scenario_levels <- scenario_order[scenario_order %in% all_scenarios]

zoo_change <- zoo_change %>%
  mutate(Scenario = factor(Scenario, levels = scenario_levels))

# Calculate mean and SD across models for each scenario, variable, and year
zoo_summary <- zoo_change %>%
  group_by(Variable, Scenario, scenario_type, Year) %>%
  summarise(
    mean_change = mean(pct_change, na.rm = TRUE),
    sd_change = sd(pct_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    upper = mean_change + sd_change,
    lower = mean_change - sd_change,
    # Ensure Scenario is a factor with same levels as zoo_change
    Scenario = factor(Scenario, levels = scenario_levels)
  )

# Define colors for scenarios (distinctive progression)
scenario_colors <- c(
  "historical" = "#8B7355",      # Brown/tan for historical
  "ssp119" = "#2E7D32",          # Dark green (lowest emissions)
  "ssp126" = "#66BB6A",          # Light green
  "ssp245" = "#1976D2",          # Blue (moderate)
  "ssp370" = "#FFA726",          # Orange (high)
  "ssp434" = "#FF6F00",          # Dark orange
  "ssp460" = "#EF5350",          # Red-orange
  "ssp534-over" = "#D32F2F",     # Red
  "ssp585" = "#6A1B9A"           # Dark purple (highest emissions)
)

# ============================================
# INDIVIDUAL MODEL PLOTS (3 plots, one per variable)
# ============================================

# Function to create plot for a single model
plot_single_model <- function(data, model_name) {
  plot_data <- data %>% filter(Model == model_name)

  # Create a complete grid of all scenarios and years to ensure all legend items appear
  # even if a model doesn't have data for certain scenarios
  all_years <- unique(data$Year)
  scenario_grid <- expand.grid(
    Year = all_years,
    Scenario = scenario_levels,
    stringsAsFactors = FALSE
  ) %>%
    mutate(Scenario = factor(Scenario, levels = scenario_levels))

  # Left join to preserve all scenario levels in the data
  plot_data <- scenario_grid %>%
    left_join(plot_data, by = c("Year", "Scenario"))

  ggplot(plot_data, aes(x = Year, y = pct_change,
                        color = Scenario, fill = Scenario)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(
      values = scenario_colors,
      limits = scenario_levels,
      breaks = scenario_levels,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = scenario_colors,
      limits = scenario_levels,
      breaks = scenario_levels,
      drop = FALSE
    ) +
    guides(color = guide_legend(title = NULL, override.aes = list(linewidth = 2)),
           fill = guide_legend(title = NULL)) +
    labs(
      x = "Year",
      y = "Change in biomass (%)",
      title = model_name
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.title = element_text(size = 9),
      legend.position = "bottom"
    )
}

# Variables to process
variables <- c("zooc", "zmicro", "zmeso")

# Create 3 separate plots, one for each variable
individual_plots <- list()

for (var in variables) {
  var_label <- case_when(
    var == "zooc" ~ "All Zooplankton",
    var == "zmicro" ~ "Microzooplankton",
    var == "zmeso" ~ "Mesozooplankton"
  )

  # Filter data for this variable and get models with data
  var_data <- zoo_change %>% filter(Variable == var)
  models_with_data <- var_data %>%
    filter(!is.na(pct_change)) %>%
    pull(Model) %>%
    unique()

  # Create list of plots for each model that has data for this variable
  model_plots <- map(models_with_data, ~plot_single_model(var_data, .x))

  # Combine into grid with 3 columns
  n_models <- length(models_with_data)
  individual_plots[[var]] <- wrap_plots(model_plots, ncol = 3) +
    plot_layout(guides = "collect") +
    plot_annotation(title = var_label,
                    theme = theme(plot.title = element_text(hjust = 0.5,
                                                            face = "bold",
                                                            size = 14))) &
    theme(legend.position = "bottom")

  # Store number of models for height calculation
  individual_plots[[paste0(var, "_n")]] <- n_models
}

# Display and save individual model plots
for (var in variables) {
  n_models <- individual_plots[[paste0(var, "_n")]]
  print(individual_plots[[var]])
  ggsave(file.path("Figures", paste0("zooplankton_individual_models_", var, ".png")),
         individual_plots[[var]],
         width = 12, height = 3 * ceiling(n_models / 3), dpi = 300)
}

# ============================================
# MULTI-MODEL MEAN PLOTS
# ============================================

# Function to create plot for a single variable
plot_variable <- function(data, var_name) {
  var_label <- case_when(
    var_name == "zooc" ~ "All Zooplankton",
    var_name == "zmicro" ~ "Microzooplankton",
    var_name == "zmeso" ~ "Mesozooplankton"
  )

  plot_data <- data %>% filter(Variable == var_name)

  ggplot(plot_data, aes(x = Year, y = mean_change,
                        color = Scenario, fill = Scenario)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(
      values = scenario_colors,
      limits = scenario_levels,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = scenario_colors,
      limits = scenario_levels,
      drop = FALSE
    ) +
    guides(color = guide_legend(title = NULL),
           fill = guide_legend(title = NULL)) +
    labs(
      x = "Year",
      y = "Change in biomass (%)",
      title = var_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

# Create individual plots
p_zooc <- plot_variable(zoo_summary, "zooc")
p_zmicro <- plot_variable(zoo_summary, "zmicro")
p_zmeso <- plot_variable(zoo_summary, "zmeso")

# Combine plots using patchwork
combined_plot <- p_zooc / p_zmicro / p_zmeso +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Display the plot
print(combined_plot)

# Save the plot
ggsave(file.path("Figures", "zooplankton_temporal_change.png"), combined_plot,
       width = 10, height = 12, dpi = 300)
