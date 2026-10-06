library(tidyverse)

# ---------------------------------------------------------------------------
# NEW PLOT: Single panel — SSPs on x-axis, mean % change on y-axis,
# grouped by variable (zooc / zmicro / zmeso) within each SSP.
# Summary layer: mean ± SD errorbar + mean point per variable × SSP.
# Background layer: individual model points (jittered, semi-transparent).
# ---------------------------------------------------------------------------

vars <- c("zooc", "zmicro", "zmeso")

# Read the single merged CSV (written by 06_Climate_Review.R).
# Rows with Model == "All models" are the pre-computed multimodel mean ± SD
# summary rows; all other rows are individual model results.
all_dat <- readr::read_csv(
  file.path("Data", "zooplankton_change_summary.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Variable = factor(Variable, levels = vars),
    Scenario = factor(Scenario),
    Model = factor(Model)
  )

# Only keep the four "headline" SSPs that appear in the sketch
headline_ssps <- c("ssp126", "ssp245", "ssp370", "ssp585")

# Individual model points — exclude the "All models" summary rows
plot_dat <- all_dat %>%
  filter(Scenario %in% headline_ssps, Model != "All models") %>%
  mutate(
    Scenario = factor(Scenario, levels = headline_ssps),
    Variable = factor(Variable, levels = vars)
  )

# Summary layer: pull the pre-computed multimodel mean ± SD directly from
# the "All models" rows rather than recomputing them here
summary_dat <- all_dat %>%
  filter(Scenario %in% headline_ssps, Model == "All models") %>%
  mutate(
    Scenario = factor(Scenario, levels = headline_ssps),
    Variable = factor(Variable, levels = vars)
  ) %>%
  rename(mean_pct = Mean_percent_change, sd_pct = SD_percent_change)

# Dodge position shared by all layers so individual points, errorbars and
# mean points all align within each SSP group
dodge <- position_dodge(width = 0.6)

# Colour palette for the three variables
var_colours <- c(
  zooc = "#1b7837", # dark green  — total zooplankton
  zmicro = "#762a83", # purple      — microzooplankton
  zmeso = "#d6604d" # red-orange  — mesozooplankton
)

p_new <- ggplot(
  plot_dat,
  aes(x = Scenario, y = Mean_percent_change, colour = Variable, fill = Variable)
) +

  # --- background: individual model points (small, semi-transparent) -------
  geom_point(
    aes(group = Variable),
    position = position_jitterdodge(
      dodge.width = 0.6,
      jitter.width = 0.25,
      seed = 42
    ),
    shape = 16,
    size = 1.5,
    alpha = 0.35
  ) +

  # --- foreground: SD errorbars --------------------------------------------
  geom_errorbar(
    data = summary_dat,
    aes(
      y = mean_pct,
      ymin = mean_pct - sd_pct,
      ymax = mean_pct + sd_pct,
      group = Variable
    ),
    position = dodge,
    width = 0.25,
    linewidth = 0.8,
    alpha = 0.9
  ) +

  # --- foreground: mean point ----------------------------------------------
  geom_point(
    data = summary_dat,
    aes(y = mean_pct, group = Variable),
    position = dodge,
    shape = 21,
    size = 3.5,
    stroke = 0.8,
    colour = "white",
    fill = NA # overridden below via scale_fill_manual
  ) +
  # Filled mean point on top (separate layer so fill shows through white stroke)
  geom_point(
    data = summary_dat,
    aes(y = mean_pct, group = Variable),
    position = dodge,
    shape = 21,
    size = 3.5,
    stroke = 0.8
  ) +

  # --- reference line at zero ----------------------------------------------
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.5) +

  # --- scales & labels -----------------------------------------------------
  scale_colour_manual(
    values = var_colours,
    labels = c(zooc = "Total", zmicro = "Small", zmeso = "Large")
  ) +
  scale_fill_manual(
    values = var_colours,
    labels = c(zooc = "Total", zmicro = "Small", zmeso = "Large")
  ) +
  scale_x_discrete(
    labels = c(
      ssp126 = "Low (SSP1-2.6)",
      ssp245 = "Medium (SSP2-4.5)",
      ssp370 = "High (SSP3-7.0)",
      ssp585 = "Very High (SSP5-8.5)"
    )
  ) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5),
    colour = "grey80",
    linewidth = 0.4,
    linetype = "solid"
  ) +
  labs(
    x = NULL,
    y = "Change in zooplankton biomass (%)",
    fill = NULL,
    colour = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = c(0.1, 0.1),
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave("Figures/PercentChange.pdf", plot = p_new, width = 10, height = 6)
ggsave(
  "Figures/PercentChange.png",
  dpi = 600,
  plot = p_new,
  width = 10,
  height = 6
)
