
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
    Model    = factor(Model)
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
  zooc   = "#1b7837",   # dark green  — total zooplankton
  zmicro = "#762a83",   # purple      — microzooplankton
  zmeso  = "#d6604d"    # red-orange  — mesozooplankton
)

p_new <- ggplot(plot_dat,
  aes(x = Scenario, y = Mean_percent_change, colour = Variable, fill = Variable)) +

  # --- background: individual model points (small, semi-transparent) -------
  geom_point(
    aes(group = Variable),
    position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.25, seed = 42),
    shape  = 16,
    size   = 1.5,
    alpha  = 0.35
  ) +

  # --- foreground: SD errorbars --------------------------------------------
  geom_errorbar(
    data = summary_dat,
    aes(y    = mean_pct,
        ymin = mean_pct - sd_pct,
        ymax = mean_pct + sd_pct,
        group = Variable),
    position = dodge,
    width    = 0.25,
    linewidth = 0.8,
    alpha    = 0.9
  ) +

  # --- foreground: mean point ----------------------------------------------
  geom_point(
    data = summary_dat,
    aes(y = mean_pct, group = Variable),
    position = dodge,
    shape    = 21,
    size     = 3.5,
    stroke   = 0.8,
    colour   = "white",
    fill     = NA          # overridden below via scale_fill_manual
  ) +
  # Filled mean point on top (separate layer so fill shows through white stroke)
  geom_point(
    data = summary_dat,
    aes(y = mean_pct, group = Variable),
    position = dodge,
    shape    = 21,
    size     = 3.5,
    stroke   = 0.8
  ) +

  # --- reference line at zero ----------------------------------------------
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +

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
  labs(
    x       = NULL,
    y       = "Change in zooplankton biomass (%)",
    fill = NULL,
    colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = c(0.1, 0.1),
    legend.direction = "vertical",
    panel.grid.major.x = element_blank()
  )

ggsave("Figures/PercentChange.pdf", plot = p_new, width = 10, height = 6)
ggsave("Figures/PercentChange.png", dpi = 600, plot = p_new, width = 10, height = 6)

# ===========================================================================
# ORIGINAL CODE (commented out — preserved for reference)
# ===========================================================================

# library(tidyverse)
# library(patchwork)
# p <- list()
# 
# # Keep each variable's change-summary tibble around (Model/Scenario/
# # Output_file etc.) for re-use in the spatial-map section below, rather than
# # re-reading the same three CSVs a second time.
# all_dat <- list()
# 
# vars <- c("zooc", "zmicro", "zmeso")
# # vars <- "zooc"
# 
# for (var in vars) {
#   
#   
#   dat <- readr::read_csv(
#     file.path("Data", paste0(var, "_change_summary.csv"))) %>%
#     filter(!str_detect(Model, "CanESM5")) %>% 
#     mutate(Scenario = factor(Scenario),
#     Model = factor(Model))
#     
#     # BUG FIX: aes() uses lazy (delayed) evaluation — mapping a variable
#     # name (e.g. `aes(y = label_y)`) does NOT capture that variable's
#     # current *value* at the moment the layer is built. It only stores a
#     # reference to the name `label_y`, which gets looked up again whenever
#     # the plot is actually drawn/rendered. Since this script loops over
#     # `vars` and reassigns `label_y` (and `y_min`/`y_max`) on every
#     # iteration, all three plots stored in the `p[[var]]` list end up
#     # sharing the SAME final value of `label_y` — whatever it was set to on
#     # the LAST loop iteration — once you print/render them after the loop
#     # finishes. This is why disabling other layers made the mismatch
#     # visible: labels appeared to "disappear" (rendered far off-panel)
#     # because they were using zmeso's y-range, not zooc's/zmicro's own.
#     #
#     # Fix: compute label_y from each var's own data ONCE, then store it as
#     # a column IN `dat` (not a loose variable) before building the plot.
#     # `ggplot(dat, ...)` evaluates and snapshots `dat` immediately, so each
#     # plot object keeps its own independent copy of label_y baked in as
#     # data — completely decoupled from later loop iterations.
#     
#     # Each point needs its own vertical "axis-style" label (the Model name)
#     # sitting below the x axis. Two things are important here:
#     #
#     #   1. Spacing: `position_jitter()` draws a random uniform offset
#     #      independently per point, so with only a dozen or so models per
#     #      Scenario, random chance can leave some points clumped together and
#     #      others spread out (this is expected — jitter is NOT evenly spaced).
#     #      `position_dodge2()` instead deterministically spreads all points
#     #      sharing an x value (Scenario) evenly across the available width,
#     #      based on their `group` (here, Model), giving consistent spacing
#     #      every time.
#     #   2. Alignment: geom_point() and geom_text() must share the *same*
#     #      position object so each label lines up under its point.
#     #
#     #   3. The label itself is anchored at y = -Inf (bottom of the panel) and
#     #      rotated 90 degrees so it reads vertically like an axis tick label.
#     #      `coord_cartesian(clip = "off")` stops it being clipped at the panel
#     #      border, and extra bottom plot margin gives it room to render.
#     dodge_pos <- position_dodge2(width = 0.9, padding = 0)
#     
#     # Vertical dividers between Scenario groups: since Scenario is a discrete
#     # x-axis factor, its levels sit at integer positions 1, 2, 3, ... To draw
#     # a line *between* each pair of adjacent scenarios (rather than through
#     # them), place `geom_vline()` at the half-integer positions 1.5, 2.5, ...
#     n_scenarios <- nlevels(dat$Scenario)
#     scenario_breaks <- seq(1.5, n_scenarios - 0.5, by = 1)
#     
#     # Anchor all Model labels at a single, fixed y position below the axis,
#     # expressed as a real data-coordinate offset (not `y = -Inf` combined
#     # with hjust/vjust > 1). Pushing labels via hjust/vjust beyond 1 offsets
#     # each one along its own rotated text-bounding-box, which is scaled to
#     # that label's individual string length — since Model names differ in
#     # length, this produces inconsistent, "jittered"-looking vertical
#     # offsets. Using a shared numeric y value with standard hjust = 1,
#     # vjust = 1 (the normal alignment for rotated axis-style labels) instead
#     # anchors every label's top-right corner at the exact same y, giving
#     # perfectly consistent alignment regardless of text length.
#     #
#     # The visible y range now needs to cover the full errorbar extents
#     # (Mean +/- SD), not just the point means, otherwise the blue error
#     # bars would be clipped at the panel border.
#     # y_min <- min(dat$Mean_percent_change - dat$SD_percent_change, na.rm = TRUE)
#     # y_max <- max(dat$Mean_percent_change + dat$SD_percent_change, na.rm = TRUE)
#     
#     y_min <- min(dat$Mean_percent_change, na.rm = TRUE)
#     y_max <- max(dat$Mean_percent_change, na.rm = TRUE)
#     
#     y_range <- y_max - y_min
#     label_y <- min(dat$Mean_percent_change, na.rm = TRUE) - 0.3 * y_range
#     
#     # Bake label_y into `dat` as a real column, rather than leaving it as a
#     # loose variable referenced by aes(). See the note above the loop for
#     # why this matters: aes() mappings are evaluated lazily against the
#     # data/environment at render time, and a loop body reuses the same
#     # environment on every iteration, so a loose variable like `label_y`
#     # would resolve to whatever the LAST loop iteration set it to once all
#     # three plots are rendered after the loop finishes. Storing it as a
#     # column makes it part of the `dat` snapshot that `ggplot(dat, ...)`
#     # captures immediately, so each plot keeps its own correct value.
#     dat <- dat %>% mutate(label_y = label_y)
#     all_dat[[var]] <- dat
#     # Per-scenario mean of Mean_percent_change, kept on the same discrete
#     # Scenario x-axis used everywhere else in this plot (rather than
#     # converting Scenario to a numeric position, which would conflict with
#     # the discrete x scale already established by aes(x = Scenario)).
#     scenario_means <- dat %>%
#     group_by(Scenario) %>%
#     summarise(mean_change = mean(Mean_percent_change, na.rm = TRUE), .groups = "drop")
#     
#     p[[var]] <- ggplot(dat, 
#       aes(x = Scenario, 
#         y = Mean_percent_change, 
#         label = Model, 
#         group = Model)) +
#         geom_vline(xintercept = scenario_breaks, color = "grey70") +
#         # geom_crossbar with ymin = ymax = mean_change collapses the box to a
#         # single horizontal line at the mean, spanning the full Scenario
#         # category width (width = 0.9 matches the dodge width used by the
#         # points, so the line lines up with the spread of points above it).
#         geom_crossbar(
#           data = scenario_means,
#           aes(x = Scenario, y = mean_change, ymin = mean_change, ymax = mean_change),
#           inherit.aes = FALSE,
#           width = 0.9,
#           middle.linewidth = 0.8,
#           color = "firebrick"
#         ) +
#         # Error bars showing +/- 1 SD around each point's mean, drawn in blue.
#         # Sharing `dodge_pos` with geom_point() keeps each error bar centered
#         # exactly on its corresponding point.
#         # geom_errorbar(
#         #   aes(ymin = Mean_percent_change - SD_percent_change,
#         #       ymax = Mean_percent_change + SD_percent_change),
#         #   position = dodge_pos,
#         #   width = 0.3,
#         #   color = "blue"
#         # ) +
#         geom_point(position = dodge_pos) +
#         geom_text(
#           aes(y = label_y),
#           position = dodge_pos,
#           angle = 45,
#           hjust = 1,
#           vjust = 1,
#           size = 3
#         ) +
#         labs(x = "Scenario", y = "Mean % change") +
#         ggtitle(var) +
#         # Without pinning ylim, adding a geom_text() below the data range would
#         # cause ggplot to auto-expand the y scale to include the label position,
#         # shifting the whole panel (and the axis line itself) down — which is
#         # what caused the axis to appear misaligned. `coord_cartesian(ylim = ...)`
#         # keeps the visible axis anchored to the true data range, while
#         # `clip = "off"` still allows the labels to render below it.
#         coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
#         theme_bw() +
#         theme(
#           axis.text.x = element_text(margin = margin(t = 5)),
#           axis.title.x = element_blank(),
#           plot.margin = margin(t = 5.5, r = 5.5, b = 80, l = 5.5)
#         )
#         
#       }
#       
#       out <- wrap_plots(p, ncol = 1) 
#       
#       ggsave("Figures/PercentChange.pdf", plot = out, width = 12, height = 8)
