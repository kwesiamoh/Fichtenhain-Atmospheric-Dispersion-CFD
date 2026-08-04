# ============================================================
# Figure 5: Directional receptor concentration comparison
# ============================================================

library(ggplot2)
library(scales)


# ============================================================
# File paths
# ============================================================

base_dir <- paste0(
  "C:/acfd2/Fichtenhain_Dispersion/",
  "Results/quantitative/windDirection"
)

input_csv <- file.path(
  base_dir,
  "Figure_5_receptor_values.csv"
)

output_png <- file.path(
  base_dir,
  "Figure_5.png"
)

output_pdf <- file.path(
  base_dir,
  "Figure_5.pdf"
)


# ============================================================
# Read and validate data
# ============================================================

if (!file.exists(input_csv)) {
  stop(
    paste(
      "Input CSV was not found:",
      input_csv
    )
  )
}

data <- read.csv(
  input_csv,
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8",
  check.names = FALSE
)

required_columns <- c(
  "scenario",
  "scenario_order",
  "sample_time_s",
  "receptor",
  "receptor_order",
  "radius_m",
  "height_m",
  "T"
)

missing_columns <- setdiff(
  required_columns,
  names(data)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (nrow(data) != 32) {
  warning(
    paste(
      "Expected 32 rows but found",
      nrow(data)
    )
  )
}

data$T <- as.numeric(data$T)

if (any(!is.finite(data$T))) {
  stop(
    "The concentration column contains missing or non-finite values."
  )
}


# ============================================================
# Scenario and receptor ordering
# ============================================================

scenario_levels <- c(
  "Wind from 0°",
  "Wind from 90°",
  "Wind from 180°",
  "Prevailing wind from 154°"
)

scenario_labels <- c(
  "Wind from 0°" =
    "Wind from 0° (N to S)",

  "Wind from 90°" =
    "Wind from 90° (E to W)",

  "Wind from 180°" =
    "Wind from 180° (S to N)",

  "Prevailing wind from 154°" =
    "Prevailing wind from 154° (SSE to NNW)"
)

receptor_levels <- c(
  "N",
  "NE",
  "E",
  "SE",
  "S",
  "SW",
  "W",
  "NW"
)

missing_scenarios <- setdiff(
  scenario_levels,
  unique(data$scenario)
)

if (length(missing_scenarios) > 0) {
  stop(
    paste(
      "The CSV is missing these scenarios:",
      paste(missing_scenarios, collapse = ", ")
    )
  )
}

data$scenario_display <- scenario_labels[
  data$scenario
]

# Reverse the factor levels so Wind from 0° appears at the top.
data$scenario_display <- factor(
  data$scenario_display,
  levels = rev(
    unname(scenario_labels)
  )
)

data$receptor <- factor(
  data$receptor,
  levels = receptor_levels
)


# ============================================================
# Reporting threshold
# ============================================================

reporting_threshold <- 1.0e-7

data$plot_T <- ifelse(
  data$T >= reporting_threshold,
  data$T,
  NA_real_
)

detected_data <- data[
  !is.na(data$plot_T),
]

non_detected_data <- data[
  is.na(data$plot_T),
]

if (nrow(detected_data) == 0) {
  stop(
    "No receptor values exceeded the reporting threshold."
  )
}


# ============================================================
# Logarithmic colour limits
# ============================================================

colour_minimum <- reporting_threshold

colour_maximum <- max(
  detected_data$plot_T
) * 1.05

minimum_exponent <- ceiling(
  log10(colour_minimum)
)

maximum_exponent <- floor(
  log10(colour_maximum)
)

legend_breaks <- 10^seq(
  minimum_exponent,
  maximum_exponent,
  by = 1
)


# ============================================================
# Scientific notation labels
# ============================================================

make_scientific_label <- function(value) {

  exponent <- floor(
    log10(value)
  )

  coefficient <- value / (
    10^exponent
  )

  sprintf(
    "%.1f %%*%% 10^{%d}",
    coefficient,
    exponent
  )
}

detected_data$cell_label <- vapply(
  detected_data$T,
  make_scientific_label,
  character(1)
)


# ============================================================
# Annotation text colour
# ============================================================

detected_data$relative_log_position <- (
  log10(detected_data$plot_T) -
    log10(colour_minimum)
) / (
  log10(colour_maximum) -
    log10(colour_minimum)
)

detected_data$text_colour <- ifelse(
  detected_data$relative_log_position >= 0.58,
  "#111111",
  "#FFFFFF"
)


# ============================================================
# Console summary
# ============================================================

cat(
  "\nWind-direction receptor summary\n"
)

cat(
  paste0(
    strrep("=", 76),
    "\n"
  )
)

for (scenario_name in scenario_levels) {

  scenario_data <- data[
    data$scenario == scenario_name,
  ]

  maximum_row <- scenario_data[
    which.max(scenario_data$T),
  ]

  cat(
    sprintf(
      "%-29s | receptor = %-2s | T = %.6e\n",
      scenario_name,
      as.character(maximum_row$receptor),
      maximum_row$T
    )
  )
}


# ============================================================
# Create figure
# ============================================================

figure_5 <- ggplot(
  data,
  aes(
    x = receptor,
    y = scenario_display
  )
) +

  # Light-grey background for values below the threshold
  geom_tile(
    fill = "#EFEFEF",
    colour = "#FFFFFF",
    linewidth = 0.8
  ) +

  # Colour only detected concentrations
  geom_tile(
    data = detected_data,
    aes(
      fill = plot_T
    ),
    colour = "#FFFFFF",
    linewidth = 0.8
  ) +

  # Scientific notation labels
  geom_text(
    data = detected_data,
    aes(
      label = cell_label,
      colour = text_colour
    ),
    parse = TRUE,
    size = 3.25,
    fontface = "plain",
    show.legend = FALSE
  ) +

  # Subtle non-detect labels
  geom_text(
    data = non_detected_data,
    label = "ND",
    colour = "#7A7A7A",
    size = 3.0,
    fontface = "plain"
  ) +

  scale_colour_identity() +

  # Cividis is perceptually uniform and colourblind safe
  scale_fill_viridis_c(
    option = "C",
    direction = 1,
    trans = "log10",
    limits = c(
      colour_minimum,
      colour_maximum
    ),
    breaks = legend_breaks,
    labels = trans_format(
      "log10",
      math_format(10^.x)
    ),
    oob = squish,
    name = expression(
      "Normalized gaseous contaminant concentration, " *
        italic(T) * " (-)"
    )
  ) +

  labs(
    x = expression(
      "Receptor direction at 200 m radius and " *
        italic(z) * " = 70 m"
    ),

    y = "Wind scenario",

    caption = paste0(
      "ND denotes ",
      "T < 1 × 10⁻⁷. ",
      "Values were sampled at t = 1800 s."
    )
  ) +

  coord_fixed(
    ratio = 0.72,
    clip = "off"
  ) +

  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = grid::unit(
        10.0,
        "cm"
      ),
      barheight = grid::unit(
        0.38,
        "cm"
      ),
      ticks = TRUE,
      frame.colour = "#444444"
    )
  ) +

  theme_classic(
    base_family = "sans",
    base_size = 10.5
  ) +

  theme(
    axis.title.x = element_text(
      size = 11,
      margin = margin(
        t = 10
      )
    ),

    axis.title.y = element_text(
      size = 11,
      margin = margin(
        r = 12
      )
    ),

    axis.text.x = element_text(
      size = 10,
      colour = "#111111",
      margin = margin(
        t = 5
      )
    ),

    axis.text.y = element_text(
      size = 9.5,
      colour = "#111111",
      margin = margin(
        r = 5
      )
    ),

    axis.ticks = element_blank(),

    axis.line = element_blank(),

    panel.border = element_rect(
      colour = "#222222",
      fill = NA,
      linewidth = 0.75
    ),

    legend.position = "bottom",

    legend.title = element_text(
      size = 9.5,
      margin = margin(
        b = 5
      )
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.caption = element_text(
      size = 8.5,
      colour = "#4D4D4D",
      hjust = 0,
      margin = margin(
        t = 8
      )
    ),

    plot.margin = margin(
      t = 10,
      r = 12,
      b = 8,
      l = 10
    )
  )


# ============================================================
# Save publication outputs
# ============================================================

ggsave(
  filename = output_png,
  plot = figure_5,
  width = 7.6,
  height = 4.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = figure_5,
  width = 7.6,
  height = 4.5,
  units = "in",
  device = "pdf",
  bg = "white"
)


# ============================================================
# Display figure
# ============================================================

print(
  figure_5
)

cat(
  "\nSaved outputs:\n"
)

cat(
  paste(
    output_png,
    "\n"
  )
)

cat(
  paste(
    output_pdf,
    "\n"
  )
)