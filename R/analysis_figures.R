library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(recipes)
library(rsample)
library(xgboost)
library(ggplot2)
library(ggpubr)
library(lubridate)

source("R/functions_analysis.R")

analysis_dir <- "analysis"
model_dir <- "MLs"
path <- file.path(analysis_dir, model_dir)

rds_dens <- read_rds(file.path(path, "med_density_xgBoostAnalysis.rds"))
rds_var <- read_rds(file.path(path, "var_density_xgBoostAnalysis.rds"))

# ylab <- c("Bias (class)", "Bias (reg)", "Percent error", "NRMSE")
ylab <- c("Median density", "Variance (density)")

model_dir <- "betaSurvival_uniqueAreaTrapSnare"
path <- file.path(analysis_dir, model_dir)
data <- read_rds(file.path(path, "abundanceScoresByPrimaryPeriod.rds")) |>
  ungroup() |>
  filter(med_density > 0) |>
  filter(var_density > 0) |>
  mutate(
    ppID = paste0(property_id, "-", PPNum),
    methods_used = as.character(methods_used)
  )

df_model_med <- subset_rename(data, "med_density", 1100)
df_model_var <- subset_rename(data, "var_density", 1100)

property_info <- data |>
  select(
    PPNum,
    property,
    property_area,
    property_id,
    density,
    med_density,
    var_density
  )

split_test_dens <- df_model_med$test |>
  rename(log_med_density = y) |>
  left_join(property_info)

split_test_var <- df_model_var$test |>
  rename(log_var_density = y) |>
  left_join(property_info)

known_density <- left_join(split_test_dens, split_test_var) |>
  mutate(
    log_pred_density = rds_dens$pred$test$pred,
    pred_density = exp(log_pred_density),
    log_pred_var = rds_var$pred$test$pred,
    pred_var = exp(log_pred_var)
  ) |>
  rename(known_density = density) |>
  select(
    PPNum,
    property_id,
    property_area,
    known_density,
    med_density,
    var_density,
    log_med_density,
    log_var_density,
    log_pred_density,
    pred_density,
    log_pred_var,
    pred_var
  )

my_theme <- function() {
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    title = element_text(size = 14)
  )
}

my_plot <- function(x, y) {
  r2 <- cor(
    pull(known_density, all_of(x)),
    pull(known_density, all_of(y))
  )^2

  known_density |>
    ggplot() +
    aes(x = .data[[x]], y = .data[[y]]) +
    geom_point(alpha = 1 / 10) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(subtitle = paste("R2 =", round(r2, 2))) +
    theme_bw() +
    my_theme()
}

g1 <- my_plot("log_med_density", "log_pred_density") +
  labs(
    x = "Log median density (from Bayes model)",
    y = "Log predicted density (from ML model)",
    title = "Log density"
  )

g2 <- my_plot("med_density", "pred_density") +
  labs(
    x = "Median density (from Bayes model)",
    y = "Predicted density (from ML model)",
    title = "Real density"
  )


gg <- ggarrange(g1, g2, nrow = 1)
gg
ggsave("density_bayes_ml_pred_obs.jpeg", path = "Plots/MLs", dpi = "print")


g3 <- my_plot("log_var_density", "log_pred_var") +
  labs(
    x = "Log posterior variance (from Bayes model)",
    y = "Log predicted variance (from ML model)",
    title = "Log variance"
  )

g4 <- my_plot("var_density", "pred_var") +
  labs(
    x = "Posterior variance (from Bayes model)",
    y = "Predicted variance (from ML model)",
    title = "Real variance"
  )

gg <- ggarrange(g3, g4, nrow = 1)
gg
ggsave("variance_bayes_ml_pred_obs.jpeg", path = "Plots/MLs", dpi = "print")


pred_truth |>
  mutate(resid = log(pred) - log(obs)) |>
  ggplot() +
  aes(x = obs, y = resid) +
  geom_point(alpha = 1 / 20) +
  geom_hline(yintercept = 0, color = "red") +
  labs(
    x = "Observed",
    y = "Residuals",
    title = ylab
    # subtitle = paste("R2 =", round(rds$pred$r2, 2))
  ) +
  theme_bw() +
  my_theme()


range(known_density$PPNum)
start_dates <- seq(ymd("2021-05-01"), by = paste(4, "week"), length.out = 40)
end_dates <- c(start_dates[-1] - 1, start_dates[length(start_dates)] + 28)

timestep_df <- tibble(end_dates) |>
  mutate(PPNum = 1:n())

density_time <- known_density |>
  left_join(timestep_df)

stat_by_window <- function(i, dates, window, value, stat) {
  startdate <- dates[i] - window
  enddate <- dates[i]
  interval <- seq(startdate, enddate, 1)

  tmp <- value[dates %in% interval]
  if (stat == "average") {
    return(mean(tmp))
  }
  if (stat == "sd") return(sd(tmp))
}

apply_sbw <- function(dates, window, value) {
  avg <- sapply(1:length(dates), function(x) {
    stat_by_window(x, dates, window, value, stat = "average")
  })
  avg[is.nan(avg)] <- NA
  sd <- sapply(1:length(dates), function(x) {
    stat_by_window(x, dates, window, value, stat = "sd")
  })
  sd[is.nan(sd)] <- NA

  return(data.frame(avg = avg, sd = sd))
}

properties <- unique(density_time$property_id)

# standardized density index
data_standard_index <- tibble()
pb <- txtProgressBar(min = 1, max = length(properties), style = 3)
for (i in seq_along(properties)) {
  property_df <- density_time |>
    filter(property_id == properties[i])

  if (nrow(property_df) < 15) {
    next
  }

  end_dates <- unique(property_df$end_dates)

  values <- property_df$pred_density
  property_df$ml_year0.5_avg <- apply_sbw(end_dates, 183, values)$avg
  property_df$ml_year0.5_sd <- apply_sbw(end_dates, 183, values)$sd
  property_df$ml_year1_avg <- apply_sbw(end_dates, 365, values)$avg
  property_df$ml_year1_sd <- apply_sbw(end_dates, 365, values)$sd
  property_df$ml_year3_avg <- apply_sbw(end_dates, 365 * 3, values)$avg
  property_df$ml_year3_sd <- apply_sbw(end_dates, 365 * 3, values)$sd

  values <- property_df$med_density
  property_df$bayes_year0.5_avg <- apply_sbw(end_dates, 183, values)$avg
  property_df$bayes_year0.5_sd <- apply_sbw(end_dates, 183, values)$sd
  property_df$bayes_year1_avg <- apply_sbw(end_dates, 365, values)$avg
  property_df$bayes_year1_sd <- apply_sbw(end_dates, 365, values)$sd
  property_df$bayes_year3_avg <- apply_sbw(end_dates, 365 * 3, values)$avg
  property_df$bayes_year3_sd <- apply_sbw(end_dates, 365 * 3, values)$sd

  property_df <- property_df |>
    mutate(
      ml_sdi_year0.5 = (pred_density - ml_year0.5_avg) / ml_year0.5_sd,
      ml_sdi_year1 = (pred_density - ml_year1_avg) / ml_year1_sd,
      ml_sdi_year3 = (pred_density - ml_year3_avg) / ml_year3_sd,
      bayes_sdi_year0.5 = (med_density - bayes_year0.5_avg) / bayes_year0.5_sd,
      bayes_sdi_year1 = (med_density - bayes_year1_avg) / bayes_year1_sd,
      bayes_sdi_year3 = (med_density - bayes_year3_avg) / bayes_year3_sd
    )

  data_standard_index <- bind_rows(data_standard_index, property_df)
  setTxtProgressBar(pb, i)
}
close(pb)

color1 <- "darkgreen"
color2 <- "darkblue"

filter_id <- data_standard_index |>
  group_by(property_id) |>
  count() |>
  ungroup() |>
  # filter(n == max(n)) |>
  pull(property_id)

filter_id <- c(
  "1.475-78-92",
  "5-10-1"
)

for (i in seq_along(filter_id)) {
  gg1 <- data_standard_index |>
    filter(property_id == filter_id[i]) |>
    mutate(
      sd = sqrt(var_density),
      ymin = pmax(0, med_density - (1.96 * sd)),
      ymax = med_density + (1.96 * sd)
    ) |>
    ggplot() +
    aes(x = end_dates, y = med_density, ymin = ymin, ymax = ymax) +
    geom_linerange(color = color1, alpha = 0.4, linewidth = 4) +
    geom_point(aes(shape = "Median density"), color = color1, size = 2) +
    geom_point(
      aes(y = known_density, shape = "Simulated (true) density"),
      size = 2
    ) +
    geom_line(aes(y = known_density), linetype = "dashed") +
    labs(
      title = "Bayes density estimates",
      x = "Date",
      y = "Density (pigs / sq. km)",
      shape = ""
    ) +
    # coord_cartesian(ylim = c(0, 3)) +
    theme_bw() +
    my_theme()

  ylim <- layer_scales(gg1)$y$get_limits()

  gg2 <- data_standard_index |>
    filter(property_id == filter_id[i]) |>
    mutate(
      sd = sqrt(pred_var),
      ymin = pmax(0, pred_density - (1.96 * sd)),
      ymax = pred_density + (1.96 * sd)
    ) |>
    ggplot() +
    aes(x = end_dates, y = pred_density, ymin = ymin, ymax = ymax) +
    geom_linerange(color = color2, alpha = 0.4, linewidth = 4) +
    geom_point(aes(shape = "Median density"), color = color2, size = 2) +
    geom_point(
      aes(y = known_density, shape = "Simulated (true) density"),
      size = 2
    ) +
    geom_line(aes(y = known_density), linetype = "dashed") +
    labs(title = "ML density estiamtes", x = "Date", y = "", shape = "") +
    coord_cartesian(ylim = ylim) +
    theme_bw() +
    my_theme()
  ggg <- ggarrange(gg1, gg2, nrow = 1, common.legend = TRUE, legend = "bottom")
  ggsave(
    paste0("timeseries_", filter_id[i], ".jpeg"),
    path = "Plots/MLs",
    dpi = "print"
  )
}


cols <- c(
  "#f03b20",
  "#feb24c",
  "#ffeda0",
  "white",
  "#edf8b1",
  "#7fcdbb",
  "#2c7fb8"
)
cats <- c(
  "Very High",
  "High",
  "Above average",
  "Average",
  "Below average",
  "Low",
  "Very low"
)

cats <- tibble(
  ymax = c(2.5, 2, 1.5, 1, -1, -1.5, -2),
  ymin = c(2, 1.5, 1, -1, -1.5, -2, -2.5),
  xmax = max(plot_df$end_dates),
  xmin = min(plot_df$end_dates),
  cat = factor(
    cats,
    levels = c(
      "Very High",
      "High",
      "Above average",
      "Average",
      "Below average",
      "Low",
      "Very low"
    )
  ),
  cols = cols
)

for (i in seq_along(filter_id)) {
  plot_df <- data_standard_index |>
    filter(property_id == filter_id[i])

  plot_df |>
    pivot_longer(
      cols = contains("_sdi_"),
      names_to = "method_time",
      values_to = "SDI"
    ) |>
    filter(!is.na(SDI)) |>
    mutate(
      method = if_else(grepl("bayes", method_time), "Bayes", "ML"),
      window = if_else(grepl("year0.5", method_time), "6 month", "x"),
      window = if_else(grepl("year1", method_time), "1 year", window),
      window = if_else(grepl("year3", method_time), "3 year", window)
    ) |>
    mutate(
      window = factor(window, levels = c("6 month", "1 year", "3 year"))
    ) |>
    ggplot() +
    aes(x = end_dates, y = SDI, linetype = method) +
    geom_rect(
      cats,
      inherit.aes = FALSE,
      mapping = aes(
        ymin = ymin,
        ymax = ymax,
        xmin = xmin,
        xmax = xmax,
        fill = cat,
        color = cat
      ),
      alpha = 0.5
    ) +
    geom_line(linewidth = 1.5) +
    # geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    facet_wrap(~window) +
    labs(
      x = "Date",
      y = "Standardized Density Index",
      title = "Standardized Density Index",
      linetype = "Estimation\nmethod",
      fill = "Density\ncategory",
      color = "Density\ncategory"
    ) +
    theme_bw() +
    my_theme() +
    theme(strip.text = element_text(size = 12))

  ggsave(
    paste0("SDI_", filter_id[i], ".jpeg"),
    path = "Plots/MLs",
    dpi = "print"
  )
}
