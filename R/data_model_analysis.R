library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(purrr)

analysis_dir <- "analysis"
model_dir <- "betaSurvival_uniqueAreaTrapSnare"
path <- file.path(analysis_dir, model_dir)
density_dirs <- paste0("density_", c(0.3, 1.475, 2.65, 3.825, 5))

map_files2 <- function(dirs_vec, file_name) {
  get_files <- function(density_dir, file_name, node) {
    sim_results <- file.path(path, density_dir)
    ls <- read_rds(file.path(sim_results, file_name))
  }

  dirs_vec |>
    map(\(x) get_files(x, file_name)) |>
    list_rbind() |>
    mutate(start_density = as.factor(start_density))
}

scores_rds <- read_rds(file.path(path, "abundanceScoresByPrimaryPeriod.rds")) |>
  select(property_id, PPNum, norm_bias_density, nm_rmse_density)

f_name <- "predicted_p_summaries.rds"
df <- map_files2(density_dirs, f_name) |>
  mutate(
    density = N / property_area,
    density_category = if_else(density <= 2, "Low", "Medium"),
    density_category = if_else(density >= 6, "High", density_category),
    density_category = if_else(density == 0, "Extinct", density_category),
    density_category = factor(
      density_category,
      levels = c("Extinct", "Low", "Medium", "High")
    ),
    recovered = if_else(p >= low_p & p <= high_p, 1, 0),
    property_id = paste(start_density, simulation, property, sep = "-")
  )


p_df <- left_join(df, scores_rds)

p_df |>
  group_by(density_category) |>
  summarise(prop_recovered = sum(recovered) / n())

my_summary <- function(df) {
  df |>
    summarise(
      low = quantile(value, 0.05),
      q1 = quantile(value, 0.25),
      med = quantile(value, 0.5),
      q3 = quantile(value, 0.75),
      high = quantile(value, 0.95),
      n = n()
    )
}

# p_df |>
#   filter(N > 0) |>
#   ggplot() +
#   aes(x = effort_per, y = p) +
#   geom_point(shape = ".") +
#   geom_smooth(method = "lm") +
#   facet_wrap(~ method, scales = "free_x") +
#   theme_bw()

# df |>
#   mutate(value = mbias_p) |>
#   group_by(density_category) |>
#   my_summary()
#   ggplot() +
#   aes(x = density_category, y = mbias_p) +
#   geom_boxplot() +
#   theme_bw()

quintiles <- p_df |>
  filter(N > 0) |>
  group_by(method) |>
  summarise(
    `20%` = quantile(effort_per, 0.2),
    `40%` = quantile(effort_per, 0.4),
    `60%` = quantile(effort_per, 0.6),
    `80%` = quantile(effort_per, 0.8)
  )

get_quintile <- function(m, q) {
  quintiles |>
    filter(method == m) |>
    pull(q)
}

create_effort_category <- function(m) {
  p_df |>
    filter(N > 0, method == m) |>
    mutate(
      effort_category = if_else(
        effort_per <= get_quintile(m, "20%"),
        "0-20%",
        "80-100%"
      ),
      effort_category = if_else(
        effort_per > get_quintile(m, "20%") &
          effort_per <= get_quintile(m, "40%"),
        "20-40%",
        effort_category
      ),
      effort_category = if_else(
        effort_per > get_quintile(m, "40%") &
          effort_per <= get_quintile(m, "60%"),
        "40-60%",
        effort_category
      ),
      effort_category = if_else(
        effort_per > get_quintile(m, "60%") &
          effort_per <= get_quintile(m, "80%"),
        "60-80%",
        effort_category
      ),
    )
}

firearms <- create_effort_category("Firearms")
fixed <- create_effort_category("Fixed wing")
helicopter <- create_effort_category("Helicopter")
snares <- create_effort_category("Snares")
traps <- create_effort_category("Traps")

cats <- bind_rows(firearms, fixed, helicopter, snares, traps) |>
  mutate(
    effort_category = factor(
      effort_category,
      levels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%")
    )
  ) |>
  mutate(
    p_category = if_else(p <= 0.2, "0-0.2", "0.8-1"),
    p_category = if_else(p > 0.2 & p <= 0.4, "0.2-0.4", p_category),
    p_category = if_else(p > 0.4 & p <= 0.6, "0.4-0.6", p_category),
    p_category = if_else(p > 0.6 & p <= 0.8, "0.6-0.8", p_category),
  ) |>
  select(
    method,
    density_category,
    effort_category,
    p_category,
    norm_bias_density,
    nm_rmse_density
  ) |>
  pivot_longer(
    cols = c(norm_bias_density, nm_rmse_density),
    names_to = "metric"
  )


my_tile <- function(df, m, x, y) {
  tmp <- df |>
    filter(metric == m) |>
    group_by(method, .data[[x]], .data[[y]]) |>
    summarise(v = median(value), n = n()) |>
    ungroup()

  l1 <- abs(min(tmp$v))
  l2 <- abs(max(tmp$v))
  ll <- max(c(l1, l2)) * 1.05

  if (grepl("bias", m)) {
    limits <- round(c(-1 * ll, ll), 2)
    breaks <- c(limits[1], 0, limits[2])
    colors = c("#e9a3c9", "#f5f5f5", "#a1d76a")
  } else {
    limits <- round(c(0, ll), 2)
    breaks <- round(seq(0, ll, length.out = 3), 2)
    colors = c("navyblue", "darkmagenta", "darkorange1")
  }

  ggplot(tmp) +
    aes(x = .data[[x]], y = .data[[y]], fill = v) +
    geom_tile() +
    facet_grid(method ~ .) +
    scale_fill_gradientn(colors = colors, limits = limits, breaks = breaks) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14, vjust = -1.5)
    )
}

x <- "effort_category"
xlab <- "Density category"

y <- "p_category"
ylab <- "Detection rate"

m <- "norm_bias_density"
nb_lab <- "nBias"
g1 <- my_tile(cats, m, x, y) +
  labs(x = xlab, y = ylab, fill = nb_lab) +
  theme(strip.background.y = element_blank(), strip.text = element_blank())

m <- "nm_rmse_density"
nrmse_lab <- "NRMSE"
g2 <- my_tile(cats, m, x, y) +
  labs(x = xlab, y = "", fill = nrmse_lab) +
  theme(axis.text.y = element_blank())

ggarrange(g1, g2, ncol = 2, labels = "auto")


f_name <- "predicted_theta_summaries.rds"
df <- map_files2(density_dirs, f_name) |>
  mutate(
    density = N / property_area,
    density_category = if_else(density <= 2, "Low", "Medium"),
    density_category = if_else(density >= 6, "High", density_category),
    density_category = if_else(density == 0, "Extinct", density_category),
    density_category = factor(
      density_category,
      levels = c("Extinct", "Low", "Medium", "High")
    ),
    recovered = if_else(theta >= low_theta & theta <= high_theta, 1, 0),
    property_id = paste(start_density, simulation, property, sep = "-")
  )

theta_df <- left_join(df, scores_rds)

theta_df |>
  group_by(method) |>
  summarise(prop_recovered = sum(recovered) / n())

theta_df |>
  filter(N > 0) |>
  ggplot() +
  aes(x = theta, y = high_theta) +
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(density_category ~ method) +
  theme_bw()


X_p <- rnorm(3)
beta1 <- rnorm(5000)
beta_p <- matrix(rnorm(15000), 5000, 3)

xx1 <- log(
  nimble::ilogit(
    beta1 +
      X_p[1] * beta_p[, 1] +
      X_p[2] * beta_p[, 2] +
      X_p[3] * beta_p[, 3]
  )
)

summary(xx1)

xx2 <- rep(NA, length(beta1))
for (i in 1:length(beta1)) {
  xx2[i] <- log(
    nimble::ilogit(
      beta1[i] + nimble::inprod(X_p, beta_p[i, ])
      # X_p[1] * beta_p[i, 1] +
      # X_p[2] * beta_p[i, 2] +
      # X_p[3] * beta_p[i, 3]
    )
  )
}
summary(xx2)


f_name <- "predicted_take_summaries.rds"
df <- map_files2(density_dirs, f_name) |>
  mutate(
    density = N / property_area,
    density_category = if_else(density <= 2, "Low", "Medium"),
    density_category = if_else(density >= 6, "High", density_category),
    density_category = if_else(density == 0, "Extinct", density_category),
    density_category = factor(
      density_category,
      levels = c("Extinct", "Low", "Medium", "High")
    ),
    recovered = if_else(take >= low_take & take <= high_take, 1, 0),
    property_id = paste(start_density, simulation, property, sep = "-")
  )

take_df <- left_join(df, scores_rds)

take_df |>
  group_by(density_category) |>
  summarise(prop_recovered = sum(recovered) / n())


f_name <- "predicted_area_summaries.rds"
df <- map_files2(density_dirs, f_name) |>
  mutate(
    density = N / property_area,
    density_category = if_else(density <= 2, "Low", "Medium"),
    density_category = if_else(density >= 6, "High", density_category),
    density_category = if_else(density == 0, "Extinct", density_category),
    density_category = factor(
      density_category,
      levels = c("Extinct", "Low", "Medium", "High")
    ),
    recovered = if_else(
      potential_area >= low_area & potential_area <= high_area,
      1,
      0
    ),
    property_id = paste(start_density, simulation, property, sep = "-")
  )

area_df <- left_join(df, scores_rds)

area_df |>
  group_by(method) |>
  summarise(prop_recovered = sum(recovered) / n())
