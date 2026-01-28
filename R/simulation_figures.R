library(tidyverse)
library(ggpubr)

source("R/functions_collate.R")

analysis_dir <- "analysis"
model_dir <- "betaSurvival_uniqueAreaTrapSnare"
density_dirs <- paste0("density_", c(0.3, 1.475, 2.65, 3.825, 5))

map_files2 <- function(dirs_vec, file_name){

  get_files <- function(density_dir, file_name, node){
    sim_results <- file.path(analysis_dir, model_dir, density_dir)
    ls <- read_rds(file.path(sim_results, file_name))
  }

  dirs_vec |>
    map(\(x) get_files(x, file_name)) |>
    list_rbind() |>
    mutate(start_density = as.factor(start_density))

}

f_name <- "all_param_samples.rds"
samples <- map_files2(density_dirs, f_name)

f_name <- "method_parameter_lookup.rds"
known_params <- map_files2(density_dirs, f_name)

f_name <- "all_beta_p.rds"
beta <- map_files2(density_dirs, f_name)

## join params to known

beta1_long <- samples |>
  select_pivot_longer("beta1") |>
  mutate(method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         position = 1)

beta_1_known <- left_join(beta1_long, beta) |>
  mutate(value = value - actual)

my_summary <- function(df){
  df |>
    summarise(low = quantile(value, 0.05),
              q1 = quantile(value, 0.25),
              med = quantile(value, 0.5),
              q3 = quantile(value, 0.75),
              high = quantile(value, 0.95),
              mu = mean(value),
              sd = sd(value))
}

method_table <- known_params |>
  select(idx, method) |>
  distinct() |>
  rename(method_idx = idx)

my_theme <- function(){
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))
}

my_linerange <- function(df){
  ggplot(df) +
    aes(x = med, xmin = low, xmax = high, y = method) +
    geom_linerange(linewidth = 2) +
    geom_linerange(aes(xmin = q1, xmax = q3), linewidth = 4) +
    geom_point(size = 7) +
    geom_point(size = 5, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(x = "Residual",
         y = "Method") +
    theme_bw() +
    my_theme()
}

# do not need to group by start density - all the same
b1_summary <- beta_1_known |>
  group_by(method_idx) |>
  my_summary()

b1_summary |>
  ungroup() |>
  left_join(method_table) |>
  my_linerange() +
  labs(title = "Intercept")

ggsave("Plots/beta1_residual.jpeg", dpi = "retina", units = "in", width = 6, height = 4)


beta_p_long <- samples |>
  select_pivot_longer("beta_p") |>
  mutate(method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         position = as.numeric(str_extract(node, "(?<=\\, )\\d")) + 1)

beta_p_known <- left_join(beta_p_long, beta) |>
  mutate(value = value - actual)

f_name <- "land_cover_lookup.rds"
land_cover <- map_files2(density_dirs, f_name)

land_cover_table <- tibble(
  position = 1:4,
  landCover = c("Intercept", "Road density", "Ruggedness", "Canopy cover")
)

# do not need to group by start density - all the same
beta_p_known |>
  group_by(method_idx, position) |>
  my_summary() |>
  ungroup() |>
  bind_rows(mutate(b1_summary, position = 1)) |>
  left_join(method_table) |>
  left_join(land_cover_table) |>
  mutate(landCover = factor(landCover, levels = c("Intercept", "Road density", "Ruggedness", "Canopy cover"))) |>
  my_linerange() +
  facet_wrap(~ landCover)

ggsave("Plots/beta_all_residual.jpeg", dpi = "print")


gamma_long <- samples |>
  select_pivot_longer("log_gamma[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         idx = idx + 3)

gamma_known <- left_join(gamma_long, known_params) |>
  select(-rho, -p_unique) |>
  mutate(log_gamma = log(gamma),
         residual = value - log_gamma)

g1 <- gamma_known |>
  group_by(method) |>
  mutate(value = residual) |>
  my_summary() |>
  ungroup() |>
  my_linerange() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 12))

ggsave("Plots/gamma_residual.jpeg", dpi = "print")

gH <- gamma_known |>
  select(simulation, start_density, method, gamma, log_gamma) |>
  distinct()

g2 <- gamma_known |>
  group_by(simulation, start_density, method) |>
  mutate(value = exp(value)) |>
  my_summary() |>
  ungroup() |>
  left_join(gH) |>
  ggplot() +
  aes(x = gamma, y = med) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(method ~ ., scales = "free_y") +
  labs(x = "Known parameter value",
       y = "Posterior median") +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 12))

ggsave("Plots/gamma_medianVsKnown.jpeg", dpi = "print")

ggarrange(g1, g2, nrow = 1, labels = "AUTO")
ggsave("Plots/gamma_medianVsKnown.jpeg", dpi = "retina", units = "in", width = 6, height = 2)


rho_long <- samples |>
  select_pivot_longer("log_rho[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d")))

rho_known <- left_join(rho_long, known_params) |>
  select(-gamma, -p_unique) |>
  mutate(log_rho = log(rho),
         residual = value - log_rho)

r1 <- rho_known |>
  group_by(method) |>
  mutate(value = residual) |>
  my_summary() |>
  ungroup() |>
  my_linerange()

ggsave("Plots/rho_residual.jpeg", dpi = "print")

rH <- rho_known |>
  select(simulation, start_density, method, rho, log_rho) |>
  distinct()

r2 <- rho_known |>
  group_by(simulation, start_density, method) |>
  mutate(value = exp(value)) |>
  my_summary() |>
  ungroup() |>
  left_join(rH) |>
  ggplot() +
  aes(x = rho, y = med) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ method, scales = "free") +
  labs(x = "Known parameter value",
       y = "Posterior median") +
  theme_bw()

ggarrange(r1, r2, nrow = 1, labels = "AUTO")
ggsave("Plots/rho_medianVsKnown.jpeg", dpi = "retina", units = "in", width = 6, height = 2)

p_long <- samples |>
  select_pivot_longer("p_mu[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         idx = idx + 3)

p_known <- left_join(p_long, known_params) |>
  select(-gamma, -rho) |>
  mutate(logit_p = boot::logit(p_unique),
         residual = value - logit_p)

p1 <- p_known |>
  group_by(method) |>
  mutate(value = residual) |>
  my_summary() |>
  ungroup() |>
  my_linerange()

ggsave("Plots/p_residual.jpeg", dpi = "print")


pH <- p_known |>
  select(simulation, start_density, method, p_unique, logit_p) |>
  distinct()

p2 <- p_known |>
  group_by(simulation, start_density, method) |>
  mutate(value = boot::inv.logit(value)) |>
  my_summary() |>
  ungroup() |>
  left_join(pH) |>
  ggplot() +
  aes(x = p_unique, y = med) +
  geom_point(size = 0.25) +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(method ~ start_density, scales = "free") +
  labs(x = "Known parameter value",
       y = "Posterior median") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggarrange(p1, p2, nrow = 1, labels = "AUTO", widths = c(1, 2))
ggsave("Plots/p_medianVsKnown.jpeg", dpi = "retina", units = "in", width = 6, height = 2)


phi_long <- samples |>
  select_pivot_longer("phi_mu") |>
  mutate(actual = 0.78,
         value = value - actual)

g1 <- phi_long |>
  my_summary() |>
  mutate(method = "Survival") |>
  my_linerange() +
  labs(y = "",
       title = "Survival") +
  theme(axis.text.y = element_blank())

psi_long <- samples |>
  select_pivot_longer("psi_phi") |>
  mutate(actual = 5,
         value = value - actual)

g2 <- psi_long |>
  group_by(start_density) |>
  my_summary() |>
  ungroup() |>
  mutate(method = start_density) |>
  my_linerange() +
  labs(y = "Starting density",
       title = "Shrinkage")

nu_long <- samples |>
  select_pivot_longer("log_nu") |>
  mutate(actual = log(5.290323),
         value = value - actual)

g3 <- nu_long |>
  my_summary() |>
  mutate(method = "Fecundity") |>
  my_linerange() +
  labs(y = "",
       title = "Fecundity") +
  theme(axis.text.y = element_blank())

gg1 <- ggarrange(g1, g3, nrow = 1, ncol = 2, labels = c("A", "B"))
g4 <- blank <- ggplot() + geom_blank() + theme_void()
gg2 <- ggarrange(g2, g4, nrow = 1, ncol = 2, widths = c(100, 1), labels = "C")

ggarrange(gg1, gg2, nrow = 2)
ggsave("Plots/vitalRates.jpeg", dpi = "retina", units = "in", width = 6, height = 4)


recover_summary <- function(df){
  df |>
    summarise(percent_recovered = round(sum(recovered) / n() * 100, 2),
              n = n())
}

percent_recovered_beta <- bind_rows(beta1_long, beta_p_long) |>
  group_by(method_idx, start_density, simulation) |>
  my_summary() |>
  ungroup() |>
  left_join(beta) |>
  mutate(recovered = if_else(actual >= low & actual <= high, 1, 0)) |>
  group_by(method_idx, position) |>
  recover_summary() |>
  ungroup() |>
  left_join(method_table) |>
  mutate(Parameter = "Beta")

percent_recovered_gamma <- gamma_long |>
  group_by(idx, start_density, simulation) |>
  my_summary() |>
  ungroup() |>
  left_join(known_params) |>
  mutate(actual = log(gamma)) |>
  mutate(recovered = if_else(actual >= low & actual <= high, 1, 0)) |>
  group_by(idx) |>
  recover_summary() |>
  ungroup() |>
  rename(method_idx = idx) |>
  left_join(method_table) |>
  mutate(Parameter = "Gamma")

percent_recovered_rho <- rho_long |>
  group_by(idx, start_density, simulation) |>
  my_summary() |>
  ungroup() |>
  left_join(known_params) |>
  mutate(actual = log(rho)) |>
  mutate(recovered = if_else(actual >= low & actual <= high, 1, 0)) |>
  group_by(idx) |>
  recover_summary() |>
  ungroup() |>
  rename(method_idx = idx) |>
  left_join(method_table) |>
  mutate(Parameter = "Rho")

percent_recovered_p <- p_long |>
  group_by(idx, start_density, simulation) |>
  my_summary() |>
  ungroup() |>
  left_join(known_params) |>
  mutate(actual = boot::logit(p_unique)) |>
  mutate(recovered = if_else(actual >= low & actual <= high, 1, 0)) |>
  group_by(idx) |>
  recover_summary() |>
  ungroup() |>
  rename(method_idx = idx) |>
  left_join(method_table) |>
  mutate(Parameter = "P")

percent_recovered_phi <- samples |>
  select_pivot_longer("phi_mu") |>
  group_by(start_density, simulation) |>
  my_summary() |>
  ungroup() |>
  mutate(actual = 0.78) |>
  mutate(recovered = if_else(actual >= low & actual <= high, 1, 0)) |>
  recover_summary() |>
  ungroup() |>
  mutate(Parameter = "Phi")

percent_recovered_psi <- samples |>
  select_pivot_longer("psi_phi") |>
  group_by(start_density, simulation) |>
  my_summary() |>
  ungroup() |>
  mutate(actual = 5) |>
  mutate(recovered = if_else(actual >= low & actual <= high, 1, 0)) |>
  recover_summary() |>
  ungroup() |>
  mutate(Parameter = "Psi")

percent_recovered_nu <- samples |>
  select_pivot_longer("log_nu") |>
  group_by(start_density, simulation) |>
  my_summary() |>
  ungroup() |>
  mutate(actual = log(5.290323)) |>
  mutate(recovered = if_else(actual >= low & actual <= high, 1, 0)) |>
  recover_summary() |>
  ungroup() |>
  mutate(Parameter = "Nu")

bind_rows(
  percent_recovered_beta,
  percent_recovered_gamma,
  percent_recovered_rho,
  percent_recovered_p,
  percent_recovered_phi,
  percent_recovered_psi,
  percent_recovered_nu
) |> View()


abundance_file <- "abundance_summaries.rds"
f_name <- abundance_file
df <- map_files2(density_dirs, f_name)

df |>
  ggplot() +
  aes(x = density, y = med_density) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ start_density) +
  theme_bw() +
  my_theme()

df |>
  filter(abundance > 0) |>
  mutate(recovered = if_else(abundance >= low_abundance &
                               abundance <= high_abundance,
                             1, 0)) |>
  group_by(start_density) |>
  recover_summary()

df |>
  filter(abundance == 0) |>
  mutate(recovered = if_else(abundance >= low_abundance &
                               abundance <= high_abundance,
                             1, 0)) |>
  group_by(start_density) |>
  recover_summary()

abundance_file <- "abundance_error_by_property.rds"
f_name <- abundance_file
df <- map_files2(density_dirs, f_name)

# create property ID for easier joining
property_ids <- df |>
  select(start_density, simulation, property) |>
  distinct() |>
  mutate(property_id = paste(start_density, simulation, property, sep = "-"))


df_property <- df |>
  left_join(property_ids) |>
  select(property_id, simulation, property, start_density, contains("density"), -density) |>
  distinct()


abundance_file <- "all_take.rds"
f_name <- abundance_file
df <- map_files2(density_dirs, f_name)

n_methods <- df |>
  left_join(property_ids) |>
  select(property_id, method) |>
  distinct() |>
  group_by(property_id) |>
  count() |>
  rename(n_methods_used = n)

n_return <- df |>
  left_join(property_ids) |>
  select(property_id, method) |>
  distinct() |>
  pivot_wider(names_from = method,
              values_from = method) |>
  unite(method,
        -c(property_id),
        sep = ", ",
        na.rm = TRUE) |>
  rename(methods_used = method) |>
  left_join(n_methods)

observed_pp <- df |>
  left_join(property_ids) |>
  select(property_id, PPNum) |>
  distinct() |>
  count(property_id) |>
  rename(n_observed_pp = n)

ts_length <- df |>
  left_join(property_ids) |>
  select(property_id, PPNum) |>
  distinct() |>
  group_by(property_id) |>
  filter(PPNum == min(PPNum) |
           PPNum == max(PPNum)) |>
  mutate(ts_length = c(0, diff(PPNum) + 1)) |>
  ungroup() |>
  filter(ts_length != 0) |>
  select(-PPNum)

take_by_property <- df |>
  left_join(property_ids) |>
  group_by(property_id, property_area, start_density) |>
  summarise(take = sum(take),
         effort = sum(effort_per),
         unit_count = sum(trap_count),
         n_total_events = n()) |>
  ungroup() |>
  left_join(n_return) |>
  left_join(observed_pp) |>
  left_join(ts_length) |>
  mutate(proportion_observed = n_observed_pp / ts_length)


property_error <- left_join(take_by_property, df_property)
glimpse(property_error)


# trends
# x = proportion_observed, y = nm_rmse_density [NONE]
# x = proportion_observed, y = rmse_density [NONE]
# x = proportion_observed, y = mbias_density [NONE]
# x = proportion_observed, y = mpe_density [NONE]

# x = ts_length, y = nm_rmse_density [NONE]
# x = ts_length, y = rmse_density [NONE]
# x = ts_length, y = mbias_density [NONE]
# x = ts_length, y = mpe_density [NONE]

# x = n_observed_pp, y = nm_rmse_density [on average no trend, but more properties with high error at low n_observed_pp]
# x = n_observed_pp, y = rmse_density [on average no trend, but more properties with high error at low n_observed_pp]
# x = n_observed_pp, y = mbias_density [on average no trend, but more properties with high error at low n_observed_pp]
# x = n_observed_pp, y = mpe_density [on average no trend, but more properties with high error at low n_observed_pp]

# x = n_total_events, y = nm_rmse_density [slightly higher error at low total events]
# x = n_total_events, y = rmse_density [slightly higher error at low total events]
# x = n_total_events, y = mbias_density [slightly higher error at low total events]
# x = n_total_events, y = mpe_density [slightly higher error at low total events]

# x = unit_count, y = nm_rmse_density [more data points at low unit count but otherwise very small trend]
# x = unit_count, y = rmse_density [more data points at low unit count but otherwise very small trend]
# x = unit_count, y = mbias_density [more data points at low unit count but otherwise very small trend]
# x = unit_count, y = mpe_density [more data points at low unit count but otherwise very small trend]

# x = effort, y = nm_rmse_density [more data points at low unit count but otherwise very small trend]
# x = effort, y = rmse_density [more data points at low unit count but otherwise very small trend]
# x = effort, y = mbias_density [more data points at low unit count but otherwise very small trend]
# x = effort, y = mpe_density [more data points at low unit count but otherwise very small trend]

# x = take, y = nm_rmse_density [very low absolute take has very high error]
# x = take, y = rmse_density [very low absolute take has very high error]
# x = take, y = mbias_density [very low absolute take has very high error]
# x = take, y = mpe_density [very low absolute take has very high error]

# x = property_area, y = nm_rmse_density [very low property area has very high error]
# x = property_area, y = rmse_density [very low property area has very high error]
# x = property_area, y = mbias_density [very low property area has very high error]
# x = property_area, y = mpe_density [very low property area has very high error]

property_error |>
  # filter(take < 162) |>
  pivot_longer(cols = c(nm_rmse_density, rmse_density, mbias_density, mpe_density),
               names_to = "metric",
               values_to = "value") |>
  filter(metric %in% c("mpe_density", "nm_rmse_density")) |>
  ggplot() +
  aes(x = property_area, y = value) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ metric, scales = "free") +
  # coord_cartesian(xlim = c(200, 2000)) +
  theme_bw()

# need to to determine best properties and what is acceptable

# best 5% in metric
q <- 0.05
qb <- 0.10

quantile(property_error$mpe_density, q)

best_properties <- property_error |>
  filter(mpe_density <= quantile(mpe_density, q),
         nm_rmse_density <= quantile(nm_rmse_density, q),
         rmse_density <= quantile(rmse_density, q),
         mbias_density >= quantile(mbias_density, qb),
         mbias_density <= quantile(mbias_density, 1 - qb))

nrow(best_properties)
nrow(best_properties) / nrow(property_error) * 100

create_range_df <- function(df ,x){
  tibble(
    metric = x,
    min = df |> pull(x) |> min(),
    max = df |> pull(x) |> max(),
    qmin = df |> pull(x) |> quantile(0.05),
    qmax = df |> pull(x) |> quantile(0.95),
    median = df |> pull(x) |> median(),
    mean = df |> pull(x) |> mean(),
    sd = df |> pull(x) |> sd()
  )
}

best_metrics <- bind_rows(
  create_range_df(best_properties, "property_area"),
  create_range_df(best_properties, "take"),
  create_range_df(best_properties, "n_total_events"),
  create_range_df(best_properties, "n_observed_pp"),
  create_range_df(best_properties, "ts_length"),
  create_range_df(best_properties, "proportion_observed"),
  create_range_df(best_properties, "effort"),
  create_range_df(best_properties, "unit_count")
)
best_metrics

write_csv(best_metrics, "analysis/betaSurvival_uniqueAreaTrapSnare/best_metrics.csv")

table(best_properties$n_methods_used)
table(best_properties$methods_used)

summary(property_error$mpe_density)
summary(property_error$nm_rmse_density)
summary(property_error$rmse_density)
summary(property_error$mbias_density)

thresh_mpe <- 25
thresh_nrmse <- 0.5
thresh_rmse <- 1
thresh_bias <- 0.5

acceptable_properties <- property_error |>
  filter(mpe_density <= thresh_mpe,
         nm_rmse_density <= thresh_nrmse,
         # rmse_density <= thresh_rmse,
         mbias_density <= thresh_bias,
         mbias_density >= -1 * thresh_bias)

nrow(acceptable_properties)
nrow(acceptable_properties) / nrow(property_error) * 100


acceptable <- property_error |>
  mutate(acceptable = mpe_density <= thresh_mpe &
           nm_rmse_density <= thresh_nrmse &
           mbias_density <= thresh_bias &
           mbias_density >= -1 * thresh_bias)


abundance_file <- "abundance_summaries.rds"
f_name <- abundance_file
df <- map_files2(density_dirs, f_name)

property_means <- df |>
  # filter(abundance > 0) |>
  left_join(property_ids) |>
  group_by(property_id, start_density) |>
  summarise(mu_known = mean(density),
            mu_est = mean(med_density))

left_join(acceptable, property_means) |>
  select(acceptable, start_density) |>
  filter(acceptable) |>
  count(start_density)

df_property_means <- left_join(acceptable, property_means)

ggpoint <- function(df){
  ggplot(df) +
    aes(x = mu_known, y = mu_est) +
    geom_point(size = 0.5) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "Mean known density",
         y = "Mean predicted density") +
    coord_cartesian(xlim = c(0, 20),
                    ylim = c(0, 30)) +
    theme_bw() +
    my_theme()
}

g1 <- df_property_means |>
  ggpoint()
g2 <- df_property_means |>
  filter(acceptable) |>
  ggpoint()

ggarrange(g1, g2, ncol = 2, labels = "AUTO")

acceptable |>
  select(take, property_area, mpe_density, nm_rmse_density, mbias_density, acceptable) |>
  # pivot_longer(cols = c(mpe_density, nm_rmse_density, mbias_density),
  #              names_to = "error_metric",
  #              values_to = "error_value") |>
  # filter(error_metric == "nm_rmse_density") |>
  ggplot() +
  aes(x = nm_rmse_density, y = take, color = acceptable) +
  geom_point() +
  facet_wrap(~ acceptable) +
  theme_bw()

acceptable |>
  select(property_area, take, n_total_events, n_observed_pp,
         ts_length, proportion_observed, effort, unit_count, acceptable) |>
  pivot_longer(cols = -acceptable,
               names_to = "property_metric",
               values_to = "value") |>
  group_by(acceptable, property_metric) |>
  my_summary() |>
  filter(property_metric %in% c("property_area", "take", "proportion_observed")) |>
  ggplot() +
  aes(x = med, xmin = low, xmax = high, y = acceptable) +
  geom_linerange(linewidth = 2) +
  geom_linerange(aes(xmin = q1, xmax = q3), linewidth = 4) +
  geom_point(size = 7) +
  geom_point(size = 5, color = "white") +
  facet_wrap(~ property_metric, scales = "free") +
  labs(x = "Value",
       y = "Acceptable property") +
  theme_bw() +
  my_theme()





acceptable_metrics <- bind_rows(
  create_range_df(acceptable_properties, "property_area"),
  create_range_df(acceptable_properties, "take"),
  create_range_df(acceptable_properties, "n_total_events"),
  create_range_df(acceptable_properties, "n_observed_pp"),
  create_range_df(acceptable_properties, "ts_length"),
  create_range_df(acceptable_properties, "proportion_observed"),
  create_range_df(acceptable_properties, "effort"),
  create_range_df(acceptable_properties, "unit_count")
)

acceptable_metrics

write_csv(acceptable_metrics, "analysis/betaSurvival_uniqueAreaTrapSnare/acceptable_metrics.csv")

table(acceptable_properties$n_methods_used)
table(acceptable_properties$methods_used)

# x = n_methods_used, y = nm_rmse_density [more outliers at low methods used but means are about the same]
# x = n_methods_used, y = rmse_density [more outliers at low methods used but means are about the same]
# x = n_methods_used, y = mbias_density [more outliers at low methods used but means are about the same]
# x = n_methods_used, y = mpe_density [more outliers at low methods used but means are about the same]

# x = methods_used, y = nm_rmse_density [more outliers at low methods used but means are about the same]
# x = methods_used, y = rmse_density [more outliers at low methods used but means are about the same]
# x = methods_used, y = mbias_density [more outliers at low methods used but means are about the same]
# x = methods_used, y = mpe_density [more outliers at low methods used but means are about the same]

property_error |>
  pivot_longer(cols = c(nm_rmse_density, rmse_density, mbias_density, mpe_density),
               names_to = "metric",
               values_to = "value") |>
  ggplot() +
  aes(y = n_methods_used, x = value, group = n_methods_used) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free") +
  coord_cartesian(xlim = c(0, 5)) +
  theme_bw()



