library(nimble)
library(parallel)
library(coda)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

path <- "out/hpc/density_1.4"
task <- "20"


rds <- read_rds(file.path(path, task, "simulation_data.rds"))
rds$psrf

burn <- "parameters_burnin.rds"
params <- read_rds(file.path(path, task, burn))
params <- params[[1]]

plot(params[,"psi_phi"])

rds$posterior_samples |>
  pull("psi_phi") |>
  # exp() |>
  hist(breaks = 50)

rds$N$property |> max()

phi_mu <- tibble(
  phi = rds$posterior_samples |> pull("phi_mu"),
  psi = rds$posterior_samples |> pull("psi_phi")
) |>
  mutate(a = phi * psi,
         b = (1-phi) * psi,
         phi_mu = rbeta(n(), a, b))

hist(phi_mu$phi_mu, breaks = 50)
quantile(phi_mu$phi_mu, c(0.025, 0.5, 0.975))


rds$N |>
  mutate(property = as.character(property)) |>
  ggplot() +
  aes(x = PPNum, y = N, color = property) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

str(rds)


n_samps <- rds$posterior_samples |>
  select(contains("xn"))

n_median <- apply(n_samps, 2, median)
n_low <- apply(n_samps, 2, quantile, 0.025)
n_high <- apply(n_samps, 2, quantile, 0.975)

tibble(
  observed = rds$N$N,
  predicted = n_median,
  ymin = n_low,
  ymax = n_high
) |>
  ggplot() +
  aes(x = observed, y = predicted, ymin = ymin, ymax = ymax) +
  geom_point() +
  # geom_linerange() +
  geom_abline(intercept = 0, slope = 1) +
  # coord_cartesian(ylim = c(0, 250)) +
  theme_bw()

N <- abundance$abundance_summaries |>
  arrange(simulation, property, PPNum) |>
  mutate(obs = if_else(obs_flag == 1, "Removal effort", "No removal effort"))

N_density <- N |>
  select(-contains("abundance")) |>
  rename(known = density,
         ymin = low_density,
         y = med_density,
         ymax = high_density)

N_abundance <- N |>
  select(-contains("density")) |>
  rename(known = abundance,
         ymin = low_abundance,
         y = med_abundance,
         ymax = high_abundance)

n_pp <- N |>
  select(simulation, property, PPNum) |>
  distinct() |>
  group_by(simulation, property) |>
  mutate(n = 1:n()) |>
  filter(n == max(n)) |>
  ungroup()


n_pp |>
  filter(n == median(n))

plot_N <- function(df, sim, p, density){
  ylab <- ifelse(density, "Density", "Abundance")
  pe <- "Posterior estimate"
  df |>
    filter(simulation == sim,
           property %in% p) |>
    ggplot() +
    aes(x = PPNum, y = known) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = "95% CI"), alpha = 0.5) +
    geom_line(aes(y = y, linetype = "Median"), linewidth = 1) +
    geom_point(aes(color = obs), size = 3) +
    scale_color_manual(name = "Known abundance", values = c("#d95f02", "#7570b3")) +
    scale_linetype_manual(pe, values = 1) +
    scale_fill_manual(pe, values = "#1b9e77") +
    labs(x = "Primary period",
         y = ylab) +
    facet_wrap(~ property) +
    theme_bw()
}

sim <- 4
p <- c(44, 84)

plot_N(N_density, sim, p, TRUE)
plot_N(N_abundance, sim, p, FALSE)

N |>
  select(simulation, property, PPNum, abundance) |>
  arrange(simulation, property, PPNum) |>
  group_by(simulation, property) |>
  filter(PPNum == min(PPNum) | PPNum == max(PPNum)) |>
  mutate(t = 1:n()) |>
  ungroup() |>
  select(-PPNum) |>
  pivot_wider(names_from = t,
              values_from = abundance) |>
  mutate(delta = `2` - `1`) |>
  filter(delta == min(delta))

sim <- 4
p <- 49
plot_N(N_density, sim, p, TRUE)
plot_N(N_abundance, sim, p, FALSE)

N |>
  select(simulation, property, PPNum) |>
  arrange(simulation, property, PPNum) |>
  group_by(simulation, property) |>
  filter(PPNum == min(PPNum) | PPNum == max(PPNum)) |>
  mutate(t = 1:n()) |>
  ungroup() |>
  pivot_wider(names_from = t,
              values_from = PPNum) |>
  mutate(delta = `2` - `1`) |>
  filter(delta == median(delta))

sim <- 4
p <- 49
plot_N(N_density, sim, p, TRUE)
plot_N(N_abundance, sim, p, FALSE)


effort_data <- read_csv("../pigs-statistical/data/insitu/effort_data.csv")
