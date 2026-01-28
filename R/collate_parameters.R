# --------------------------------------------------------------------
#
# Script for collating parameters from simulation results
#
# John Foster
#
# --------------------------------------------------------------------

message("\n=== PARAMETERS ===\n")

# will get config_name from here
source("R/functions_collate.R")
config <- config::get(config = config_name)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

read_path <- get_path("read", config_name, task_id)

density_tasks <- list.files(read_path)
message("Tasks to collate ", length(density_tasks))

nodes <- "parameters"
tasks_ls <- get_tasks(density_tasks, read_path, nodes)
all_samples <- tasks_ls$all_samples
all_beta_p <- tasks_ls$all_beta_p
all_methods <- tasks_ls$all_methods
all_area <- tasks_ls$all_area
all_theta <- tasks_ls$all_theta
all_p <- tasks_ls$all_p
all_psrf <- tasks_ls$all_psrf
all_land_cover <- tasks_ls$all_land_cover

path <- get_path("write", config_name, task_id)
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_rds(all_methods, file.path(path, "method_parameter_lookup.rds"))
write_rds(all_land_cover, file.path(path, "land_cover_lookup.rds"))

data_model_summaries <- list(
  theta = all_theta,
  potential_area = all_area,
  p = all_p
)

write_rds(data_model_summaries, file.path(path, "data_model_summaries.rds"))
message("\ndata model summaries done\n")

rm(list = c("all_theta", "all_area", "all_p", "data_model_summaries"))
gc()

all_param_samples <- all_samples |>
  select(
    simulation,
    start_density,
    contains("beta1"),
    contains("beta_p"),
    contains("log_gamma["),
    contains("log_rho["),
    contains("p_mu["),
    contains("log_nu"),
    contains("phi_mu"),
    contains("psi_phi")
  )

write_rds(all_param_samples, file.path(path, "all_param_samples.rds"))
write_rds(all_beta_p, file.path(path, "all_beta_p.rds"))

recovery_list <- list()
residual_list <- list()

## capture probability intercepts ------
beta1_long <- all_samples |>
  select_pivot_longer("beta1") |>
  mutate(
    method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
    position = 1
  )

recovery_list$beta1 <- recov_beta1(beta1_long, all_psrf)
residual_list$beta1 <- resid_beta1(beta1_long)

message("\ncapture intercepts done\n")

## capture probability covariates ------
beta_p_long <- all_samples |>
  select_pivot_longer("beta_p") |>
  mutate(
    method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
    position = as.numeric(str_extract(node, "(?<=\\, )\\d")) + 1
  )

recovery_list$beta_p <- recov_beta_p(beta_p_long, all_psrf)
residual_list$beta_p <- resid_beta_p(beta_p_long)

message("\ncapture covariates done\n")

## gamma ------

gH <- all_methods |>
  select(idx, gamma, method, simulation) |>
  filter(method %in% c("Snares", "Traps")) |>
  mutate(idx = idx - 3) |>
  rename(actual = gamma)

gamma_long <- all_samples |>
  select_pivot_longer("log_gamma[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = exp(value))

recovery_list$gamma <- recov_gamma(gamma_long, gH, all_psrf)
residual_list$gamma <- resid_gamma(gamma_long, gH)

message("\nsaturation constant done\n")

## rho ------
rH <- all_methods |>
  select(idx, rho, method, simulation) |>
  rename(actual = rho)

rho_long <- all_samples |>
  select_pivot_longer("log_rho[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = exp(value))

recovery_list$rho <- recov_gamma(rho_long, rH, all_psrf)
residual_list$rho <- resid_gamma(rho_long, rH)

message("\nsearch area done\n")

## unique area ------
pH <- all_methods |>
  filter(idx %in% c(4, 5)) |>
  select(idx, p_unique, method, simulation) |>
  rename(actual = p_unique) |>
  mutate(idx = idx - 3)

p_mu_long <- all_samples |>
  select_pivot_longer("p_mu[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = ilogit(value))

recovery_list$p_mu <- recov_gamma(p_mu_long, pH, all_psrf)
residual_list$p_mu <- resid_gamma(p_mu_long, pH)

message("\nunique area done\n")

## litter size ------
actual <- 5.290323
ls_long <- all_samples |>
  select_pivot_longer("log_nu") |>
  mutate(value = exp(value))

ls_recovery <- ls_long |>
  group_by(simulation, node, start_density) |>
  my_summary() |>
  mutate(actual = actual) |>
  ungroup() |>
  recovered() |>
  left_join(all_psrf)

ls_residual <- ls_long |>
  mutate(actual = actual) |>
  mutate(value = value - actual) |>
  group_by(node, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$litter_size <- ls_recovery
residual_list$litter_size <- ls_residual

message("\nliter size done\n")

## survival ------
actual <- config$phi_mu
phi_long <- all_samples |>
  select_pivot_longer("phi_mu")

recovery_list$phi_mu <- recov_phi(phi_long, actual, all_psrf)
residual_list$phi_mu <- resid_phi(phi_long, actual)

actual <- config$psi_phi
psi_phi_long <- all_samples |>
  select_pivot_longer("psi_phi")

recovery_list$psi_phi <- recov_phi(psi_phi_long, actual, all_psrf)
residual_list$psi_phi <- resid_phi(psi_phi_long, actual)

message("\nsurvival done\n")

write_rds(recovery_list, file.path(path, "parameterRecovery.rds"))
write_rds(residual_list, file.path(path, "parameterResidual.rds"))


message("\n=== PARAMETERS DONE ===\n")
