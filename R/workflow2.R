library(nimble)
library(coda)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

top_dir <- config$top_dir
out_dir <- config$out_dir
dev_dir <- config$dev_dir
model_dir <- config$model_dir
project_dir <- config$project_dir
start_density <- config$start_density
density_dir <- paste0("density_", start_density)

path <- file.path(
  top_dir,
  project_dir,
  out_dir,
  dev_dir,
  model_dir,
  density_dir
)
message("Simulations will be written to\n   ", path)

task_id <- Sys.getpid()
message("  Task ID: ", task_id)

dest <- file.path(path, task_id)
if (!dir.exists(dest)) {
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
}

# need a lookup table for property IDs and how may methods they employ ----
n_method_lookup <- df |>
  select(property, method) |>
  distinct() |>
  count(property)

# get the proportion of properties that use n methods
# for determining the number of properties to simulate
n_rel <- n_method_lookup |>
  count(n, name = "n_sum") |>
  mutate(rel_prop = n_sum / sum(n_sum), n_simulate = ceiling(rel_prop * 100))

# -----------------------------------------------------------------
# 1-method properties ----
# -----------------------------------------------------------------
message("Create 1-method properties")
source("R/one_method_properties.R")
n_pp <- config$n_pp # the number of primary periods to simulate
method_1 <- one_method_properties(df, n_rel$n_simulate[1], n_pp)

# -----------------------------------------------------------------
# n-method properties ----
# -----------------------------------------------------------------

source("R/n_method_properties.R")

message("Create 2-method properties")
method_2 <- n_method_properties(df, n_rel$n_simulate[2], 2, n_pp)
message("Create 3-method properties")
method_3 <- n_method_properties(df, n_rel$n_simulate[3], 3, n_pp)
message("Create 4-method properties")
method_4 <- n_method_properties(df, n_rel$n_simulate[4], 4, n_pp)
message("Create 5-method properties")
method_5 <- n_method_properties(df, n_rel$n_simulate[5], 5, n_pp)

properties <- c(
  method_1,
  method_2,
  method_3,
  method_4,
  method_5
)

n_properties <- length(properties)
assertthat::are_equal(n_properties, sum(n_rel$n_simulate))

# -----------------------------------------------------------------
# Simulate eco/take dynamics ----
# -----------------------------------------------------------------

## the first thing we need to do is assign properties to counties
max_counties <- config$max_counties
weights <- runif(max_counties)
norm_weights <- weights / sum(weights)
properties_per_county <- rmulti(1, n_properties, norm_weights)
properties_per_county <- properties_per_county[properties_per_county > 0]
n_county <- length(properties_per_county)

## randomly assign which county each property goes in
order_county <- tibble(
  order_property = sample.int(n_properties),
  order_county = rep(1:n_county, times = properties_per_county)
) |>
  arrange(order_property) |>
  pull(order_county)

## we need covariate data
land_cover <- matrix(rnorm(n_county * 3), n_county, 3)

## we need the effect of X on detection probability p
beta_p <- matrix(rnorm(20), 5, 4)

# method lookup table (data model parameters)
method_lookup <- tibble(
  idx = 1:5,
  method = c("Firearms", "Fixed wing", "Helicopter", "Snares", "Traps"),
  p_unique = c(runif(1), 0, 0, runif(2)),
  rho = c(
    runif(1, 0.01, 5), # firearms; p_mu[1]
    runif(1, 1, 200), # fixed wing
    runif(1, 1, 200), # helicopter
    runif(1, 0.1, 15), # snare; gamma[1], p_mu[2]
    runif(1, 0.1, 15)
  ), # traps; gamma[2], p_mu[3]
  gamma = c(
    0,
    0,
    0,
    rgamma(1, 7.704547, 4.41925),
    rgamma(1, 3.613148, 3.507449)
  )
)

message("Data model parameters")
method_lookup

phi_mu <- config$phi_mu
psi_phi <- config$psi_phi
data_dir <- config$data_dir

message("Simulate swine dynamics")
source("R/eco_dynamics.R")
all_dynamics <- 1:n_properties |>
  map(
    \(x) {
      simulate_dm(
        properties[[x]],
        x,
        order_county[x],
        phi_mu,
        psi_phi,
        land_cover[order_county[x], ],
        beta_p,
        start_density,
        method_lookup,
        file.path(top_dir, data_dir, "insitu/effort_data.csv")
      )
    }
  ) |>
  compact() |> # some properties will be empty because they go extinct during spin-up
  list_rbind()

N <- all_dynamics |>
  select(PPNum, N, property, county, property_area) |>
  distinct() |>
  mutate(
    property = as.numeric(as.factor(property)),
    density = N / property_area,
    n_id = 1:n()
  )

take <- all_dynamics |>
  filter(!is.na(take)) |>
  mutate(property = as.numeric(as.factor(property)))

known_values <- list(
  N = N,
  take = take,
  phi_mu = phi_mu,
  psi_phi = psi_phi,
  method_lookup = method_lookup,
  beta_p = beta_p,
  land_cover = land_cover
)

write_rds(known_values, file.path(dest, "knownValues.rds"))

# -----------------------------------------------------------------
# Fit MCMC ----
# -----------------------------------------------------------------
message("Fit MCMC")
source("R/prep_nimble.R")
nimble_data <- prep_nimble(N, take, land_cover)
constants <- nimble_data$constants
data <- nimble_data$data

source("R/inits.R")
inits_test <- inits(data, constants)

custom_samplers <- tribble(
  ~node     , ~type      ,
  "log_nu"  , "slice"    ,
  "phi_mu"  , "slice"    ,
  "psi_phi" , "slice"    ,
  "log_rho" , "AF_slice"
)

source("R/model_removal_dm.R")

monitors_add <- c("xn", "p", "log_theta")

n_iter <- config$n_iter
n_chains <- config$n_chains

message(
  "Fitting MCMC with ",
  n_iter,
  " iterations across ",
  n_chains,
  " chains"
)

source("R/fit_mcmc.R")
samples <- fit_mcmc(
  modelCode,
  data,
  constants,
  n_iter,
  n_chains,
  custom_samplers,
  monitors_add
)

# -----------------------------------------------------------------
# Check MCMC ----
# -----------------------------------------------------------------
message("Check MCMC")
params_check <- c(
  "beta_p",
  "beta1",
  "log_gamma",
  "log_rho",
  "phi_mu",
  "psi_phi",
  "log_nu",
  "p_mu"
)

# message("MCMC warnings")
# warnings()

message("Checking MCMC")
source("R/check_mcmc.R")
n_mcmc <- config$n_mcmc
check <- check_mcmc(samples, params_check, n_mcmc, dest)

samples_draw <- check$posterior_samples

source("R/functions_predict.R")
post <- data_posteriors(samples_draw, constants, data)

# -----------------------------------------------------------------
# Write to disk ----
# -----------------------------------------------------------------

out_list <- list(
  task_id = task_id,
  posterior_samples = as_tibble(samples_draw),
  posterior_take = post$y,
  posterior_p = post$p,
  posterior_potential_area = post$potential_area,
  take = take,
  N = N,
  method_lookup = method_lookup,
  beta_p = beta_p,
  constants = constants,
  start_density = start_density,
  data = data,
  psrf = check$psrf,
  effective_samples = check$effective_samples,
  burnin = check$burnin,
  converged = check$converged,
  bad_mcmc = check$bad_mcmc
)

message("Write to disk")
write_rds(
  out_list,
  file.path(dest, "simulation_data.rds")
)
return(out_list)
