library(stringr)
library(dplyr)
library(tidyr)
library(nimble)
library(readr)
library(coda)

config_name <- "hpc_production"

get_path <- function(type, config_name, task_id) {
  config <- config::get(config = config_name)
  top_dir <- config$top_dir
  out_dir <- config$out_dir
  analysis_dir <- config$analysis_dir
  dev_dir <- config$dev_dir
  model_dir <- config$model_dir
  project_dir <- config$project_dir

  start_vec <-
    c("0.3", "1.475", "2.65", "3.825", "5")
  start_density <- start_vec[task_id]
  density_dir <- paste0("density_", start_density)

  if (type == "read") {
    path <- file.path(
      top_dir,
      project_dir,
      out_dir,
      dev_dir,
      model_dir,
      density_dir
    )
  }
  if (type == "write") {
    path <- file.path(
      top_dir,
      project_dir,
      analysis_dir,
      dev_dir,
      model_dir,
      density_dir
    )
  }
  return(path)
}

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])


add_ids <- function(df, task_id, s_density) {
  df |>
    mutate(simulation = task_id, start_density = s_density)
}

bind_samples <- function(all_bind, ls, t_id, dens) {
  samples <- ls$posterior_samples |>
    add_ids(t_id, dens)
  bind_rows(all_bind, samples)
}

get_y <- function(ls, t_id, dens) {
  y_pred <- ls$posterior_take
  colnames(y_pred) <- 1:ncol(y_pred)
  y_pred |>
    as_tibble() |>
    pivot_longer(
      cols = everything(),
      names_to = "p_id",
      values_to = "value"
    ) |>
    add_ids(t_id, dens)
}

take_summary <- function(y_long, t_id, dens) {
  y_long |>
    group_by(p_id) |>
    my_summary() |>
    ungroup() |>
    add_ids(t_id, dens) |>
    mutate(p_id = as.numeric(p_id))
}

take_error_by_observation <- function(y_long, rds_take, t_id, dens) {
  y_long |>
    group_by(p_id) |>
    take_calc() |>
    ungroup() |>
    select(-nm_rmse, -sd_ratio) |>
    left_join(rds_take) |>
    add_ids(t_id, dens)
}

take_effort_summary <- function(ebo) {
  ebo |>
    group_by(start_density, simulation, property, PPNum, method) |>
    summarise(
      low_effort = quantile(effort_per, 0.025),
      q1_effort = quantile(effort_per, 0.25),
      median_effort = quantile(effort_per, 0.5),
      mean_effort = mean(effort_per),
      sum_effort = sum(effort_per),
      q3_effort = quantile(effort_per, 0.75),
      high_effort = quantile(effort_per, 0.975),
      low_trap_count = quantile(trap_count, 0.025),
      q1_trap_count = quantile(trap_count, 0.25),
      median_trap_count = quantile(trap_count, 0.5),
      mean_trap_count = mean(trap_count),
      sum_trap_count = sum(trap_count),
      q3_trap_count = quantile(trap_count, 0.75),
      high_trap_count = quantile(trap_count, 0.975),
      n_reps = n()
    )
}

take_error_by_simulation <- function(y_long) {
  y_long |>
    group_by(start_density, simulation) |>
    take_calc() |>
    ungroup()
}

take_error_by_simulation_method <- function(y_long) {
  y_long |>
    group_by(simulation, start_density, method) |>
    take_calc() |>
    ungroup()
}

take_error_by_property <- function(y_long) {
  y_long |>
    group_by(simulation, start_density, property) |>
    take_calc() |>
    ungroup()
}

get_take <- function(ls, t_id, dens) {
  ls$take |>
    add_ids(t_id, dens) |>
    mutate(p_id = 1:n())
}

bind_post_summaries <- function(all_bind, node, ls, t_id, dens) {
  df <- ls[[node]]
  tb <- as_tibble(df) |>
    add_ids(t_id, dens) |>
    mutate(p_id = 1:n(), parameter = node)
  bind_rows(all_bind, tb)
}

bind_N <- function(all_bind, ls, t_id, dens) {
  obs_flag <- ls$take |>
    select(property, county, PPNum) |>
    distinct() |>
    mutate(obs_flag = 1)

  N <- ls$N |>
    add_ids(t_id, dens) |>
    left_join(obs_flag) |>
    mutate(obs_flag = if_else(is.na(obs_flag), 0, obs_flag))
  bind_rows(all_bind, N)
}

bind_beta_p <- function(all_bind, ls, t_id, dens) {
  # need a lookup table for known data model covariates
  bH <- tibble(
    method_idx = rep(1:nrow(ls$beta_p), ncol(ls$beta_p)),
    position = rep(1:ncol(ls$beta_p), each = nrow(ls$beta_p)),
    actual = as.numeric(ls$beta_p)
  ) |>
    add_ids(t_id, dens)
  bind_rows(all_bind, bH)
}

bind_methods <- function(all_bind, ls, t_id, dens) {
  method_lookup <- ls$method_lookup |>
    add_ids(t_id, dens)
  bind_rows(all_bind, method_lookup)
}

bind_land_cover <- function(all_bind, ls, t_id, dens) {
  land_cover <- ls$data$X_p
  colnames(land_cover) <- c("c_road_den", "c_rugged", "c_canopy")
  land_cover <- land_cover |>
    as_tibble() |>
    mutate(county = 1:n()) |>
    add_ids(t_id, dens)
  bind_rows(all_bind, land_cover)
}

bind_psrf <- function(all_bind, ls, t_id, dens) {
  names <- rownames(ls$psrf)
  psrf <- ls$psrf |>
    as_tibble() |>
    mutate(node = names) |>
    add_ids(t_id, dens)
  bind_rows(all_bind, psrf)
}

select_pivot_longer <- function(df, node) {
  df |>
    select(contains(node), simulation, start_density) |>
    pivot_longer(cols = -c(simulation, start_density), names_to = "node")
}

recovered <- function(df) {
  df |>
    mutate(parameter_recovered = if_else(actual >= low & actual <= high, 1, 0))
}

my_summary <- function(df) {
  df |>
    summarise(
      low = quantile(value, 0.05),
      q1 = quantile(value, 0.25),
      med = quantile(value, 0.5),
      q3 = quantile(value, 0.75),
      high = quantile(value, 0.95),
      mu = mean(value),
      sd = sd(value)
    )
}

recov_beta1 <- function(df, psrf) {
  df |>
    group_by(simulation, node, method_idx, position, start_density) |>
    my_summary() |>
    left_join(all_beta_p) |>
    ungroup() |>
    recovered() |>
    left_join(psrf)
}

resid_beta1 <- function(df) {
  df |>
    left_join(all_beta_p) |>
    mutate(value = value - actual) |>
    group_by(node, position, method_idx, start_density) |>
    my_summary() |>
    ungroup()
}

recov_beta_p <- function(df, psrf) {
  df |>
    group_by(simulation, node, method_idx, position, start_density) |>
    my_summary() |>
    left_join(all_beta_p) |>
    ungroup() |>
    recovered() |>
    left_join(psrf)
}

resid_beta_p <- function(df) {
  df |>
    left_join(all_beta_p) |>
    mutate(value = value - actual) |>
    group_by(node, position, method_idx, start_density) |>
    my_summary() |>
    ungroup()
}

recov_gamma <- function(df, H, psrf) {
  df |>
    group_by(simulation, node, idx, start_density) |>
    my_summary() |>
    left_join(H) |>
    ungroup() |>
    recovered() |>
    left_join(psrf)
}

resid_gamma <- function(df, H) {
  df |>
    left_join(H) |>
    mutate(value = value - actual) |>
    group_by(node, idx, start_density) |>
    my_summary() |>
    ungroup()
}

recov_phi <- function(df, known, psrf) {
  df |>
    group_by(simulation, node, start_density) |>
    my_summary() |>
    mutate(actual = known) |>
    ungroup() |>
    recovered() |>
    left_join(psrf)
}

resid_phi <- function(df, known) {
  phi_residual <- df |>
    mutate(actual = known) |>
    mutate(value = value - actual) |>
    group_by(node, start_density) |>
    my_summary() |>
    ungroup()
}

get_xn <- function(df, H) {
  df |>
    select_pivot_longer("N[") |>
    filter(!is.na(value)) |>
    mutate(n_id = as.numeric(str_extract(node, "(?<=\\[)\\d*"))) |>
    left_join(H) |>
    filter(!is.na(abundance)) |>
    mutate(estimated_density = value / property_area)
}

get_post_take <- function(df, H) {
  df |>
    pivot_longer(cols = -c(simulation, start_density), names_to = "p_id") |>
    filter(!is.na(value)) |>
    mutate(p_id = as.numeric(p_id)) |>
    left_join(H)
}

take_calc <- function(df) {
  df |>
    summarise(
      mae = mean(abs(value - take)),
      mpe = mean(abs((value + 1) - (take + 1)) / (take + 1)) * 100,
      mbias = mean(value - take),
      rmse = sqrt(mean((value - take)^2)),
      rmsle = sqrt(mean((log(value + 1) - log(take + 1))^2)),
      nm_rmse = rmse / mean(take),
      sd_ratio = sd(value) / sd(take)
    )
}

get_tasks <- function(density_tasks, path, nodes) {
  if (nodes == "parameters") {
    all_samples <- tibble()
    all_beta_p <- tibble()
    all_methods <- tibble()
    all_area <- tibble()
    all_theta <- tibble()
    all_p <- tibble()
    all_psrf <- tibble()
    all_land_cover <- tibble()
  }

  if (nodes == "abundance") {
    all_samples <- tibble()
    all_N <- tibble()
  }

  if (nodes == "take") {
    all_y <- tibble()
    all_take <- tibble()
    all_by_observation <- tibble()
    all_effort_summary <- tibble()
    all_by_simulation <- tibble()
    all_by_simulation_method <- tibble()
    all_by_property <- tibble()
  }

  message("Loop through tasks...")

  pb <- txtProgressBar(max = length(density_tasks), style = 1)
  for (i in seq_along(density_tasks)) {
    task_id <- density_tasks[i]

    rds_file <- file.path(path, task_id, "simulation_data.rds")

    if (file.exists(rds_file)) {
      rds <- read_rds(rds_file)
    } else {
      next
    }

    psrf <- rds$psrf |>
      as_tibble() |>
      mutate(node_names = rownames(rds$psrf)) |>
      filter(node_names != "psi_phi")

    bad_mcmc <- rds$bad_mcmc | any(psrf$`Point est.` > 1.1)

    if (bad_mcmc) {
      next
    }

    start_density <- rds$start_density

    if (nodes == "parameters") {
      all_samples <- bind_samples(all_samples, rds, task_id, start_density)
      all_beta_p <- bind_beta_p(all_beta_p, rds, task_id, start_density)
      all_methods <- bind_methods(all_methods, rds, task_id, start_density)
      all_area <- bind_post_summaries(
        all_area,
        "posterior_potential_area",
        rds,
        task_id,
        start_density
      )
      all_theta <- bind_post_summaries(
        all_theta,
        "posterior_theta",
        rds,
        task_id,
        start_density
      )
      all_p <- bind_post_summaries(
        all_p,
        "posterior_p",
        rds,
        task_id,
        start_density
      )
      all_psrf <- bind_psrf(all_psrf, rds, task_id, start_density)
      all_land_cover <- bind_land_cover(
        all_land_cover,
        rds,
        task_id,
        start_density
      )
    }

    if (nodes == "abundance") {
      all_samples <- bind_samples(all_samples, rds, task_id, start_density)
      all_N <- bind_N(all_N, rds, task_id, start_density)
    }

    if (nodes == "take") {
      y_long <- get_y(rds, task_id, start_density)
      rds_take <- get_take(rds, task_id, start_density)

      y_summary <- take_summary(y_long, task_id, start_density)
      yy <- y_summary |> left_join(rds_take)
      all_y <- bind_rows(all_y, yy)
      all_take <- bind_rows(all_take, rds_take)

      y_long_take <- y_long |>
        mutate(p_id = as.integer(p_id)) |>
        left_join(rds_take)
      ebo <- take_error_by_observation(
        y_long_take,
        rds_take,
        task_id,
        start_density
      )
      all_by_observation <- bind_rows(all_by_observation, ebo)

      tes <- take_effort_summary(ebo)
      all_effort_summary <- bind_rows(all_effort_summary, tes)

      by_simulation <- take_error_by_simulation(y_long_take)
      all_by_simulation <- bind_rows(all_by_simulation, by_simulation)

      by_simulation_method <- take_error_by_simulation_method(y_long_take)
      all_by_simulation_method <- bind_rows(
        all_by_simulation_method,
        by_simulation_method
      )

      by_property <- take_error_by_property(y_long_take)
      all_by_property <- bind_rows(all_by_property, by_property)
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  ls <- list()
  if (nodes == "parameters") {
    ls$all_samples <- all_samples
    ls$all_beta_p <- all_beta_p
    ls$all_methods <- all_methods
    ls$all_area <- all_area
    ls$all_theta <- all_theta
    ls$all_p <- all_p
    ls$all_psrf <- all_psrf
    ls$all_land_cover <- all_land_cover
  }

  if (nodes == "abundance") {
    ls$all_samples <- all_samples
    ls$all_N <- all_N
  }

  if (nodes == "take") {
    ls$all_y <- all_y
    ls$all_take <- all_take
    ls$all_by_observation <- all_by_observation
    ls$all_effort_summary <- all_effort_summary
    ls$all_by_simulation <- all_by_simulation
    ls$all_by_simulation_method <- all_by_simulation_method
    ls$all_by_property <- all_by_property
  }

  return(ls)
}
