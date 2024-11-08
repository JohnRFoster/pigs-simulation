# --------------------------------------------------------------------
#
# Script for collating predicted values
#
# John Foster
#
# --------------------------------------------------------------------

message("\n=== PARAMETERS ===\n")

# will get config_name from here
source("R/functions_collate.R")
source("R/functions_predict.R")
config <- config::get(config = config_name)

args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])

read_path <- get_path("read", config_name, array_id)

density_tasks <- list.files(read_path)
message("Tasks to collate ", length(density_tasks))

all_take <- tibble()
all_p <- tibble()
all_area <- tibble()
all_theta <- tibble()

pb <- txtProgressBar(max = length(density_tasks), style = 1)
for(i in seq_along(density_tasks)){
# for(i in 1:1){

  task_id <- density_tasks[i]

  rds_file <- file.path(read_path, task_id, "simulation_data.rds")

  if(file.exists(rds_file)){
    rds <- read_rds(rds_file)
  } else {
    next
  }


  psrf <- rds$psrf |>
    as_tibble() |>
    mutate(node_names = rownames(rds$psrf)) |>
    filter(node_names != "psi_phi")

  bad_mcmc <- rds$bad_mcmc | any(psrf$`Upper C.I.` > 1.1)

  if(bad_mcmc) next

  start_density <- rds$start_density

  samples <- rds$posterior_samples
  constants <- rds$constants
  constants$samples <- as.matrix(samples)
  data <- rds$data

  message("NIMBLE")
  print(str(constants))
  print(str(data))

  ls <- data_posteriors(samples, constants, data)

  N <- rds$N
  take <- rds$take


  message("Take")

  print(str(ls))
  print(str(take))


  extract_ls <- function(ls, dfy, dft){

    df <- purrr::pluck(ls, dfy)

    id <- paste0("c", 1:ncol(df))

    y <- as.matrix(df)
    colnames(y) <- id

    take <- dft |>
      mutate(id = id)

    y |>
      as_tibble() |>
      mutate(iter = 1:nrow(y)) |>
      pivot_longer(cols = -iter,
                   names_to = "id") |>
      left_join(take)

  }

  y_iter <- extract_ls(ls, "y_pred", take)
  p_iter <- extract_ls(ls, "p_pred", take)
  area_iter <- extract_ls(ls, "potential_area_pred", take)
  theta_iter <- extract_ls(ls, "theta_pred", take)

  take_summaries <- y_iter |>
    group_by(PPNum, N, take, method, effort_per, trap_count, theta, p,
             potential_area, property, county, property_area) |>
    summarise(low_take = quantile(value, 0.05),
              med_take = quantile(value, 0.5),
              high_take = quantile(value, 0.95),
              var_take = var(value),
              mbias_take = mean(value - take),
              rmse_take = sqrt(mean((value - take)^2)),
              rmsle_take = sqrt(mean((log(value + 1) - log(take + 1))^2)))

  p_summaries <- p_iter |>
    group_by(PPNum, N, take, method, effort_per, trap_count, theta, p,
             potential_area, property, county, property_area) |>
    summarise(low_p = quantile(value, 0.05),
              med_p = quantile(value, 0.5),
              high_p = quantile(value, 0.95),
              var_p = var(value),
              mbias_p = mean(value - p),
              rmse_p = sqrt(mean((value - p)^2)),
              norm_rmse_p = rmse_p / mean(p),
              rmsle_p = sqrt(mean((log(value + 1) - log(p + 1))^2)))

  area_summaries <- area_iter |>
    group_by(PPNum, N, take, method, effort_per, trap_count, theta, p,
             potential_area, property, county, property_area) |>
    summarise(low_area = quantile(value, 0.05),
              med_area = quantile(value, 0.5),
              high_area = quantile(value, 0.95),
              var_area = var(value),
              mbias_area = mean(value - potential_area),
              rmse_area = sqrt(mean((value - potential_area)^2)),
              norm_rmse_area = rmse_area / mean(potential_area),
              rmsle_area = sqrt(mean((log(value + 1) - log(potential_area + 1))^2)))

  theta_summaries <- theta_iter |>
    group_by(PPNum, N, take, method, effort_per, trap_count, theta, p,
             potential_area, property, county, property_area) |>
    summarise(low_theta = quantile(value, 0.05),
              med_theta = quantile(value, 0.5),
              high_theta = quantile(value, 0.95),
              var_theta = var(value),
              mbias_theta = mean(value - theta),
              rmse_theta = sqrt(mean((value - theta)^2)),
              norm_rmse_theta = rmse_theta / mean(theta),
              rmsle_theta = sqrt(mean((log(value + 1) - log(theta + 1))^2)))

  take_summaries <- take_summaries |> add_ids(task_id, start_density)
  p_summaries <- p_summaries |> add_ids(task_id, start_density)
  area_summaries <- area_summaries |> add_ids(task_id, start_density)
  theta_summaries <- theta_summaries |> add_ids(task_id, start_density)

  all_take <- bind_rows(all_take, take_summaries)
  all_p <- bind_rows(all_p, p_summaries)
  all_area <- bind_rows(all_area, area_summaries)
  all_theta <- bind_rows(all_theta, theta_summaries)


  setTxtProgressBar(pb, i)
}
close(pb)

path <- get_path("write", config_name, array_id)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

write_rds(all_take, file.path(path, "predicted_take_summaries.rds"))
write_rds(all_p, file.path(path, "predicted_p_summaries.rds"))
write_rds(all_area, file.path(path, "predicted_area_summaries.rds"))
write_rds(all_theta, file.path(path, "predicted_theta_summaries.rds"))







