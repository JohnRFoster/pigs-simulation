## workflow for fitting gradient boosting on nbaf hpc

start_time <- Sys.time()

Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = FALSE)
renv::load("/home/john.foster/pigs-simulation/")

library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(recipes)
library(rsample)
library(xgboost)
library(doParallel)
library(foreach)
source("R/functions_analysis.R")

analysis_dir <- "analysis"
model_dir <- "betaSurvival_uniqueAreaTrapSnare"
path <- file.path(analysis_dir, model_dir)

config_name <- "hpc_production"
config <- config::get(config = config_name)
top_dir <- config$top_dir
analysis_dir <- config$analysis_dir
dev_dir <- config$dev_dir
model_dir <- config$model_dir
project_dir <- config$project_dir
path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, model_dir)

data <- read_rds(file.path(path, "abundanceScoresByPrimaryPeriod.rds")) |>
  group_by(property_id) |>
  mutate(delta = c(NA, diff(PPNum))) |>
  ungroup() |>
  # filter(med_density > 0,
  #        var_density > 0) |>
  mutate(mbias_density_class = as.numeric(mbias_density > 0)) |>
  rename(mbias_density_reg = mbias_density)

responses <- c(
  "nm_rmse_density",
  "mpe_density",
  "mbias_density_reg",
  "mbias_density_class",
  "med_density",
  "var_density"
)

eta_grid <- tibble(
  responses = responses,
  eta = 0.1,
  task = 1:length(responses)
)

# hyperparameter grid
hyper_grid <- expand_grid(
  max_depth = 3:8,
  min_child_weight = 0.5,
  subsample = 0.5,
  colsample_bytree = 0.5,
  gamma = c(0, 1, 10, 100, 1000),
  lambda = c(0, 1e-2, 0.1, 1, 100, 1000),
  alpha = c(0, 1e-2, 0.1, 1, 100, 1000),
  rmse = 0,
  trees = 0
)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
if (is.na(task_id)) {
  task_id <- 5
}
message("task id: ", task_id)

task_grid <- eta_grid |>
  filter(task == task_id)

y <- task_grid |>
  pull(responses)

message("\ny: ", y)

eta <- task_grid |>
  pull(eta)
message("eta: ", eta)
# message("Number of threads: ", Sys.getenv("OMP_NUM_THREADS"))

array_grid <- bind_cols(hyper_grid, task_grid)

total_cores <- 96
n_threads <- 2
n_models_per_loop <- total_cores / n_threads
n_loops <- ceiling(nrow(array_grid) / n_models_per_loop)

array_grid <- array_grid |>
  mutate(task = rep(1:n_loops, each = n_models_per_loop)[1:nrow(array_grid)])

df_model <- subset_rename(data, y, 1100)

baked_data <- my_recipe(
  df_model$train |> select(-PPNum, -property, -property_id, -simulation_id),
  df_model$test |> select(-PPNum, -property, -property_id, -simulation_id)
)

train_data <- baked_data$df_train
X <- train_data |>
  select(-y) |>
  as.matrix()
Y <- train_data |> pull(y)

objective <- if_else(
  y == "mbias_density_class",
  "binary:logistic",
  "reg:squarederror"
)
objective <- "reg:tweedie"

fit_xgBoost <- function(i, array_grid, n_threads) {
  set.seed(123)
  out_grid <- array_grid[i, ]
  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 5000,
    objective = objective,
    metrics = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    nthread = n_threads,
    verbose = 0,
    params = list(
      eta = array_grid$eta[i],
      max_depth = array_grid$max_depth[i],
      min_child_weight = array_grid$min_child_weight[i],
      subsample = array_grid$subsample[i],
      colsample_bytree = array_grid$colsample_bytree[i],
      gamma = array_grid$gamma[i],
      lambda = array_grid$lambda[i],
      alpha = array_grid$alpha[i]
    )
  )

  out_grid$rmse <- min(m$evaluation_log$test_rmse_mean)
  out_grid$trees <- m$best_iteration
  out_grid
}

path <- file.path(
  top_dir,
  project_dir,
  analysis_dir,
  dev_dir,
  "gradientBoosting"
)

message("Begin grid search...")
for (j in seq_len(n_loops)) {
  model_time <- Sys.time()

  filename <- file.path(path, paste0(j, "_", y, "_xgbTree_tweedie.rds"))

  if (file.exists(filename)) {
    next
  }

  J <- array_grid |> filter(task == j)

  cl <- makeCluster(nrow(J))
  registerDoParallel(cl)

  out <- foreach::foreach(
    i = 1:nrow(J),
    .combine = rbind,
    .inorder = FALSE,
    .packages = c("xgboost")
  ) %dopar%
    fit_xgBoost(i, J, n_threads) |>
    suppressPackageStartupMessages() |>
    as_tibble()

  stopCluster(cl)

  write_rds(out, filename)

  total_time <- Sys.time() - model_time
  message("\n[", j, "/", n_loops, "]")
  print(round(total_time, 2))
}

message("Grid seach complete!")

out_files <- file.path(path, list.files(path))
all_out <- out_files |>
  purrr::map_dfr(readRDS) |>
  distinct()

out <- all_out |>
  as_tibble() |>
  arrange(rmse)

message("All fits")
print(out)

filename <- file.path(path, paste0("xgbTree_tweedie_", y, ".rds"))
write_rds(out, filename)

message("=== DONE ===\n\n")

total_time <- Sys.time() - start_time
message("Elapsed time: ")
print(total_time)
