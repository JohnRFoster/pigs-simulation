library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(recipes)
library(rsample)
library(xgboost)
source("R/functions_analysis.R")

ylab <- "var_density"
message("Response: ", ylab)
analysis_dir <- "analysis"
model_dir <- "MLs"
path <- file.path(analysis_dir, model_dir)

best_model <- read_rds(file.path(path, paste0("xgbTree_", ylab, ".rds"))) |>
  group_by(responses) |>
  filter(rmse == min(rmse)) |>
  filter(min_child_weight == min(min_child_weight)) |>
  ungroup()

model_dir <- "betaSurvival_uniqueAreaTrapSnare"
path <- file.path(analysis_dir, model_dir)
# path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, config$model_dir)
data <- read_rds(file.path(path, "abundanceScoresByPrimaryPeriod.rds")) |>
  group_by(property_id) |>
  mutate(delta = c(NA, diff(PPNum))) |>
  ungroup() |>
  filter(med_density > 0) |>
  filter(var_density > 0)

df_model <- subset_rename(data, ylab, 1100)

split_train <- df_model$train |>
  select(-PPNum, -property, -property_id, -simulation_id)
split_test <- df_model$test |>
  select(-PPNum, -property, -property_id, -simulation_id)

# split_test <- data |>
#   rename(y = all_of(ylab),
#          sum_take_d = sum_take_density) |>
#   group_by(property_id) |>
#   mutate(delta = c(NA, diff(PPNum))) |>
#   ungroup() |>
#   select(all_of(colnames(split_train)))

# predict to all data not in training set!

train_best_pred <- function(data_train, data_test, y, best_model) {
  baked_data <- my_recipe(data_train, data_test)
  train <- baked_data$df_train
  test <- baked_data$df_test

  tune_grid <- best_model |>
    filter(responses == y)

  params <- list(
    eta = pull(tune_grid, eta),
    max_depth = pull(tune_grid, max_depth),
    min_child_weight = pull(tune_grid, min_child_weight),
    subsample = pull(tune_grid, subsample),
    colsample_bytree = pull(tune_grid, colsample_bytree),
    gamma = pull(tune_grid, gamma),
    lambda = pull(tune_grid, lambda),
    alpha = pull(tune_grid, alpha)
  )

  X <- train |>
    select(-y) |>
    as.matrix()
  Y <- train |> pull(y)

  objective <- if_else(
    y == "mbias_density_class",
    "binary:logistic",
    "reg:squarederror"
  )
  # objective <- "reg:tweedie"

  # train final model
  fit <- xgboost(
    params = params,
    data = X,
    label = Y,
    nrounds = pull(tune_grid, trees),
    objective = objective,
    verbose = 0
  )

  newdata <- test |>
    select(-y) |>
    as.matrix()
  pred <- predict(fit, newdata)

  test <- test |>
    mutate(pred = pred)

  rmse <- sqrt(mean((test$y - pred)^2))
  r2 <- cor(test$y, pred)^2
  vi <- xgb.importance(model = fit)
  vi$gainRelative <- vi$Gain / max(vi$Gain)

  return(
    list(
      fit = fit,
      train = X,
      Y = Y,
      test = test,
      rmse = rmse,
      r2 = r2,
      vi = vi
    )
  )
}

partial_dependence <- function(col, train_dat, fit) {
  df <- pdp::partial(
    fit,
    col,
    train = train_dat,
    plot = FALSE,
    parallel = FALSE
  )
  as_tibble(df)
}

best_pred <- train_best_pred(split_train, split_test, ylab, best_model)
best_pred$vi

df_test <- best_pred$test
df_train <- best_pred$train
fit <- best_pred$fit

message("Single dependence...")
col <- best_pred$vi$Feature[best_pred$vi$gainRelative >= 0.1]
single_dependence <- list()
pb <- txtProgressBar(max = length(col), style = 3)
for (i in seq_along(col)) {
  pdp_s <- partial_dependence(col[i], df_train, fit)
  single_dependence[[col[i]]] <- pdp_s
  setTxtProgressBar(pb, i)
}
close(pb)

out <- list(
  pred = best_pred,
  single_dependence = single_dependence
)

message("Write rds")
# path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, model_dir)
analysis_dir <- "analysis"
model_dir <- "MLs"
path <- file.path(analysis_dir, model_dir)
dest <- file.path(path, paste0(ylab, "_xgBoostAnalysis.rds"))
write_rds(out, dest)

message("Done!")
