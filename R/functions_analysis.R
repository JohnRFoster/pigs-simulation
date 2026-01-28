my_recipe <- function(df_train, df_test) {
  require(rsample)
  require(recipes)

  blueprint <- recipe(y ~ ., data = df_train) |>
    step_dummy(all_nominal_predictors()) |>
    step_nzv(all_predictors()) |>
    step_center(all_numeric_predictors()) |>
    step_scale(all_numeric_predictors())

  prepare <- prep(blueprint, training = df_train)

  baked_train <- bake(prepare, new_data = df_train)
  baked_test <- bake(prepare, new_data = df_test)

  return(list(df_train = baked_train, df_test = baked_test))
}

subset_rename <- function(df, y, n_sample = 500) {
  require(dplyr)
  require(rsample)
  set.seed(5)

  dat <- df |>
    ungroup() |>
    mutate(
      simulation_id = stringr::str_extract(
        property_id,
        "[[:graph:]]*(?=-[[:digit:]]*$)"
      ),
      simulation_id = as.numeric(as.factor(simulation_id)),
      methods_used = as.factor(methods_used)
    ) |>
    rename(y = all_of(y), sum_take_d = sum_take_density) |>
    select(
      -contains("density"),
      -contains("abundance"),
      -extinct,
      -recovered,
      -obs_flag
    )

  if (
    y %in% c("nm_rmse_density", "mpe_density", "med_density", "var_density")
  ) {
    dat <- dat |> mutate(y = log(y))
  }

  # create strata by decile
  # each property will belong to a decile of each land cover variable
  df_strata <- dat |>
    mutate(
      simulation_strata = make_strata(simulation_id, breaks = 10),
      canopy_strata = make_strata(c_canopy, breaks = 10),
      rugged_strata = make_strata(c_rugged, breaks = 10),
      road_den_strata = make_strata(c_road_den, breaks = 10)
    ) |>
    select(property_id, contains("strata")) |>
    distinct()

  col_sample <- function(dfs, col) {
    min_sample <- dfs |>
      pull(all_of(col)) |>
      table() |>
      min()

    dfs |>
      group_by(.data[[col]]) |>
      slice_sample(n = min(min_sample, n_sample)) |>
      ungroup() |>
      pull(property_id)
  }

  simulation_sample <- col_sample(df_strata, "simulation_strata")
  canopy_sample <- col_sample(df_strata, "canopy_strata")
  rugged_sample <- col_sample(df_strata, "rugged_strata")
  road_den_sample <- col_sample(df_strata, "road_den_strata")

  props <- c(simulation_sample, canopy_sample, rugged_sample, road_den_sample)
  df_sample <- df_strata |>
    filter(property_id %in% unique(props))

  message("Number of properties in each simulation strata:")
  print(table(df_sample$simulation_strata))
  message("Number of properties in each canopy strata:")
  print(table(df_sample$canopy_strata))
  message("Number of properties in each rugged strata:")
  print(table(df_sample$rugged_strata))
  message("Number of properties in each road density strata:")
  print(table(df_sample$road_den_strata))

  train <- dat |>
    filter(property_id %in% props)

  test <- dat |>
    filter(!property_id %in% props)

  test_per <- round(nrow(test) / nrow(dat), 2)
  train_per <- round(nrow(train) / nrow(dat), 2)
  message(test_per, " / ", train_per, " [test / train]")

  list(train = train, test = test)
}
