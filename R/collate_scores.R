library(readr)
library(dplyr)
library(tidyr)
library(purrr)

config <- config::get(config = "hpc_production")
top_dir <- config$top_dir
out_dir <- config$out_dir
analysis_dir <- config$analysis_dir
dev_dir <- config$dev_dir
model_dir <- config$model_dir
project_dir <- config$project_dir
path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, model_dir)
density_dirs <- list.files(path)
density_dirs <- grep("density_", density_dirs, value = TRUE)

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

message("Get abundance summaries")
abundance_file <- "abundance_summaries.rds"
f_name <- abundance_file
df <- map_files2(density_dirs, f_name)

# create property ID for easier joining
property_ids <- df |>
  select(start_density, simulation, property) |>
  distinct() |>
  mutate(property_id = paste(start_density, simulation, property, sep = "-"))

message("\nSample size by start density:")
sample_size <- df |>
  select(start_density, simulation) |>
  distinct() |>
  group_by(start_density) |>
  count() |>
  rename(`Start density` = start_density, `n simulations` = n)
print(sample_size)

abundance_file <- "abundance_error_by_observation.rds"
f_name <- abundance_file
abundance_error_by_observation <- map_files2(density_dirs, f_name) |>
  left_join(property_ids) |>
  mutate(
    nm_rmse_abundance = if_else(rmse_abundance == 0, 0, nm_rmse_abundance),
    nm_rmse_density = if_else(rmse_density == 0, 0, nm_rmse_density)
  )

abundance_summary <- map_files2(density_dirs, "abundance_summaries.rds") |>
  left_join(property_ids)

f_name <- "land_cover_lookup.rds"
df <- map_files2(density_dirs, f_name)

a_join <- left_join(abundance_summary, df)
b_join <- left_join(a_join, abundance_error_by_observation)

density <- b_join |>
  select(
    start_density,
    PPNum,
    contains("property"),
    contains("abundance"),
    contains("density"),
    contains("c_"),
    obs_flag
  ) |>
  mutate(
    recovered = if_else(
      abundance >= low_abundance & density <= high_abundance,
      "Recovered",
      "Not Recovered"
    ),
    extinct = if_else(abundance == 0, "Extinct", "Extant")
  )

message("\nProportion recovered by extinct vs extant:")
density |>
  group_by(extinct) |>
  count(recovered) |>
  mutate(proportion = round(n / sum(n), 3))


f_name <- "all_take.rds"
all_take <- map_files2(density_dirs, f_name) |>
  left_join(property_ids) |>
  select(-theta, -p, -property, -obs_flag, -simulation, -p_id) |>
  mutate(effort = effort_per * trap_count)

sum_take <- all_take |>
  select(property_id, PPNum, sum_take) |>
  distinct()

density_obs <- density |>
  # filter(obs_flag == 1) |>
  left_join(sum_take) |>
  mutate(sum_take_density = sum_take / property_area) |>
  group_by(property_id) |>
  mutate(delta = c(0, diff(PPNum)))

n_methods_pp <- all_take |>
  select(property_id, PPNum, method) |>
  distinct() |>
  group_by(property_id, PPNum) |>
  count() |>
  rename(n_methods_used = n)

n_return <- all_take |>
  select(property_id, PPNum, method) |>
  distinct() |>
  pivot_wider(names_from = method, values_from = method) |>
  unite(method, -c(property_id, PPNum), sep = ", ", na.rm = TRUE) |>
  group_by(property_id, method) |>
  mutate(return_interval = c(0, diff(PPNum))) |>
  ungroup() |>
  rename(methods_used = method) |>
  left_join(n_methods_pp)

take_joint_return <- left_join(all_take, n_return)

take_effort_summary <- take_joint_return |>
  group_by(property_id, PPNum, methods_used) |>
  summarise(
    mean_effort_per_unit = mean(effort_per),
    sum_effort_per_unit = sum(effort_per),
    mean_effort = mean(effort),
    sum_effort = sum(effort),
    mean_unit_count = mean(trap_count),
    sum_unit_count = sum(trap_count),
    n_reps_pp = n()
  ) |>
  ungroup()

data_final_join <- left_join(density_obs, take_effort_summary)

path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, model_dir)
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_rds(
  data_final_join,
  file.path(path, "abundanceScoresByPrimaryPeriod.rds")
)
