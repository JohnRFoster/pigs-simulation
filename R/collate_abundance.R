# --------------------------------------------------------------------
#
# Script for collating abundance from simulation results
#
# John Foster
#
# --------------------------------------------------------------------

message("\n=== ABUNDANCE ===\n")

# will get config_name from here
source("R/functions_collate.R")
config <- config::get(config = config_name)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

read_path <- get_path("read", config_name, task_id)

density_tasks <- list.files(read_path)
message("Tasks to collate ", length(density_tasks))

nodes <- "abundance"
tasks_ls <- get_tasks(density_tasks, read_path, nodes)
all_samples <- tasks_ls$all_samples
all_N <- tasks_ls$all_N

path <- get_path("write", config_name, task_id)
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

## abundance ------
abundance <- all_N |>
  rename(abundance = N)

rm(all_N)
gc()

xn <- all_samples |>
  get_xn(abundance)

vals <- c("abundance", "density", "value", "estimated_density")

n_attributes <- xn |>
  select(-value, -estimated_density) |>
  distinct()

abundance_summaries <- xn |>
  select(simulation, property, PPNum, all_of(vals)) |>
  group_by(simulation, property, PPNum) |>
  summarise(
    low_abundance = quantile(value, 0.05),
    med_abundance = quantile(value, 0.5),
    high_abundance = quantile(value, 0.95),
    var_abundance = var(value),
    cv_abundance = sd(value) / mean(value),
    low_density = quantile(estimated_density, 0.05),
    med_density = quantile(estimated_density, 0.5),
    high_density = quantile(estimated_density, 0.95),
    var_density = var(estimated_density),
    cv_density = sd(estimated_density) / mean(estimated_density)
  ) |>
  ungroup() |>
  left_join(n_attributes)

write_rds(abundance_summaries, file.path(path, "abundance_summaries.rds"))
rm(abundance_summaries)
gc()
message("\nposterior abundance done\n")

error_by_observation <- xn |>
  select(simulation, property, PPNum, all_of(vals)) |>
  group_by(simulation, property, PPNum) |>
  summarise(
    mpe_abundance = mean(abs((value + 1) - (abundance + 1)) / (abundance + 1)) *
      100,
    mpe_density = mean(
      abs((estimated_density + 0.1) - (density + 0.1)) / (density + 0.1)
    ) *
      100,
    mae_abundance = mean(abs(value - abundance)),
    mae_density = mean(abs(estimated_density - density)),
    mbias_abundance = mean(value - abundance),
    mbias_density = mean(estimated_density - density),
    norm_bias_abundance = mean((value - abundance) / abundance),
    norm_bias_density = mean((estimated_density - density) / density),
    rmse_abundance = sqrt(mean((value - abundance)^2)),
    rmse_density = sqrt(mean((estimated_density - density)^2)),
    rmsle_abundance = sqrt(mean((log(value + 1) - log(abundance + 1))^2)),
    rmsle_density = sqrt(mean(
      (log(estimated_density + 1) - log(density + 1))^2
    ))
  ) |>
  ungroup() |>
  arrange(simulation, property, PPNum) |>
  left_join(n_attributes) |>
  mutate(
    nm_rmse_abundance = if_else(
      abundance == 0,
      rmse_abundance,
      rmse_abundance / abundance
    ),
    nm_rmse_density = if_else(
      density == 0,
      rmse_density,
      rmse_density / density
    )
  )

write_rds(
  error_by_observation,
  file.path(path, "abundance_error_by_observation.rds")
)
rm(error_by_observation)
gc()
message("\nabundance error by observation done\n")

error_by_property <- xn |>
  select(simulation, property, PPNum, all_of(vals)) |>
  group_by(simulation, property) |>
  summarise(
    mpe_abundance = mean(abs((value + 1) - (abundance + 1)) / (abundance + 1)) *
      100,
    mpe_density = mean(
      abs((estimated_density + 0.1) - (density + 0.1)) / (density + 0.1)
    ) *
      100,
    mae_abundance = mean(abs(value - abundance)),
    mae_density = mean(abs(estimated_density - density)),
    mbias_abundance = mean(value - abundance),
    mbias_density = mean(estimated_density - density),
    rmse_abundance = sqrt(mean((value - abundance)^2)),
    rmse_density = sqrt(mean((estimated_density - density)^2)),
    rmsle_abundance = sqrt(mean((log(value + 1) - log(abundance + 1))^2)),
    rmsle_density = sqrt(mean(
      (log(estimated_density + 1) - log(density + 1))^2
    )),
    nm_rmse_abundance = rmse_abundance / mean(abundance),
    nm_rmse_density = rmse_density / mean(density),
    nr_rmse_abundance = rmse_abundance / (max(abundance) - min(abundance)),
    nr_rmse_density = rmse_density / (max(density) - min(density))
  ) |>
  ungroup() |>
  arrange(simulation, property) |>
  left_join(n_attributes)

write_rds(error_by_property, file.path(path, "abundance_error_by_property.rds"))
rm(error_by_property)
gc()
message("\nabundance error by property done\n")

error_by_simulation <- xn |>
  select(simulation, start_density, all_of(vals)) |>
  group_by(simulation, start_density) |>
  summarise(
    mpe_abundance = mean(abs((value + 1) - (abundance + 1)) / (abundance + 1)) *
      100,
    mpe_density = mean(
      abs((estimated_density + 0.1) - (density + 0.1)) / (density + 0.1)
    ) *
      100,
    mae_abundance = mean(abs(value - abundance)),
    mae_density = mean(abs(estimated_density - density)),
    mbias_abundance = mean(value - abundance),
    mbias_density = mean(estimated_density - density),
    rmse_abundance = sqrt(mean((value - abundance)^2)),
    rmse_density = sqrt(mean((estimated_density - density)^2)),
    rmsle_abundance = sqrt(mean((log(value + 1) - log(abundance + 1))^2)),
    rmsle_density = sqrt(mean(
      (log(estimated_density + 1) - log(density + 1))^2
    )),
    nm_rmse_abundance = rmse_abundance / mean(abundance),
    nm_rmse_density = rmse_density / mean(density),
    ns_rmse_abundance = rmse_abundance / sd(abundance),
    ns_rmse_density = rmse_density / sd(density),
    nr_rmse_abundance = rmse_abundance / (max(abundance) - min(abundance)),
    nr_rmse_density = rmse_density / (max(density) - min(density)),
    sd_ratio_abundance = sd(value) / sd(abundance),
    sd_ratio_density = sd(estimated_density) / sd(density)
  ) |>
  ungroup()


write_rds(
  error_by_simulation,
  file.path(path, "abundance_error_by_simulation.rds")
)
rm(error_by_simulation)
rm(xn)
gc()

message("\nabundance error by simulation done\n")

message("=== ABUNDANCE DONE ===")
