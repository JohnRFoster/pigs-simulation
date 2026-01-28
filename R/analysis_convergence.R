library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(lubridate)


# script to collate convergence diagnostics across all known parameter and landcover nodes

# will get config_name from here
source("R/functions_collate.R")
source("R/functions_predict.R")
config <- config::get(config = config_name)

args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])

read_path <- get_path("read", config_name, array_id)

density_tasks <- list.files(read_path)
message("Tasks to collate ", length(density_tasks))

check <- numeric(length(density_tasks))

all_beta_p <- tibble()
all_land_cover <- tibble()
all_methods <- tibble()
all_psrf <- tibble()
all_sim_converged <- tibble()

method_h <- tibble(
	method_idx = 1:5,
	method = c("Firearms", "Fixed wing", "Helicopter", "Snares", "Traps")
)

pb <- txtProgressBar(max = length(density_tasks), style = 1)
for (i in seq_along(density_tasks)) {
	task_id <- density_tasks[i]

	rds_file <- file.path(read_path, task_id, "simulation_data.rds")
	# rds_file <- "out/hpc/density_5/1/simulation_data.rds"

	if (file.exists(rds_file)) {
		rds <- read_rds(rds_file)
		rds2 <- read_rds(file.path(read_path, task_id, "knownValues.rds"))

		start_density <- rds$start_density

		psrf <- rds$psrf |>
			as_tibble() |>
			mutate(
				node_names = rownames(rds$psrf),
				node_converged = if_else(`Point est.` > 1.1, FALSE, TRUE)
			) |>
			add_ids(task_id, start_density)

		all_psrf <- bind_rows(all_psrf, psrf)

		sim_converged <- tibble(
			simulation = task_id,
			start_density = start_density,
			converged = all(psrf$`Point est.` <= 1.1),
			bad_mcmc = rds$bad_mcmc,
			sim_died = FALSE
		)

		all_beta_p <- bind_beta_p(all_beta_p, rds, task_id, start_density) |>
			mutate(node_names = paste0("beta_p[", method_idx, ", ", position, "]")) |>
			left_join(method_h, by = "method_idx")

		all_methods <- bind_methods(all_methods, rds, task_id, start_density) |>
			pivot_longer(
				cols = c(p_unique, rho, gamma),
				names_to = "parameter",
				values_to = "actual",
				values_drop_na = TRUE
			) |>
			mutate(
				node_names = case_when(
					parameter == "p_unique" ~ paste0("p_mu[", idx - 3, "]"),
					parameter == "rho" ~ paste0("log_rho[", idx, "]"),
					parameter == "gamma" ~ paste0("log_gamma[", idx - 3, "]")
				)
			) |>
			rename(method_idx = idx) |>
			select(-parameter)

		lc <- rds2$land_cover
		colnames(lc) <- 1:3
		land_cover <- lc |>
			as_tibble() |>
			mutate(county = 1:n()) |>
			pivot_longer(
				cols = -county,
				names_to = "position",
				values_to = "value"
			) |>
			add_ids(task_id, start_density)

		all_land_cover <- bind_rows(all_land_cover, land_cover)
	} else {
		sim_converged <- tibble(
			simulation = task_id,
			start_density = c(0.3, 1.475, 2.65, 3.825, 5)[array_id],
			converged = NA,
			bad_mcmc = NA,
			sim_died = TRUE
		)
	}

	all_sim_converged <- bind_rows(all_sim_converged, sim_converged)

	setTxtProgressBar(pb, i)
}

path <- get_path("write", config_name, array_id)
if (!dir.exists(path)) {
	dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

all_params <- bind_rows(all_beta_p, all_methods)
write_rds(all_methods, file.path(all_params, "all_known_parameters.rds"))

write_rds(all_land_cover, file.path(path, "all_land_cover_lookup.rds"))
write_rds(all_psrf, file.path(path, "all_psrf.rds"))
write_rds(all_sim_converged, file.path(path, "all_sim_converged.rds"))
