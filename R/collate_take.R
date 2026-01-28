# --------------------------------------------------------------------
#
# Script for collating take from simulation results
#
# John Foster
#
# --------------------------------------------------------------------

message("=== TAKE ===")

# will get config_name from here
source("R/functions_collate.R")
config <- config::get(config = config_name)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

read_path <- get_path("read", config_name, task_id)

density_tasks <- list.files(read_path)
message("Tasks to collate ", length(density_tasks))

nodes <- "take"
tasks_ls <- get_tasks(density_tasks, read_path, nodes)
all_y <- tasks_ls$all_y
all_take <- tasks_ls$all_take

path <- get_path("write", config_name, task_id)
if (!dir.exists(path)) {
	dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_rds(all_take, file.path(path, "all_take.rds"))

take_summaries <- all_y

write_rds(take_summaries, file.path(path, "take_summaries.rds"))
rm(take_summaries)
gc()
message("\nposterior take summaries done\n")

error_by_observation <- tasks_ls$all_by_observation

write_rds(
	error_by_observation,
	file.path(path, "take_error_by_observation.rds")
)
message("\nposterior take error by observation done\n")

take_effort_summary <- tasks_ls$all_effort_summary

write_rds(take_effort_summary, file.path(path, "take_effort_summary.rds"))
rm(take_effort_summary)
rm(error_by_observation)
gc()
message("\nposterior take error summary done\n")

error_by_simulation <- tasks_ls$all_by_simulation

write_rds(error_by_simulation, file.path(path, "take_error_by_simulation.rds"))
rm(error_by_simulation)
gc()
message("\nposterior take error by observation done\n")

error_by_simulation_method <- tasks_ls$all_by_simulation_method

write_rds(
	error_by_simulation_method,
	file.path(path, "take_error_by_simulation_method.rds")
)
rm(error_by_simulation_method)
gc()
message("\nposterior take error by simulation method done\n")

error_by_property <- tasks_ls$all_by_property

write_rds(error_by_property, file.path(path, "take_error_by_property.rds"))
rm(error_by_property)
gc()
message("\nposterior take error by property done\n")


message("=== TAKE DONE ===")
