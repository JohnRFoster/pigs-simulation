N <- all_dynamics |>
  select(PPNum, N, property, county, property_area) |>
  distinct() |>
  mutate(
    property = as.numeric(as.factor(property)),
    density = N / property_area,
    n_id = 1:n()
  )

for (c in unique(N$county)) {
  data_county <- N |>
    select(-n_id) |>
    filter(county == c)

  property_areas <- data_county |>
    select(property, property_area) |>
    distinct()

  all_property_areas <- sum(property_areas$property_area)

  proportion_county_sampled <- 0.2
  county_area <- all_property_areas / proportion_county_sampled

  effort <- data_county |>
    group_by(PPNum) |>
    summarise(area_sampled = sum(property_area), total_pigs = sum(N)) |>
    ungroup() |>
    mutate(proportion_area_sampled = area_sampled / county_area)

  process_model <- function(N, zeta, a_phi, b_phi) {
    phi <- rbeta(1, a_phi, b_phi)
    lambda <- N * zeta / 2 + N * phi
    rpois(1, lambda)
  }
  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi
  zeta <- 0.405833

  pp_start <- min(data_county$PPNum)
  pp_end <- max(data_county$PPNum)
  pp_seq <- pp_start:pp_end
  Nc <- rep(0, length(pp_seq))

  N_spin <- round(config$start_density * county_area) # initial abundance
  for (i in 1:6) {
    N_spin <- process_model(N_spin, zeta, a_phi, b_phi)
  }

  Nc[1] <- N_spin
  for (i in 2:length(pp_seq)) {
    Nc[i] <- process_model(Nc[i - 1], zeta, a_phi, b_phi)
  }

  latent_population <- tibble(
    PPNum = pp_seq,
    N = Nc,
    property = 1000 + c,
    property_area = county_area - all_property_areas,
    density = N / property_area
  )

  N_county <- bind_rows(data_county, latent_population)

  county_population <- N_county |>
    group_by(PPNum) |>
    summarise(M = sum(N)) |>
    ungroup()

  M <- N_county |>
    select(PPNum, N, property) |>
    pivot_wider(names_from = property, values_from = N) |>
    right_join(county_population) |>
    arrange(PPNum) |>
    left_join(effort) |>
    mutate(proportion_pigs_sample = total_pigs / M) |>
    select(-area_sampled)

  print(c)
  print(cor(M$proportion_pigs_sample, M$proportion_area_sampled))
}
