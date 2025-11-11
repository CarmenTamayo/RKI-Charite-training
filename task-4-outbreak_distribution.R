## Simulate Outbreak Size Distributions Across Offspring Models
#
# This script simulates transmission chains to explore how outbreak size and
# length vary when changing the reproduction number (R) and the assumed offspring
# distribution. It builds on the concepts outlined in these vignettes:
# https://epiverse-trace.github.io/epichains/articles/epichains.html
# https://epiverse-trace.github.io/epichains/articles/projecting_incidence.html
# currently maintained in
# how-to guide: https://epiverse-trace.github.io/howto/analyses/simulate_transmission/epichains-outbreak-size.html
# Load required packages --------------------------------------------------
# epichains provides the chain simulations, ggplot2 handles plotting, and
# epiparameter delivers curated epidemiological parameters for empirical runs.

library(epichains)
library(ggplot2)
library(epiparameter)

# Choose a seed so results are reproducible -------------------------------

set.seed(1)

# Define outbreak parameter space -----------------------------------------
# Build the scenario grid exploring how different reproduction numbers (R)
# shape the distribution of outbreak statistics under a Poisson offspring model.

statistic <- "size"
offspring_dist <- "rpois"
R <- seq(0.1, 1.0, 0.1)

scenarios <- expand.grid(
  offspring_dist = offspring_dist,
  statistic = statistic,
  R = R,
  stringsAsFactors = FALSE
)

scenarios
# number of simulations to run
n_chains <- 1000

# outbreak size groupings
# Predefine interpretable bins to categorise outbreak sizes/lengths.
breaks <- c(0, 2, 5, 10, 20, 50, 100, Inf)

# Simulate outbreak size distribution (Poisson) ---------------------------
# Generate Monte Carlo realisations of outbreak sizes for each R under Poisson
# offspring, truncating chains that grow beyond the largest bin.

outbreak_list <- vector(mode = "list", length = nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  offspring_dist_fun <- match.fun(scenarios[i, "offspring_dist"])
  
  outbreak_list[[i]] <- epichains::simulate_chain_stats(
    n_chains = n_chains,
    statistic = scenarios[i, "statistic"],
    offspring_dist = offspring_dist_fun,
    lambda = scenarios[i, "R"],
    stat_threshold = breaks[length(breaks) - 1] + 1
  )
}

# Group outbreak sizes ----------------------------------------------------
# Convert raw simulations into proportions per size bin for each scenario.

intervals <- lapply(outbreak_list, cut, breaks = breaks)

prop <- lapply(intervals, function(interval) {
  table(interval) / sum(table(interval))
})

outbreak_size_list <- lapply(prop, as.data.frame)

for (i in seq_len(nrow(scenarios))) {
  outbreak_size_list[[i]]$R <- scenarios[i, "R"]
  outbreak_size_list[[i]]$offspring_dist <- scenarios[i, "offspring_dist"]
  outbreak_size_list[[i]]$statistic <- scenarios[i, "statistic"]
}

outbreak_size1 <- do.call(rbind, outbreak_size_list)

head(outbreak_size1)

# Plot outbreak size distribution -----------------------------------------

ggplot2::ggplot(data = outbreak_size1) +
  ggplot2::geom_col(
    mapping = ggplot2::aes(x = as.factor(R), y = Freq, fill = interval)
  ) +
  ggplot2::scale_x_discrete(name = "Reproduction number (R)") +
  ggplot2::scale_y_continuous(name = "Proportion of outbreaks") +
  ggplot2::scale_fill_brewer(
    name = "Outbreak size",
    palette = "Spectral"
  ) +
  ggplot2::theme_bw()

# Change transmission chain statistic to length ---------------------------
# Reuse the scenario grid but switch to measuring the number of generations.

scenarios$statistic <- "length"

# Simulate outbreak length distribution -----------------------------------
# Run the same Poisson-based simulations to capture how long chains persist.

outbreak_list <- vector(mode = "list", length = nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  offspring_dist_fun <- match.fun(scenarios[i, "offspring_dist"])
  
  outbreak_list[[i]] <- epichains::simulate_chain_stats(
    n_chains = n_chains,
    statistic = scenarios[i, "statistic"],
    offspring_dist = offspring_dist_fun,
    lambda = scenarios[i, "R"],
    stat_threshold = breaks[length(breaks) - 1] + 1
  )
}

# Group outbreak lengths --------------------------------------------------
# Summarise the simulated chain lengths into the predefined bins.

intervals <- lapply(outbreak_list, cut, breaks = breaks)

prop <- lapply(intervals, function(interval) {
  table(interval) / sum(table(interval))
})

outbreak_length_list <- lapply(prop, as.data.frame)

for (i in seq_len(nrow(scenarios))) {
  outbreak_length_list[[i]]$R <- scenarios[i, "R"]
  outbreak_length_list[[i]]$offspring_dist <- scenarios[i, "offspring_dist"]
  outbreak_length_list[[i]]$statistic <- scenarios[i, "statistic"]
}

outbreak_length <- do.call(rbind, outbreak_length_list)

head(outbreak_length)

# Plot outbreak length distribution ---------------------------------------
# Assess how higher R elevates the chance of longer-lasting chains.

ggplot2::ggplot(data = outbreak_length) +
  ggplot2::geom_col(
    mapping = ggplot2::aes(x = as.factor(R), y = Freq, fill = interval)
  ) +
  ggplot2::scale_x_discrete(name = "Reproduction number (R)") +
  ggplot2::scale_y_continuous(name = "Proportion of outbreaks") +
  ggplot2::scale_fill_brewer(
    name = "Outbreak length",
    palette = "Spectral"
  ) +
  ggplot2::theme_bw()


# Change offspring distribution to Negative binomial ----------------------
# Introduce dispersion (k) to model superspreading while keeping R variation.

statistic <- "size"
offspring_dist <- "rnbinom"
R <- seq(0.1, 1.0, 0.1)
k <- c(0.1, 5, 10, 1000)

scenarios <- expand.grid(
  offspring_dist = offspring_dist,
  statistic = statistic,
  R = R,
  k = k,
  stringsAsFactors = FALSE
)

# Simulate outbreak size distribution (Negative binomial) -----------------
# Allow overdispersion so that even subcritical R can yield large outbreaks.

outbreak_list <- vector(mode = "list", length = nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  offspring_dist_fun <- match.fun(scenarios[i, "offspring_dist"])
  
  outbreak_list[[i]] <- epichains::simulate_chain_stats(
    n_chains = n_chains,
    statistic = scenarios[i, "statistic"],
    offspring_dist = offspring_dist_fun,
    mu = scenarios[i, "R"],
    size = scenarios[i, "k"],
    stat_threshold = breaks[length(breaks) - 1] + 1
  )
}

# Group outbreak sizes ----------------------------------------------------
# Calculate per-scenario proportions across size bins for plotting.

intervals <- lapply(outbreak_list, cut, breaks = breaks)

prop <- lapply(intervals, function(interval) {
  table(interval) / sum(table(interval))
})

outbreak_size_list <- lapply(prop, as.data.frame)

for (i in seq_len(nrow(scenarios))) {
  outbreak_size_list[[i]]$R <- scenarios[i, "R"]
  outbreak_size_list[[i]]$k <- scenarios[i, "k"]
  outbreak_size_list[[i]]$offspring_dist <- scenarios[i, "offspring_dist"]
  outbreak_size_list[[i]]$statistic <- scenarios[i, "statistic"]
}

outbreak_size2 <- do.call(rbind, outbreak_size_list)

head(outbreak_size2)

# Plot outbreak size distribution (Negative binomial) ---------------------
# Compare how dispersion alters outbreak-size probabilities across R values.

ggplot2::ggplot(data = outbreak_size2) +
  ggplot2::geom_col(
    mapping = ggplot2::aes(x = as.factor(R), y = Freq, fill = interval)
  ) +
  ggplot2::scale_x_discrete(name = "Reproduction number (R)") +
  ggplot2::scale_y_continuous(name = "Proportion of outbreaks") +
  ggplot2::scale_fill_brewer(
    name = "Outbreak size",
    palette = "Spectral"
  ) +
  ggplot2::facet_wrap(
    facets = c("k"),
    labeller = ggplot2::label_both
  ) +
  ggplot2::theme_bw()


# Simulating outbreak sizes for past outbreaks ----------------------------
# Retrieve empirical offspring distributions to illustrate pathogen-specific
# transmission patterns captured in the epiparameter database.

offspring_dists <- epiparameter::epiparameter_db(
  epi_name = "offspring distribution"
)

# Check all empirical distributions are the same --------------------------
# Confirm a single family (e.g., negative binomial) to streamline sampling.

length(unique(vapply(offspring_dists, family, FUN.VALUE = character(1)))) == 1

# Simulate empirical outbreak size distributions --------------------------
# Recreate outbreak scenarios for each disease using its fitted mean (R) and
# dispersion (k) estimates.

outbreak_list <- vector(mode = "list", length = length(offspring_dists))

for (i in seq_along(offspring_dists)) {
  offspring_dist_fun <- match.fun(paste0("r", family(offspring_dists[[i]])))
  
  outbreak_list[[i]] <- epichains::simulate_chain_stats(
    n_chains = n_chains,
    statistic = "size",
    offspring_dist = offspring_dist_fun,
    mu = epiparameter::get_parameters(offspring_dists[[i]])[["mean"]],
    size = epiparameter::get_parameters(offspring_dists[[i]])[["dispersion"]],
    stat_threshold = breaks[length(breaks) - 1] + 1
  )
}

# Group outbreak sizes ----------------------------------------------------

# Summarise each empirical scenario into outbreak-size proportions.
# paste suffix as some diseases have multiple offspring distributions
diseases <- make.unique(
  vapply(offspring_dists, `[[`, FUN.VALUE = character(1), "disease")
)

intervals <- lapply(outbreak_list, cut, breaks = breaks)

prop <- lapply(intervals, function(interval) {
  table(interval) / sum(table(interval))
})

outbreak_size_list <- lapply(prop, as.data.frame)

for (i in seq_along(offspring_dists)) {
  outbreak_size_list[[i]]$R <- epiparameter::get_parameters(offspring_dists[[i]])[["mean"]]
  outbreak_size_list[[i]]$k <- epiparameter::get_parameters(offspring_dists[[i]])[[
    "dispersion"
  ]]
  outbreak_size_list[[i]]$offspring_dist <- family(offspring_dists[[i]])
  outbreak_size_list[[i]]$disease <- diseases[i]
  outbreak_size_list[[i]]$statistic <- "size"
}

outbreak_size3 <- do.call(rbind, outbreak_size_list)

outbreak_size3 <- outbreak_size3 |>
  dplyr::mutate(
    disease = stringr::str_remove(disease, "\\.\\d+$"),
    interval = factor(interval, levels = unique(interval)), # preserve order
    # Create a label combining R and k for each disease
    disease_label = paste0(disease, "\n(R=", R, ", k=", k, ")")
  )

head(outbreak_size3)

# Plot empirical outbreak sizes -------------------------------------------
# Visualise differences between diseases by presenting outbreak-size mixtures.

ggplot2::ggplot(data = outbreak_size3) +
  ggplot2::geom_col(
    mapping = ggplot2::aes(
      x = as.factor(disease_label),
      y = Freq,
      fill = interval
    )
  ) +
  ggplot2::scale_x_discrete(
    name = "Disease"
  ) +
  ggplot2::scale_y_continuous(name = "Proportion of outbreaks") +
  ggplot2::scale_fill_brewer(
    name = "Outbreak size",
    palette = "Spectral"
  ) +
  coord_flip() +
  ggplot2::theme_bw()

