# Transmission clusters and superspreading events
#
# This script simulates an overdispersed outbreak, visualises the resulting
# contact and transmission networks, fits alternative offspring distributions,
# and summarises key superspreading indicators such as cluster sizes,
# transmission concentration, and extinction probabilities.
#
# Recommended reading for more details:
# https://epiverse-trace.github.io/superspreading/articles/epidemic_risk.html
# https://epiverse-trace.github.io/superspreading/articles/proportion_transmission.html
# https://epiverse-trace.github.io/howto/analyses/simulate_transmission/superspreading-probability-extintion.html

# Load required R packages ------------------------------------------------

library(simulist) # To simulate outbreaks
library(epicontacts) # To create and plot contact networks
library(superspreading) # To estimate superspreading parameters
library(fitdistrplus) # To fit distribution models

# Choose a seed that results in suitable and reproducible outbreak --------

set.seed(1)

# Simulate outbreak -------------------------------------------------------

## Contact distribution R = 1.2 and k = 0.5
## (i.e. overdispersed contact and transmission).

outbreak <- simulist::sim_outbreak(
  contact_distribution = function(x) dnbinom(x = x, mu = 1.2, size = 0.4),
  prob_infection = 0.5,
  outbreak_size = c(50, 100)
)

# Plot contact network ----------------------------------------------------

contact_net <- epicontacts::make_epicontacts(
  linelist = outbreak$linelist,
  contacts = outbreak$contacts,
  id = "case_name",
  from = "from",
  to = "to",
  directed = TRUE
)

plot(contact_net)


# Plot transmission network -----------------------------------------------

transmission_net <- outbreak$contacts[outbreak$contacts$was_case == TRUE, ]

transmission_net <- epicontacts::make_epicontacts(
  linelist = outbreak$linelist,
  contacts = transmission_net,
  id = "case_name",
  from = "from",
  to = "to",
  directed = TRUE
)

transmission_net

plot(transmission_net)

# Extract secondary case data from outbreak -------------------------------

contacts <- outbreak$contacts

# subset to contacts that caused transmission
infections <- contacts[contacts$was_case == TRUE, ]

# Tabulate number of infections from each infector
secondary_cases <- table(infections$from)

# Calculate number of infections that did not cause any secondary cases
num_no_transmit <- sum(!infections$to %in% infections$from)

# Bind number of secondary cases with those that had no onward transmission
all_cases <- sort(c(rep(0, num_no_transmit), secondary_cases))

# Fit and compare offspring distributions ---------------------------------

# fit a set of offspring distribution models:
# - Poisson
# - Geometric
# - Negative Binomial
pois_fit <- fitdistrplus::fitdist(data = all_cases, distr = "pois")
geom_fit <- fitdistrplus::fitdist(data = all_cases, distr = "geom")
nbinom_fit <- fitdistrplus::fitdist(data = all_cases, distr = "nbinom")

# compare model fits
model_tbl <- superspreading::ic_tbl(
  pois_fit, geom_fit, nbinom_fit
)
model_tbl

# Extract parameters from best fit model ----------------------------------

R <- nbinom_fit$estimate[["mu"]]
k <- nbinom_fit$estimate[["size"]]

# print estimates
message("R = ", signif(R, digits = 4), "\n", "k = ", signif(k, digits = 4))

# Estimate proportions of cases occur in clusters of >= a given si --------

superspreading::proportion_cluster_size(
  R = R, 
  k = k, 
  cluster_size = c(2, 5, 10)
)


# Estimate proportion of cases causing 80% transmission -------------------

superspreading::proportion_transmission(
  R = R, 
  k = k,
  prop_transmission = 0.8
)

# Estimate probability of outbreak extinction -----------------------------

superspreading::probability_extinct(
  R = R,
  k = k,
  num_init_infect = 1
)

# calculate upper confidence interval for R estimate to explore extinction probability
R_upper_bound <- nbinom_fit$estimate[["mu"]] + (qnorm(0.975) * nbinom_fit$sd[["mu"]])

superspreading::probability_extinct(
  R = R_upper_bound, 
  k = k, 
  num_init_infect = 1
)

# increase number of initial introductions seeding outbreaks to assess risk
superspreading::probability_extinct(
  R = R_upper_bound,
  k = k,
  num_init_infect = 10
)

# apply small control measure on transmission to see affect on extinction probability
superspreading::probability_extinct(
  R = R,
  k = k,
  num_init_infect = 10,
  ind_control = 0.1
)
