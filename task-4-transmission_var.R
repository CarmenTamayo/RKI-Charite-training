# Transmission clusters and superspreading events
#
# This script simulates an overdispersed outbreak, visualises the resulting
# contact and transmission networks, fits alternative offspring distributions to outbreak data,
# and explores different metrics influenced by superspreading, such as the
# proportion of secondary cases from different cluster sizes, heterogeneity of transmission,
# and probability of outbreak extinction. 

#
# Recommended reading for more details:
# https://epiverse-trace.github.io/superspreading/articles/epidemic_risk.html
# https://epiverse-trace.github.io/superspreading/articles/proportion_transmission.html
# https://epiverse-trace.github.io/howto/analyses/simulate_transmission/superspreading-probability-extintion.html

# Load required R packages ------------------------------------------------

library(simulist)       # To simulate stochastic outbreaks with overdispersion
library(epicontacts)    # To build and visualise contact networks
library(superspreading) # To quantify superspreading metrics and probabilities
library(fitdistrplus)   # To fit offspring and transmission distributions

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
# Including all contacts from cases, even if these didn't result in disease transmission

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
# Including only contacts that resulted in onward transmission

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
# We create a vector that includes all cases for the complete offspring distribution-
# zeros for non-transmitters and the observed counts for transmitters.
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
# These parameters correspond to the reproduction number and the dispersion parameter, k

R <- nbinom_fit$estimate[["mu"]]
k <- nbinom_fit$estimate[["size"]]

# print estimates
message("R = ", signif(R, digits = 4), "\n", "k = ", signif(k, digits = 4))

# Estimate proportions of cases in different cluster sizes --------
# The function `proportion_cluster_size` tells us the proportion of cases among
# all transmission events that originate within a transmission cluster of a 
# given size, based on specific offspring distribution parameters. 
# In this case, we look at clusters of 2, 5, and 10 cases, given our R and k values.
# Despite the low R, given the high heterogeneity in transmission, over 46% of cases
# happened in clusters of at least 5 people.

superspreading::proportion_cluster_size(
  R = R, 
  k = k, 
  cluster_size = c(2, 5, 10)
)


# Estimate proportion of cases causing 80% transmission -------------------
# ~18% of cases are responsible for 80% of secondary infections in this simulated outbreak; 
# the remaining ~82% of cases together generate only 20% of the spread.

superspreading::proportion_transmission(
  R = R, 
  k = k,
  prop_transmission = 0.8
)

# Estimate probability of outbreak extinction -----------------------------
# Using a branching process with the previously estimated R and k, we calculate
# the probability of outbreak extinction in a population with 1 initial infected individual.
# In this case, given the low R and high dispersion, the probability = 1.

superspreading::probability_extinct(
  R = R,
  k = k,
  num_init_infect = 1
)

# Calculate upper confidence interval for R estimate to explore extinction probability
# We take the higher end of R values to increase the probability of a sustained epidemic,
# and evaluate the impact of modifying the initial number of infected individuals and
# the introduction of control measures on the probability of outbreak extinction.

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
  R = R_upper_bound,
  k = k,
  num_init_infect = 10,
  ind_control = 0.1
)
