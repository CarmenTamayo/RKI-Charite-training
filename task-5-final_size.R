# Estimate epidemic size with heterogeneous contacts
#
# This script explores how varying contact patterns reshape epidemic final
# sizes, using social mixing matrices to demonstrate heterogeneity effects.
#
# Recommended reading for more details:
# https://epiverse-trace.github.io/finalsize/articles/varying_contacts.html
# https://epiverse-trace.github.io/howto/analyses/simulate_transmission/finalsize-attack-rate-heterogeneity.html

# Load packages --------------------------------------------------------------------
library(finalsize)  # To compute epidemic attack rates under heterogeneity
library(socialmixr) # To obtain age-structured contact matrices
library(tidyverse)  # To wrangle data and visualise outputs


# Simple quick calculation with homogenous mixing ----------------------------------
# This quick call uses the analytic SIR final-size relation under homogeneous mixing.
# It returns the overall attack rate when everyone mixes uniformly and R0 = 2.
# Treat it as a baseline before layering in age structure and heterogeneity.
r0_input <- 2
finalsize::final_size(r0 = r0_input)

# Set up the transmission model -------------------------------------------
# We now move to an age-structured version of the same calculation.
# The steps below mirror the contact-matrix setup in task 6 so the objects can be reused.

# load contact and population data from socialmixr::polymod
# `socialmixr::polymod` stores diary-based contact matrices for several countries.
# We grab the UK data, restricting to five age bands to keep the example readable.
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 5, 18, 40, 65),
  symmetric = TRUE
)

# prepare contact matrix and demography vector for use in model
# The finalsize package expects the next-generation matrix in column-major order,
# hence we transpose the socialmixr output so that each column represents contacts
# by individuals in age group j.
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population
names(demography_vector) <- rownames(contact_matrix)

# scale the contact matrix so the largest eigenvalue is 1.0
# this is to ensure that the overall epidemic dynamics correctly reflect
# the assumed value of R0
# Scaling by the dominant eigenvalue normalises the matrix so that multiplying
# it by `r0_input` gives a next-generation matrix with spectral radius R0.
contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$values))

# divide each row of the contact matrix by the corresponding demography
# this reflects the assumption that each individual in group {j} make contacts
# at random with individuals in group {i}
# This converts absolute contact counts into per-capita contact rates,
# matching the finalsize model’s assumption about encounter probabilities.
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)

# all individuals are equally and highly susceptible
# Here we set up susceptibility parameters. For simplicity everyone is equally susceptible,
# but the objects are matrices so you can plug in age-specific susceptibility profiles later.
n_susc_groups <- 1L
susc_guess <- 1.0

susc_uniform <- matrix(
  data = susc_guess,
  nrow = n_demo_grps,
  ncol = n_susc_groups
)

# Final size calculations also need to know the proportion of each demographic group {𝑖} 
# that falls into the susceptibility group {𝑗}. This distribution of age groups into 
# susceptibility groups can be represented by the demography-susceptibility distribution matrix.
p_susc_uniform <- matrix(
  data = 1.0,
  nrow = n_demo_grps,
  ncol = n_susc_groups
)

# Using final_size to estimate the proportion of infected individuals per demographic group
output <- finalsize::final_size(
  r0 = r0_input,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susc_uniform,
  p_susceptibility = p_susc_uniform
)

output

# Visualising the results 
output %>% 
  mutate(demo_grp = as_factor(demo_grp)) %>% 
  ggplot(aes(x = demo_grp, y = p_infected)) +
  geom_col() +
  ylim(0,1) +
  labs(
    x = "Age group",
    y = "Proportion infected",
    title = "Final size of an SIR epidemic",
    subtitle = "Fully susceptible population"
  ) +
  theme_bw()

# The bar chart highlights how heterogenous contacts skew attack rates:
# age groups with higher per-capita mixing face larger epidemic burdens
# even though the overall R0 is fixed at 2.
