## Simulate, clean, validate and plot outbreak data

# Load required R packages ------------------------------------------------
library(simulist)   # To simulate outbreak line lists and contacts
library(cleanepi)   # To validate and clean linelist fields
library(numberize)  # To convert messy text fields into numeric values
library(linelist)   # To tag and manage linelist metadata
library(incidence2) # To summarise incidence over time
library(tibble)     # To organise data as tibbles
library(tidyr)      # To reshape linelist data tidily
library(dplyr)      # To manipulate and summarise data frames

# Choose a seed that results in suitable and reproducible outbreak --------
set.seed(1)

# Simulate outbreak -------------------------------------------------------
line_list <- simulist::sim_linelist() %>% 
  # to tibble for tidier printing
  tibble::as_tibble()

head(line_list)

# Create messy line list data ---------------------------------------------
line_list <- simulist::messy_linelist(line_list, inconsistent_dates = TRUE)
head(line_list)

# Tag line list of data validation ----------------------------------------
# see what tags are available
linelist::tags_names()

# In this case the tags have the same name as the simulated line list's columns
# However, tags can be named differently from column names
line_list <- linelist::make_linelist(
  x = line_list,
  date_onset = "date_onset",
  date_admission = "date_admission",
  date_outcome = "date_outcome",
  date_reporting = "date_reporting",
  id = "id",
  age = "age",
  outcome = "outcome",
  gender = "sex"
  
)
head(line_list)

# Scan line list data for issues ------------------------------------------
# see {cleanepi} website: https://epiverse-trace.github.io/cleanepi/
cleanepi::scan_data(line_list)

# Clean line list ---------------------------------------------------------
# Using {numberize} to convert relevant columns to digits
line_list$age <- numberize::numberize(line_list$age)
line_list$id <- numberize::numberize(line_list$id)

# Check id column
line_list <- cleanepi::check_subject_ids(line_list, target_columns = "id", range = c(1, nrow(line_list)))
print_report(data = line_list, "incorrect_subject_id")

# Tidy column names and remove duplicated rows
line_list <- line_list %>%
  cleanepi::standardize_column_names() %>%
  cleanepi::remove_constants() %>%
  cleanepi::remove_duplicates()
print_report(data = line_list, "found_duplicates")

# Standardise date format
date_columns <- colnames(line_list)[startsWith(colnames(line_list), "date_")]

line_list <- line_list %>%
  cleanepi::standardize_dates(target_columns = date_columns)
print_report(data = line_list, "date_standardization")

# Clean inconsistent sex and case_type columns using a data dictionary 
# Find inconsistencies
line_list %>% count(sex)
line_list %>% count(case_type)
line_list %>% count(outcome)

# Define dictionary
dat_dictionary <- tibble::tribble(
  ~options,     ~values,     ~grp,
  # ---- SEX ----
  "1",          "male",      "sex",      
  "2",          "female",    "sex",      
  "M",          "male",      "sex",      
  "F",          "female",    "sex",      
  "m",          "male",      "sex",      
  "f",          "female",    "sex",
  "Male",       "male",      "sex",
  "male",       "male",      "sex",
  "malb",       "male",      "sex",
  "mald",       "male",      "sex",
  "mmle",       "male",      "sex",
  "myle",       "male",      "sex",
  "Female",     "female",    "sex",
  "female",     "female",    "sex",
  "femvle",     "female",    "sex",
  "femyle",     "female",    "sex",
  
  # ---- CASE TYPE ----
  "suspected",  "suspected", "case_type",
  "Auspected",  "suspected", "case_type",
  
  "confirmed",  "confirmed", "case_type",
  "Ionfirmed",  "confirmed", "case_type",
  "Oonfirmed",  "confirmed", "case_type",
  "confirjed",  "confirmed", "case_type",
  "confirued",  "confirmed", "case_type",
  "confurmed",  "confirmed", "case_type",
  "conhirmed",  "confirmed", "case_type",
  "conyirmed",  "confirmed", "case_type",
  
  "probable",   "probable",  "case_type",
  "Trobable",   "probable",  "case_type",
  "prgbable",   "probable",  "case_type",
  "probabve",   "probable",  "case_type",
  "probbble",   "probable",  "case_type",
  "probhble",   "probable",  "case_type",
  "proiable",   "probable",  "case_type",
  "pvobable",   "probable",  "case_type",
  
  # ---- OUTCOME ----
  "died",        "died",       "outcome",
  "Died",        "died",       "outcome",
  "Mied",        "died",       "outcome",
  "Kied",        "died",       "outcome",
  "diod",        "died",       "outcome",
  
  "recovered",   "recovered",  "outcome",
  "Recovered",   "recovered",  "outcome",
  "Lecovered",   "recovered",  "outcome",
  "recovehed",   "recovered",  "outcome",
  "recoverad",   "recovered",  "outcome",
  "recoverpd",   "recovered",  "outcome",
  "recwvered",   "recovered",  "outcome",
  "rkcovered",   "recovered",  "outcome",
  "rpcovered",   "recovered",  "outcome",
  "rucovered",   "recovered",  "outcome"
)

# Apply dictionary
line_list <- line_list %>% 
  cleanepi::clean_using_dictionary(
    dictionary = dat_dictionary
  )

# Validate clean line list ------------------------------------------------
# Line list is now valid after cleaning
line_list_validated <- linelist::validate_linelist(line_list)

# Now, get data frame with tagged columns only
line_list_validated_tags <- linelist::tags_df(line_list_validated)
head(line_list_validated_tags)

# Aggregate and visualise data --------------------------------------------
# See visualising line list data vignette: https://epiverse-trace.github.io/simulist/articles/vis-linelist.html
# Aggregate to daily incidence data
daily <- incidence2::incidence(
  x = line_list_validated_tags,
  date_index = "date_onset",
  interval = "daily",
  complete_dates = TRUE
)
plot(daily)

# Aggregate to epiweek incidence data
weekly <- incidence2::incidence(
  x = line_list_validated_tags,
  date_index = "date_onset",
  interval = "epiweek",
  complete_dates = TRUE
)
plot(weekly)

# Aggregate and plot onset, hospital admission and death
weekly_by_type <- line_list_validated_tags %>% 
  incidence2::incidence(
    date_index = c("date_onset","date_admission","date_outcome"),
    interval = "epiweek",
    complete_dates = TRUE
  )
plot(weekly_by_type)
