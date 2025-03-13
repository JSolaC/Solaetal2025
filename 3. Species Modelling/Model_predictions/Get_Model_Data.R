
# ============================================================================
# Supplementary R Code: Global Biomass and Migration Model Predictions
# ============================================================================
# Authors: Sola et al.
# Date: February 2025
# Journal: Nature Portfolio
# ============================================================================

# Overview:
# This script processes model outputs  to predict biomass trends, and integrates manual historical data to reconstruct 
# biomass estimates for key migratory species.

# Key Considerations:
# - Bayesian modeling: Predictions are generated using previously trained Bayesian models 
#   for different taxonomic groups.
# - Log Ratio Computation: Population size changes are quantified using log ratios.
# - Biomass estimation: Predictions are backtransformed and scaled based on 
#   known biomass values.
# - Dataset integration: Model predictions and manual estimates are combined for a 
#   comprehensive view of biomass fluctuations.

# Outputs:
# - Predicted biomass trends for migratory species.
# - Processed datasets containing Bayesian model predictions with credible intervals.
# - Integrated dataset including both model-based and manually curated biomass estimates.
# - Exported CSV files for downstream ecological analyses and visualization.

# ============================================================================


# 0- Set up the Workspace -------

# Define the required packages
packages <- c("tidyverse", "brms", "tools")

# Install missing packages
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.rstudio.com/")
})

# Load all packages
lapply(packages, library, character.only = TRUE)

# Load functions
'%ni%' <- Negate('%in%') # conditional for instances within a vector NOT included in another vector
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
} # mutate only at selected rows

# Set up the working directory
setwd("/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/KAUST_Projects/Migrations_Projects/Data_wrangling/Publication/")

#

# 1- Load initial datasets ####

Species_list <- read.csv('1-SpeciesList/Output/Species_list_all_out.csv') %>%
  dplyr::select(-X)

Species_migratory_biomass <- read.csv("2-SpeciesBiomass/Outputs/Biomass_all_spp.csv") %>%
  filter(Migratory=="Y")

LPRAM_data_mass <- read.csv("3-SpeciesModelling/Model_data/Output_datasets/LPRAM_data_filtered.csv") %>%
  filter(!is.na(Pop_size_LR))

Migrations_manual_data <- read.csv("2-SpeciesBiomass/Outputs/Mass_Migrations_data.csv") %>%
  
  # Keep unique records based on specified columns
  distinct("Species"=Binomial, "Organisation_level2"=Taxa.group, "Year"=Historical.Year, "Biomass_predicted"=Hitorical.Biomass, "Biomass_manual"=Biomass) %>%
  filter(!is.na(Biomass_predicted) & Biomass_predicted!="" & (Organisation_level2=="Terrestrial mammals" | Organisation_level2=="Marine mammal" | # Add Terrestrial and Marine mammals
                                        Species=="Ectopistes migratorius" | Species=="Melanoplus spretus" | Species=="Pinguinus impennis")) %>% # Add EX species
  # Make sure variables are OK
  mutate(
    Species = word(Species, 1, 2),
    Biomass_predicted = as.numeric(Biomass_predicted),
  ) %>%
  
  # check if biomass predicted needs to be changed due to wrong manual data input
  left_join(Species_migratory_biomass %>% dplyr::select("Species"=Binomial, Biomass)) %>%
  
  # Keep only essential variables
  distinct(Species, Organisation_level2, Year, Biomass_predicted, Biomass_manual, Biomass)

#
# 2- Get predictions from models (before running, make sure to include datasets in Step 9 from the modelling preparation R script 'Migratory_Species_Model_Data_Pub') ------------

# A) Load previously saved RDS models for analysis
TM_model <- readRDS("3-SpeciesModelling/Modelling/Modelling_output/TM_timeseries.RDS")
MM_model <- readRDS("3-SpeciesModelling/Modelling/Modelling_output/MM_timeseries.RDS")
AF_model <- readRDS("3-SpeciesModelling/Modelling/Modelling_output/AF_timeseries.RDS")
SB_model <- readRDS("3-SpeciesModelling/Modelling/Modelling_output/SB_timeseries.RDS")
ST_model <- readRDS("3-SpeciesModelling/Modelling/Modelling_output/ST_timeseries.RDS")
MI_model <- readRDS("3-SpeciesModelling/Modelling/Modelling_output/MI_timeseries.RDS")
TB_model <- readRDS("3-SpeciesModelling/Modelling/Modelling_output/TB_timeseries.RDS")
MF_model <- readRDS("3-SpeciesModelling/Modelling/Modelling_output/MF_timeseries.RDS")

# B) Obtain predictions
  # This function generates predictions for a given model and dataset.
  # Arguments:
  # - model: a trained model object.
  # - data: a data frame containing species and year data.
  # - model_name: a string specifying the model name.

obtain_predictions <- function(model, data, model_name) {
 
   # Create a unique species-year dataset
  species_year_data <- data %>%
    dplyr::select(Binomial, Year) %>%
    # Add placeholders for 'Site' and 'TimeSeriesID' based on existing columns
    mutate(
      Site = Binomial,
      TimeSeriesID = Binomial
    ) %>%
    distinct()
  
  # Predict values, allowing for new levels in Site
  predictions <- posterior_epred(
    model, newdata = species_year_data, 
    allow_new_levels = TRUE, summary = TRUE
  )
  
  # Create a structured data frame of predictions with credible intervals
  pred_df <- data.frame(
    Species = species_year_data$Binomial,
    Year = species_year_data$Year,
    Prediction = apply(predictions, 2, mean),
    CI_Lower = apply(predictions, 2, quantile, probs = 0.025),
    CI_Upper = apply(predictions, 2, quantile, probs = 0.975),
    Model = model_name
  )
  
  return(pred_df)
}

# Obtain predictions
Predictions <- rbind(
  obtain_predictions(TM_model, TM_time_data_model, "Terrestrial mammal"),
  obtain_predictions(MM_model, MM_time_data_model, "Marine mammal"),
  obtain_predictions(AF_model, AF_time_data_model, "Anadromous fish"),
  obtain_predictions(SB_model, SB_time_data_model, "Seabirds"),
  obtain_predictions(ST_model, ST_time_data_model, "Sea turtles"),
  #obtain_predictions(MI_model, MI_time_data_model, "Marine invertebrates"),
  obtain_predictions(TB_model, TB_time_data_model, "Terrestrial birds"),
  obtain_predictions(MF_model, MF_time_data_model, "Marine fish")
)

# C) Backtransform data and estimate biomass comparing estimates to max Year

Predictions2 <- Predictions %>%
  
  left_join(Species_migratory_biomass %>% dplyr::select("Species"=Binomial, Organisation_level2, Biomass)) %>%
  filter(!is.na(Biomass)) %>%  # Keep records with biomass data

  mutate(
    Prediction_exp = exp(Prediction),
    CI_Lower_exp = exp(CI_Lower),
    CI_Upper_exp = exp(CI_Upper)
  ) %>%
  
  group_by(Species) %>%
  # since we have data for current date and not starting date, we need to use max Year value as 0
  mutate(
    Prediction2 = Prediction_exp / max(Prediction_exp[Year == max(Year)], na.rm = TRUE),
    Biomass_predicted = Biomass * Prediction2,
    Biomass_Lower_predicted = Biomass * (CI_Lower_exp/Prediction_exp),
    Biomass_Upper_predicted = Biomass * (CI_Upper_exp/Prediction_exp),
    max_biomass = max(Biomass_predicted) # this allows to select all species that at some point could have been classified as mass migrations - as per definition with species at present
  ) %>%

  # This is to make sure that each species has 1 data point per species and year
  group_by(Organisation_level2, Species, Year) %>%
  mutate(Count = n())


# D) Add manual data for marine and terrestrial mammals

Missing_data <- Species_migratory_biomass %>%
  filter(Binomial %in% Migrations_manual_data$Species) %>%
  mutate(
    Year=2020,
    Biomass_Lower_predicted = min(Biomass_97, Biomass_2),
    Biomass_Upper_predicted = max(Biomass_97, Biomass_2)
    ) %>%
  mutate(Biomass = ifelse(Binomial=="Ectopistes migratorius" | Binomial=="Melanoplus spretus" | Binomial=="Pinguinus impennis", 0, Biomass)) %>%
  mutate(Biomass = ifelse(Binomial=="Oryx dammah" & Year == 2020, 21900000, Biomass)) %>%
  dplyr::select("Species"=Binomial, Organisation_level2, Year, "Biomass_predicted"=Biomass, Biomass_Lower_predicted, Biomass_Upper_predicted) %>%
  bind_rows(Migrations_manual_data)

Predictions3 <- Predictions2 %>%
  
  filter(Species!="Oryx dammah") %>% # It is computed based on the biomass before disturbance, instead of current biomass
  
  dplyr::select(Species, Organisation_level2, Year, Biomass_predicted, Biomass_Lower_predicted, Biomass_Upper_predicted) %>%
  
  bind_rows(Missing_data) %>%
  
  # for cases where manual data overlaps with predicted data, keep max values only!
  group_by(Species, Year) %>%
  filter(Biomass_predicted == max(Biomass_predicted)) # should we add CI for manual data?

#write.csv(Predictions3, "3-SpeciesModelling/Model_predictions/LP_Predictions.csv")
#

