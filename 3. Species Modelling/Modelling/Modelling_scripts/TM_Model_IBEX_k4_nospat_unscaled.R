# SET UP THE WORKING SPACE ####

# Load necessary libraries
packages <- c("tidyverse", "cmdstanr", "brms", "rstan")
lapply(packages, library, character.only = TRUE)

set_cmdstan_path("/home/codinajs/.cmdstan/cmdstan-2.36.0")
cmdstanr::cmdstan_path()

# Function to normalize the spatial distance matrix
norm_range <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Load data using fread for better performance
model_data <- read.csv("TM_Model_data.csv") %>% # [1:100,] %>%
  dplyr::select(TimeSeriesID, Binomial, Latitude, Longitude, Site, Region, Year, Pop_size_LR) %>%
  mutate(
    Binomial = factor(Binomial),
    Site = factor(Site),
    TimeSeriesID = factor(TimeSeriesID)
  )

# -----------------------------------------------
# Specify Informative Priors
# -----------------------------------------------
priors <- c(
  # prior(normal(0, 1), class = "b"),              # not needed when using tensors
  prior(exponential(1), class = "sd"),            
  prior(normal(0, 0.5), class = "ar")            
)

# -----------------------------------------------
# Define Model Formula
# -----------------------------------------------
model_formula <- bf(
  Pop_size_LR ~                                        
    t2(Year, Binomial, bs = c("cr", "re"), k = c(4, 4)) +                                  
    (-1 + Year | TimeSeriesID) +                                     
    (1|Site/Region),                       
  sigma ~ s(Year, k=4, bs="cr"),                                  
  autocor = ~ ar(time = Year, gr = TimeSeriesID, p = 1, cov = TRUE)  
)

# -----------------------------------------------
# Function to Fit the Model with Combined Methods
# -----------------------------------------------
fit_single_model <- function(dataset, priors, model_formula, iterations = 2000) {
  
  message("Fitting model using brms...")
  fit <- brm(
    formula = model_formula, 
    data = dataset,
    family = gaussian(link = "identity"), 
    iter = iterations,  
    chains = 4,  
    cores = 4,  
    prior = priors,
    backend = "cmdstanr",  
    control = list(adapt_delta = 0.9, max_treedepth = 12),
    init = 0.1,
    refresh = 500  
  )
  
  # Return the fitted model
  return(fit)
}

# Run the model

gc()

Bayesian_model <- fit_single_model(dataset = model_data, priors = priors, model_formula = model_formula)

gc()

# Save the result
saveRDS(Bayesian_model, "TM_timeseries.RDS")
