# SET UP THE WORKING SPACE ####

# Load necessary libraries
packages <- c("tidyverse", "cmdstanr", "brms", "rstan")
lapply(packages, library, character.only = TRUE)

set_cmdstan_path("/home/codinajs/.cmdstan/cmdstan-2.35.0")
cmdstanr::cmdstan_path()

# Function to normalize the spatial distance matrix
norm_range <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

cluster_coordinates <- function(data, distance_threshold = 100, min_points = 2) {
  
  # Step 1: Convert lat/lon to X, Y coordinates (X = Longitude, Y = Latitude) in kilometers
  # Explanation: We project latitude and longitude to X and Y in km
  # - 1 degree of latitude ≈ 111.32 km
  # - Longitude varies based on the latitude, so we adjust by cos(latitude)
  coords <- data.frame(
    Latitude = data$Latitude, 
    Longitude = data$Longitude
  )
  
  # Convert lat/lon to X and Y coordinates (X, Y) in km
  coords$X <- coords$Longitude * 111.32 * cos(mean(coords$Latitude) * pi / 180)  # Adjust longitude to km
  coords$Y <- coords$Latitude * 111.32  # Convert latitude to km
  
  # Step 2: Apply DBSCAN clustering on the X, Y coordinates
  # Explanation: Instead of a large geodesic distance matrix, DBSCAN works directly on X, Y (in km)
  cluster_result <- dbscan::dbscan(coords[, c("X", "Y")], eps = distance_threshold, minPts = min_points)
  
  # Step 3: Attach cluster labels to the original data
  # Explanation: Each point is assigned a "Cluster" label, with 0 meaning it's an outlier (not part of any cluster)
  coords$Cluster <- cluster_result$cluster
  
  # Step 4: Calculate the new "central" coordinates (mean latitude, mean longitude) for each cluster
  # Explanation: Calculate the mean latitude and longitude for each cluster (centroid of the cluster)
  cluster_centroids <- coords %>%
    filter(Cluster > 0) %>%  # Exclude outliers (Cluster = 0)
    group_by(Cluster) %>%
    summarize(
      New_Latitude = mean(Latitude, na.rm = TRUE),
      New_Longitude = mean(Longitude, na.rm = TRUE)
    )
  
  # Step 5: Merge the new centroid coordinates back to the original data
  # Explanation: Each point will now have the new central coordinates for the cluster it belongs to
  coords <- coords %>%
    left_join(cluster_centroids, by = "Cluster")
  
  # Return the data with cluster assignments and central coordinates
  return(coords)
}

# Load data using fread for better performance
model_data <- read.csv("TB_Model_data.csv") %>% # [1:100,]
  dplyr::select(TimeSeriesID, Binomial, Latitude, Longitude, Site, Region, Year, Pop_size_LR) %>%
  mutate(Pop_size_LR_scale = scale(Pop_size_LR) + 1e-6) %>%
  mutate(
    Binomial = factor(Binomial),
    Site = factor(Site),
    TimeSeriesID = factor(TimeSeriesID)
  )

model_data_clustered <- cluster_coordinates(model_data) 

model_data <- model_data %>%
  mutate(
    Latitude = model_data_clustered$New_Latitude,
    Longitude = model_data_clustered$New_Longitude,
    Site = paste(model_data_clustered$New_Latitude, model_data_clustered$New_Longitude, sep = "_")
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
    control = list(adapt_delta = 0.8, max_treedepth = 10),
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
saveRDS(Bayesian_model, "TB_timeseries.RDS")
