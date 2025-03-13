
# ============================================================================
# Supplementary R Code: Global Biomass and Migration Data Processing
# ============================================================================
# Authors: Sola et al.
# Date: February 2025
# Journal: Nature Portfolio
# ============================================================================

# Overview:
# This script processes species biomass, migration patterns, and population trends 
# by integrating multiple datasets, including the Living Planet and RAM Legacy 
# Stock Assessment databases. The goal is to generate a unified dataset for analyzing 
# biomass fluctuations and species movement across various taxa.

# Key Considerations:
# - Species selection: A predefined list of migratory species is used to filter datasets.
# - Log Ratio Computation: Population size changes are quantified using log ratios, 
#   with adjustments for time gaps.
# - Spatial structuring: Latitude and longitude are refined to account for 
#   spatial autocorrelation issues.
# - Dataset integration: Living Planet and RAM data are combined into a 
#   single dataset for further modeling.
# - Mass migrations: Species with significant biomass contributions are 
#   identified and the overall dataset is filtered to only include these.
# - Missing data evaluation: Patterns of missing data are explored to assess 
#   potential biases in the datasets.
# - Modeling datasets: Processed data is split by taxonomic groups for 
#   subsequent modeling.

# Outputs:
# - A unified dataset containing biomass estimates, taxonomic classifications, 
#   and migration status.
# - Separate filtered datasets for different taxa, focusing on mass migrations.
# - Processed species-level information with log-transformed population trends 
#   and spatial attributes.
# - Diagnostic outputs, including correlation matrices and missing data visualizations.
# - Exported CSV files ready for ecological modeling and statistical analyses.

# ============================================================================


# 0- Set up the Workspace -------

# Define the required packages
packages <- c("tidyverse", "car", "patchwork", "ggcorrplot")

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
# 1- Load Species list ------------------

Migratory_synonyms <- read.csv('1-SpeciesList/Output/Species_list_all_out.csv') %>%
 
 filter(Migratory=="Y") %>%
  
 # Remove unnecessary columns and rename final variables
 dplyr::select(-X, Original, Binomial)

#
# 2- Get Living Planet data --------

#After a manual review of all units in the LP dataset, we categorised each unit type into a broader category
LP_units <- read.csv("3-SpeciesModelling/Model_data/Initial_datasets/LPD2022_metrics.csv")

# PART 1: Load data, clean and obtain variables

LP_data <- read.csv("3-SpeciesModelling/Model_data/Initial_datasets/LPD2022_public.csv") %>%
  
 # A - FROM WIDE TO LONG DATASET 
 
 # to transform population size values, convert dataset from wide to long
 pivot_longer(
  cols=35:105,
  names_to = "Year",
  values_to = "Pop_size"
 ) %>%
 
 # B - DATA CLEANING
 
 filter(
  # remove empty rows - to reduce computation times in the following steps
  Pop_size!="NULL" &
   # remove instances where no Latitude or Longitude is provided - as we will need that information for the modelling and it's only 42 data points that are missing that information
   Latitude!="NULL"
 ) %>%
 
 mutate(
  # convert Pop_size to a numeric variable
  Pop_size = as.numeric(Pop_size),
  # replace _ for blank space - so that we can later join this dataset with the migratory species list
  Binomial = str_trim(gsub("_"," ",Binomial)),
  # convert Year to a numeric variable (so that it quantifies time)
  Year=as.numeric(gsub("X","",Year)),
 ) %>%
 
 mutate(
  
  # Backtransform the log data: since we will need to comput the log ratio for all time series, any instances that are already transform may introduce further noise into the data
  Pop_size2 = ifelse(Units.class=="Log5 Abundance", Pop_size^5, 
            ifelse(Units.class=="Log Abundance estimate" | Units.class=="log Abundance" | Units.class=="Log Abundance" | Units.class=="Log Biomass" | Units.class=="Ln Abundance estimate" | 
                Units.class=="Ln CPUE" | Units.class=="log CPUE" | Units.class=="Log CPUE" | Units.class=="Log Breeding", Pop_size^exp(1),
               ifelse(Units.class=="Log10 Abundance" | Units.class=="Log10 CPUE", Pop_size^10,
                   ifelse(Units.class=="log Abundance + 1", (Pop_size^exp(1))-1, Pop_size
                   )))),
  # Classify response variables into general categories: Biomass, Biomass estimate, Abundance, Abundance estimate (incl. here CPUE and Breeding estimates) and Unknown
  Units.class2 = ifelse(Units.class=="Density" | Units.class=="Log5 Abundance" | Units.class=="Log10 Abundance" | Units.class=="Log Abundance" | Units.class=="log Abundance" | Units.class=="Log Abundance + 1", "Abundance",
             ifelse(Units.class=="Density estimate" | Units.class=="Abundannce estimate" | Units.class=="Log Abundance estimate" | Units.class=="Ln Abundance estimate", "Abundance estimate",
                 ifelse(Units.class=="Log Biomass" | Units.class=="log Biomass" | Units.class=="Biomass density", "Biomass",
                    ifelse(Units.class=="?", "Unknown",
                        ifelse(Units.class=="CPUE" | Units.class=="Log CPUE" | Units.class=="Log10 CPUE" | Units.class=="log CPUE" | Units.class=="Ln CPUE" | Units.class=="CPUE estimate", "Abundance estimate", 
                           ifelse(Units.class=="Breeding" | Units.class=="Log Breeding", "Abundance estimate", Units.class))))))
  
 ) %>%
  
  left_join(LP_units %>% dplyr::select(-ID, "Units"=Metric, "Units_Category"=Category)) %>%
 
 # C - EXTRACT INFORMATION PER TIME SERIES
 
 # Here, we group the data per time sires (Time Series ID = ID)
 group_by(ID) %>%
 
 # to get lagged Pop size and Year, we need to make sure that the dataset is arranged by Year
 arrange(Year) %>%
 
 mutate(
  # For sensitivity analysis and potential factors influencing the relationship, we extract the first and last year, as well as the duration of the time series
  First_Year = min(Year),
  Last_Year = max(Year),
  Duration_Years=n(),
  # we will need this lagged pop size later for sequential LR
  Pop_size_lag = lag(Pop_size2, n = 1),
  # how many years there is between a data point and its previous data point? We may need that later for sensitivity analyses or to compute alternative LR variables
  Year_gap = Year - lag(Year, n = 1)
 ) %>%
 
 # D - OBTAIN THE RESPONSE VARIABLE: log RATIO (LR)
 
 # to compute LR, we need to consider each time series separately, as LR is computed by comparing each time point to the first time point in the time series
 group_by(ID) %>%
 mutate(
  # two considerations when computing the log ratio: 
  # both the numerator and denominator present a '+1'. That is for two reasons: a) to avoid '-Inf' values when denominator is 0; and b) to avoid '+Inf' when numerator is 0
  # the LR calculated here quantified the total rate of change (which ultimately will depend on 'Duration_Years')
  Pop_size_LR = log((Pop_size2+1)/(Pop_size2[Year==min(Year)]+1))
  # we do not consider Pop_size_LR_std since standardising assumes a linear relationship between Pop size and Year, which is not necessarily representative of the data
 ) %>%
 
 # to compute sequential LR, we compare every time point to the previous time point
 # for LR_seq, remember to remove rows with Year_gap==NA - since these ones will be the first years in the dataset, and therefore lack any previous year to be compared to
 rowwise() %>%
 dplyr::mutate(
  # the LR calculated here quantified the rate of change in comparison to the previous time points (regardless of how many time units are in between the two values)
  Pop_size_LR_seq = log((Pop_size2+1)/(Pop_size_lag+1)),
  # finally, the sequential LR can be further standardised per year gap (so that when comparisons are across years further appart, that potential increase in change is accounted for)
  Pop_size_LR_seq_std = log(((Pop_size2+1)/(Pop_size_lag+1))/Year_gap)
 ) %>%
 
 # E - OBTAIN SPATIAL VARIABLES
 
 ungroup() %>%
 mutate(
  # convert Location to a factor variable, as this may be needed in the subsequent modelling
  Location = as.factor(Location),
  # obtain a variable indicating an intermediate spatial structure between location (where sampling was done) and region (roughly equivalent to continent or biogeographic region)
  Site = as.character(paste(round(as.numeric(Latitude), digits = 2),
               round(as.numeric(Longitude), digits = 2),
               sep=";")),
  
  # for potential complications in the spatial autocorrelation matrices, add a small random number to Latitude and Longitude values that are added in the last decimal position
  # the reason why we are doing this is that lme CorSpatial cannot process perfect correlations
  Latitude_ca = as.numeric(Latitude) + runif(1, 0.00001, 0.0001),
  Longitude_ca = as.numeric(Longitude) + runif(1, 0.00001, 0.0001)
 ) %>%
 
 # check whether there are >1 data points with the same new spatial coordinates
 group_by(Latitude_ca, Longitude_ca) %>%
 mutate(n_spat=n()) 

# PART 2: Remove non-migratory species

LP_filtered <- LP_data %>%
 
 # A - JOIN DATASET WITH MIGRATORY SPECIES LIST
 
 # this will allow to account for the multiple synonyms used to refer to the same species, and thus allow for the correct matching of species with other datasets in this study
 left_join(Migratory_synonyms) %>%
 # we need to compare cases where the species name in the LP dataset is not the main species name used in the Migratory species list
 mutate(Binomial_old = Binomial) %>%
 # indicate cases where the main species name differs between the migratory species list and the LP dataset
 mutate(
  # keep the name as shown in the migratory species list
  Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial),
  # indicate which species name differ - to check if there are any inconsistencies
  Binomial_check = ifelse(Binomial_old!=Original & !is.na(Original), "different", "equal")
 ) %>%
 # Remove cases not matched by the migratory species list
 filter(!is.na(Original))


#write.csv(LP_filtered, "3-SpeciesModelling/Model_data/Output_datasets/LP_data_model.csv")

#
# 3- Get RAM data -------------------------

# load data for the RAM dataset
load("3-SpeciesModelling/Model_data/Initial_datasets/DBdata[asmt][v4.64].RData")

# A- Get dataset, select variables and obtain response variables
RAM_data <- stock %>% 
 dplyr::select(stockid, scientificname, stocklong, areaid, primary_FAOarea, GRSF_uuid) %>% 
 left_join(metadata %>% dplyr::select(stockid, areaname, region)) %>%
 left_join(timeseries) %>%
 filter(!is.na(tsvalue)) %>%
  # obtain log ratio
  group_by(stockid, tsid) %>%
  dplyr::mutate(
    # two considerations when computing the log ratio: 
    # both the numerator and denominator present a '+1'. That is for two reasons: a) to avoid '-Inf' values when denominator is 0; and b) to avoid '+Inf' when numerator is 0
    # the LR calculated here quantified the total rate of change (which ultimately will depend on 'Duration_Years')
    Pop_size_LR = log((tsvalue+1)/(tsvalue[tsyear==min(tsyear)]+1)),
    Pop_size = ((tsvalue+1)/(tsvalue[tsyear==min(tsyear)]+1)),
    # we do not consider Pop_size_LR_std since standardising assumes a linear relationship between Pop size and Year, which is not necessarily representative of the data
    Years_n = n()
  ) %>% # similar to before, there's a warning message that originates from the for WHAKE4T and COD2J3KL. In both cases, it's only one year out of >25 year time series, and it's because they're negative numbers.
  separate(tsid, c("Units2"), sep="-", remove=FALSE)

# B- LOAD DATASET THAT QUANTIFIES THE USEFULNESS OF EACH METRIC TO TRACK TEMPORAL DYNAMICS FOR POPULATIONS
RAM_units <- read.csv("3-SpeciesModelling/Model_data/Initial_datasets/RAM2_Units_scores.csv") %>%
  # remove metrics tracking pop change with mortality - difficult to interpret using LnRR!
  mutate(Score=ifelse(Metric=="Z"|Metric=="F"|Metric=="F/Z"|Metric=="FdivFmgt"|Metric=="FdivFmsy"|Metric=="tFmsy", 0,Score)) %>%
  dplyr::select("Units2"=Metric, Score)
  
# C- FILTER DATASET TO OBTAIN METRICS THAT ALLOW TO TRACK POP CHANGES OVER TIME
RAM_filtered <- RAM_data %>%
 
 # A - JOIN METRIC SCORE AND REMOVE LOW QUALITY SCORES
  
  left_join(RAM_units) %>%
  filter(Score>=7) %>% # remove metrics with a score below 7
  
 rename(Binomial = scientificname) %>%
 
 # B - JOIN DATASET WITH MIGRATORY SPECIES LIST
 
 # this will allow to account for the multiple synonyms used to refer to the same species, and thus allow for the correct matching of species with other datasets in this study
 left_join(Migratory_synonyms) %>%
 # we need to compare cases where the species name in the LP dataset is not the main species name used in the Migratory species list
 mutate(Binomial_old = Binomial) %>%
 # indicate cases where the main species name differs between the migratory species list and the LP dataset
 mutate(
  # keep the name as shown in the migratory species list
  Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial),
  # indicate which species name differ - to check if there are any inconsistencies
  Binomial_check = ifelse(Binomial_old!=Original & !is.na(Original), "different", "equal")
 ) %>%
 # Remove cases not matched by the migratory species list
 filter(!is.na(Original)) %>%
 filter(!is.na(Pop_size_LR)) # remove the two cases with NAs (see warning message at the end of the previous pipeline)

# D- INCLUDE SPATIAL COORDINATES TO THE DATASET - using a custom-made coordinate dataset made by manually checking the location of each stock considered in the dataset
RAM_spat <- RAM_filtered %>%
  left_join(read.csv("3-SpeciesModelling/Model_data/Initial_datasets/RAM_area_coords.csv") %>% 
              ungroup() %>%
              dplyr::select(region, areaname, areaid, Latitude, Longitude) %>%
              distinct(region, areaname, areaid, .keep_all = TRUE))

#write.csv(RAM_spat, "3-SpeciesModelling/Model_data/Output_datasets/RAM_data_model.csv")

#
# 4- Join the two datasets ####

# A- Join the two datasets
LPRAM_data <- bind_rows(
  # load the LP dataset, select variables and add/transform variables to enable comparison
  read.csv("3-SpeciesModelling/Model_data/Output_datasets/LP_data_model.csv") %>%
    dplyr::select(ID, Binomial, Year, "Pop_size"=Pop_size2, Pop_size_LR, "Units2"=Units_Category, "Units"=Units.class2, Region, Site, Location, Latitude, Longitude) %>%
    mutate(
      ID=as.character(ID),
      Dataset = "LivingPlanet"
    ),
  # load the RAM dataset, select variables and add/transform variables to enable comparison
  read.csv("3-SpeciesModelling/Model_data/Output_datasets/RAM_data_model.csv") %>%
    mutate(
      ID=paste(assessid,tsid), # this helps solve the issue that dimensionless timeseries are repeated - and can only be told apart because of their different assessment id
      Dataset = "RAML"
      ) %>%
    dplyr::select(ID, Binomial, "Year"=tsyear, Pop_size, Pop_size_LR, "Units"=Units2, Dataset, "Region"=region, "Site"=areaname, "Location"=areaid, Latitude, Longitude)
)

# B- Get migratory species information - only for species for which we have biomass
Migratory_species_taxonomy <- read.csv('2-SpeciesBiomass/Outputs/Biomass_all_spp.csv') %>%
  filter(Migratory=="Y") %>%
  distinct(Binomial, Organisation_level2, Migration, Migration2)

# C- Add migratory data to global dataset - taxonomy information
LPRAM_data_migratory <- LPRAM_data %>%
  # Keep only species that are migratory
  filter(Binomial %in% Migratory_species_taxonomy$Binomial) %>%
  
  #Add taxonomic information
  left_join(Migratory_species_taxonomy) %>%
  mutate(TimeSeriesID = paste(Dataset,ID,Units)) %>%
  
  #obtain min and max year per time series
  group_by(TimeSeriesID) %>%
  dplyr::mutate(
    Year_min = min(Year),
    Year_max = max(Year)
    )

# D- Correct for the removal of some rows in the previous step and account for missing data
LPRAM_data_migratory2 <- LPRAM_data_migratory %>%
  
  # as a blank canvas, create a dataset with all combinations of years and time series
  ungroup() %>%
  tidyr::expand(TimeSeriesID,Year) %>%
  
  # Add the response variable for each Year x TimeSeriesID combination (and leave NAs where applicable)
  left_join(LPRAM_data_migratory %>% dplyr::select(TimeSeriesID, Year, Pop_size, Pop_size_LR)) %>%
  # Similarly, add all covariates for each TimeSeriesID (as some years may or may not be present?)
  left_join(LPRAM_data_migratory %>% distinct(TimeSeriesID, Binomial, Year_min,Year_max,Region, Site, Location, Latitude,Longitude, Organisation_level2, Units, Units2)) %>%
  # Acotate each time series according to their min and max year registered - this way we also have NAs corresponding to the gaps in years (used later to evaluate missing data)
  filter(Year>=Year_min & Year<=Year_max)

# E- 
LPRAM_data_migratory3 <- LPRAM_data_migratory2 %>%
  
  # Obtain a simplified Units column
  mutate(col_split = str_split(Units, "-")) %>%  
  unnest_wider(col_split, names_sep = "_") %>%
  
  #Select rows and prune variables
  dplyr::select(TimeSeriesID, Binomial, Year, Pop_size, Pop_size_LR, Year_min, Year_max, Region, Site, Location, Latitude, Longitude, Organisation_level2, Units, Units2, "Units3"=col_split_1) %>%
  mutate(Units2=str_trim(Units2)) %>%
  
  # Get information used for missing data assessments later
  group_by(TimeSeriesID) %>%
  dplyr::mutate(
    total_data_points = n(),
    missing_count = sum(is.na(Pop_size_LR)),
    year_count = (Year_max-Year_min+1) # we need to add one year to the difference to count the starting year as well
  ) %>%
  dplyr::mutate(
    missing_percentage = (missing_count / total_data_points) * 100
    ) %>%
  filter(total_data_points > 2) # since we only want to know trends, we include all time series with at least 3 time points

#write.csv(LPRAM_data_migratory3, "3-SpeciesModelling/Model_data/Output_datasets/LPRAM_data_model.csv")

#
# 5- Filter data within the LPRAM dataset to target mass migrations -------

LPRAM_data_mass <- LPRAM_data_migratory3 %>%

  # A- Compute multiple levels of geographical location in the model
  dplyr::mutate(
    Region_Site=Site,
    Location=iconv(Location, "UTF-8", "ASCII", sub = ""),
    Site=paste(Latitude,Longitude,sep = "_")
  ) %>%
  
  # B- Remove empty rows signaling missing data - do not run if wanting to check missing data influence on data!
  #filter(!is.na(Pop_size_LR)) %>%
  
  # C- Add biomass data to the dataset
  left_join(read.csv("2-SpeciesBiomass/Outputs/Biomass_all_spp.csv") %>% 
              dplyr::select(Binomial, Biomass) %>% distinct(), by="Binomial") %>%
  
  # D- filter species by biomass - keep only species which at some point in their time series they may have been > 10e11
  group_by(Binomial) %>%
  mutate(Biomass_max = exp(max(Pop_size_LR, na.rm = T)) * Biomass) %>%
  # Only count species with estimated max biomass > 10e11
  filter(Biomass_max > 100000000000 | Organisation_level2=="Marine invertebrates" | Organisation_level2=="Sea turtles") %>% 
  
  # E- Establish a hierarchical order in terms of quality across variables - so that when multiple variables co-occur, only keep the best ones
  ungroup() %>%
  mutate(
    Units_class = case_when(
      Units3=="Abundance" | Units3=="Biomass" | Units3=="TN" | Units3=="TB" ~ 1,
      Units3=="Abundance estimate" | Units3=="Biomass estimate" | Units3=="TBbest" ~ 2,
      Units3=="SSB" ~ 3,
      Units3=="survB_absolute" | Units3=="STB" ~ 4,
      Units3=="SSBm" | Units3=="SSBf" ~ 5,
      Units3=="R" ~ 6,
      Units3=="TBdivTBmsy" | Units3=="TNdivTNmsy" | Units3=="TBdivTBmgt" ~ 7,
      Units3=="SSBdivSSBmsy" | Units3=="SSBdivSSBmgt" | Units3=="tSSBmsy" ~ 8,
      Units3=="Unknown" ~ 9
    )
  ) %>%
  
  # F- For each species and location, keep only the time series with the best units: given that there's high correlation among units (see section below), that helps reduce collinearity, increases accuracy of estiamtes (based on better information) and added noise in the data
  group_by(Binomial, Site) %>%
  filter(Units_class == min(Units_class)) %>%
  
  # After checking duplicated rows, these have all the same response and covariate values, so simply removing them should be ok.
  distinct()

write.csv(LPRAM_data_mass, "LPRAM_data_filtered.csv")

#
# 6- Check correlations across the different units (run only after excluding step F above) ------------

Check_data <- LPRAM_data_mass %>%
  
  group_by(Binomial, Site) %>%
  mutate(n_Units = length(unique(Units3))) %>%
  
  ungroup() %>%
  filter(n_Units > 1) %>%
  
  group_by(Binomial, Site, Year, Units3, Units_class) %>%
  summarise(Pop_size_LR = mean(Pop_size_LR))


Check_data2 <- Check_data %>%
  dplyr::select(-Units_class) %>%
  pivot_wider(names_from = Units3, values_from = Pop_size_LR) %>%
  ungroup() %>%
  dplyr::select(-Binomial, -Site, -Year) %>%
  mutate(across(everything(), as.numeric))

# Compute correlation matrix (now numeric)
cor_matrix <- cor(Check_data2, use = "pairwise.complete.obs")

print(cor_matrix)

ggcorrplot(cor_matrix, method = "circle", type = "lower", 
           lab = TRUE, lab_size = 3, colors = c("blue", "white", "red"), 
           title = "Correlation Matrix of Units3 Levels (Pop_size_LR)")

#
# 7- Check multicollineality in the data ---------------

Check_multicollinearity <- function(input_dataset){

# Fit a linear model with both numerical and categorical variables
model_mixed <- lm(Pop_size_LR ~ total_data_points + Year + Units2 + Latitude + Longitude, data = input_dataset)

#alias(model_mixed) # If we include binomial, there's perfect multicollinearity among many species

# Calculate the generalized VIF (GVIF) for the model
output <- vif(model_mixed)

return(output)

}
Check_multicollinearity2 <- function(input_dataset){
  
  # Fit a linear model with both numerical and categorical variables
  model_mixed <- lm(Pop_size_LR ~ total_data_points + Year, data = input_dataset)
  
  #alias(model_mixed) # If we include binomial, there's perfect multicollinearity among many species
  
  # Calculate the generalized VIF (GVIF) for the model
  output <- vif(model_mixed)
  
  return(output)
  
}

Check_multicollinearity(LPRAM_data_mass)

Check_multicollinearity(LPRAM_data_mass %>% filter(Organisation_level2=="Terrestrial mammal"))
Check_multicollinearity(LPRAM_data_mass %>% filter(Organisation_level2=="Marine mammal"))

Check_multicollinearity(LPRAM_data_mass %>% filter(Organisation_level2=="Marine fish"))
Check_multicollinearity(LPRAM_data_mass %>% filter(Organisation_level2=="Anadromous fish"))

Check_multicollinearity(LPRAM_data_mass %>% filter(Organisation_level2=="Terrestrial birds"))
Check_multicollinearity(LPRAM_data_mass %>% filter(Organisation_level2=="Seabirds"))

Check_multicollinearity(LPRAM_data_mass %>% filter(Organisation_level2=="Sea turtles"))

Check_multicollinearity2(LPRAM_data_mass %>% filter(Organisation_level2=="Marine invertebrates"))

#
# 8- Investigate whether missing data is linked to covariates in the data -----------

Plot_missing_data <- function(dataset_input){
  
# HOW MUCH DATA ARE MISSING PER TIME SERIES?

# Calculate the proportion of missing values per time series
All_data_migratory_missing_summary <- dataset_input

# Visualizing the distribution of missing percentages  - if any values < 0 it means that there is more data points than years (that is, >1 timeseries clustered together!)
Plot_A <- ggplot(All_data_migratory_missing_summary %>%
         distinct(TimeSeriesID, missing_percentage)
       , aes(x = missing_percentage)) +
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept=100,col="red", alpha=0.2) +
  geom_vline(xintercept=0,col="red", alpha=0.2) +
  labs(title = "a) Distribution of Missing Data Percentages Across Time Series", x = "% Missing Data", y = "Number of Time Series")

Plot_A_b <- ggplot(All_data_migratory_missing_summary %>%
                     distinct(TimeSeriesID, total_data_points, missing_percentage)
                   , aes(y = missing_percentage, x = total_data_points)) +
  geom_point(alpha=0.1)+
  #xlim(0,200)+
  labs(title = "b) Percentage of missing data by number of data points", x = "Number of Data Points", y = "% Missing Data")

# WHEN IS DATA MISSING (After the start of each time series?)

# Explore missing data over time (by Year or Time Period)
All_data_migratory_missing_by_year <- dataset_input %>%
  group_by(Year) %>%
  summarise(
    missing_count = sum(is.na(Pop_size_LR)),
    total_data_points = n(),
    missing_percentage = (missing_count / total_data_points) * 100
  )

# Visualize missing data trend over time
Plot_B <- ggplot(All_data_migratory_missing_by_year, aes(x = Year, y = missing_percentage)) +
  geom_line() +
  labs(title = "c) Percentage of Missing Data Over Time", x = "Year", y = "% Missing Data")

# WHERE IS DATA MISSING

# Explore missing data by location (e.g., Latitude, Longitude)
All_data_migratory_missing_by_location <- dataset_input %>%
  filter(!is.na(Latitude) & !is.na(Longitude)) %>%
  group_by(Latitude, Longitude) %>%
  summarise(
    missing_count = sum(is.na(Pop_size_LR)),
    total_data_points = n(),
    missing_percentage = (missing_count / total_data_points) * 100
  )

# Visualize missing data distribution by geographical location
Plot_C <- ggplot(All_data_migratory_missing_by_location, aes(x = Longitude, y = Latitude, color = missing_percentage)) +
  geom_point() +
  scale_color_viridis_c() +
  ylim(-90,90) + xlim(-180,180)+
  labs(title = "d) Spatial Distribution of Missing Data", x = "Longitude", y = "Latitude", color = "% Missing Data")

# IS MISSING DATA LINKED TO UNITS OR TAXONOMIC GROUP?

# Investigating missing data across species or categorical variables
All_data_migratory_missing_by_covariate <- dataset_input %>%
  group_by(Units3) %>%
  summarise(
    missing_count = sum(is.na(Pop_size_LR)),
    total_data_points = n(),
    missing_percentage = (missing_count / total_data_points) * 100
  )

# Visualizing missing data across species
Plot_D <- ggplot(All_data_migratory_missing_by_covariate, aes(x = reorder(Units3, -missing_percentage), y = missing_percentage)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "e) Percentage of Missing Data by Units class", x = "Units class", y = "% Missing Data") +
  scale_x_discrete(labels = function(x) substr(x, 1, 15)) +
  ylim(0,100) +
  theme(axis.text.y = element_text(size = 5))

# Investigating missing data across species or categorical variables
All_data_migratory_missing_by_covariate2 <- dataset_input %>%
  group_by(Organisation_level2) %>%
  summarise(
    missing_count = sum(is.na(Pop_size_LR)),
    total_data_points = n(),
    missing_percentage = (missing_count / total_data_points) * 100
  )

# Visualizing missing data across species
Plot_E <- ggplot(All_data_migratory_missing_by_covariate2, aes(x = reorder(Organisation_level2, -missing_percentage), y = missing_percentage)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "f) Percentage of Missing Data by Taxonomic group", x = "Taxonomic group", y = "% Missing Data") +
  ylim(0,100) +
  theme()

Plot_list <- Plot_A+Plot_A_b+Plot_B+Plot_C+Plot_D+Plot_E

return(Plot_list)

}

Plot_missing_data(LPRAM_data_mass)

Plot_missing_data(LPRAM_data_mass %>% filter(Organisation_level2=="Terrestrial mammal"))
Plot_missing_data(LPRAM_data_mass %>% filter(Organisation_level2=="Marine mammal"))

Plot_missing_data(LPRAM_data_mass %>% filter(Organisation_level2=="Marine fish"))
Plot_missing_data(LPRAM_data_mass %>% filter(Organisation_level2=="Anadromous fish"))

Plot_missing_data(LPRAM_data_mass %>% filter(Organisation_level2=="Terrestrial birds"))
Plot_missing_data(LPRAM_data_mass %>% filter(Organisation_level2=="Seabirds"))

Plot_missing_data(LPRAM_data_mass %>% filter(Organisation_level2=="Sea turtles"))

Plot_missing_data(LPRAM_data_mass %>% filter(Organisation_level2=="Marine invertebrates"))

#
# 9- Split and save datasets for modelling (make sure to include step F in stage 5 of this script) ---------

TM_time_data_model <- LPRAM_data_mass %>%
  filter(Organisation_level2=="Terrestrial mammal") %>%
  filter(!is.na(Pop_size_LR))

MM_time_data_model <- LPRAM_data_mass %>%
  filter(Organisation_level2=="Marine mammal") %>%
  filter(!is.na(Pop_size_LR))

AF_time_data_model <- LPRAM_data_mass %>%
  filter(Organisation_level2=="Anadromous fish" | Organisation_level2=="Freshwater fish") %>%
  filter(!is.na(Pop_size_LR))

MF_time_data_model <- LPRAM_data_mass %>%
  filter(Organisation_level2=="Marine fish") %>%
  filter(!is.na(Pop_size_LR))

TB_time_data_model <- LPRAM_data_mass %>%
  filter(Organisation_level2=="Terrestrial birds") %>%
  filter(!is.na(Pop_size_LR))

TB_SBWB_time_data_model <- TB_time_data_model %>%
  left_join(read.csv("2-SpeciesBiomass/Outputs/Biomass_all_spp.csv") %>% 
              dplyr::select(Binomial, Organisation_level2, Organisation_level3)) %>%
  filter(Organisation_level3=="Shorebirds" | Organisation_level3=="Waterbirds")

TB_LB_time_data_model <- TB_time_data_model %>%
  left_join(read.csv("2-SpeciesBiomass/Outputs/Biomass_all_spp.csv") %>% 
              dplyr::select(Binomial, Organisation_level2, Organisation_level3)) %>%
  filter(Organisation_level3!="Shorebirds" & Organisation_level3!="Waterbirds")

SB_time_data_model <- LPRAM_data_mass %>%
  filter(Organisation_level2=="Seabirds") %>%
  filter(!is.na(Pop_size_LR))

ST_time_data_model <- LPRAM_data_mass %>%
  filter(Organisation_level2=="Sea turtles") %>%
  filter(!is.na(Pop_size_LR))

MI_time_data_model <- LPRAM_data_mass %>%
  filter(Organisation_level2=="Marine invertebrates") %>%
  filter(!is.na(Pop_size_LR))

write.csv(MM_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/MM_Model_data.csv")  
write.csv(TM_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/TM_Model_data.csv")  
write.csv(AF_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/AF_Model_data.csv")  
write.csv(MF_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/MF_Model_data.csv")  
write.csv(TB_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/TB_Model_data.csv")  
write.csv(TB_SBWB_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/TB_SBWB_Model_data.csv")  
write.csv(TB_LB_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/TB_LB_Model_data.csv")  
write.csv(ST_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/ST_Model_data.csv")  
write.csv(MI_time_data_model, "3-SpeciesModelling/Modelling/Modelling_datasets/MI_Model_data.csv")  

#