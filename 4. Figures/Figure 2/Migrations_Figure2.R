
# ============================================================================
# Supplementary R Code: Species temporal Trends in Biomass Change
# ============================================================================
# Authors: Sola et al.
# Date: February 2025
# Journal: Nature Portfolio
# ============================================================================

# Overview:
# This script processes temporal trends in species biomass estimates, computes 
# biomass change rates, and generates visualizations illustrating changes in 
# species biomass over time.

# Key Considerations:
# - Species selection: A predefined dataset of migratory species is used.
# - Data processing:
#   - Identifies the years corresponding to minimum, maximum, and current biomass.
#   - Computes biomass change, change rate, and proportional change.
#   - Adjusts missing or inconsistent biomass values.
# - Visualization:
#   - Species temporal changes in biomass.
#   - Biomass change rate plots: Comparison of biomass change rates for different taxa.
#   - Species-specific trends: Identification of species showing the most extreme changes.
# - Statistical adjustments:
#   - Log transformations are applied for better visualization.
#   - Credible intervals are included to account for uncertainty.
#   - Scaling techniques are used to balance variations in biomass data.

# Outputs:
# - A structured dataset containing biomass change values per species, rates, and proportional changes.
# - Processed trends in biomass change across migratory species.
# - Visualizations of biomass changes over time, highlighting key taxa and species.
# - Exported figures for use in publication and supplementary materials.

# ============================================================================


# 0- Set up the Workspace -------

# Define the required packages
packages <- c("tidyverse", "scales", "ggbreak", "ggh4x", "rnaturalearth", "sf", "raster", "PupillometryR", "cowplot")

# Install missing packages
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.rstudio.com/")
})

# Load all packages
lapply(packages, library, character.only = TRUE)

# Load functions
z_trans <- function(myVar, na.rm){(myVar - mean(myVar, na.rm = T)) / sd(myVar, na.rm = T)}
'%ni%' <- Negate('%in%') # conditional for instances within a vector NOT included in another vector
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
} # mutate only at selected rows

# Set up the working directory
setwd("/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/KAUST_Projects/Migrations_Projects/Data_wrangling/Publication/")

#
# 1 - Load datasets ####

Species_biomass_data <- read.csv("3-SpeciesModelling/Model_predictions/LP_Predictions.csv")  %>%
  
  # Not needed anymore as groups work - but leave here just in case!
  mutate(Organisation_level2 = ifelse(Organisation_level2 == "Terrestrial mammal", "Terrestrial mammals",Organisation_level2)) %>%
  
  group_by(Species, Organisation_level2) %>%
  mutate(
    Data = "Predicted",
    Year_before = min(Year),                                # Earliest year in the data
    Year_after = ifelse(Species=="Oryx dammah", 2006, max(Year))                              # Latest year in the data - issue with O.dammah
  )

#
# FIGURE 2A ######
# 1) Obtain datasets for the figure -----------

# need to recalculate this! I need to get the max biomass BEFORE min!!!
# add max before min point - if max is before min, then max - but if not use max before!

Species_biomass_temporal1 <- Species_biomass_data %>%
  
  # Step 1: Filter and preprocess data
  #bind_rows(Missing_species_biomass) %>% # Should I add data for EX species at present????
  #group_by(Species) %>%
  #filter(n() > 1) %>% # Retain species with more than one record.
  
  # Step 2: Calculate biomass statistics to compute biomass changes and rates
  group_by(Species, Organisation_level2) %>%
  dplyr::summarise(
    
    # Years corresponding to min, max, current and last year before current
    min_year = min(Year[which.min(Biomass_predicted)], Year[which.max(Biomass_predicted)]),  # First date corresponding to either max or min biomass
    max_year = max(Year[which.min(Biomass_predicted)], Year[which.max(Biomass_predicted)]),  # Second date corresponding to either max or min biomass
    # In cases where the year with the min value occurs before the year with the max value, and min_year is not the min of the time series
    prev_year = ifelse(Year[which.min(Biomass_predicted)] < Year[which.max(Biomass_predicted)] & Year[which.min(Biomass_predicted)] > min(Year), 
                       Year[Year < min_year & Biomass_predicted == max(Biomass_predicted[Year < min_year])], # select the year with max biomass predicted that occurs before min_year
                       NA), # if not, return NA
    current_year = max(Year),  # Most recent year in the dataset.
    
    # Mean biomass value
    min_biomass = Biomass_predicted[Year == min_year],  
    max_biomass = Biomass_predicted[Year == max_year],  
    prev_biomass = Biomass_predicted[Year == prev_year],
    current_biomass = Biomass_predicted[Year == current_year],  
    
    # Lower CI biomass value
    min_biomass_Lower = Biomass_Lower_predicted[Year == min_year],  
    max_biomass_Lower = Biomass_Lower_predicted[Year == max_year],  
    prev_biomass_Lower = Biomass_Lower_predicted[Year == prev_year],
    current_biomass_Lower = Biomass_Lower_predicted[Year == current_year],  
    
    # Higher CI biomass value
    min_biomass_Upper = Biomass_Upper_predicted[Year == min_year],  
    max_biomass_Upper = Biomass_Upper_predicted[Year == max_year],  
    prev_biomass_Upper = Biomass_Upper_predicted[Year == prev_year],
    current_biomass_Upper = Biomass_Upper_predicted[Year == current_year]  
    
  ) %>%
  
  distinct()

Species_biomass_temporal2 <- rbind( 
  
  # Step 3: Split dataset into Change 1 and Change 2
  Species_biomass_temporal1 %>% 
    ungroup() %>%
    mutate(Time="1") %>% # Change 1: the earliest in the dataset
    rowwise() %>%
    mutate(
      # values corresponding to current go to max
      max_year = ifelse(!is.na(prev_year), min_year, max_year),
      max_biomass = ifelse(!is.na(prev_year), min_biomass, max_biomass),
      max_biomass_Lower = ifelse(!is.na(prev_year), min_biomass_Lower, max_biomass_Lower),
      max_biomass_Upper = ifelse(!is.na(prev_year), min_biomass_Upper, max_biomass_Upper),
      # values corresponding to last go to min
      min_year = ifelse(!is.na(prev_year), prev_year, min_year),
      min_biomass = ifelse(!is.na(prev_year), prev_biomass, min_biomass),
      min_biomass_Lower = ifelse(!is.na(prev_year), prev_biomass_Lower, min_biomass_Lower),
      min_biomass_Upper = ifelse(!is.na(prev_year), prev_biomass_Upper, min_biomass_Upper)
    ) %>%
    dplyr::select(Species, Organisation_level2, Time, min_year, max_year, min_biomass, max_biomass, min_biomass_Lower, min_biomass_Upper, max_biomass_Lower, max_biomass_Upper)
  , 
  
  Species_biomass_temporal1 %>% 
    mutate(
      Time="2", # Change 2: the most recent one in the dataset - but what happens if there is only one change (that is, if current biomass is also either max or min?)
      # values corresponding to last go to min
      min_year = ifelse(!is.na(prev_year), min_year, max_year),
      min_biomass = ifelse(!is.na(prev_year), min_biomass, max_biomass),
      min_biomass_Lower = ifelse(!is.na(prev_year), min_biomass_Lower, max_biomass_Lower),
      min_biomass_Upper = ifelse(!is.na(prev_year), min_biomass_Upper, max_biomass_Upper),
      # values corresponding to current go to max
      max_year = current_year,
      max_biomass = current_biomass,
      max_biomass_Lower = current_biomass_Lower,
      max_biomass_Upper = current_biomass_Upper
    ) %>%
    dplyr::select(Species, Organisation_level2, Time, min_year, max_year, min_biomass, max_biomass, min_biomass_Lower, min_biomass_Upper, max_biomass_Lower, max_biomass_Upper)
  
) %>%
  
  group_by(Species) %>%
  filter(
    (!(duplicated(paste(min_year, max_year))) | Time == 1) &
      min_year != max_year
  ) %>%
  
  
  # Step 4: Compute biomass changes and rates
  group_by(Species, Organisation_level2, Time, min_year, max_year) %>%
  
  dplyr::summarise(
    
    # Step 4.1: Compute Biomass change values
    Biomass_change =  max_biomass - min_biomass,  # Mean biomass change
    Biomass_Lower_change = max_biomass_Lower - min_biomass_Lower,  # Lower biomass CI change
    Biomass_Upper_change = max_biomass_Upper - min_biomass_Upper, # Upper biomass CI change
    Biomass_CI_change = ifelse(Biomass_change>0, 
                               abs(Biomass_change-max(Biomass_Upper_change,Biomass_Lower_change))+Biomass_change, 
                               -abs(Biomass_change-min(Biomass_Upper_change,Biomass_Lower_change))+Biomass_change),
    
    # Step 4.2: Compute Biomass change rate values
    Biomass_change_rate = Biomass_change / (max_year - min_year),
    Biomass_Lower_change_rate = Biomass_Lower_change / (max_year - min_year),
    Biomass_Upper_change_rate = Biomass_Upper_change / (max_year - min_year),
    Biomass_CI_change_rate = ifelse(Biomass_change_rate>0, 
                                    abs(Biomass_change_rate-max(Biomass_Upper_change_rate,Biomass_Lower_change_rate))+Biomass_change_rate, 
                                    -abs(Biomass_change_rate-min(Biomass_Upper_change_rate,Biomass_Lower_change_rate))+Biomass_change_rate),
    
    # Step 4.3: Compute Proportional Biomass change values
    Biomass_change_proportional = (max_biomass - min_biomass) / min_biomass,
    Biomass_Lower_change_proportional = (max_biomass_Lower - min_biomass_Lower) / min_biomass_Lower,
    Biomass_Upper_change_proportional = (max_biomass_Upper - min_biomass_Upper) / min_biomass_Upper,
    Biomass_CI_change_proportional = ifelse(Biomass_change_proportional>0, 
                                            abs(Biomass_change_proportional-max(Biomass_Upper_change_proportional,Biomass_Lower_change_proportional))+Biomass_change_proportional, 
                                            -abs(Biomass_change_proportional-min(Biomass_Upper_change_proportional,Biomass_Lower_change_proportional))+Biomass_change_proportional),
    
    # Step 4.4: Compute Proportional Biomass change rate values
    Biomass_change_proportional_rate = Biomass_change_proportional / (max_year - min_year),
    Biomass_Lower_change_proportional_rate = Biomass_Lower_change_proportional / (max_year - min_year),
    Biomass_Upper_change_proportional_rate = Biomass_Upper_change_proportional / (max_year - min_year),
    Biomass_CI_change_proportional_rate = ifelse(Biomass_change_proportional_rate>0, 
                                                 abs(Biomass_change_proportional_rate-max(Biomass_Upper_change_proportional_rate,Biomass_Lower_change_proportional_rate))+Biomass_change_proportional_rate, 
                                                 -abs(Biomass_change_proportional_rate-min(Biomass_Upper_change_proportional_rate,Biomass_Lower_change_proportional_rate))+Biomass_change_proportional_rate)
    
  ) %>%
  
  # Step 5: Compute variables to aid with graphic representation
  ungroup() %>%
  mutate(
    
    Change = ifelse(Biomass_change > 0, "Positive", "Negative"),  # Classify change.
    
    Selected_year = ifelse(Biomass_change>0, max_year, min_year), # Year selection based on biomass change direction.
    Selected_year2 = ifelse(Selected_year < 1800, 1800, Selected_year),  # Adjust selected year for historical consistency.
    Year_to_depletion = ifelse(Biomass_change > 0, max_year - min_year, min_year - max_year),
    Year_to_depletion2 = ifelse(Year_to_depletion < -200, -200, 
                                ifelse(Year_to_depletion > 100, 100, Year_to_depletion)),
    
    Biomass_change_propotional_limited = ifelse(Biomass_change_proportional > 1, 1, Biomass_change_proportional),  # Limit proportional changes.
    Biomass_change_log = log(abs(Biomass_change), base = 10) * (Biomass_change / abs(Biomass_change)),  # Log of absolute biomass change.
    Biomass_change_rate_log = log(abs(Biomass_change_rate), base = 10),  # Log of biomass change rate.
    Biomass_change_proportional_log = ifelse(Biomass_change_proportional > 1, log(Biomass_change_proportional, base=10)+1, Biomass_change_proportional),
    Biomass_change_proportional_rate_log = ifelse(Biomass_change_proportional_rate > 1, log(Biomass_change_proportional_rate), Biomass_change_proportional_rate)  # Cap rate changes.
    
  ) %>%
  
  group_by(Species) %>%
  mutate(Biomass_arrange = sum(Biomass_change_proportional_rate_log)) %>%
  filter(!is.na(Biomass_change_proportional_log))

Species_biomass_temporal_selection <- Species_biomass_temporal2 %>%
  
  # Step 1: Select relevant columns for analysis
  dplyr::select(
    Species, Organisation_level2, Change, Time,
    Biomass_change_log,
    Biomass_change_proportional, Biomass_Lower_change_proportional, Biomass_Upper_change_proportional, Biomass_CI_change_proportional,
    Biomass_change_proportional_rate, Biomass_Lower_change_proportional_rate, Biomass_Upper_change_proportional_rate, Biomass_CI_change_proportional_rate
  ) %>%
  
  # Step 2: Handle NaN and Inf values across all numeric columns
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.) | is.infinite(.), NA, .))) %>%
  filter(!is.na(Biomass_change_proportional)) %>% # Filter out rows where Biomass_change_proportional is NA
  
  # Step 3: Compute CI and log values
  rowwise() %>%
  mutate(
    # Biomass change proportional
    Biomass_change_proportional_log = ifelse(Biomass_change_proportional > 1, log(Biomass_change_proportional, base=10)+1, Biomass_change_proportional),
    Biomass_change_proportional_log_abs = abs(Biomass_change_proportional_log),
    Biomass_change_proportional_CI_log = ifelse(Biomass_CI_change_proportional > 1, log(Biomass_CI_change_proportional, base=10)+1, Biomass_CI_change_proportional),
    # Biomass change proportional rate
    Biomass_change_proportional_rate_log = ifelse(Biomass_change_proportional_rate > 1, log(Biomass_change_proportional_rate, base=10)+1, Biomass_change_proportional_rate),
    Biomass_change_proportional_rate_log_abs = abs(Biomass_change_proportional_rate_log),
    Biomass_change_proportional_rate_CI_log = ifelse(Biomass_CI_change_proportional_rate > 1, log(Biomass_CI_change_proportional_rate, base=10)+1, Biomass_CI_change_proportional_rate)
  ) %>%
  
  group_by(Species) %>%
  mutate(
    Biomass_change_proportional_net = sum(Biomass_change_proportional, na.rm=T),
    Biomass_change_proportional_rate_net = sum(Biomass_change_proportional_rate, na.rm=T)
  )

Species_biomass_rate_selection <-  rbind(
  Species_biomass_temporal_selection %>%
    ungroup() %>%
    filter(Time==2) %>%
    slice_max(Biomass_change_proportional, n = 20),
  Species_biomass_temporal_selection %>%
    ungroup() %>%
    filter(Time==2) %>%
    slice_min(Biomass_change_proportional, n = 20)
) #%>%
#dplyr::select(Species, Biomass_change_log, Biomass_change_proportional_rate, Biomass_change)

#
# 2) Plot species trends using ggplot ------------

# Create a common scale for all plots

combined_range <- range(c(min(Species_biomass_temporal2$Biomass_change_log, na.rm=T), max(Species_biomass_temporal2$Biomass_change_log, na.rm=T),
                          min(Species_biomass_temporal_selection$Biomass_change_log, na.rm=T), max(Species_biomass_temporal_selection$Biomass_change_log, na.rm=T),
                          min(Species_biomass_rate_selection$Biomass_change_log, na.rm=T), max(Species_biomass_rate_selection$Biomass_change_log, na.rm=T)))

common_limits <- c(combined_range[1], combined_range[2])

# Produce plots

plot2a <- ggplot() +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = -1, linetype = "dashed", alpha = 0.2) +
  geom_point(Species_biomass_temporal2, mapping=aes(x = Selected_year2, y = Biomass_change_proportional_log, col = Biomass_change_log), alpha=0.9) +
  #geom_rug(Species_biomass_temporal_benchmarks3, mapping = aes(x = min_year), sides = "t", alpha=0.2, length = unit(0.01, "npc")) +  
  #geom_rug(Species_biomass_depletion, mapping=aes(y = Biomass_change_limited), sides = "r", alpha=0.2, length = unit(0.01, "npc")) +
  #scale_colour_manual(values = c("Positive" = "green4", "Negative" = "red3")) +
  scale_color_gradientn(
    limits = common_limits,
    colors = c("red3", "white", "green4"),
    values = scales::rescale(c(-15, -10, 10, 15)),
    oob = scales::squish, # Ensures values outside range are handled
    name = "Biomass Change (log g)"
  ) +
  ggbreak::scale_y_break(c(1, 1.15), scales = c(0.09, 2)) +
  scale_size(range = c(0.01, 10)) +
  labs(x = "Year", y = "Biomass change", color = "Biomass Change (log g)") +
  theme_classic() +
  theme(
    #legend.position = "none",
    axis.text.y = element_text(size = 6)
    #legend.position = c(0.15, 0.15)#,
    #legend.background = element_blank() 
  ) 

ggsave("Plot2a.pdf", plot = plot2a, dpi = 300)

plot2b <- ggplot() +
  #geom_errorbar(Species_biomass_temporal_selection, 
  #              mapping=aes(x = reorder(Species, -Biomass_change_proportional_rate_net), 
  #                          ymin = Biomass_change_proportional_rate_log, ymax=Biomass_change_proportional_rate_CI_log),
  #              width=0.1, alpha=0.5) +
  geom_bar(Species_biomass_temporal_selection, 
           mapping=aes(x = reorder(Species, -Biomass_change_proportional_rate_net), 
                       y = Biomass_change_proportional_rate, fill = Biomass_change_log),
           stat = "identity", alpha=1) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = -1, linetype = "dashed", alpha = 0.2) +
  #scale_fill_manual(values = c("Positive" = "green3", "Negative" = "red3")) +
  scale_fill_gradientn(
    limits = common_limits,
    colors = c("red3", "white", "green4"),
    values = scales::rescale(c(-15, -10, 10, 15)),
    oob = scales::squish, # Ensures values outside range are handled
    name = "Biomass Change (g)"
  ) +
  ggbreak::scale_y_break(c(1,1.1), scales = c(0.15, 2)) +
  labs(y = "Biomass change rate (year -1)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y = element_text(size=7)
  )

ggsave("Plot2b.pdf", plot = plot2b, dpi = 300)

plot2c <- ggplot() +
  geom_bar(Species_biomass_rate_selection, 
           mapping=aes(x = reorder(Species, -Biomass_change_proportional), 
                       y = Biomass_change_proportional_rate, fill = Biomass_change_log),
           stat = "identity", alpha=1) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_hline(yintercept = -1, linetype = "dashed", alpha = 0.2) +
  #scale_fill_manual(values = c("Positive" = "green3", "Negative" = "red3")) +
  scale_fill_gradientn(
    limits = common_limits,
    colors = c("red3", "white", "green4"),
    values = scales::rescale(c(-15, -10, 10, 15)),
    oob = scales::squish, # Ensures values outside range are handled
    name = "Biomass Change (g)"
  ) +
  ggbreak::scale_y_break(c(1,1.1), scales = c(0.15, 5)) +
  labs(y = "Biomass change rate (year -1)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.x = element_blank(), #element_text(angle = 45, margin = margin(t = 45)),
    axis.text.y = element_text(size=7)
  )

ggsave("Plot2c.pdf", plot = plot2c, dpi = 300)

#