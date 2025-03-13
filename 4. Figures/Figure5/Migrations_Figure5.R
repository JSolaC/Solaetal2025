
# ============================================================================
# Supplementary R Code: Produce migration distance and biomass shift values for Mapping
# ============================================================================
# Authors: Sola et al.
# Date: February 2025
# Journal: Nature Portfolio
# ============================================================================

# Overview:
# This script analyzes global migration patterns and biomass changes across 
# taxonomic groups. It integrates data on species movement, historical biomass, 
# and depletion trends to visualize large-scale shifts in migration systems.

# Key Components:
# - Migration distance calculations: 
#   - Combines manual data, biomass estimates, and spatial modeling.
# - Biomass weighting:
#   - Adjusts migration distances based on species-specific biomass proportions.
#   - Accounts for regional variability in biomass shifts.
# - Data transformations:
#   - Normalizes migration distances for visualization.
#   - Applies log transformations for proportional biomass change analysis.
# - Color Mapping:
#   - Custom non-linear color scale to represent biomass depletion trends.
#   - Colors assigned based on species-specific biomass decline or recovery.

# Outputs:
# - Processed datasets with weighted migration distances and biomass shifts.
# - Global-scale visualizations of species movement patterns.
# - Publication-ready maps showing biomass declines over time.

# ============================================================================


# 0- Set up the Workspace -------

# Define the required packages
packages <- c("tidyverse", "scales", "ggbreak", "ggh4x", "rnaturalearth", "sf", "raster", "PupillometryR", "cowplot", "splines","ggspatial")

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
# 1- Load datasets ---------

Species_list <- read.csv('1-SpeciesList/Output/Species_list_all_out.csv') %>%
  dplyr::select(-X)

Mass_migrations_manual <- read.csv("2-SpeciesBiomass/Outputs/Mass_Migrations_data.csv") 

Fish_pops <- read.csv("4-Figures/Figure3/Fish_pop_migrations_biomass.csv") %>%
  
  # standardise names
  left_join(Species_list, by="Binomial") %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial)) %>%
  # select variables and remove duplicates
  group_by(Binomial, Region) %>%
  dplyr::summarise(Biomass_proportion=sum(Biomass_proportion))

Species_migratory_biomass <- read.csv("2-SpeciesBiomass/Outputs/Biomass_all_spp.csv") %>%
  filter(Migratory=="Y")

#
# Figure 5 (need to run Figure 1 steps 1 to 4, and run Figure 4) ---------------------

# Add color to the dataset
assign_colors <- function(x) {
  
  # Define color stops (corresponding to values_map)
  color_stops <- c("black", "darkred", "red", "white", "darkgreen", "#47D45A")
  
  # Define custom value mapping (non-linear scale)
  values_map <- c(-1, -0.99, -0.3, 0, 0.3, 1)
  
  # Ensure x is within the expected range (-1 to 1)
  x <- pmax(pmin(x, 1), -1)
  
  # Generate high-resolution scale (every 0.01)
  high_res_x <- seq(-1, 1, by = 0.01)  # 201 values
  
  # Interpolate positions for high-res scale
  normalized_x <- approx(values_map, seq(0, 1, length.out = length(values_map)), xout = high_res_x)$y
  
  # Create a high-resolution color gradient (1000 colors for smoothness)
  gradient_colors <- colorRampPalette(color_stops)(1000)
  
  # Map normalized_x values to gradient colors
  color_lookup <- gradient_colors[ceiling(normalized_x * 999) + 1]
  
  # Use `findInterval()` to match x to the closest high_res_x value
  closest_indices <- findInterval(x, high_res_x, all.inside = TRUE)
  
  # Return assigned colors based on the closest match
  return(color_lookup[closest_indices])
}

# MIGRATION MAP - SCALES: distance, biomass, depletion (1950s), depletion (1800s-1900s)

Species_biomass_raw <- Species_migratory_biomass %>%
  ungroup() %>%
  dplyr::select(Binomial, Biomass)

Migration_Map_data <- Mass_migrations_manual %>%
  
  dplyr::select(Binomial, Population, Population.proportion, Taxa.group, Region, Population.proportion, 
                Movement_Type, "Migration.distance1"=Migration.distance..km., "Migration.distance2"=Migration.distance2..km2.) %>%
  
  left_join(Fish_pops) %>%
  mutate(Population.proportion = ifelse(Taxa.group=="Oceanodromous fish" | Taxa.group=="Anadromous" | Taxa.group=="Marine invertebrates", Biomass_proportion, Population.proportion)) %>%
  mutate(Population.proportion = ifelse(is.na(Population.proportion), 1, Population.proportion)) %>%
  mutate(Population.proportion = ifelse(Binomial=="Carcharhinus limbatus" & Population=="Africa - Mediterranean", 0, Population.proportion)) %>% # there's no biomass data for C. limbatus in A-Med
  
  left_join(Species_list, by="Binomial") %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial)) %>%
  
  left_join(Historical_biomass_data %>% rename(Binomial=Species)) %>%
  left_join(Species_biomass_raw) %>%
  mutate(Biomass_after=ifelse(is.na(Biomass_after), Biomass, Biomass_after)) %>%
  mutate(Biomass_before=ifelse(is.na(Biomass_before),Biomass_after, Biomass_before)) %>%
  
  # there's some species for which biomass is missiing! Check!
  mutate(
    Biomass_after_proportion = Biomass_after*Population.proportion,
    Biomass_before_proportion = Biomass_before*Population.proportion
  ) %>%
  
  bind_rows(
    .,
    filter(., Binomial == "Bison bison") %>%
      mutate(Binomial="Bison bison EX")
  ) %>%
  
  mutate_cond(Binomial=="Bison bison EX" | Binomial=="Oryx dammah" | Binomial=="Ectopistes migratorius", Biomass_after_proportion = Biomass_before) %>%
  
  dplyr::select(Binomial, Population, Region, Taxa.group, Organisation_level2, 
                Biomass_before, Biomass_after, Biomass_before_proportion, Biomass_after_proportion, 
                Migration.distance1, Migration.distance2) %>%
  
  # convert migration distance to numeric
  dplyr::mutate_at(c(6:11), as.numeric) %>%
  
  mutate(Organisation_level2=ifelse(is.na(Organisation_level2), Taxa.group, Organisation_level2)) %>%
  mutate(Organisation_level2=ifelse(Organisation_level2=="Land birds", "Terrestrial birds", Organisation_level2)) %>%
  mutate(Organisation_level2=ifelse(Organisation_level2=="Aquatic mammals", "Marine mammal", Organisation_level2)) %>%
  mutate(Biomass_change=(Biomass_after-Biomass_before)/(Biomass_before+1)) %>%
  mutate(Biomass_change_proportion_caped = ifelse(Biomass_change>1,1,Biomass_change)) %>%
  mutate(Biomass_change_colour=assign_colors(Biomass_change_proportion_caped)) %>%
  filter(Binomial!="Cervus canadensis" & Binomial!="Stercorarius lonnbergi") # need to remove this entry - it is a bad translation of C elpahus (due to taxize mixing things up) and S. lonnbergi not in list anymore!


# Use this dataset for individual maps (including >10e14) or to check individual species for collective map (excluding >10e14)
Migration_Map_data1 <- Migration_Map_data %>%
  
  filter(Biomass_after_proportion>9.99e10) %>%
  filter(Biomass_after_proportion<1e14) %>%
  
  mutate(
    Biomass_max=max(Biomass_after_proportion, na.rm=T)
  ) %>%
  
  mutate(
    Biomass_length_log = rescale(log(Biomass_after_proportion+1) / log(Biomass_max), to = c(0.02, 1)) * (6), # that is, in relation to the legend corresponding to 5 (1x10^13), what the largest value in the dataset should be
    # to standardise migration distance on the map, I used the metric distance (cm) on the map scale corresponding to 800 km. Thus, to obtain a metric (cm) representation, I need to obtain a conversion factor indicating how many cm is one km (1.88/800)
    # for cases when area rather than distance is recorded, we assume a circumfernce and estimate the diameter as the max distance 
    Distance2_length = (sqrt(Migration.distance2/3.141593)*2)*(1.88/800),
    Distance_length = Migration.distance1*1.88/800
  ) %>%
  
  dplyr::select(Binomial, Population, Region, Organisation_level2, 
                Biomass_before, Biomass_after, Biomass_change, Biomass_change_proportion_caped, Biomass_change_colour,
                Biomass_before_proportion, Biomass_after_proportion, Biomass_length_log, 
                Distance_length, Distance2_length) %>%
  
  ungroup() %>%
  mutate(
    Biomass_max=max(Biomass_after_proportion, na.rm=T),
    Biomass_length_log = rescale(log(Biomass_after_proportion) / log(Biomass_max), to = c(0.02, 1)) * (5)
  )


Migration_Map_data2 <- Migration_Map_data %>%
  left_join(Mass_migrations_manual  %>% dplyr::select(Binomial, "Migration_system"=Migration.System)) %>%
  left_join(Mass_migrations_manual  %>% dplyr::select(Binomial, Movement_Type) %>% filter(Movement_Type!="all")) %>%
  distinct(Binomial, Organisation_level2, Region, Migration_system, Movement_Type, 
           Biomass_before, Biomass_after, Biomass_before_proportion, Biomass_after_proportion) %>%
  filter(Biomass_before != Biomass_after) %>%
  #na.omit() %>%
  
  mutate_cond(Binomial=="Carcharhinus limbatus" & Migration_system=="High to Low Latitudes", Organisation_level2="Marine fish") %>%
  mutate(Organisation_level2=ifelse(Organisation_level2=="Oceanodromous fish", "Marine fish", Organisation_level2)) %>%
  mutate(Organisation_level2=ifelse(Organisation_level2=="Anadromous", "Anadromous fish", Organisation_level2)) %>%
  mutate(Organisation_level2=ifelse(Organisation_level2=="Terrestrial mammal", "Terrestrial mammals", Organisation_level2)) %>%
  
  # summarise other bird movements into three main classes (given the small amount of biomass on each movement type)
  mutate(Movement_Type=ifelse((Organisation_level2=="Terrestrial birds" | Organisation_level2=="Seabirds") & is.na(Region), NA, Movement_Type)) %>%
  mutate(Migration_system=ifelse(Migration_system=="Low to High Latitudes and Coastal/Inshore to Inland" | Migration_system =="Low to High Latitudes and Inland to Coastal" | Migration_system=="Low to High Latitudes and Low to High Altitudes",
                                 "Low to High Latitudes", Migration_system)) %>%
  mutate(Migration_system=ifelse(is.na(Movement_Type) & (Migration_system=="Low to High Latitudes" | Migration_system=="Hight to Low Latitudes" | Migration_system=="High to Low Latitudes" | Migration_system=="Transequatorial"), "Other Latitude movements", Migration_system)) %>%
  mutate(Migration_system=ifelse(is.na(Movement_Type) & (Migration_system=="Coastal/Inshore to Inland" | Migration_system=="Offshore to Inshore/Coastal"), "Offshore-inshore movements", Migration_system)) %>%
  #mutate(Migration_system=ifelse(is.na(Movement_Type) & (Migration_system=="Low to High Altitude" | Migration_system=="Unspecified Altitude Migration"), "Altitude movements", Migration_system)) %>% # altitude movements are marginal
  mutate(Migration_system=ifelse(is.na(Movement_Type) & (Migration_system=="Low to High Altitude" | Migration_system=="Unspecified Altitude Migration" | Migration_system=="Productivity Unspecified" | Migration_system=="Wet to Dry Areas" |
                                                           Migration_system=="Unspecified" | Migration_system=="East to West" | is.na(Migration_system)), "Other movements", Migration_system)) %>%
  #mutate_cond(Organisation_level2=="Seabirds" & Migration_system=="Offshore to Inshore/Coastal", Movement_Type="nomadic") %>%
  
  filter(
    (Binomial!="Carcharhinus limbatus" | Migration_system!="Offshore to Inshore/Coastal") & # only in Med & N Africa - and there's no data for this
      (Binomial!="Balaenoptera physalus" | Movement_Type!="nomadic" | Region=="Mediterranean") & # nomadic only in the Mediterranean
      (Binomial!="Balaenoptera physalus" | Movement_Type!="to-and-fro" | Region!="Mediterranean") & # to-and-fro not in the Mediterranean
      (Binomial!="Thunnus albacares" | Migration_system!="High to Low Latitudes" | Region!="E Pacific") & # ony E to W movement there
      (Binomial!="Thunnus albacares" | Migration_system!="East to West" | Region!="Temperate and Tropical Waters") # only H to L latitudes there
  ) %>%
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%
  
  group_by(Organisation_level2, Movement_Type, 
           Region, Migration_system) %>%
  dplyr::summarise(
    Binomial = paste0(Binomial, collapse = ","),
    Biomass_before_proportion=sum(Biomass_before_proportion, na.rm=T),
    Biomass_after_proportion=sum(Biomass_after_proportion, na.rm=T)
  ) %>%
  mutate(Biomass_change_proportion = (Biomass_after_proportion - Biomass_before_proportion)/(Biomass_before_proportion+1)) %>%
  mutate(Biomass_change_proportion_caped = ifelse(Biomass_change_proportion>1,1,Biomass_change_proportion)) %>%
  mutate(Biomass_change_colour=assign_colors(Biomass_change_proportion_caped)) %>%
  
  filter(Biomass_after_proportion>9.99e10) %>%
  filter(Biomass_after_proportion<1e14) %>%
  
  ungroup() %>%
  mutate(
    Biomass_max=max(Biomass_after_proportion, na.rm=T),
    Biomass_length_log = rescale(log(Biomass_after_proportion+1) / log(Biomass_max), to = c(0.02, 1)) * (5)
  ) %>%
  dplyr::select(Binomial, Organisation_level2, Region, Movement_Type, Migration_system, Biomass_length_log, Biomass_change_proportion, Biomass_change_colour, Biomass_after_proportion) %>%
  mutate(Biomass_ratio=Biomass_after_proportion/10e12)

# In the North Sea, opposite fluxes between S. scombrus (out, larger) and T. trachurus (in, smaller) -> (Sc - Tt) = ((1.056300e+12-159242728971)-(8.460492e+11-2.125784e+11))/(8.460492e+11-2.125784e+11)

# BASEMAP

world = ne_countries(scale = "medium", returnclass = "sf")
ggplot()+
  geom_sf(data = world,
          fill = alpha("grey", 1),
          color = alpha("grey", 1)) +
  annotation_scale(style="ticks",
                   pad_x = unit(13.25, "cm"),
                   pad_y = unit(.75, "cm")) +
  theme_minimal() +
  #scalebar(data = world, location = "bottomright", dist = 1000,
  #         dist_unit = "km", transform = TRUE,  model = "WGS84") +
  #addscalebar(plotunit = "m") +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank() #remove minor gridlines
  )

#