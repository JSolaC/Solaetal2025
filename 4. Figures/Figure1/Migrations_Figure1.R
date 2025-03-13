
# ============================================================================
# Supplementary R Code: Global Biomass and Migration Data Visualization
# ============================================================================
# Authors: Sola et al.
# Date: February 2025
# Journal: Nature Portfolio
# ============================================================================

# Overview:
# This script processes species biomass and migration datasets, integrates human 
# population biomass trends, and generates visual representations of biomass 
# changes across taxonomic groups over time.

# Key Considerations:
# - Species selection: A predefined list of migratory species is used to filter datasets.
# - Data integration: Biomass estimates from multiple sources are combined, including:
#   - Model predictions for species biomass.
#   - Manually curated historical biomass data.
#   - Human population biomass data.
# - Data transformation: Missing data is addressed, and biomass estimates are adjusted 
#   for consistency across datasets.
# - Visualization:
#   - **Figure 1A**: Total migratory biomass change.
#   - **Figure 1B**: Percentage biomass change per taxonomic group.
#   - **Figure 1C**: Mean biomass change across species, highlighting outliers.
# - Statistical adjustments:
#   - Log scaling is applied where necessary to account for large differences in biomass.
#   - Uncertainty is incorporated through credible intervals.

# Outputs:
# - A combined dataset containing biomass estimates, taxonomic classifications, and 
#   migration status.
# - Processed and visualized trends in biomass change across taxa.
# - Bar plots and boxplots illustrating biomass shifts at taxonomic and species levels.
# - Exported visualizations for use in publication and supplementary materials.

# ============================================================================


# 0- Set up the Workspace -------

# Define the required packages
packages <- c("tidyverse", "scales", "ggbreak", "ggh4x", "rnaturalearth", "sf", "raster", "PupillometryR", "cowplot", "Cairo")

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
    Year_after = max(Year)
  )

Human_data <- read.csv("/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/KAUST_Projects/Migrations_Projects/Data_wrangling/Final_R_scripts/Datasets/Human_migrations.csv")

# Add species that do not change over time, mantaining the same levels of biomass
Missing_species_biomass <- read.csv("2-SpeciesBiomass/Outputs/Biomass_all_spp.csv") %>%
  filter(Migratory == "Y" & Binomial %ni% Species_biomass_data$Species) %>%
  rowwise() %>%
  mutate(
    Biomass_before = Biomass,
    Biomass_after = Biomass, 
    Biomass_before_Upper = pmax(Biomass_97, Biomass_2, na.rm = TRUE),
    Biomass_before_Lower = pmin(Biomass_97, Biomass_2, na.rm = TRUE),
    Biomass_after_Upper = pmax(Biomass_97, Biomass_2, na.rm = TRUE),
    Biomass_after_Lower = pmin(Biomass_97, Biomass_2, na.rm = TRUE),
    Year_before = 2000.5,
    Year_after = 2020.5,
    Organisation_level2 = ifelse(Organisation_level2=="Aquatic mammal" | Organisation_level2=="Freshwater mammal", "Marine mammal",
                                 ifelse(Organisation_level2=="Land birds", "Terrestrial birds", 
                                        ifelse(Organisation_level2=="Terrestrial mammal", "Terrestrial mammals", Organisation_level2)))
  ) %>%
  dplyr::select(
    "Species" = Binomial, Organisation_level2, 
    Year_before, Year_after,
    Biomass_before, Biomass_before_Upper, Biomass_before_Lower,
    Biomass_after, Biomass_after_Upper, Biomass_after_Lower
  ) %>%
  distinct()

#
# Figure 1A ####

# Step 1: Process human biomass data, excluding the COVID years --------
Human_biomass_data <- Human_data %>%
  mutate(
    Species = "Homo sapiens",
    Organisation_level2 = "Humans",
    Biomass_before = Biomass[Year==1500],
    Biomass_after = Biomass[Year==2019],
    Year_before = 1500,
    Year_after = 2019
  ) %>%
  distinct(Species, Organisation_level2, Biomass_before, Biomass_after, Year_before, Year_after) %>%
  mutate(
    Biomass_before_Upper = Biomass_before,
    Biomass_after_Upper = Biomass_after,
    Biomass_before_Lower = Biomass_before,
    Biomass_after_Lower = Biomass_after
  )

#
# Step 2: Process species biomass data ----------

# What about species with manual historical biomass - but current biomass not inlcuded in the modelling (e.g. bison bison)?

Historical_biomass_data <- Species_biomass_data %>%
  
  group_by(Species, Organisation_level2, Year_before, Year_after) %>%
  
  dplyr::summarise(
    
    Biomass_before = ifelse(Data!="Missing_EX", Biomass_predicted[Year == min(Year)], Biomass_predicted),  # Biomass at the start
    Biomass_after = ifelse(Data!="Missing_EX", Biomass_predicted[Year == max(Year)], 0),   # Biomass at the end
    
    # I need to include values that do not present Upper!
    Biomass_before_Upper = ifelse(Data!="Missing_EX", Biomass_Upper_predicted[Year == min(Year)], Biomass_Upper_predicted),  # Biomass at the start
    Biomass_after_Upper = ifelse(Data!="Missing_EX", Biomass_Upper_predicted[Year == max(Year)],0),   # Biomass at the end
    
    # I need to include values that do not present Upper!
    Biomass_before_Lower = ifelse(Data!="Missing_EX", Biomass_Lower_predicted[Year == min(Year)], Biomass_Lower_predicted),  # Biomass at the start
    Biomass_after_Lower = ifelse(Data!="Missing_EX", Biomass_Lower_predicted[Year == max(Year)],0),   # Biomass at the end
    
  ) %>%
  
  distinct() %>%
  # Combine human and species biomass data
  bind_rows(Human_biomass_data) %>%
  
  # fix missing rows
  mutate(
    Biomass_before_Upper = ifelse(!is.na(Biomass_before_Upper), Biomass_before_Upper, Biomass_before), # replace cases with Upper = NA for biomass
    Biomass_after_Upper = ifelse(!is.na(Biomass_after_Upper),Biomass_after_Upper,Biomass_after), # replace cases with Upper = NA for biomass
    Biomass_before_Lower = ifelse(!is.na(Biomass_before_Lower), Biomass_before_Lower, Biomass_before), # replace cases with Upper = NA for biomass
    Biomass_after_Lower = ifelse(!is.na(Biomass_after_Lower),Biomass_after_Lower,Biomass_after) # replace cases with Upper = NA for biomass
  )

# To represent the total amount of migratory biomass accounting for species for which there is no historical biomass estimate
# Replace Historical_biomass_data for Historical_biomass_data2 in Steps 3 and 4 - in Supplementary Materials
Historical_biomass_data2 <- Historical_biomass_data %>%
  bind_rows(Missing_species_biomass)

Historical_biomass_data3 <- Historical_biomass_data %>%
  group_by(Organisation_level2) %>%
  summarise(Species_n = n())

#
# Step 3: Create a combined data table for plotting ------------
Migration_table_combined <- Historical_biomass_data %>%
  group_by(Organisation_level2) %>%
  dplyr::summarise(
    Biomass_before = sum(Biomass_before, na.rm = TRUE),
    Biomass_after = sum(Biomass_after, na.rm = TRUE),
    Biomass_before_Upper = sum(Biomass_before_Upper, na.rm = TRUE),
    Biomass_after_Upper = sum(Biomass_after_Upper, na.rm = TRUE)
  ) %>%
  tibble::column_to_rownames('Organisation_level2')

Migration_table_before <- Migration_table_combined %>%
  rownames_to_column(var = "Organisation_level2") %>%  
  arrange(Biomass_before) %>%
  mutate(Biomass_before = Biomass_before + 1) %>%
  dplyr::select(Organisation_level2, Biomass_before) %>%  # Explicitly specify dplyr::select
  column_to_rownames(var = "Organisation_level2") %>%
  as.matrix()

Migration_table_after <- Migration_table_combined %>%
  rownames_to_column(var = "Organisation_level2") %>%  
  arrange(Biomass_after) %>%
  mutate(Biomass_after = Biomass_after + 1) %>%
  dplyr::select(Organisation_level2, Biomass_after) %>%  # Explicitly specify dplyr::select
  column_to_rownames(var = "Organisation_level2") %>%
  as.matrix()

#
# Step 4: Plot the combined data as a barplot -------------
# Define color palette for the bars
groups <- rownames(Migration_table_after)  # Extract group names
color_mapping <- setNames(viridis::viridis(length(groups)), groups)  # Assign a unique color to each group

# Get current plot size from RStudio graphics window
#plot_size <- dev.size("in")

#CairoSVG("biomass_before_barplot.svg", width = plot_size[1], height = plot_size[2])

# Plot bar chart showing biomass before and after for each group
bar_centers_before <- barplot(Migration_table_before,
                              log = "y",            # Log scale for better visualization
                              col = color_mapping[rownames(Migration_table_before)],    # Use color palette for better distinction
                              border = NA,          # Remove border around bars
                              ylim = c(1000000000,10000000000000000),
                              names.arg = NA,
                              args.legend = list(x = "topright")  # Optional: add a legend
)

# Calculate the total Biomass and Upper Bound for 'Before'
total_biomass_before <- sum(Migration_table_combined$Biomass_before, na.rm = TRUE)
total_biomass_before_upper <- sum(Migration_table_combined$Biomass_before_Upper, na.rm = TRUE)

# Add a single error bar at the top of the plot
arrows(x0 = mean(bar_centers_before),  # Center the error bar horizontally
       y0 = total_biomass_before,
       x1 = mean(bar_centers_before),
       y1 = total_biomass_before_upper,
       angle = 90,
       code = 2,
       length = 0.1,
       col = "black")

# Close the SVG file
#dev.off()

# Get current plot size from RStudio graphics window
#plot_size <- dev.size("in")

#CairoSVG("biomass_after_barplot.svg", width = plot_size[1], height = plot_size[2])

bar_centers_after <- barplot(Migration_table_after,
                             log = "y",            # Log scale for better visualization
                             col = color_mapping[rownames(Migration_table_after)],    # Use color palette for better distinction
                             border = NA,          # Remove border around bars
                             ylim = c(1000000000,10000000000000000),
                             names.arg = NA,
                             args.legend = list(x = "topright")  # Optional: add a legend
)

# Calculate the total Biomass and Upper Bound for 'After'
total_biomass_after <- sum(Migration_table_combined$Biomass_after, na.rm = TRUE)
total_biomass_after_upper <- sum(Migration_table_combined$Biomass_after_Upper, na.rm = TRUE)

# Add a single error bar at the top of the plot
arrows(x0 = mean(bar_centers_after),  # Center the error bar horizontally
       y0 = total_biomass_after,
       x1 = mean(bar_centers_after),
       y1 = total_biomass_after_upper,
       angle = 90,
       code = 2,
       length = 0.1,
       col = "black")

# Close the SVG file
#dev.off()

# Function to display color-label correspondence
show_color_mapping <- function(data_table, color_scale = viridis::viridis) {
  # Extract labels from the table's row names
  labels <- rownames(data_table)
  
  # Generate the color palette based on the number of labels
  colors <- color_scale(length(labels))
  
  # Create a data frame to map labels to colors
  color_mapping <- data.frame(Label = labels, Color = colors)
  
  # Print the mapping table for inspection
  print(color_mapping)
  
  # Plot to show correspondence between labels and colors
  plot(1:length(labels), rep(1, length(labels)), pch = 15, cex = 2,
       col = colors, xlab = "Label Index", ylab = "Color", main = "Color Correspondence")
  text(1:length(labels), rep(1.2, length(labels)), labels, srt = 90, adj = 1, cex = 0.8)
}

# Example usage with Migration_table_combined
show_color_mapping(Migration_table_after)

#
# FIGURE 1B ####

# 1) Sum biomass per Taxonomic group and compute % Biomass change --------------
# BIOMASS CHANGE 

#PART 1: Calculate biomass change (%) for historical biomass and weight it on historical biomass

Historical_biomass_change_group <- Historical_biomass_data %>%
  
  # Filter species that cannot be included in this assessment
  filter(
    Organisation_level2!= "Terrestrial invertebrates" & # remove as only historical spp is the extinct Melanoplus spretus
      Year_before != Year_after
  ) %>%
  
  # Fill gaps for historical manual estimates without CI data
  mutate(
    Sign = ifelse(Biomass_after >= Biomass_before, "Positive", "Negative"),
    Biomass_before_Upper = ifelse(is.na(Biomass_before_Upper), Biomass_before, Biomass_before_Upper),
    Biomass_before_Lower = ifelse(is.na(Biomass_before_Lower), Biomass_before, Biomass_before_Lower),
    Biomass_after_Upper = ifelse(is.na(Biomass_after_Upper), Biomass_after, Biomass_after_Upper),
    Biomass_after_Lower = ifelse(is.na(Biomass_after_Lower), Biomass_after, Biomass_after_Lower)
  ) %>%
  
  # Obtain Biomass change at the taxonomic group level - Gross changes to find winers and losers
  group_by(Organisation_level2, Sign) %>%
  mutate(
    Data_n=n(),
    Species_n = length(unique(Species)), # check the number of species included for each taxonomic group
    Biomass_before_Sign = sum(Biomass_before),
    Biomass_change_Sign = (sum(Biomass_after) - sum(Biomass_before)) / sum(Biomass_before)
  ) %>%
  
  # Obtain Biomass change overall at the taxonomic group level
  group_by(Organisation_level2) %>%
  mutate(
    Biomass_change_all = (sum(Biomass_after) - sum(Biomass_before)) / sum(Biomass_before),
    Sign_all = ifelse(Biomass_change_all>0, "Positive", "Negative"),
    Biomass_change_all_Upper = (sum(Biomass_after_Upper) - sum(Biomass_before_Upper)) / sum(Biomass_before_Upper),
    Biomass_change_all_Lower = (sum(Biomass_after_Lower) - sum(Biomass_before_Lower)) / sum(Biomass_before_Lower),
    Biomass_change_all_CI = ifelse(Sign_all == "Positive", 
                                   abs(Biomass_change_all-max(Biomass_change_all_Upper,Biomass_change_all_Lower))+Biomass_change_all, 
                                   -abs(Biomass_change_all-min(Biomass_change_all_Upper,Biomass_change_all_Lower))+Biomass_change_all
    )
  ) %>%
  
  distinct(Organisation_level2, Data_n, Species_n, Sign, Biomass_change_Sign, Biomass_before_Sign, Sign_all, Biomass_change_all, Biomass_change_all_Upper, Biomass_change_all_Lower, Biomass_change_all_CI) %>%
  
  # Weight Biomass change for winers and losers based on their relative biomass contribution - so that changes are proportional to the overall biomass change
  group_by(Organisation_level2) %>%
  mutate(Biomass_change_Sign_w = Biomass_change_Sign * (Biomass_before_Sign/sum(Biomass_before_Sign)))


plot1b <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", alpha=0.4)+
  geom_hline(yintercept=0, alpha=0.4)+
  geom_hline(yintercept=-1, linetype="dashed", alpha=0.4)+
  geom_bar(Historical_biomass_change_group %>% distinct(Organisation_level2, Biomass_change_all, Biomass_change_all_CI, Sign_all)
           , mapping=aes(x = reorder(Organisation_level2, Biomass_change_all), y = Biomass_change_all, fill = Sign_all), stat = "identity", alpha=0.6) +
  geom_bar(Historical_biomass_change_group
           , mapping=aes(x = reorder(Organisation_level2, Biomass_change_all), y = Biomass_change_Sign_w, fill = Sign), stat = "identity", alpha=0.6) +
  geom_errorbar(Historical_biomass_change_group %>% distinct(Organisation_level2, Biomass_change_all, Biomass_change_all_CI)
                , mapping=aes(x = reorder(Organisation_level2, Biomass_change_all), ymax = Biomass_change_all_CI, ymin = Biomass_change_all), stat = "identity", alpha=0.8, width=0.001) +
  scale_fill_manual(values = c("Positive" = "green4", "Negative" = "red3")) +
  scale_y_break(c(1, 1.1, 2.5, 35000), scales = c(0.09, 0.01, 2)) +
  labs(x = "Species", y = "Biomass Change") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Plot1b.pdf", plot = plot1b, dpi = 300)

#
# FIGURE 1C --------
#
# 1) mean change per species ---------

# Here, species >10e12 are still not added since the point that we are trying to make is that biomass patterns are driven by a few species rather than being diversity-driven

Historical_biomass_change_species <- Historical_biomass_data %>%
  
  # Filter species that cannot be included in this assessment
  filter(
    #Species %ni% unique(Missing_species_biomass[Missing_species_biomass$Data=="Missing",]$Species) & # remove non-extinct species with low biomass
    Organisation_level2!= "Terrestrial invertebrates" & # remove as only historical spp is the extinct Melanoplus spretus
      Organisation_level2 != "Humans"
  ) %>%
  
  # Obtain Biomass change at the species level
  mutate(
    Biomass_change_spp = (Biomass_after - Biomass_before) / Biomass_before,
    Sign = ifelse(Biomass_after >= Biomass_before, "Positive", "Negative")
  ) %>%
  
  filter(Biomass_change_spp!=Inf)

# same colours as viridis!

Plot1c <- ggplot(Historical_biomass_change_species, mapping=aes(x=reorder(Organisation_level2, -Biomass_change_spp, mean), y=Biomass_change_spp,
                                                      fill = Organisation_level2, col = Organisation_level2))+
  geom_hline(yintercept=1, linetype="dashed", alpha=0.4)+
  geom_hline(yintercept=0, alpha=0.4)+
  geom_hline(yintercept=-1, linetype="dashed", alpha=0.4)+
  #geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, alpha=0.5, col=NA)+
  geom_point(aes(size=abs(Biomass_after-Biomass_before)),position = position_jitter(width = .15), size = .25)+
  #note that here we need to set the x-variable to a numeric variable and bump it to get the boxplots to line up with the rainclouds. 
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = .5, colour = "BLACK") + 
  guides(fill=FALSE, colour = FALSE) +
  scale_y_break(c(1,2,20, 30), scales = c(0.09, 2)) +
  scale_fill_manual(values = color_mapping) +    # Apply viridis colors to fill
  scale_color_manual(values = color_mapping) +   # Apply viridis colors to color
  # ggtitle("Raincloud Plot")+
  labs(x = "Species", y = "Biomass Change") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Plot1c.pdf", plot = Plot1c, dpi = 300)

#

# Summary figures for text --------

Historical_biomass_data_wildlife <- Historical_biomass_data %>% filter(Species!="Homo sapiens")

Biomass_change <- (sum(Historical_biomass_data_wildlife$Biomass_after)-sum(Historical_biomass_data_wildlife$Biomass_before))/sum(Historical_biomass_data_wildlife$Biomass_before)

Biomass_difference <- abs(sum(Historical_biomass_data_wildlife$Biomass_after)-sum(Historical_biomass_data_wildlife$Biomass_before))

Biomass_difference_elephants <- (Biomass_difference / 6146000) / 1000000 # biomass change in # million elephans
Biomass_difference_bluewhales <- (Biomass_difference / 93968437.5) / 1000000 # biomass change in # million blue whales

# Comparison with Bar-On et al. 2018 estimates for global biosphere biomass

# Define parameters
carbon_mass_gtc <- 2  # Gt C
migratory_biomass_g <- 2.25 * 10^14  # Migratory biomass in grams

# Convert Gt C to grams
carbon_mass_g <- carbon_mass_gtc * 10^15  

# Define carbon fraction range
carbon_fractions <- c(0.15, 0.30)  # 15% to 30% C content

# Perform calculations using dplyr
results <- tibble(carbon_fraction = carbon_fractions) %>%
  mutate(
    biomass_g = carbon_mass_g / carbon_fraction,  # Convert to biomass in grams
    biomass_pg = biomass_g / 10^15,  # Convert to petagrams (Pg)
    migratory_percentage = (migratory_biomass_g / biomass_g) * 100  # Percentage of migratory biomass
  )


#