
# ============================================================================
# Supplementary R Code: Migration System Shifts and Biomass Change Analysis
# ============================================================================
# Authors: Sola et al.
# Date: February 2025
# Journal: Nature Portfolio
# ============================================================================

# Overview:
# This script analyzes biomass changes across different migration systems, 
# quantifying historical shifts in species movement. 
# It integrates migration type, biomass estimates, and depletion trends.

# Key Components:
# - Migration System Classification:
#   - Categorizes migrations by latitude, altitude, coastal, depth, and productivity.
#   - Groups less significant or unclear movements into "Other Movements".
# - Biomass Change Estimations:
#   - Computes total biomass before and after per migration system.
#   - Includes confidence intervals for uncertainty analysis.
# - Visualization:
#   - Figure 5.1: Barplot comparing biomass before and after migration.
#   - Figure 5.2: Percentage change in biomass by migration system.
# - Color Mapping:
#   - Uses custom color scales for positive (green) and negative (red) biomass shifts.
# - Data Transformations:
#   - Log-transformed biomass changes to compare across orders of magnitude.
#   - Normalization techniques applied to scaling error bars and plots.

# Outputs:
# - Processed datasets summarizing biomass change per migration system.
# - Barplots illustrating historical migration shifts.
# - Credible interval estimates for robust interpretation.
# - Figures for submission illustrating migration-related biomass losses.

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
# FIGURE 4 (run Figure 1 steps 1-4 before running the script) -----------------

Migration_map_data_unfiltered <- read.csv("2-SpeciesBiomass/Outputs/Mass_Migrations_data.csv")  %>%
  dplyr::select("Species" = Binomial, "Taxa_group"=Taxa.group, "Migration_system"=Migration.System) %>%
  mutate(Species = word(Species, 1,2)) %>%
  
  left_join(Historical_biomass_data %>% dplyr::select(Species, Organisation_level2, 
                                                      Biomass_before, Biomass_before_Upper, Biomass_before_Lower, 
                                                      Biomass_after, Biomass_after_Upper, Biomass_after_Lower
  )) %>%
  
  mutate(
    Migration_system2 = case_when(
      
      Migration_system == "Coastal/Inshore to Inland" ~ "Land - Inland",
      Migration_system == "Inland to Coastal/Inshore" ~ "Land - Inshore",
      Migration_system == "Inland" ~ "Land - Unspecified",
      
      Migration_system == "Deep to Shallow" ~ "Depth - Shallower",
      Migration_system == "Shallow to Deep" ~ "Depth - Deeper",
      
      Migration_system == "East to West" ~ "Longitude - West",
      Migration_system == "West to East" ~ "Longitude - East",
      
      Migration_system == "Low to High Altitude" ~ "Altitude - Higher",
      Migration_system == "High to Low Altitude" | Migration_system == "High to Low Altitudes" ~ "Altitude - Lower",
      Migration_system == "Unspecified Altitude Migration" ~ "Altitude - Unspecified",
      
      Migration_system == "Low to High Latitudes" | Migration_system == "Low to High Latitudes_EX"  ~ "Latitude - Higher", 
      Migration_system == "High to Low Latitudes" ~ "Latitude - Lower",
      Migration_system == "Transequatorial" ~ "Latitude - transequatorial",
      Migration_system == "Latitude Unspecified" | Migration_system == "Latitudinal unspecified" ~ "Latitude - Unspecified",
      
      Migration_system == "Inshore/Coastal to Offshore" ~ "Coastal - Offshore",
      Migration_system == "Offshore to Coastal/Inshore" | Migration_system == "Offshore to Inshore/Coastal" ~ "Coastal - Inshore",
      Migration_system == "Unclear inshore-offshore" ~ "Coastal - Unclear",
      
      Migration_system == "Wet to Dry Areas" ~ "Productivity - Lower",
      Migration_system == "Productivity Unspecified" ~ "Productivity - Unspecified",
      
      Migration_system == "Low to High Latitudes and Coastal/Inshore to Inland" | Migration_system == "Low to High Latitudes and Inland to Coastal" ~ "Mixed Latitude and Land",
      Migration_system == "Low to High Latitudes and Low to High Altitudes" ~ "Mixed Latitude and Altitude",
      Migration_system == "Unspecified" ~ "Mixed, unclear or negligible"
      
    )
  )

Migration_map_data <- Migration_map_data_unfiltered %>%
  
  filter(
    !is.na(Species) &
      !is.na(Migration_system) & Migration_system!="DVM" &
      !is.na(Biomass_after) & !is.na(Biomass_before)
  ) %>%
  
  distinct()

Migration_map_data2 <- Migration_map_data %>%
  group_by(Migration_system2) %>%
  summarise(Species_n=n())

Migration_map_summary <- Migration_map_data %>%
  
  group_by(Organisation_level2, Migration_system2) %>%
  dplyr::summarise(
    Biomass_before = sum(Biomass_before),
    Biomass_before_Upper = sum(Biomass_before_Upper),
    Biomass_before_Lower = sum(Biomass_before_Lower),
    Biomass_after = sum(Biomass_after),
    Biomass_after_Upper = sum(Biomass_after_Upper),
    Biomass_after_Lower = sum(Biomass_after_Lower)
  ) %>%
  
  group_by(Migration_system2) %>%
  mutate(
    Biomass_after_all=sum(Biomass_after),
    Biomass_before_all=sum(Biomass_before)
  ) %>%
  
  ungroup() %>%
  mutate(
    Migration_system3 = ifelse((Biomass_before_all<100000000000 & Biomass_after_all<100000000000) | 
                                 Migration_system2=="Mixed, unclear or negligible" | Migration_system2=="Productivity - Unspecified" | Migration_system2=="Coastal - Unclear", 
                               "Other_Movements", Migration_system2)
  ) %>%
  
  group_by(Organisation_level2, "Migration_system2"=Migration_system3) %>%
  dplyr::summarise(
    Biomass_before = sum(Biomass_before),
    Biomass_before_Upper = sum(Biomass_before_Upper),
    Biomass_before_Lower = sum(Biomass_before_Lower),
    Biomass_after = sum(Biomass_after),
    Biomass_after_Upper = sum(Biomass_after_Upper),
    Biomass_after_Lower = sum(Biomass_after_Lower)
  )


# Figure 4.1: Biomass before and after per migration system and organismal group -------

# Step 1: Data preparation
Migration_table_combined <- Migration_map_summary %>%
  
  dplyr::select(-Biomass_before_Lower, -Biomass_after_Lower) %>%
  
  pivot_longer(
    cols = c(Biomass_before, Biomass_after, Biomass_before_Upper, Biomass_after_Upper),
    names_to = "Biomass_type",
    values_to = "Biomass"
  ) %>%
  pivot_wider(
    names_from = Organisation_level2,
    values_from = Biomass,
    values_fill = 0
  ) 

# Step 1: Add a total biomass column for sorting
Migration_table_combined_sorted <- Migration_table_combined %>%
  filter(!str_detect(Biomass_type, "Upper")) %>%
  rowwise() %>%
  mutate(Total_Biomass = sum(c_across(`Anadromous fish`:`Terrestrial mammals`), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Migration_system2) %>%
  mutate(Group_Total_Biomass = sum(Total_Biomass)) %>%
  mutate(Group_Total_Biomass = ifelse(Migration_system2=="Other_Movements", 0, Group_Total_Biomass)) %>%
  ungroup() %>%
  arrange(desc(Group_Total_Biomass), Migration_system2) %>%
  dplyr::select(-Total_Biomass, -Group_Total_Biomass) %>%
  unite("RowID", Migration_system2, Biomass_type, sep = "_", remove = FALSE) %>%
  tibble::column_to_rownames("RowID")

# Step 2: Prepare the corresponding upper bounds
Migration_table_combined_sorted_upper <- Migration_table_combined %>%
  filter(str_detect(Biomass_type, "Upper")) %>%
  mutate(
    Biomass_type = str_remove(Biomass_type, "_Upper"),  # Remove '_Upper' from Biomass_type
    RowID = paste(Migration_system2, Biomass_type, sep = "_")  # Create RowID
  ) %>%  # Create RowID
  arrange(match(RowID, rownames(Migration_table_combined_sorted))) %>%  # Align with the sorted table
  tibble::column_to_rownames("RowID")

# Step 3: Convert the table to a matrix for plotting
Migration_table <- t(as.matrix(Migration_table_combined_sorted[, -c(1, 2)]))  # Exclude 'Migration_system2' and 'Biomass_type' columns

# Step 4: Define color palette for the bars

groups <- rownames(Migration_table_after)  # Extract group names
color_mapping <- setNames(viridis::viridis(length(groups)), groups)  # Assign a unique color to each group

Migration_table[Migration_table <= 0] <- 1  # Replace zero or negative values with a small positive value

# Adjust the bottom margin
par(mar = c(20, 4, 4, 2) + 0.1)  # Increase the bottom margin (first value in `mar`)

# Get current plot size from RStudio graphics window
plot_size <- dev.size("in")

Cairo::CairoSVG("biomass_before_barplot.svg", width = plot_size[1], height = plot_size[2])

# Step 5: Adjust x-axis grouping for Biomass_before and Biomass_after
bar_centers <- barplot(
  Migration_table,
  log = "y",            # Log scale for the y-axis
  col = color_mapping[rownames(Migration_table)],  # Use color mapping for each group
  border = NA,          # Remove borders around bars
  ylim = c(1e11, 1e15),  # Set y-axis limits for better visualization
  #legend.text = TRUE,   # Add legend to the chart
  args.legend = list(x = "topright", title = "Organisation Level 2"),  # Position and title of the legend
  cex.axis = 0.8,  # Scale down the axis label text
  #xlab = "Migration System and Biomass Type",  # X-axis label
  ylab = "Biomass (g)", # Y-axis label
  las = 2                       # Rotate x-axis labels to vertical
)

# Step 5: Add error bars for each bar
# Total biomass and upper bounds must align with Migration_table order
total_biomass <- rowSums(Migration_table_combined_sorted[, -c(1, 2)], na.rm = TRUE)
total_biomass_upper <- rowSums(Migration_table_combined_sorted_upper[, -c(1, 2)], na.rm = TRUE)

arrows(
  x0 = bar_centers,               # X positions of bars
  y0 = total_biomass,             # Lower ends of error bars
  x1 = bar_centers,               # Same as x0
  y1 = total_biomass_upper,       # Upper ends of error bars
  angle = 90,                     # Right angle for arrowheads
  code = 2,                       # Both ends have arrowheads
  length = 0.01,                   # Length of arrowheads
  col = "black"                   # Color of error bars
)

# Close the SVG file
dev.off()

# Figure 4.2: percentage change per migratory system ---------

Migration_map_change <- Migration_map_summary %>%
  group_by(Migration_system2) %>%
  dplyr::summarise(
    Biomass_before = sum(Biomass_before),
    Biomass_after = sum(Biomass_after),
    Biomass_before_Upper = sum(Biomass_before_Upper),
    Biomass_after_Upper = sum(Biomass_after_Upper),
    Biomass_before_Lower = sum(Biomass_before_Lower),
    Biomass_after_Lower = sum(Biomass_after_Lower)
  ) %>%
  rowwise() %>%
  mutate(
    Biomass_change_log_abs = log(abs((Biomass_after - Biomass_before))) * ((Biomass_after - Biomass_before)/abs((Biomass_after - Biomass_before))),
    Biomass_change_sign = ifelse(Biomass_after > Biomass_before, "Positive", "Negative"),
    Biomass_change = (Biomass_after - Biomass_before) / Biomass_before,
    Biomass_change_Upper = (Biomass_after_Upper - Biomass_before_Upper) / Biomass_before_Upper,
    Biomass_change_Lower = (Biomass_after_Lower - Biomass_before_Lower) / Biomass_before_Lower,
    Biomass_change_CI = abs(max(Biomass_change-Biomass_change_Upper, Biomass_change-Biomass_change_Lower))
  ) %>%
  mutate(
    Biomass_change_CI = ifelse(Biomass_change>0, Biomass_change_CI, -Biomass_change_CI)
  )

plot5b <- ggplot(Migration_map_change, aes(x = reorder(Migration_system2, Biomass_change), 
                                 y = Biomass_change, 
                                 fill = Biomass_change_sign)) +
  geom_hline(yintercept=1, linetype="dashed", alpha=0.4)+
  geom_hline(yintercept=0, alpha=0.4)+
  geom_hline(yintercept=-1, linetype="dashed", alpha=0.4)+
  geom_bar(Migration_map_change, mapping=aes(x = reorder(Migration_system2, Biomass_change), 
                                             y = Biomass_change, 
                                             fill = Biomass_change_sign), 
           stat = "identity", alpha=0.8) +
  geom_errorbar(Migration_map_change, mapping=aes(x = reorder(Migration_system2, Biomass_change), 
                                                  ymax = Biomass_change+Biomass_change_CI, ymin = Biomass_change,
                                                  fill = Biomass_change_sign),
                stat = "identity", alpha=0.8, width=0.001) +
  scale_fill_manual(values = c("red3", "green4")) +  # Color scheme: red for negative, green for positive
  scale_y_break(c(1,1.1), scales = c(0.09, 2)) +
  labs(y = "") +
  theme_classic() +  # Apply theme_bw
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )

ggsave("Plot5b.pdf", plot = plot5b, dpi = 300)

#