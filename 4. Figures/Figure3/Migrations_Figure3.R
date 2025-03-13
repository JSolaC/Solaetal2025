
# ============================================================================
# Supplementary R Code: Migration Distances and Trophic Shifts
# ============================================================================
# Authors: Sola et al.
# Date: February 2025
# Journal: Nature Portfolio
# ============================================================================

# Overview:
# This script analyzes migration distances and biomass changes over time, integrating 
# species-level movement data with biomass estimates to assess historical and 
# contemporary shifts in mass migrations and trophic distributions.

# Key Considerations:
# - Migration distance estimation: 
#   - Integrates multiple sources to estimate migration distances.
#   - Accounts for population proportion when estimating migration distance at 
#     the species level.
#   - Standardizes movement types (e.g., nomadic, to-and-fro, local).
# - Biomass weighting:
#   - Migration distances are weighted based on biomass proportions across species.
#   - Historical and current biomass estimates are included to track shifts.
#
# Outputs:
# - A structured dataset containing weighted migration distances and biomass changes.
# - Processed trends in species movement distances before and after the reference period.
# - Visual representations of migration shifts across different taxonomic groups.
# - Exported figures for publication and supplementary materials.

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

# 1- Load datasets ---------

Species_list <- read.csv('1-SpeciesList/Output/Species_list_all_out.csv') %>%
  dplyr::select(-X)

Mass_migrations_manual <- read.csv("2-SpeciesBiomass/Outputs/Mass_Migrations_data.csv") 

Distance_pops <- read.csv("4-Figures/Figure3/Fish_pop_migrations_biomass.csv") %>%
  
  # standardise names
  left_join(Species_list, by="Binomial") %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial)) %>%
  # select variables and remove duplicates
  group_by(Binomial, Region) %>%
  dplyr::summarise(Biomass_proportion=sum(Biomass_proportion))

Mass_migrations_DietTrophic <- read.csv("4-Figures/Figure3/MassMigrations_DietTrophic.csv")

#
# FIGURE 3C (run after running code for Figure 1 Steps 1 to 4) ####

# A - Obtain Migration distance data per species
Migration_distance_data_unfiltered <- Mass_migrations_manual %>%
  
  dplyr::select(Binomial, Taxa.group, Region, Population.proportion, Movement_Type, "Migration.distance1"=Migration.distance..km., "Migration.distance2"=Migration.distance2..km2.) %>%
  left_join(Distance_pops) %>%
  mutate(Population.proportion2 = ifelse(!is.na(Biomass_proportion), Biomass_proportion, Population.proportion)) %>%
  mutate(Population.proportion2 = ifelse(is.na(Population.proportion2),1,Population.proportion2)) %>%
  distinct(Binomial, Taxa.group, Population.proportion2, Movement_Type, Migration.distance1, Migration.distance2) %>%
  
  # 1.1 summarise information per movement type and binomial
  group_by(Binomial, Population.proportion2, Movement_Type) %>%
  dplyr::summarise(
    Migration.distance1=mean(as.numeric(Migration.distance1), na.rm=T),
    Migration.distance2=mean(as.numeric(Migration.distance2), na.rm=T)
  ) %>%
  
  # 1.2 for area-based distances, find linear distance as the diameter of the area
  mutate(Migration.distance2 = sqrt(Migration.distance2/pi) * 2) %>%
  mutate(Migration.distance = ifelse(!is.na(Migration.distance1), Migration.distance1, Migration.distance2)) #%>% # replace NAs for area-based estimates


Migration_distance_data <- Migration_distance_data_unfiltered %>%
  # 1.3 remove rows that contain NAs or are DVM movement
  filter(
    !is.na(Binomial) &
      !is.na(Movement_Type) &
      Movement_Type!="DVM" & # remove DVM
      Movement_Type != "" & # remove rows not indicating movement type
      (Migration.distance1!="NaN" | Migration.distance2!="NaN")
  ) %>%
  
  # 1.4 join with historical biomass dataset and filter for rows containing current and historical biomass
  left_join(Historical_biomass_data %>% rename(Binomial=Species)) %>%
  filter(!is.na(Biomass_before) & !is.na(Biomass_after)) %>%
  
  mutate(
    Biomass_before = Biomass_before * Population.proportion2,
    Biomass_after = Biomass_after * Population.proportion2,
    Biomass_before_Upper = Biomass_before_Upper * Population.proportion2,
    Biomass_after_Upper = Biomass_after_Upper * Population.proportion2,
    Biomass_before_Lower = Biomass_before_Lower * Population.proportion2,
    Biomass_after_Lower = Biomass_after_Lower * Population.proportion2
  ) %>%
  
  # 1.5 Combine data for all movement types
  bind_rows(
    mutate(., Movement_Type = "all")
  ) %>%
  mutate(Movement_Type = factor(Movement_Type, levels=c("all", "nomadic", "to-and-fro", "local")))

Migration_distance_spp <- Migration_distance_data %>%
  filter(Movement_Type=="all") %>%
  group_by(Organisation_level2) %>%
  summarise(Species_n = n())

# B - Compute Migration distance metrics per taxa group for the historical and current time points
Migration_distance_summary <- Migration_distance_data %>%
  
  # 2.1 Group by organismal group and movement type
  group_by(Organisation_level2, Movement_Type) %>%
  
  # 2.2 Calculate weighted migration distances
  mutate(across(
    c(Biomass_before, Biomass_after, Biomass_before_Upper, Biomass_after_Upper, Biomass_before_Lower, Biomass_after_Lower),
    ~ Migration.distance * (.x / sum(.x, na.rm = TRUE)),
    .names = "Migration_distance_weighted_{.col}"
  )) %>%
  
  # 2.4 Pivot longer to handle "before" and "after" columns dynamically
  pivot_longer(
    cols = contains("_before") | contains("_after"),
    names_to = c(".value", "Time", "Suffix"),
    names_pattern = "^(.*)_(before|after)(_Upper|_Lower)?$"
  ) %>%
  
  group_by(Binomial) %>%
  mutate(
    Suffix = ifelse(Suffix=="_Upper", "CI_Upper", 
                    ifelse(Suffix=="_Lower", "CI_Lower", "Biomass_weighted")),
    Year=ifelse(Time=="before", min(Year, na.rm=T), max(Year, na.rm = T))
  ) %>%
  
  dplyr::select(-Biomass) %>%
  
  # 2.5 Separate Upper and Mean into their respective columns
  pivot_wider(
    names_from = Suffix,
    values_from = Migration_distance_weighted_Biomass
  ) %>%
  
  # 2.6 Reformat and reorder the Time column
  mutate(Time = factor(Time, levels = c("before", "after"), labels = c("Before", "After"))) %>%
  # 2.7 There's an issue with CI, which is that when considering larger biomass, that does not imply that weighted distance is larger - as weighted distance depends on the relative biomass which may result in lower overall migration distance
  
  group_by(Organisation_level2, Time, Movement_Type) %>%
  dplyr::summarise(
    Biomass_weighted = sum(Biomass_weighted),
    CI_Upper = sum(CI_Upper),
    CI_Lower = sum(CI_Lower)
  ) %>%
  mutate(Weighted_CI = abs(CI_Upper-CI_Lower)/2+Biomass_weighted) %>%
  
  group_by(Organisation_level2) %>%
  mutate(Biomass_weighted_all = mean(Biomass_weighted[Movement_Type=="all"])) %>%
  mutate(Movement_Type = factor(Movement_Type, levels=c("all", "nomadic", "to-and-fro", "local")))

# C - Plotting
plot3c <- ggplot() +
  geom_bar(Migration_distance_summary, 
           mapping=aes(x=reorder(paste(Organisation_level2),-as.numeric(Biomass_weighted_all), mean), y=Biomass_weighted, alpha=Time, fill=Organisation_level2),
           stat="identity", position = position_dodge2(preserve = "single", padding = 0), col="black", size = 0.2) +
  geom_errorbar(Migration_distance_summary, 
                mapping=aes(x=reorder(paste(Organisation_level2),-as.numeric(Biomass_weighted_all), mean), 
                            ymax=Weighted_CI, ymin=Biomass_weighted, alpha=Time),
                col="black",
                position = position_dodge(width = 0.9),  # Ensure error bars align with bars
                width = 0.05) +
  geom_point(Migration_distance_data, mapping=aes(x=Organisation_level2, y=Migration.distance1/1.5, group=Organisation_level2), 
             alpha=0.1, position=position_jitter(h=0.15,w=0.15)) +
  ylim(0,525)+
  scale_y_continuous(sec.axis = sec_axis(~ .*1.5, "Migration distance (km)")) +
  facet_wrap(~Movement_Type, nrow = 1, drop=FALSE)+
  scale_alpha_manual(values=c(1,0.75,0.5))+
  scale_fill_manual(values = color_mapping) +    # Apply viridis colors to fill
  scale_color_manual(values = color_mapping) +   # Apply viridis colors to color
  #scale_fill_manual(values=c("purple","blue4","red3","green4")) +
  #scale_colour_manual(values=c("purple","blue4","red3","green4")) +
  ylab("Weighted migration distance (km)") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust = 1),
    #axis.text.y = element_text(angle=90, vjust = 1),
    legend.position = "none"
  )

ggsave("Plot3c.pdf", plot = plot3c, dpi = 300)

#
# FIGURE 3A -----------

Migration_distance_change <- Migration_distance_summary %>%
  
  # remove all
  filter(Movement_Type=="all") %>%
  
  group_by(Organisation_level2) %>%
  
  dplyr::summarise(
    Migration_change = (Biomass_weighted[Time=="After"] - Biomass_weighted[Time=="Before"]) / Biomass_weighted[Time=="Before"],
    Migration_change_Upper = (CI_Upper[Time=="After"] - CI_Upper[Time=="Before"]) / CI_Upper[Time=="Before"],
    Migration_change_Lower = (CI_Lower[Time=="After"] - CI_Lower[Time=="Before"]) / CI_Lower[Time=="Before"]
  ) %>%
  
  mutate(
    Sign = ifelse(Migration_change>0, "Positive", "Negative"),
    Migration_change_CI = ifelse(Sign=="Positive",
                                 (abs(Migration_change_Upper-Migration_change_Lower)/2)+Migration_change,
                                 -(abs(Migration_change_Upper-Migration_change_Lower)/2)+Migration_change)
  )

plot3a <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", alpha=0.4)+
  geom_hline(yintercept=0, alpha=0.4)+
  geom_hline(yintercept=-1, linetype="dashed", alpha=0.4)+
  geom_bar(Migration_distance_change
           , mapping=aes(x = reorder(Organisation_level2, Migration_change), y = Migration_change, fill=Sign), stat = "identity", alpha=0.8) +
  geom_errorbar(Migration_distance_change
                , mapping=aes(x = reorder(Organisation_level2, Migration_change), ymax = Migration_change_CI, ymin = Migration_change, fill=Sign), stat = "identity", alpha=0.8, width=0.001) +
  scale_fill_manual(values = c("Positive" = "green4", "Negative" = "red3")) +
  #scale_y_break(c(1, 1.1), scales = c(0.09, 0.01, 2)) +
  labs(x = "Species", y = "Biomass Change") +
  theme_classic() +
  theme(
    legend.position = "none",
    #axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Plot3a.pdf", plot = plot3a, dpi = 300)

#
# FIGURE 3 (discarded) -----------

Migration_distance_change_Mov <- Migration_distance_summary %>%
  
  # remove all
  #filter(Movement_Type=="all") %>%
  ungroup() %>%
  group_by(Movement_Type, Time) %>%
  
  dplyr::summarise(
    Biomass_weighted = sum(Biomass_weighted),
    CI_Upper = sum(CI_Upper),
    CI_Lower = sum(CI_Lower)
  ) %>%
  
  dplyr::summarise(
    Migration_change = (Biomass_weighted[Time=="After"] - Biomass_weighted[Time=="Before"]) / Biomass_weighted[Time=="Before"],
    Migration_change_Upper = (CI_Upper[Time=="After"] - CI_Upper[Time=="Before"]) / CI_Upper[Time=="Before"],
    Migration_change_Lower = (CI_Lower[Time=="After"] - CI_Lower[Time=="Before"]) / CI_Lower[Time=="Before"]
  ) %>%
  
  mutate(
    Sign = ifelse(Migration_change>0, "Positive", "Negative"),
    Migration_change_CI = ifelse(Sign=="Positive",
                                 (abs(Migration_change_Upper-Migration_change_Lower)/2)+Migration_change,
                                 -(abs(Migration_change_Upper-Migration_change_Lower)/2)+Migration_change)
  ) %>%
  mutate(Movement_order = case_when(
    Movement_Type =="local" ~ 1,
    Movement_Type == "to-and-fro" ~ 2,
    Movement_Type == "nomadic" ~ 3,
    Movement_Type == "all" ~ 4
  ))

plot3b <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", alpha=0.4)+
  geom_hline(yintercept=0, alpha=0.4)+
  geom_hline(yintercept=-1, linetype="dashed", alpha=0.4)+
  geom_bar(Migration_distance_change_Mov
           , mapping=aes(x = reorder(Movement_Type, Movement_order), y = Migration_change, fill=Sign), stat = "identity", alpha=0.8) +
  geom_errorbar(Migration_distance_change_Mov
                , mapping=aes(x = reorder(Movement_Type, Movement_order), ymax = Migration_change_CI, ymin = Migration_change, fill=Sign), stat = "identity", alpha=0.8, width=0.001) +
  scale_fill_manual(values = c("Positive" = "green4", "Negative" = "red3")) +
  #scale_y_break(c(1, 1.1), scales = c(0.09, 0.01, 2)) +
  labs(x = "Species", y = "Biomass Change") +
  theme_classic() +
  theme(
    legend.position = "none"
    #axis.text.y = element_text(size = 6),
    #axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Plot3b.pdf", plot = plot3b, dpi = 300)

#
# FIGURE 3D ------------

Migration_trophic_diet <- Mass_migrations_DietTrophic %>%
  
  mutate(Binomial = str_trim(word(Binomial, 1,2, sep=" "))) %>%
  left_join(Historical_biomass_data %>% rename(Binomial=Species)) %>%
  
  # need to remove, otherwise returns rows with 0s
  filter(!is.na(Biomass_before) & !is.na(Biomass_after)) %>%
  filter(Organisation_level2 != "Freshwater fish" & Organisation_level2!="Terrestrial invertebrates") %>%
  
  # 2.4 Pivot longer to handle "before" and "after" columns dynamically
  pivot_longer(
    cols = contains("Biomass_"),
    names_to = c(".value", "Time", "Suffix"),
    names_pattern = "^(.*)_(before|after)(_Upper|_Lower)?$"
  ) %>%
  mutate(Suffix = ifelse(Suffix=="_Upper", "CI_Upper", 
                         ifelse(Suffix=="_Lower", "CI_Lower", "Biomass_mean"))) %>%
  
  # 2.5 Separate Upper and Mean into their respective columns
  pivot_wider(
    names_from = Suffix,
    values_from = Biomass
  ) %>% 
  
  pivot_longer(
    cols=Herbivore:Omnivore,
    names_to="Trophic_group3",
    values_to="Presence"
  ) %>%
  
  filter(Presence=="x") %>%
  
  mutate(
    Trophic_group = ifelse(Trophic_group3=="Primary_Consumer" | Trophic_group3=="Secondary_Consumer" | Trophic_group3=="Tertiary_Consumer",
                           "Consumer", Trophic_group3),
    Time=factor(Time,levels=c("before","after"))
  )

Migration_trophic_diet_group <- Migration_trophic_diet %>%
  
  group_by(Organisation_level2, Time, Trophic_group) %>%
  dplyr::summarise(
    Species_n = n(),
    Biomass_group = sum(Biomass_mean),
    Biomass_Upper = sum(CI_Upper),
    Biomass_Lower = sum(CI_Lower)
  )

Migration_trophic_diet2 <- Migration_trophic_diet %>%
  group_by(Trophic_group) %>%
  summarise(Species_n = n())

plot3e <- ggplot(Migration_trophic_diet_group) +
  geom_bar(
    mapping=aes(x=reorder(paste(Organisation_level2),-as.numeric(Biomass_group), mean), 
                y=Biomass_group, alpha=Time, fill=Organisation_level2),
    stat="identity", position = position_dodge2(preserve = "single", padding = 0.2), col="black", size = 0.2) +
  geom_errorbar(
    mapping=aes(x=reorder(paste(Organisation_level2),-as.numeric(Biomass_group), mean), 
                ymax=Biomass_Upper, ymin=Biomass_group, alpha=Time),
    col="black",
    position = position_dodge(width = 0.9),  # Ensure error bars align with bars
    width = 0.05) +
  facet_wrap(~ Trophic_group, nrow = 1)+
  scale_y_log10(limits=c(5000000000, 1000000000000000), oob = scales::squish) +
  scale_alpha_manual(values=c(1,0.75,0.5))+
  scale_fill_manual(values = color_mapping) +  # Apply viridis colors to fill
  scale_color_manual(values = color_mapping) +  # Apply viridis colors to error bar colors
  ylab("Biomass (g)") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust = 1),
    axis.text.y = element_text(angle=90),
    legend.position = "none"
  )

ggsave("Plot3e.pdf", plot = plot3e, dpi = 300)

#
# FIGURE 3B --------------

Migration_trophic_change <-  Migration_trophic_diet_group %>%
  
  group_by(Trophic_group, Time) %>%
  dplyr::summarise(
    Biomass_group = sum(Biomass_group),
    Biomass_Upper = sum(Biomass_Upper),
    Biomass_Lower = sum(Biomass_Lower)
  ) %>%
  
  group_by(Trophic_group) %>%
  dplyr::summarise(
    Biomass_change = (Biomass_group[Time=="after"] - Biomass_group[Time=="before"]) / Biomass_group[Time=="before"],
    Biomass_change_Upper = (Biomass_Upper[Time=="after"] - Biomass_Upper[Time=="before"]) / Biomass_Upper[Time=="before"],
    Biomass_change_Lower = (Biomass_Lower[Time=="after"] - Biomass_Lower[Time=="before"]) / Biomass_Lower[Time=="before"]
  ) %>%
  
  mutate(
    Sign = ifelse(Biomass_change>0, "Positive", "Negative"),
    Biomass_change_CI = ifelse(Sign=="Positive",
                               (abs(Biomass_change_Upper-Biomass_change_Lower)/2)+Biomass_change,
                               -(abs(Biomass_change_Upper-Biomass_change_Lower)/2)+Biomass_change)
  )


plot3d <- ggplot() +
  geom_hline(yintercept=1, linetype="dashed", alpha=0.4)+
  geom_hline(yintercept=0, alpha=0.4)+
  geom_hline(yintercept=-1, linetype="dashed", alpha=0.4)+
  geom_bar(Migration_trophic_change
           , mapping=aes(x = reorder(Trophic_group, Biomass_change), y = Biomass_change, fill=Sign), stat = "identity", alpha=0.8) +
  geom_errorbar(Migration_trophic_change
                , mapping=aes(x = reorder(Trophic_group, Biomass_change), ymax = Biomass_change_CI, ymin = Biomass_change, fill=Sign), stat = "identity", alpha=0.8, width=0.001) +
  scale_fill_manual(values = c("Positive" = "green4", "Negative" = "red3")) +
  #scale_y_break(c(1, 1.1), scales = c(0.09, 0.01, 2)) +
  labs(x = "Species", y = "Biomass Change") +
  theme_classic() +
  theme(
    legend.position="none",
    axis.title.x = element_blank(),
    #axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("Plot3d.pdf", plot = plot3d, dpi = 300)

#