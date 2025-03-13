# =======================================================================
# Supplementary R Code: Global Biomass and Migration Data Processing
# =======================================================================
# Authors: Sola et al.
# Date: February 2025
# Journal: Nature Portfolio
# =======================================================================
# Overview:
# This script integrates species taxonomy, biomass estimates, and human migration data 
# to create a unified dataset for global biomass and species movement analyses.
#
# Key Considerations:
# - Species biomass is compiled from multiple sources, including marine, terrestrial, 
#   and extinct species datasets.
# - Migratory species require synonym resolution to account for naming inconsistencies (using the list developed in the previous script).
# - Human migration data is incorporated from international and domestic travel statistics, as well as estimates on the number of nomadic individuals.
# - Mass migrations of various taxa are identified and processed separately.
#
# Outputs:
# - A final dataset containing biomass estimates, taxonomic classifications, and 
#   migration status for all species.
# - A human migration dataset incorporating international and domestic movements.
# - A processed list of species involved in mass migrations, categorised by taxa.
# =======================================================================

# 0- Set up Workspace -------------------

# Define the required packages
packages <- c("tidyverse", "brms", "countrycode", "stringi")

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
# 1- Load Initial datasets ####

Species_list <- read.csv('1-SpeciesList/Output/Species_list_all_out.csv') %>%
  dplyr::select(-X)

Species_taxonomy <- bind_rows(
  # migratory species
  read.csv('1-SpeciesList/Output/Species_list_taxonomy_migratory.csv')  %>%
  dplyr::select(-X)
  , 
  # non-migratory species
  read.csv('1-SpeciesList/Output/NonMigratory_species_list_taxa.csv')
  ) %>%
  dplyr::select(-X, -FID) %>%
  distinct(Binomial, .keep_all = TRUE)

#
# 2- Load Biomass datasets ####

global_tmammals <- read.csv("2-SpeciesBiomass/Datasets_Biomass/theglobalbiomassofwildanimals.csv") %>%
  mutate(Binomial = str_trim(binomial)) %>%
  dplyr::select(Binomial, 
                "Biomass"=biomass_g, "Biomass_97"=biomass_g_97pt5, "Biomass_2"=biomass_g_2pt5, 
                "Abundance"=estimated_population, "Abundance_97"=estimated_population_97pt5, "Abundance_2"=estimated_population_2pt5) %>%
  left_join(Species_list) %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial))

global_marine <-  read.csv("2-SpeciesBiomass/Datasets_Biomass/wild_marine_biomass.csv") %>%
  mutate(Binomial = str_trim(Binomial)) %>%
  dplyr::select(Binomial, "Biomass"=biomass_g)%>%
  left_join(Species_list) %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial))

global_turtles <- read.csv('2-SpeciesBiomass/Datasets_Biomass/global_turtles.csv') %>%
  mutate(Binomial = str_trim(word(Binomial, 1,2, sep=" "))) %>%
  left_join(Species_list) %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial)) %>%
  dplyr::select(Binomial, Biomass, Abundance, Body_Mass)

global_extinct <- read.csv('2-SpeciesBiomass/Datasets_Biomass/global_EX_species.csv') %>%
  mutate(Binomial = str_trim(word(Binomial, 1,2, sep=" "))) %>%
  #mutate(Migratory="Y") %>% # regardless, the ones listed are all migratory
  left_join(Species_list) %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial)) %>%
  dplyr::select(Binomial, "Organisation_level_extinct"=Organisation_level2, Biomass, Abundance, Body_Mass, Migratory)

global_fish <- read.csv("2-SpeciesBiomass/Datasets_Biomass/global_fish_data.csv") %>%
  dplyr::select(-X) %>%
  mutate(Binomial = str_trim(Binomial)) %>%
  left_join(Species_list) %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial))

global_marine_inv <- read.csv('2-SpeciesBiomass/Datasets_Biomass/global_marine_inv.csv') %>%
  mutate(Binomial = str_trim(word(Binomial, 1,2, sep=" "))) %>%
  left_join(Species_list) %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial)) %>%
  dplyr::select(Binomial, "Biomass_inv"=Biomass, "Abundance_inv"=Abundance, Body_Mass, Migratory) %>%
  left_join(global_fish %>% dplyr::select(Binomial, Biomass, Abundance)) %>%
  mutate_cond(Binomial=="Callinectes sapidus", Biomass=Body_Mass*Abundance)

global_birds_traits <-  read.csv("2-SpeciesBiomass/Datasets_Biomass/AVONET_data.csv") %>%
  mutate(Binomial = str_trim(Species1)) %>%
  dplyr::select(Binomial, "Body_Mass"=Mass) %>%
  left_join(Species_list) %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial))

global_birds <-  read.csv("2-SpeciesBiomass/Datasets_Biomass/Bird_PNAS_dataset.csv") %>%
  mutate(Binomial = str_trim(Scientific.name)) %>%
  dplyr::select(Binomial, "Abundance"=Abundance.estimate, "Abundance_97"=X95..Upper.CI, "Abundance_2"=X95..Lower.CI) %>%
  left_join(Species_list) %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial)) %>%
  left_join(global_birds_traits %>% dplyr::select(Binomial, Body_Mass)) %>%
  mutate(
    Biomass=ifelse(!is.na(Body_Mass), Abundance*as.numeric(Body_Mass), NA),
    Biomass_97=ifelse(!is.na(Body_Mass), Abundance_97*as.numeric(Body_Mass), NA),
    Biomass_2=ifelse(!is.na(Body_Mass), Abundance_2*as.numeric(Body_Mass), NA)
  ) %>%
  filter(!is.na(Body_Mass))

# OBTAIN GLOBAL DATASET

global_all_spp <- bind_rows(
  data.frame(global_tmammals, Tax_group = "Terrestrial mammal"), 
  data.frame(global_marine, Tax_group = "Marine mammal"), 
  data.frame(global_fish, Tax_group = "Fish") %>% filter(Binomial %ni% global_marine_inv$Binomial), 
  data.frame(global_birds, Tax_group = "Birds"),
  data.frame(global_turtles, Tax_group = "Seaturtles"),
  data.frame(global_marine_inv, Tax_group = "Marine invertebrates"),
  data.frame(global_extinct, Tax_group = "Extinct species")
  ) %>%
  dplyr::select(-Original) %>%
  left_join(Species_taxonomy) %>%
  mutate(
    Migratory = ifelse(Tax_group!="Marine invertebrates" & Tax_group!="Terrestrial invertebrates" & Tax_group!="Extinct species", ifelse(!is.na(Migration2), "Y", "N"), Migratory),
    Biomass = ifelse(is.na(Biomass_97) & Biomass==0, NA, Biomass)
    ) %>%
  mutate(Organisation_level2 = ifelse(is.na(Organisation_level2), Tax_group, Organisation_level2)) %>%
  mutate(Organisation_level2 = ifelse(Tax_group=="Extinct species", Organisation_level_extinct, Organisation_level2)) %>%
  mutate(
    Organisation_level2 = ifelse(Organisation_level2=="Terrestrial invertebrates" & Tax_group=="Marine invertebrates", "Marine invertebrates", Organisation_level2),
    Organisation_level3 = ifelse(Organisation_level3=="Terrestrial invertebrates" & Tax_group=="Marine invertebrates", "Marine invertebrates", Organisation_level3)
    ) %>%
  # Remove repeated instances
  distinct(Binomial, Biomass, .keep_all=TRUE) %>%   
  # For species with the same name (due to the use of synonyms), sum overall biomass
  group_by(Binomial) %>%
  summarise(
    Biomass = sum(Biomass, na.rm = TRUE),  # Summing biomass
    across(everything(), first, .names = "{.col}")  # Keeping the first occurrence of all other columns
  ) %>%
  ungroup()

write.csv(global_all_spp, "2-SpeciesBiomass/Outputs/Biomass_all_spp.csv")

#
# 3- Load human datasets ####

# LOAD DATASET ON NATIONAL TRIPS AND SUMMARISE DATA

# OECD data obtained from: https://www.oecd-ilibrary.org/economics/data/oecd-tourism-statistics/domestic-tourism-edition-2019_0f8b05d0-en (accessed on  24th October 2024)
National_trips_OECD <- read.csv("2-SpeciesBiomass/Datasets_Homosapiens/NationalTourismData.csv") %>%
  mutate_cond(Country=="China (People's Republic of)", Country="China") %>%
  mutate_cond(Country=="United States" , Country="United States Of America") %>%
  filter(Variable=="Total domestic trips") %>%
  group_by(Country, Variable) %>%
  dplyr::mutate(LP_value = log(as.numeric(Value)/as.numeric(Value[Year==min(Year)])))

# UNWTO data obtained from: https://www.unwto.org/tourism-statistics/key-tourism-statistics (accessed on 28th October 2024)
National_trips_UNWTO <- read.csv("2-SpeciesBiomass/Datasets_Homosapiens/UNTourism_DomesticTourism_Trips.csv") %>%
  filter(Trip_Type=="Total trips") %>%
  pivot_longer(
    cols=X1995:X2022,
    names_to="Year",
    values_to="Value"
  ) %>%
  filter(Value!="") %>%
  mutate(
    Country = str_to_title(Country),
    Year=as.numeric(gsub("X","",Year)),
    Value=as.numeric(gsub(',','',Value))*1000
  ) %>%
  dplyr::select(-Trip_Type, -Units)

# When joining datasets, keep UNWTO data since it is more complete and extensive and add countries for which there is no data
National_trips_all <- National_trips_UNWTO %>% ungroup() %>% dplyr::select(Country, Year, Value) %>% mutate(Dataset="UNWTO") %>%
  # add only countries that are not already in the UNWTO dataset
  bind_rows(
  National_trips_OECD %>% ungroup() %>% dplyr::select(Country, Year, Value) %>% mutate(Dataset="OECD") %>%
    filter(Country %ni% unique(National_trips_UNWTO$Country))
  ) %>%
  group_by(Country, Dataset) %>%
  dplyr::mutate(LP_value = log(as.numeric(Value)/as.numeric(Value[Year==min(Year)])))

ggplot() +
  geom_point(National_trips_all, mapping=aes(x=Year, y=LP_value, col=Country, size=Value), alpha=0.2) +
  geom_line(National_trips_all, mapping=aes(x=Year, y=LP_value, col=Country), alpha=0.4) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() + 
  theme(legend.position = "none")

# 
National_trips_all_summary <- National_trips_all %>%
  group_by(Country, Dataset) %>%
  filter(Year==min(Year) | Year==max(Year))

# GET TEMPORAL DATA FOR INTERNATIONAL TRAVELLING

# UNWTO data obtained from: https://www.unwto.org/tourism-statistics/key-tourism-statistics (accessed on 28th October 2024)
International_outbound_trips_UNWTO <- read.csv("2-SpeciesBiomass/Datasets_Homosapiens/UNTourism_OutboundTourism_Transport.csv") %>%
  filter(Transport_mode=="Total departures") %>%
  pivot_longer(
    cols=X1995:X2022,
    names_to="Year",
    values_to="Value"
  ) %>%
  filter(Value!="") %>%
  mutate(
    Country = str_to_title(Country),
    Year=as.numeric(gsub("X","",Year)),
    Value=as.numeric(gsub(',','',Value))*1000
  ) %>%
  dplyr::select(-Transport_mode, -Units) %>%
  group_by(Country) %>%
  dplyr::mutate(LP_value = log(as.numeric(Value)/as.numeric(Value[Year==min(Year)])))

ggplot() +
  geom_point(International_outbound_trips_UNWTO, mapping=aes(x=Year, y=LP_value, col=Country, size=Value), alpha=0.2) +
  geom_line(International_outbound_trips_UNWTO, mapping=aes(x=Year, y=LP_value, col=Country), alpha=0.4) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() + 
  theme(legend.position = "none")

International_outbound_trips_all_summary <- International_outbound_trips_UNWTO %>%
  group_by(Country) %>%
  filter(Year==min(Year) | Year==max(Year))

# GET TEMPORAL DATA FOR NOMADIC PEOPLE

Nomadic_Manual <- read.csv("2-SpeciesBiomass/Datasets_Homosapiens/Nomadic_Numbers_Manual.csv") %>%
  dplyr::select("Country"=Origin, "Value"=5, 6) %>%
  mutate(Value = Value*1000000) %>%
  mutate(Year=2019)

# SUMMARISE ALL DATA (redo with new datasets)

Human_migration_data <- rbind(
  International_outbound_trips_UNWTO %>% mutate(Migration_Type="to_and_fro"),
  National_trips_all %>% ungroup() %>% dplyr::select(-Dataset) %>% mutate(Migration_Type="local"),
  Nomadic_Manual %>% dplyr::select(-Body.weight) %>% mutate(Migration_Type="nomadic")
  ) %>%
  filter(!is.na(Value)) %>%
  # to convert species characters, first we need to do this
  mutate(Country = stri_trans_general(Country, "Latin-ASCII")) %>%
  dplyr::mutate(
    Country = case_when(
      Country == "Cote D�Ivoire" ~ "Côte d'Ivoire",
      Country == "Cura�Ao" ~ "Curaçao",
      Country == "Dijbouti" ~ "Djibouti",
      Country == "Kosovo" ~ "Croatia",  # countrycode recognizes "Kosovo" directly
      Country == "Lao People�S Democratic Republic" ~ "Laos",
      Country == "Namibia and Angola" ~ "Namibia",  # Choose one; or handle separately if needed
      Country == "Saba" ~ "Cuba",  # Caribbean island
      Country == "SADR" ~ "Western Sahara",  # SADR (Sahrawi Arab Democratic Republic)
      Country == "Sint Eustatius" ~ "Cuba",  # Caribbean island
      Country == "Bonaire" ~ "Cuba",  # Caribbean island
      Country == "Tasmania" ~ "Australia",  # Tasmania is a part of Australia
      Country == "Tibet" ~ "China",  # Tibet is part of China
      Country == "CAnada" ~ "Canada",
      Country == "USA" ~ "United States",
      Country == "United States Of America" ~ "United States",
      Country == "T�Rkiye" ~ "Turkey",
      TRUE ~ Country  # Keep the original name if it doesn't need replacement
    )
  ) %>%
  mutate(Continent = countrycode(Country, origin = "country.name", destination = "continent")) %>%
  mutate(
    Continent = case_when(
      Country %in% c("Australia", "New Zealand") ~ "Australia",
      Country %in% c("Greenland") ~ "Europe",
      Country %in% c("Canada", "United States", "Puerto Rico", "Mexico") ~ "North America",
      Country %in% c("Grenada", "Turks And Caicos Islands", "Saint Lucia", "Guadeloupe", "Saint Vincent And The Grenadines", "Haiti", "Jamaica", "Martinique", "Montserrat", "Saint Kitts And Nevis", "United States Virgin Islands", "British Virgin Islands", "Cayman Islands", "Cuba", "Curaçao", "Dominica", "Dominican Republic", "Anguilla", "Antigua And Barbuda", "Aruba","Bahamas" , "Barbados", "Bermuda", "Panama", "Costa Rica", "Guatemala", "El Salvador", "Honduras", "Nicaragua", "Belize", "Trinidad And Tobago", "Sint Maarten (Dutch Part)") ~ "Central America",
      Country %in% c("French Guiana", "Venezuela, Bolivarian Republic Of", "Bolivia, Plurinational State Of", "Brazil", "Argentina", "Colombia", "Chile", "Peru", "Venezuela", "Ecuador", "Bolivia", "Paraguay", "Uruguay", "Guyana", "Suriname") ~ "South America",
      TRUE ~ Continent  # For countries outside the Americas
    )
  ) %>%
  mutate(
    Body_weight = case_when(
      Continent == "Europe" ~ 70.8,
      Continent == "Asia" ~ 57.7,
      Continent == "North America" ~ 80.7,
      Continent == "Africa" ~ 60.7,
      Continent == "Central America" ~ 67.9,
      Continent == "South America" ~ 67.9,
      Continent == "Australia" ~ 74.1,
      Continent == "Oceania" ~ 57.7,
      TRUE ~ NA_real_  # For any unmatched continent
    )
  ) %>%
  mutate(Biomass = Value*Body_weight*1000)
  
Human_migration_data_continent <- Human_migration_data %>%
  group_by(Continent, Year) %>%
  summarise(
    Abundance = sum(Value),
    Biomass = sum(Biomass)
    ) %>%
  # this information comes from Table 2 in Lucassen, J., & Lucassen, L. (2014). Measuring and quantifying cross-cultural migrations: an introduction. In Globalising Migration History (pp. 1-54). Brill.
  bind_rows(
    data.frame(
      "Continent" = c("Europe", "Europe"),
      #"Migration_Type" = c("to_and_fro", "to_and_fro"),
      "Year" = c(1500, 1850),
      "Biomass" = c(9800000 * 70.8, 102500000 * 70.8)
    )
  ) %>%
  filter(Year<2020)

Human_migration_data_continent2 <- Human_migration_data %>%
  group_by(Country, Migration_Type, Year) %>%
  summarise(
    Abundance = sum(Value),
    Biomass = sum(Biomass)
  ) %>%
  group_by(Country, Migration_Type) %>%
  summarise(
    Abundance = Abundance[which.max(Year)],
    Biomass = Biomass[which.max(Year)]
  ) %>%
  group_by(Migration_Type) %>%
  summarise(
    Abundance = sum(Abundance),
    Biomass = sum(Biomass)
  ) %>%
  ungroup() %>%
  mutate(
    Biomass_proportional = Biomass / sum(Biomass)
  )

Human_body_weight <- Human_migration_data_continent %>%
  filter(Year==2019) %>%
  ungroup() %>%
  summarise(
    Abundance=sum(Abundance),
    Biomass=sum(Biomass)
  ) %>%
  mutate(BodyWeight=Biomass/Abundance)

ggplot(Human_migration_data_continent, aes(x=Year, y=Biomass, col=Continent)) +
  geom_line() +
  #facet_wrap(~Migration_Type) +
  scale_y_log10() +
  xlim(1990,2022) +
  theme_bw()

Human_growth_ratios <- Human_migration_data_continent %>%
  filter(Continent == "Europe", Year %in% c(1500, 1850, 2019)) %>%
  summarise(
    Ratio_1500 = Biomass[Year == 2019] / Biomass[Year == 1500],
    Ratio_1850 = Biomass[Year == 2019] / Biomass[Year == 1850]
  )

Human_biomass_historical <- Human_migration_data_continent %>%
  filter(Continent != "Europe") %>%
  filter(Year == 2019) %>% # Use 2019 biomass as the base
  mutate(
    Biomass_1500 = Biomass / Human_growth_ratios$Ratio_1500,
    Biomass_1850 = Biomass / Human_growth_ratios$Ratio_1850
  )

Human_biomass_historical2 <- Human_biomass_historical %>%
  dplyr::select(Continent, Biomass_1500, Biomass_1850) %>%
  pivot_longer(cols = starts_with("Biomass"), 
               names_to = "Year", 
               values_to = "Biomass") %>%
  mutate(
    Year = case_when(
      Year == "Biomass_1500" ~ 1500,
      Year == "Biomass_1850" ~ 1850
    )
  )

Human_migration_summary <- Human_migration_data_continent %>%
  bind_rows(Human_biomass_historical2) %>%
  group_by(Year) %>%
  summarise(Biomass=sum(Biomass))

#write.csv(Human_migration_summary, "Human_migrations.csv")

#
# 4- Find mass migrations ####

Mass_Migrations_spp <- read.csv("2-SpeciesBiomass/Outputs/Biomass_all_spp.csv") %>%
  #filter(Migratory=="Y") %>%
  filter(Organisation_level2!="Birds" & Organisation_level2!="Fish" & Organisation_level2!="Terrestrial invertebrates") %>%
  mutate(Organisation_level2 = ifelse(Organisation_level2=="Bats", "Terrestrial mammal", Organisation_level2)) %>%
  arrange(desc(Biomass)) %>%
  group_by(Organisation_level2, Migratory) %>%
  mutate(
    Position = 1:n(),
    CumSum=cumsum(Biomass),
    CumPer=cumsum(Biomass)/sum(na.omit(Biomass)),
    Richness = n()
  ) %>%
  dplyr::select(Position,Binomial, Migration, Biomass, CumSum, CumPer, Richness)

Terrestrial_mammals_Mass <- Mass_Migrations_spp %>%
  filter(Migratory=="Y" & Organisation_level2=="Terrestrial mammal" & CumPer<0.9)

Marine_mammals_Mass <- Mass_Migrations_spp %>%
  filter(Migratory=="Y" & Organisation_level2=="Marine mammal" & CumPer<0.9)

Terrestrial_birds_Mass <- Mass_Migrations_spp %>%
  filter(Migratory=="Y" & Organisation_level2=="Terrestrial birds" & CumPer<0.9)

Seabirds_Mass <- Mass_Migrations_spp %>%
  filter(Migratory=="Y" & Organisation_level2=="Seabirds" & CumPer<0.9)

MarineFish_Mass <- Mass_Migrations_spp %>%
  filter(Migratory=="Y" & Organisation_level2=="Marine fish" & CumPer<0.9)

AnaFish_Mass <- Mass_Migrations_spp %>%
  filter(Migratory=="Y" & (Organisation_level2=="Anadromous fish" | Organisation_level2=="Freshwater fish"))

MarineInv_Mass <- Mass_Migrations_spp %>%
  filter(Migratory=="Y" & Organisation_level2=="Marine invertebrates")

#write.csv(Mass_Migrations_spp, "/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/KAUST_Projects/Migrations_Projects/Data_wrangling/Final_R_scripts/Datasets/Mass_Migrations.csv")

# Plots - accumulation curves for mass migrations

Plot_PercentageBiomass <- function(data){
  
  out <- ggplot()+
    geom_line(data, mapping=aes(x=as.numeric(Position), y=CumPer, col=Migratory), alpha=0.5) +
    geom_hline(yintercept=0.9, color="grey", linetype="dashed", alpha=0.5)+
    #geom_vline(xintercept=20, color="red", linetype="dashed", alpha=0.5)+
    facet_wrap(~Organisation_level2, nrow=2, scales="free") +
    ylim(0,1)+
    #xlim(0,300)+
    ylab("Proportion of total biomass") +
    xlab("Number of species") +
    theme_bw() +
    theme(legend.position = "none")
  
  return(out)
}

Plot_PercentageBiomass(Mass_Migrations_spp %>% filter(Migratory=="Y"))

# Plots - comparison between migratory and non-migratory biomass per taxonomic group

Mass_Migrations_comparison <- Mass_Migrations_spp %>%
  group_by(Organisation_level2, Migratory) %>%
  summarise(Biomass = sum(Biomass), .groups = "drop") %>%  # Drop groups to avoid issues
  mutate(Migratory=ifelse(Migratory=="Y", "Migratory", "Non-Migratory")) %>%
  filter(!is.na(Migratory) & Biomass > 0) %>%  
  ungroup() %>%  # Ensure ungrouped before complete()
  complete(Organisation_level2, Migratory, fill = list(Biomass = NA))  # Add missing combinations

ggplot(Mass_Migrations_comparison, aes(x = Organisation_level2, y = Biomass, fill = Migratory)) +
  geom_col(position = "dodge", na.rm = FALSE) +  # na.rm = FALSE ensures missing values aren't removed
  scale_y_log10() +  
  scale_x_discrete(drop = FALSE) +  # Force ggplot to respect missing x-axis values
  scale_fill_manual(values = c("Migratory" = "#E69F00", "Non-Migratory" = "#56B4E9")) +
  xlab("")+
  ylab("Biomass (g)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for clarity


#
#
# 5- Load mass migrations datasets ####

Mass_migrations <- read.csv("2-SpeciesBiomass/Datasets_Biomass/mass_migrations.csv", fileEncoding = "ISO-8859-1") %>%  
  
  # Curate and filter by Species name
  mutate(Binomial = str_trim(word(Binomial, 1,2, sep=" "))) %>%
  filter(Binomial!="" & Binomial!=" ") %>%

  #OBTAIN MIGRATION INFORMATION
  rowwise() %>%
  
  #Modify movement type classes
  mutate(Movement_Type = str_trim(tolower(Movement.type))) %>% 
  mutate(Movement_Type = gsub("partial; ", "", Movement_Type)) %>% # remove partial if it's the first option
  mutate(Movement_Type = gsub("partial;", "", Movement_Type)) %>% # remove partial if it's the first option
  mutate(Movement_Type = gsub("seasonal; |season; ", "", Movement_Type)) %>% # remove partial if it's the first option
  mutate(Movement_Type = gsub("seasonal;|season;", "", Movement_Type)) %>% # remove partial if it's the first option
  mutate(Movement_Type = gsub("seasonal", "", Movement_Type)) %>% # remove partial if it's the first option
  mutate(Movement_Type = gsub("\\?", NA, Movement_Type)) %>% # remove partial if it's the first option
  
  separate(Movement_Type, c("Movement_Type"), sep=";") %>%
  
  mutate(Movement_Type = gsub("technical migrant|estuarine-inshore|inshore-offshore|local migrant|altitude|altitudinal|estuarine|river-inshore|depth|tidal|potadromous|anadromous", "local", Movement_Type)) %>%
  mutate(Movement_Type = gsub("intracontinental-home ranging|intercontinental-nomadising|interoceanic|intraoceanic|intercontinental|intracontinental|oceanodromous", "to-and-fro", Movement_Type)) %>%
  mutate(Movement_Type = gsub("to-and-fro-local|to-and-fro-nomadic", "to-and-fro", Movement_Type)) %>%
  mutate(Movement_Type = gsub("nomadising|home ranging|amphidromous|dispersal|nomadic-local|nomadisng", "nomadic", Movement_Type)) %>%
  mutate(Movement_Type = gsub("dvm", "DVM", Movement_Type)) %>%
  
  left_join(Species_list %>% filter(Migratory=="Y"), by="Binomial") %>%
  mutate(Binomial = ifelse(Binomial!=Original & !is.na(Original), Original, Binomial))


# check species not considered in the list of migratory species - only extinct species and other artificial groupings of species
Mass_migration_na <- Mass_migrations %>% 
  filter(!is.na(Binomial) & is.na(Original))

#write.csv(Mass_migrations, "Mass_Migrations_data.csv")

#

