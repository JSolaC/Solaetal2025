# =======================================================================
# Supplementary R Code: Processing and Merging Taxonomic Species Lists
# =======================================================================
# Authors: Sola et al.
# Date: 2025 - February
# Journal: Nature Portfolio
# =======================================================================
# Overview:
# This script refines taxonomic species lists by standardising names, resolving synonyms, and merging multiple species datasets into a final unified list.
#
# Key Considerations:
# - Species can be referred to using different scientific names (synonyms).
# - Synonyms are compiled into a reference library to ensure accurate merging.
# - Data sources include online repositories and biomass datasets (see Methods and Supplementary Materials for further reference).
# - Marine invertebrates were initially excluded from migratory species lists and have been manually verified and added at the end.
#
# Computational Constraints:
# - The script involves API queries and large dataset processing, making execution time-consuming (several days). Pre-processed datasets are provided for reproducibility.
# - Due to missing taxonomic data, some species were manually assigned classifications.
#
# Issues and Limitations:
# - Potential misclassifications may arise due to automated synonym resolution.
# - Some species require manual verification to ensure accurate taxonomy.
#
# Outputs:
# - The final dataset includes all unique species, with synonyms standardized and duplicates removed.
# =======================================================================

# 0- Set up the Workspace -------

# Define the required packages
packages <- c("tidyverse", "taxize", "janitor", "splitstackshape", "readxl", "gdata")

# Install missing packages
lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cran.rstudio.com/")
})

# Load all packages
lapply(packages, library, character.only = TRUE)

# Load functions
'%ni%' <- Negate('%in%') # conditional for instances within a vector NOT included in another vector

# Set up the working directory
setwd("/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/KAUST_Projects/Migrations_Projects/Data_wrangling/Publication/")

#
# 1- Load initial datasets ------------------

# Raw species list from websites
Species_list_raw <- read.csv("1-SpeciesList/Datasets/Species_List_Raw.csv")

# Raw species list from the biomass datasets (see Methods and Supplementary Materials for references)
Species_nonmigrant <- read.csv('1-SpeciesList/Datasets/Biomass_allspp_nomigandmig.csv') %>%
  filter(Migrant_species!="Y")

#
# 2- Obtain species taxonomy (no need to run this part) ------------------

# A) Migratory species -------

# retrieve taxonomic data per species to classify them into groups
spp_list <- Species_list_raw %>%
  mutate(Binomial=str_replace_all(Scientific.name, "[[:punct:]]", " ")) %>%
  mutate(Binomial1=gsub("\\�"," ",Binomial)) %>%
  mutate(Binomial2=str_trim(Binomial1)) %>%
  mutate(Binomial3 = word(Binomial2, 1,2, sep=" ")) %>%
  mutate(n_words = str_count(Binomial3, '\\w+')) %>%
  filter(n_words>1) %>% # checked cases with only 1 word, and they don't contribute new species
  distinct(Binomial3, .keep_all = F) 

spp_list1 <- as.character(spp_list$Binomial3)[1:2916]
spp_list2 <- as.character(spp_list$Binomial3)[2917:5832]

migration_spp_info <- classification(spp_list2, db = "itis") 

species_taxonomic_data <- migration_spp_info %>%
  as.list() %>%
  lapply(function(x) {t(x)}) %>%
  lapply(function(x) {as.data.frame(x)[c(3,2,1),]}) %>%
  lapply(function(x) {
    as.data.frame(x)%>%
      row_to_names(row_number = 2)
  }) %>%
  lapply(function(x) {as.data.frame(x)}) %>%
  #discard(~ nrow(.x) == 0) %>%
  bind_rows() %>% #flag cases: 181,182,186,187. Total rows= 4839 
  filter(!is.na(species)) # 250 species not found, which matches what I got from the function classification

# The resulting species list after accounting for all
Species_list_taxonomy <- rbind(
  read.csv('1-SpeciesList/Datasets/Taxonomy/Species_list_taxsize_1.csv'),
  read.csv('1-SpeciesList/Datasets/Taxonomy/Species_list_taxsize_2.csv')
)

#
# B) Non-migratory species ------------

no_migrant_spp_list1<- as.character(unique(Species_nonmigrant$Binomial))[1:500]
no_migrant_spp_list2<- as.character(unique(Species_nonmigrant$Binomial))[501:1000]
no_migrant_spp_list3<- as.character(unique(Species_nonmigrant$Binomial))[1001:2000]
no_migrant_spp_list4<- as.character(unique(Species_nonmigrant$Binomial))[2001:3000]
no_migrant_spp_list5<- as.character(unique(Species_nonmigrant$Binomial))[3001:4000]
no_migrant_spp_list6<- as.character(unique(Species_nonmigrant$Binomial))[4001:5000]
no_migrant_spp_list7<- as.character(unique(Species_nonmigrant$Binomial))[5001:6000]
no_migrant_spp_list8<- as.character(unique(Species_nonmigrant$Binomial))[6001:7000]
no_migrant_spp_list9<- as.character(unique(Species_nonmigrant$Binomial))[7001:8000]
no_migrant_spp_list10<- as.character(unique(Species_nonmigrant$Binomial))[8001:9000]
no_migrant_spp_list11<- as.character(unique(Species_nonmigrant$Binomial))[9001:10000]
no_migrant_spp_list12<- as.character(unique(Species_nonmigrant$Binomial))[10001:11000]
no_migrant_spp_list13<- as.character(unique(Species_nonmigrant$Binomial))[11001:12000]
no_migrant_spp_list14<- as.character(unique(Species_nonmigrant$Binomial))[12001:12407]

no_migration_spp_info <- classification(no_migrant_spp_list1, db = "itis") %>%
  as.list() %>%
  lapply(function(x) {t(x)}) %>%
  lapply(function(x) {as.data.frame(x)[c(3,2,1),]}) %>%
  lapply(function(x) {
    as.data.frame(x)%>%
      row_to_names(row_number = 2)
  }) %>%
  lapply(function(x) {as.data.frame(x)}) %>%
  bind_rows() %>%
  filter(!is.na(species)) %>% 
  dplyr::select(species, phylum, class, subclass, superclass, infraphylum, order, suborder, superfamily, family, subfamily, genus) 

# The resulting dataset
no_migration_spp_info_taxa <- bind_rows(
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info1.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info2.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info3.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info4.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info5.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info6.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info7.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info8.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info9.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info10.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info11.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info12.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info13.csv"),
  read.csv("1-SpeciesList/Datasets/Taxonomy/no_migration_spp_info14.csv")
)

no_migration_spp_no_info_taxa <- Species_nonmigrant %>%
  dplyr::select("species"=Binomial) %>%
  filter(species %ni% no_migration_spp_info_taxa$species)

Non_migrant_species <- bind_rows(no_migration_spp_info_taxa, no_migration_spp_no_info_taxa)

#
# 3- Group Species by Taxonomy ---------
# A) Migratory species ------------------

# After obtaining species taxonomy, the dataset was manually processed to summarise taxonomy into three levels or organisation for species that did not present taxonomic information
Species_list_manual <- read.csv("1-SpeciesList/Datasets/Species_List_Manual.csv")

# Obtain taxonomic groups from data 
Obtain_spp_groups <- function(dataset){
  
  dataset_out1 <- dataset %>%
    
    dplyr::mutate(FID=1:n()) %>% 
    
    # Define tha main organisation levels
    mutate(Organisation_level=NA) %>%
    mutate(Organisation_level = ifelse(superclass=="Chondrichthyes" | superclass=="Actinopterygii" | superclass=="Sarcopterygii" | infraphylum=="Agnatha", 
                                       "Fish", Organisation_level)) %>% # since there's no data for FW fish, there's no need to make the distinction
    mutate(Organisation_level = ifelse(phylum=="Arthropoda" | phylum=="Mollusca" | phylum=="Echinodermata", "Arthropods", Organisation_level)) %>%
    mutate(Organisation_level = ifelse(class=="Mammalia", "Mammal", Organisation_level)) %>%
    mutate(Organisation_level = ifelse(class=="Aves", "Birds", Organisation_level)) %>%
    mutate(Organisation_level= ifelse(class=="Reptilia", "Reptiles", Organisation_level)) %>%
    mutate(Organisation_level= ifelse(class=="Amphibia", "Amphibians", Organisation_level)) %>%
    
    mutate(Organisation_level2=Organisation_level) %>%
    
    # Divide birds into terrestrial and seabirds
    mutate(Organisation_level2 = ifelse(order=="Sphenisciformes" | family=="Phaethontidae" | order=="Procellariiformes" | (order=="Suliformes" & genus!="Anhinga") | family=="Laridae" | family=="Stercorariidae" | family=="Alcidae",  
                                        "Seabirds", Organisation_level2) ) %>%
    mutate(Organisation_level2 = ifelse(Organisation_level=="Birds" & Organisation_level2!="Seabirds", 
                                        "Terrestrial birds", Organisation_level2)) %>%
    mutate(Organisation_level2 = ifelse(Organisation_level2=="Land birds", "Terrestrial birds", Organisation_level2)) %>%
    
    # Divide mammals into marine and terrestrial mammals
    mutate(Organisation_level2 = ifelse(order=="Cetacea" | order=="Sirenia" | family=="Odobenidae" | family=="Otariidae" | family=="Phocidae" | Binomial=="Ursus maritimus", "Marine mammal", Organisation_level2)) %>%
    mutate(Organisation_level2 = ifelse(Organisation_level=="Mammal" & Organisation_level2!="Marine mammal", "Terrestrial mammal", Organisation_level2))%>%
    mutate(Organisation_level2 = ifelse(Organisation_level2=="Aquatic mammal" | Organisation_level2=="Freshwater mammal", "Marine mammal", Organisation_level2)) %>%
    
    # Divide invertebrates into terrestrial and marine
    mutate(Organisation_level2 = replace(Organisation_level2, order=="Decapoda" | order=="Euphausiacea" | order=="Stomatopoda" | order=="Isopoda" | order=="Amphipoda" | order=="Leptostraca" | 
                                           subclass=="Pteriomorphia" | subclass=="Heterodonta" | subclass=="Anomalodesmata" | subclass=="Palaeoheterodonta" | subclass=="Protobranchia" | class=="Cephalopoda", "Marine invertebrates")) %>%
    mutate(Organisation_level2 = ifelse(Organisation_level=="Arthropods" & Organisation_level2!="Marine arthropods", "Terrestrial invertebrates", Organisation_level2)) %>%
    
    mutate(Organisation_level3=Organisation_level2) %>%
    mutate(Organisation_level3 = replace(Organisation_level3, order=="Chiroptera", "Bats")) %>%
    mutate(Organisation_level3 = replace(Organisation_level3, class=="Chondrichthyes", "Sharks")) %>%
    mutate(Organisation_level3 = replace(Organisation_level3, superfamily=="Chelonioidea", "Marine reptiles"))
  
  dataset_out2 <- dataset_out1 %>%
    mutate(Organisation_level2 = ifelse (Organisation_level3=="Marine reptiles", "Sea turtles", Organisation_level2)) %>%
    mutate(Organisation_level2 = ifelse (Organisation_level2=="Lizzard" | Organisation_level2=="Snakes" | Organisation_level2=="Turles & tortoises", "Reptiles", Organisation_level2)) %>%
    mutate(Organisation_level2 = ifelse (Organisation_level2=="Terrestrial mammal" | Organisation_level2=="Bats", "Terrestrial mammal", Organisation_level2)) %>%
    mutate(Organisation_level2 = ifelse (Organisation_level2=="Terrestrial arthropods" | Organisation_level2=="Butterflies & moths" | Organisation_level2=="Dragonflies" | 
                                           Organisation_level2=="Flies, aphids & others", "Terrestrial invertebrates", Organisation_level2)) %>%
    
    mutate(Organisation_level2 = ifelse(Organisation_level=="Fish", "Marine fish", Organisation_level2)) %>%
    
    mutate(Organisation_level2 = ifelse (Organisation_level=="Arthropods" & Organisation_level2!="Terrestrial invertebrates",
                                         "Marine invertebrates", Organisation_level2)) 
  
  return(dataset_out2)
  
}

Species_list_taxonomy2 <- Species_list_taxonomy %>%
  dplyr::rename(Binomial=species) %>%
  mutate(n_words = str_count(Binomial, '\\w+')) %>%
  filter(n_words<3) %>% # make sure that no subspecies are in the dataset
  distinct(Binomial, .keep_all=T) %>%
  
  dplyr::mutate(FID=1:n()) %>% 
  
  Obtain_spp_groups() %>%
  
  left_join(Species_list_manual, by="Binomial") %>%
  
  mutate(
    Organisation_level = ifelse(is.na(Organisation_level.x), Organisation_level.y, Organisation_level.x),
    Organisation_level2 = ifelse((Organisation_level2.x=="Marine or freshwater fish" | is.na(Organisation_level2.x)) & !is.na(Organisation_level2.y), Organisation_level2.y, Organisation_level2.x),
    Organisation_level3 = ifelse(is.na(Organisation_level3.y) | Organisation_level3.y=="", Organisation_level3.x, Organisation_level3.y)
  ) %>%
  
  dplyr::select(Binomial, phylum, infraphylum, superclass, class, subclass, order, superfamily, family, genus, Organisation_level, Organisation_level2, Organisation_level3, Migration)

#Add species for which taxonomy was not found (and therefore were removed from the dataset)
Species_without_taxonomy <- Species_list_manual %>%
  filter(!is.na(Binomial) & Binomial!="") %>%
  mutate(Binomial = word(Binomial, 1,2, sep=" ")) %>%
  distinct(Binomial, .keep_all=T) %>% 
  dplyr::select("Binomial"=Binomial1, Organisation_level, Organisation_level2, Organisation_level3, Migration) %>%
  filter(Binomial %ni% Species_list_taxonomy$Binomial)

# Join the taxonomy list and the missing taxonomy into one list + group into migratory classes

Species_list_taxonomy3 <- Species_list_taxonomy2 %>%
  
  bind_rows(Species_without_taxonomy) %>%
  
  #Obtain migratory information across species
  mutate(Migration2 = ifelse(Migration=="partial" | Migration=="intracontinental and partial" | Migration=="Adults sedentary" | 
                               Migration=="partly migratory" | Migration=="partial and intracontinental", 
                        "Partial migration", Migration ) ) %>%
  mutate(Migration2 = ifelse(Migration=="amphidromous" | Migration=="diadromous" | Migration=="oceano-estuarine" | Migration=="anadromous" | 
                          Migration=="amphidromous?" | Migration=="anadromous?",
                        "Ocean-River",Migration2) ) %>%
  mutate(Migration2 = ifelse( Migration=="catadromous" | Migration=="limnodromous" ,
                         "River-Ocean",Migration2) ) %>%
  mutate(Migration2 = ifelse(Migration=="local migrant" | Migration=="range extension" |  Migration=="Small scale displacemet" | 
                               Migration=="Local and seasonal movements" | Migration=="Sedentary" | 
                               Migration=="Local movements from unsuitable habitats" | Migration=="local",
                        "Local migration", Migration2) ) %>%
  mutate(Migration2 = ifelse(Migration=="technical migrant" | Migration=="Groms migrant" | Migration=="emigration", 
                        "Migratory", Migration2)) %>%
  mutate(Migration2 = ifelse(Migration=="intercontinental" | Migration=="interoceanic" | Migration=="intracontinental" | Migration=="intraoceanic",
                        "Large scale migration", Migration2) ) %>%
  mutate(Migration2 = ifelse(Migration=="altitudinal movement" | Migration=="altitudinal migrant",
                        "Altitude migration", Migration2) ) %>%
  mutate(Migration2 = ifelse(Migration=="data deficient" | Migration=="Possibly migratory" | Migration=="possibly migratory" | Migration=="", 
                        "Missing data", Migration2) ) %>%
  mutate(Migration2 = ifelse(Migration=="nomadising", "Nomad migration", Migration2)) %>%
  mutate(Migration2 = ifelse(Migration=="potamodromous", "Freshwater migration", Migration2)) %>%
  mutate(Migration2 = ifelse(Migration=="oceanodromous", "Ocean migration", Migration2)) %>%
  
  mutate(Migration2 = ifelse(is.na(Migration2) | Migration2=="", "No information", Migration2)) %>%
  mutate(Organisation_level2 = ifelse (Organisation_level=="Fish" & Migration2=="Freshwater migration", "Freshwater fish", Organisation_level2)) %>%
  mutate(Organisation_level2 = ifelse (Migration2=="Ocean-River", "Anadromous fish", Organisation_level2)) %>%
  mutate(Organisation_level3 = ifelse (Migration2=="Ocean migration" & Organisation_level=="Fish", "Marine fish", Organisation_level3)) %>%
  
  # For an unknown reason, this species is not classified as a terrestrial bird - checked through the dataset again and it seems to be an isolated case
  mutate(Organisation_level2 = ifelse(Binomial=="Alopochen aegyptiacus", "Terrestrial birds", Organisation_level2)) %>%
  
  # For each species, filter repeated cases that present NA in phyllum while the other row does not present NAs
  group_by(Binomial) %>%
  # Remove cases where Phylum is NA *only if* there are non-NA values for that Binomial
  filter(!(is.na(phylum) & any(!is.na(phylum)))) %>%
  ungroup()

#write.csv(Species_list_taxonomy3, "1-SpeciesList/Output/Species_list_taxonomy_migratory.csv")

#
# B) Non-migratory species ####

NonMigratory_species_taxa <- Non_migrant_species %>%
  
  dplyr::rename(Binomial = species) %>%
  Obtain_spp_groups() %>%
  # For each species, filter repeated cases that present NA in phyllum while the other row does not present NAs
  group_by(Binomial) %>%
  # Remove cases where Phylum is NA *only if* there are non-NA values for that Binomial
  filter(!(is.na(phylum) & any(!is.na(phylum)))) %>%
  ungroup()

# Given that we will not be using this list to get species-specific information, we resolved by not manually checking each individual species. 
#write.csv(NonMigratory_species_taxa, "NonMigratory_species_list_taxa.csv")

#
# 4- Find synonyms and remove repeated species (no need to run this part) ----------
# A) Migratory species ####

# Look for synonyms

Fill_entries <- function(x){
  if(is_empty(x)){
    out <- data.frame(
      #"Binomial"= NA,
      "sub_tsn"= NA ,
      'acc_name'= NA,
      "acc_tsn"= NA,
      "acc_author"= NA,
      "syn_author"= NA,
      "syn_name"= NA,
      "syn_tsn"= NA,
      "x"= NA
    )} else{
      out <- x
    }
  
  return(out)
  
}

Species_corrected_synonyms_list <- synonyms(c(Species_list_taxonomy3$Binomial), db="itis", rows=1)

Species_corrected_synonyms <- Species_corrected_synonyms_list %>%
  
  # put the list into a dataframe - filling in empty entries with NAs
  lapply(function(x) {as.data.frame(Fill_entries(x))}) %>%
  data.table::rbindlist( idcol = "Binomial", fill=TRUE ) %>%
  
  # transform the dataset into a long format with only two columns in order to check for reapeated names
  dplyr::select(Binomial, acc_name, syn_name) %>%
  pivot_longer(
    cols = acc_name:syn_name,
    names_to = "Type",
    values_to = "Synonyms"
  ) %>%
  
  #remove duplicated rows - after not considering subspecies
  mutate(Synonyms = word(str_trim(Synonyms), 1,2, sep=" ")) %>%
  distinct(Binomial, Synonyms, .keep_all = T) %>%
  
  #for cases where Synonyms and Binomials match, remove the synonym
  rowwise() %>%
  mutate(Synonyms = ifelse(Synonyms==Binomial, NA, Synonyms)) %>%
  
  # identify cases where the binomial is contained in synonyms and the synonyms are contained in the binomial
  ungroup() %>%
  mutate(
    Binomial_Repeated = ifelse(Binomial %in% Synonyms, 1, 0),
    Synonym_Repeated = ifelse(Synonyms %in% Binomial, 1, 0)
  )

# remove species already present in the filtered list (as we know that these ones have been manually and appropiately allocated). 
# the only repeated species is Myotis aurascens, which was removed at first instance since it was a synonym of Myotis hajastanicus, while being a different species
Species_corrected_synonyms2 <- Species_corrected_synonyms %>%
  
  mutate(Synonyms=ifelse(is.na(Synonyms), Binomial, Synonyms)) %>%
  filter(Binomial!="") %>%
  filter(Binomial %ni% Species_filtered_synonyms$Binomial & Binomial %ni% Species_filtered_synonyms$Synonyms &
           Synonyms %ni% Species_filtered_synonyms$Synonyms & Synonyms %ni% Species_filtered_synonyms$Binomial)

# Due to the presence of species complexes and some inaccuracies in the matching of synonyms with binomials, species were checked one by one to ensure that consistent criteria is applied across all cases
#write.csv(Species_corrected_synonyms, "/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/KAUST_Projects/Migrations_Projects/Data_wrangling/Final_R_scripts/Datasets/Species_corrected_synonyms.csv")
#write.csv(Species_corrected_synonyms2, "/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/KAUST_Projects/Migrations_Projects/Data_wrangling/Final_R_scripts/Datasets/Species_corrected_synonyms2.csv")

Species_filtered_synonyms <- bind_rows(
  read.csv("1-SpeciesList/Datasets/Synonyms/Species_corrected_synonyms_done.csv"),
  read.csv('1-SpeciesList/Datasets/Synonyms/Species_corrected_synonyms2.csv') %>%
    dplyr::select(Binomial, Synonyms)
  ) %>%
  mutate(Synonyms = ifelse(Synonyms=="" | Synonyms==Binomial,NA,Synonyms)) %>%
  mutate(
    Binomial = word(str_trim(Binomial), 1,2, sep=" "),
    Synonyms = word(str_trim(Synonyms), 1,2, sep=" ")
  ) %>% # modify synonyms to exclude subspecies
  distinct(Binomial, Synonyms)

# ADD REPEATED NAMES BASED ON THEIR COMMON NAME! (ANTIGONE VS GRUS!)

Species_list_raw_duplicated <- Species_list_raw %>%
  mutate(Binomial=str_replace_all(Scientific.name, "[[:punct:]]", " ")) %>%
  mutate(Binomial1=gsub("\\�"," ",Binomial)) %>%
  mutate(Binomial2=str_trim(Binomial1)) %>%
  mutate(Binomial3 = word(Binomial2, 1,2, sep=" ")) %>%
  mutate(Common.name1=str_replace_all(Common.name, "[[:punct:]]", " ")) %>%
  mutate(Common.name2=gsub("\\�"," ",Common.name1)) %>%
  mutate(Common.name3=str_to_lower(Common.name2)) %>%
  mutate(Common.name4=str_trim(Common.name3)) %>%
  distinct(Binomial3, Common.name4) %>%
  group_by(Common.name4) %>% dplyr::mutate(Count=n()) %>% filter(Count==2 & Count <10) %>%
  #remove cases with the same common name refering to different species
  filter(
      Common.name4!="white catfish" & #Different species
      Common.name4!="spotted lanternfish" & #Different species
      Common.name4!="snaggletooth" & #Different species
      Common.name4!="moon fish" & #Different species
      Common.name4!="largetooth sawfish" & #Different species
      Common.name4!="indian mottled eel" & #Different species
      Common.name4!="indian mackerel" & #Different species
      Common.name4!="haddock" & #Different species
      Common.name4!="goonch" & #Different species
      Common.name4!="freshwater moray" & #Different species
      Common.name4!="fan tailed warbler" & #Different species
        Common.name4!="elephant fish" & #Different species
        Common.name4!="croaker" & #Different species
        Common.name4!="cocco s lantern fish" & #Different species
        Common.name4!="bullhead" & #Different species
        Common.name4!="bleeker s whipray" & #Different species
        Common.name4!="atlantic sturgeon" & #Different species
        Common.name4!="black skipjack" #Different species
    ) %>%
  left_join(Species_filtered_synonyms %>% dplyr::select("Binomial3"=Binomial) %>% mutate(Type="Binomial")) %>%
  left_join(Species_filtered_synonyms %>% dplyr::select("Binomial3"=Synonyms) %>% mutate(Type="Synonyms")) %>%
  distinct(Binomial3, Common.name4, Type) %>%
  group_by(Common.name4) %>% dplyr::mutate(non_na_count = sum(!is.na(Type))) %>%
  # now that we know which species appear twice, we have more information to remove other species that are different
  filter(
    Common.name4!="yellow legged gull" & #Different species
      Common.name4!="woolly necked stork" & #Different species
      Common.name4!="queen mackerel" & #Different species
      Common.name4!="lake trout" & #Different species
      Common.name4!="chub" #Different species
  ) %>%
  # This removes species already present as synonyms (so no need to check common names)
  left_join(Species_filtered_synonyms %>% dplyr::select("Binomial3"=Binomial, Synonyms)) %>%
  filter(Binomial3 %ni% Synonyms) %>%
  
  # now, identify common names with only one binomial left (so that can be removed since they don't have any undetected synonyms)
  group_by(Common.name4) %>% dplyr::mutate(Binomial_n=length(unique(Binomial3))) %>%
  # Remove cases with only one binomial left
  filter(Binomial_n>1) %>%
  
  # also, check cases with more than 1 synonym that may have NAs
  group_by(Common.name4, Binomial3) %>% dplyr::mutate(Synonyms_n=length(unique(Synonyms))) %>%
  # cases with only one synonym, change any NAs to "none"
  mutate(Synonyms = ifelse(Synonyms_n==1 & is.na(Synonyms), "none", Synonyms)) %>%
  # remove NAs from cases with >1 synonym
  filter(!is.na(Synonyms)) %>%
  # back transform "none" synonyms to NA
  mutate(Synonyms = ifelse(Synonyms=="none", NA, Synonyms)) %>%
  
  #We are left with "true" undetected synonyms. Which we classify into priority (name that will be the Binomial) and the others will be synonyms
  group_by(Common.name4) %>% dplyr::mutate(Priority= dplyr::first(Binomial3)) %>%
  dplyr::select("Binomial"=Priority, "Synonyms"=Binomial3, Common.name4) %>%
  filter(Binomial!=Synonyms)


# The outcome of this process is two datasets
Species_synonyms_list <- Species_filtered_synonyms %>%
  bind_rows(Species_list_raw_duplicated %>% ungroup() %>% dplyr::select(Binomial, Synonyms)) %>%
  filter(Binomial %ni% Synonyms) %>%
  distinct(Binomial, Synonyms)

#write.csv(Species_synonyms_list, "Species_synonyms_list.csv")

Migratory_species_list <- read.csv("1-SpeciesList/Datasets/Species_synonyms_list.csv")

#
# B) Non-migratory species ####

NonMigrant_Species_corrected_synonyms_list <- synonyms(c(NonMigratory_species_taxa$Binomial), db="itis", rows=1)

NonMigrant_Species_corrected_synonyms <- NonMigrant_Species_corrected_synonyms_list %>%
  
  # put the list into a dataframe - filling in empty entries with NAs
  lapply(function(x) {as.data.frame(Fill_entries(x))}) %>%
  data.table::rbindlist( idcol = "Binomial", fill=TRUE ) %>%
  
  # transform the dataset into a long format with only two columns in order to check for reapeated names
  dplyr::select(Binomial, acc_name, syn_name) %>%
  pivot_longer(
    cols = acc_name:syn_name,
    names_to = "Type",
    values_to = "Synonyms"
  ) %>%
  
  #remove duplicated rows - after not considering subspecies
  mutate(Synonyms = word(str_trim(Synonyms), 1,2, sep=" ")) %>%
  distinct(Binomial, Synonyms, .keep_all = T) %>%
  
  #for cases where Synonyms and Binomials match, remove the synonym
  rowwise() %>%
  mutate(Synonyms = ifelse(Synonyms==Binomial, NA, Synonyms)) %>%
  
  # remove migratory species 
  filter(Binomial %ni% unique(Species_filtered_synonyms$Binomial) & Synonyms %ni% unique(Species_filtered_synonyms$Synonyms)) %>%
  
  # identify cases where the binomial is contained in synonyms and the synonyms are contained in the binomial
  ungroup() %>%
  mutate(
    Binomial_Repeated = ifelse(Binomial %in% Synonyms, 1, 0),
    Synonym_Repeated = ifelse(Synonyms %in% Binomial, 1, 0)
  )

# Due to the presence of species complexes and some inaccuracies in the matching of synonyms with binomials, I have resolved to check species by species to ensure consistent criteria is applied across all cases
#write.csv(NonMigrant_Species_corrected_synonyms, "NonMigrant_Species_corrected_synonyms.csv")

# remove species already present in the filtered list (as we know that these ones have been manually and appropiately allocated). 
# the only repeated species is Myotis aurascens, which was removed at first instance since it was a synonym of Myotis hajastanicus, while being a different species
NonMigrant_Species_corrected_synonyms2 <- NonMigrant_Species_corrected_synonyms %>%
  
  mutate(Synonyms=ifelse(is.na(Synonyms), Binomial, Synonyms)) %>%
  filter(Binomial!="") %>%
  filter(Binomial %ni% NonMigrant_Species_corrected_synonyms$Binomial & Binomial %ni% NonMigrant_Species_corrected_synonyms$Synonyms &
           Synonyms %ni% NonMigrant_Species_corrected_synonyms$Synonyms & Synonyms %ni% NonMigrant_Species_corrected_synonyms$Binomial)

# The resulting datasets from all this process

NonMigrant_Species_synonyms_list <- read.csv("NonMigrant_Species_corrected_synonyms_done.csv") %>%
  mutate(Synonyms = ifelse(Synonyms=="",NA,Synonyms)) %>%
  distinct(Binomial, Synonyms) %>%
  mutate(
    Binomial = word(str_trim(Binomial), 1,2, sep=" "),
    Synonyms = word(str_trim(Synonyms), 1,2, sep=" ")
  ) # modify synonyms to exclude subspecies

#write.csv(NonMigrant_Species_synonyms_list, "NonMigrant_Species_synonyms_list.csv")

Non_migratory_species_list <- read.csv("1-SpeciesList/Datasets/NonMigrant_Species_synonyms_list.csv")

#
# 5- Join the two datasets (and include missing marine invertebrates migrant species) ####

# Species synonyms list

# After running this code for all species (taking a total of 1-2 weeks), marine invertebrates were also added to the dataset. To facilitate the addition of these species, we included them at this step.
Mar_inv_mig <- read.csv('/Users/jordisola/Library/Mobile Documents/com~apple~CloudDocs/Documents/KAUST_Projects/Migrations_Projects/Data_wrangling/Final_R_scripts/Datasets/Species_list_marine_inv_migratory.csv')

Species_list_all <- bind_rows(
  Migratory_species_list %>% dplyr::select(Binomial, Synonyms) %>% mutate(Migratory="Y"),
  # If non-migrant species are found in the migratory marine invertebrates list - then they are migratory
  Non_migratory_species_list %>% dplyr::select(Binomial, Synonyms) %>% filter(Binomial %in% Mar_inv_mig$Binomial) %>% mutate(Migratory="Y"),
  # If they are not found - then they are not migratory
  Non_migratory_species_list %>% dplyr::select(Binomial, Synonyms) %>% filter(Binomial %ni% Mar_inv_mig$Binomial) %>% mutate(Migratory="N")
)

# Since Binomial is not included as a synonym, and to allow for the effective conversion of all potential synonyms into the same Binomial (and not introduce NAs in the dataset), we also include Binomial as a Synonym (even though that involves repeating Binomial across the two columns)
Species_list_all_out <- bind_rows(
  Species_list_all %>% dplyr::select("Original"=Binomial, "Binomial"=Synonyms, Migratory),
  Species_list_all %>% dplyr::select("Original"=Binomial, "Binomial"=Binomial, Migratory)
) %>% 
  distinct(Original, Binomial, .keep_all=TRUE) %>%
  filter(!is.na(Binomial) )

#write.csv(Species_list_all_out, "Species_list_all_out.csv")

#