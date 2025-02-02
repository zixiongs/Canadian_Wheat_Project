library(tidyverse)
library(janitor)
library(readxl)

# Cleaning 2022 data ------------------------------------------------------

biomass_raw_2022 <- read_xlsx("data/data-raw/Biomass Files/0_Biotron2022_Datasheet_2023.10.21.xlsx", sheet = "Dates and biomass") |> 
  clean_names()

biomass_2022 <- biomass_raw_2022 |> 
  select(year, biome, treatment_corr, pot_label, genotype, tiller_num, plant_height, ear_biomass, shoot_biomass) |> 
  mutate(pot_label = str_remove_all(pot_label, " ")) |> 
  #filter(rowSums(across(everything(), is.na)) < 4) |> 
  mutate(temperature_level = str_extract(treatment_corr, "T[0-8]"),
         co2_level = str_extract(treatment_corr, "AC|EC"), .after = treatment_corr) |> 
  mutate(above_ground_biomass = ear_biomass + shoot_biomass) |> 
  mutate(genotype = str_replace_all(genotype, " ", "_")) |> 
  mutate(genotype = str_replace_all(genotype, "\\bBrandon\\b", "AAC_Brandon")) # Rename this genotype to match with Eric

saveRDS(biomass_2022, "data/data-input/biomass_2022.rds")



  
