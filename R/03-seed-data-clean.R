library(tidyverse)
library(janitor)
library(readxl)

seed_raw <- read_xlsx("data/data-raw/Seed Data/FINAL - 2022 Biotron wheat grain (cleaned seed weights).xlsx",
                      skip = 5,
                      col_names = c("treatment", "label", "weight_g")) |> 
  janitor::clean_names() |> 
  rename(grain_yield_g = weight_g) |> 
  mutate(label = str_remove_all(label, " ")) |> 
  mutate(pot_label = paste(treatment, label, sep = "-")) |> 
  mutate(grain_yield_g = na_if(grain_yield_g, "EMPTY")) |> 
  mutate(notes = if_else(is.na(grain_yield_g), "empty bag", ""))

biomass <- readRDS("data/data-input/biomass_2022.rds")

meta_data_2022 <- biomass |> select(1:7)

seed <- left_join(seed_raw, meta_data_2022, by = "pot_label")

saveRDS(seed,"data/data-input/seed_2022.rds")
