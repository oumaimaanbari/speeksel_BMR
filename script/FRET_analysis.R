library(tidyverse)

# Inlezing

fret_raw <- read_tsv("~/Documents/School/BMR/speeksel_BMR/raw_data/FRET/FRET_2910.TSV",
                     skip = 6,
                     locale = locale(decimal_mark = ","),
                     col_names = TRUE)
names(fret_raw)[1:2] <- c("Well", "Content")
fret_raw <- fret_raw %>% select(-last_col())

# Tidy data

fret <- fret_raw %>%
  pivot_longer(cols = matches("min$"),
               names_to = "Time",
               values_to = "Fluorescence") %>%
  mutate(Time = parse_number(Time))


