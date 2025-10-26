#!/usr/bin/env
## This scripnt can  be used to automize the results of the BCA assay.
## Input = 96-well plate CSV file | Output = CSV file with total protein concentration per donor

# ---- Step 1 Preparations ----

# Load libraries
library(tidyverse)
library(stringr)

# Load necessary files

data_path <- readline(prompt = "File path to raw data (e.g. raw_data/test.csv): ")

plate_data <- read.csv(data_path)
plate_layout <- read_csv("metadata/BCA_plate_layout")
output_file <- "data/total_protein_donors.csv"

# Couple results to plate layout

data <- left_join(plate_data, plate_layout, by = "Well")

# ---- Step 2 Start calculations ----

# Correction with PBS blank

blank_mean <- data %>%
  filter(Type == "Blank") %>%
  summarise(mean_blank = mean(Abs, na.rm = TRUE)) %>%
  pull(mean_blank)

# Corrected data results

data <- data %>% mutate(Abs_corr = Abs - blank_mean)

# Mean and SD calculation per sample

summary <- data %>%
  filter(Sample != "PBS") %>%
  group_by(Sample, Condition, Conc, Type, Dilution) %>%
  summarise(
    mean_abs = mean(Abs_corr, na.rm = TRUE),
    SD_abs = sd(Abs_corr, na.rm = TRUE),
    .groups = "drop"
  ) 

# ---- Step 3 Callibration curve ----

# Standards curve lineair 

standards <- filter(summary, Type == "Standard")

fit <- lm(mean_abs ~ Conc, data = standards)
a <- coef(fit)[2]
b <- coef(fit)[1]
r2 <- summary(fit)$r.squared

callibration_curve <- ggplot(standards, aes(x = Conc, y = mean_abs)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_abs - SD_abs, ymax = mean_abs + SD_abs)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "BCA standard curve", 
       subtitle = paste0("Abs = ", round(a,6), "x Conc + ", round(b,6), "   |   R² = ", round(r2,5)),
       x = "Concentration (µg/mL)",
       y = "Absorbance (blank corrected)") +
  theme_minimal()

# User viewing check GO/NO GO

print(callibration_curve)
cat("\n Equation: Absorption = ", round(a,6), "x Concentration (µg/mL) + ", round(b,6), "   |   R² = ", round(r2,5), "\n\n")
answer <- readline(prompt = "Go through with calculations? (y/n)")
if (tolower(answer) == "n") {
  stop("\n \n Calculations halted by user.")
}


# ---- Step 4 Total protein calculation ----

# Total protein concentration corrected with dilution factor

summary <- summary %>% 
  filter(Type != "Standard") %>%
  mutate(
    Calc_conc = (mean_abs - b) / a,
    total_protein = Calc_conc * Dilution)

# Picking Donor dilution

donors <- c("Donor1", "Donor2", "Donor3")
results <- list()

for (donor in donors) {
 
  cat("\n\n", donor, "overview per dilution \n\n") 
  
  donor_data <- summary %>% 
  filter(Sample == donor) %>%
  select(Condition,Dilution,mean_abs,SD_abs,Calc_conc, total_protein)
  
  if (nrow(donor_data) == 0) {
    message("\n\n No data found for ", donor, ". \n SKIPPING.\n")
    next
  }
  
  print(donor_data)
  
  chosen_dil <- readline(prompt = paste0("Choose best dilution for ", donor, " (5/10/20):"))
  
  chosen_summary <- summary %>%
  filter(Sample == donor, Dilution == chosen_dil) %>%
  select(Sample, Condition, total_protein) %>%
  pivot_wider(names_from = Condition, values_from = total_protein,
              names_prefix = "total_protein_") %>%
  mutate(Date = Sys.Date()) %>%
  select(Sample, Date, starts_with("total_protein_"))
  
  cat("\n\n Chosen dilution:", chosen_dil, "x. \n\n")
  print(chosen_summary)
  
  results[[donor]] <- chosen_summary
}

# Combining all donor results

final_summary <- bind_rows(results)
cat("\n\n Overview of all donors: \n\n")
print(final_summary)


# ---- Step 5 Control check ----

# Control group overview

cat("\n\n Control group overview \n\n")

controls <- summary %>%
  filter(Type == "Control", Sample %in% c("PBS_Saliva", "Alpha-amylase")) %>%
  select(Sample, Conc, mean_abs, SD_abs, Calc_conc)

print(controls)

# Control groups ceck by user
answer2 <- readline(prompt = "Check control groups! Go through with calculations? (y/n)")

if (tolower(answer2) == "n") {
  stop("\n\n Calculations halted by user.")
}

# ---- Step 6 Save new donor(s)----

# Read or create file and determine n of donors present

if (file.exists(output_file)) {
  old_data <- read_csv(output_file, show_col_types = FALSE)
  start_n <- nrow(old_data) + 1
} else {
  old_data <- tibble()
  start_n <- 1
}

# Add donor IDs to summary

final_summary <- final_summary %>%
  mutate(
    Donor_ID = paste0("Donor_", seq(start_n, start_n + nrow(final_summary) - 1))
  ) %>%
  select(Donor_ID, Date, total_protein_Centrifuged, total_protein_Uncentrifuged)

# Definitive user check

cat("\n\n FINAL CHECK: \n\n Add the following donors?")
print(final_summary)
answer3 <- readline(prompt = "\n\n FINAL CHECK: Add these donors to the file? (y/n): ")

if (tolower(answer3) == "n") {
  stop("\n\n Data not added. \n\n")
}

# Combine data

write_csv(bind_rows(old_data, final_summary), output_file)

cat("\n\n Data updated in donor file: ", output_file, "\n\n")
cat("\n Donor IDs assigned from ", start_n, "to ", start_n + nrow(final_summary) - 1, ". \n")


