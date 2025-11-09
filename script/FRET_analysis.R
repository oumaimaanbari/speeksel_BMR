#!/usr/bin/env
## This script can  be used to automize the results of the FRET assay.
## Input = 96-well plate TSV file decimal "," | Output = CSV file with F/min and graphs

# ---- Step 1 Preparations ----

# Load libraries
library(tidyverse)

# Load necessary files

file_name <- readline(prompt = "File name (e.g. FRET_2910.tsv) ")

if (!file.exists(paste0("raw_data/FRET/", file_name))) {
  stop(file_name, " File not found!")
}

plate_data <- read_tsv(paste0("raw_data/FRET/", file_name),
                       skip = 6,
                       locale = locale(decimal_mark = ","),
                       col_names = TRUE)
names(plate_data)[1:2] <- c("Well", "Content")
plate_data <- plate_data %>% select(-last_col())

plate_layout <- read_csv("metadata/FRET_plate_layout.csv")

# Tidy data

data <- plate_data %>%
  pivot_longer(cols = matches("min$"),
               names_to = "Time_min",
               values_to = "Fluorescence_RFU") %>%
  mutate(Time_min = parse_number(Time_min))

data$Content <- factor(data$Content)


# Couple results to plate layout

data <- left_join(data, plate_layout, by = "Well")


# ---- Step 2 Start Calculations ----

# Average and STDEV


# Mean and SD calculation
fret_summary <- data %>%
  group_by(Type, Sample, Conc, Time_min, Condition) %>% 
  summarise(mean_RFU = mean(Fluorescence_RFU, na.rm = TRUE),
            sd_RFU = sd(Fluorescence_RFU, na.rm = TRUE),
            .groups = "drop")

# ---- Step 3 Plotting loop ----


# Combine Sample with Condition 
fret_plot <- fret_summary %>%
  mutate(Sample_label = if_else(is.na(Condition), 
                                Sample, 
                                paste(Sample, Condition)))

fret_plot <- fret_plot %>% 
  mutate(Sample_label = parse_factor(Sample_label, levels = unique(Sample_label)))

# Donor vector
donors <- levels(parse_factor(fret_plot$Sample, levels = unique(fret_plot$Sample)))
donors <- samples[str_starts(samples, "Donor")]

# Curve including PBS and Trypsin
for (donor in donors) {
  message("Plotting for:  ", donor)
  
  conc <- fret_plot %>% 
    filter(Sample == donor) %>%
    pull(Conc) %>%
    unique()
  
  for (c in conc) {
    
    plot_summary <- fret_plot %>% 
      filter(Sample %in% c(donor, "PBS", "Trypsin"),
             Conc == c | is.na(Conc))
    
    plot <- ggplot(plot_summary,
    aes(Time_min, mean_RFU, colour = Sample_label)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_RFU - sd_RFU, ymax = mean_RFU + sd_RFU), 
                  width = 0.3, alpha = 0.5) + 
    geom_ribbon(aes(ymin = mean_RFU - sd_RFU, ymax = mean_RFU + sd_RFU, fill = Sample_label),
                  alpha = 0.1, colour = NA) +
    theme_minimal() +
    labs(
      title = paste0("FRET assay ", donor, " (Protein concentration: ", c, " µg/mL)"),
      subtitle = "Centrifuged vs Uncentrifuged, incl PBS and Trypsin controls",
      x = "Time (minutes)",
      y = "Mean relative fluorescence units (RFU)"
    )
    
    ggsave(
      filename = paste0("output/figures/FRET_", donor, "_", c, ".png"), 
      plot = plot
    )
  }
}


