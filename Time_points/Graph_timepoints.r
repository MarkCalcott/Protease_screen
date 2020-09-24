################################################################################
#
# Author: Mark Calcott
# Date: 20200910
# Overview: Graphing the ratio of emission between Venus and eCFP
#           at 30 min timepoints
#
################################################################################

# Load relevant packages
library(dplyr)
library(ggplot2)
library(reshape2)

#  Color scheme
no_inhibition <- "#0072B2"
complete_inhibition <- "#D55E00"

############################# Load data and rearrange ########################## 

# Load file
read_files <- function(filename){
  file <- read.csv(filename, header = TRUE, skip = 43)
  
  # Rearrange data
  data <- file %>%
    filter(Repeat %in% c("MeasA:Result", "MeasB:Result")) %>%
    melt(id = c("Well", "Repeat")) %>%
    mutate(value = as.numeric(value), variable = (as.numeric(variable)-1)/2) %>%
    mutate(replicate = as.numeric(substr(Well, 2, 3))) %>%
    mutate(treatment = ifelse(replicate %in% c(1, 3, 5), "Protease", "No protease")) %>%
    filter(replicate %in% 1:6 & !is.na(value))
    
  # Calculate ratio
  data_A <- data %>%
    filter(Repeat == "MeasA:Result")
  data_B <- data %>%
    filter(Repeat == "MeasB:Result")
  
  data_A$ratio <- data_A$value/data_B$value
  data_A
}



############################# Rearrange data  ################################## 

# Read data
data <- read_files("Measurements_over_time.csv")

# Average data points for each replicate
replicate_data_points <- data %>%
  group_by(treatment, replicate, variable) %>%
  summarise(average = mean(ratio))

# Average data across treatment groups
treatment_average <- data %>%
  group_by(treatment, variable) %>%
  summarise(average = mean(ratio))


# Draw graph
ggplot(treatment_average, aes(variable, average, fill = treatment)) +
  geom_line() +
  geom_point(data = replicate_data_points, aes(variable, average, colour = treatment), alpha = 0.4) +
  labs(y = "Ratio of emission at 528 nm to 477 nm", x = "Time (hours)") +
  theme_classic() +
  geom_vline(aes(xintercept = 4), alpha = 0.4, linetype = "dashed") +
  scale_color_manual(values=c(complete_inhibition, no_inhibition)) + 
  theme(legend.position = "none") +
  geom_text(aes(label = "No protease"), x = 8, y = 4, color = complete_inhibition) +
  geom_text(aes(label = "Protease"), x = 8, y = 2, color = no_inhibition)
