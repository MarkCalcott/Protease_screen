################################################################################
#
# Author: Mark Calcott
# Date: 20200910
# Overview: Graphing the emission spectra with individual points shown
#
################################################################################

# Load relevant packages
library(dplyr)
library(ggplot2)
library(reshape2)

# Color scheme
no_inhibition <- "#0072B2"
complete_inhibition <- "#D55E00"


############################### Load data ###################################### 

# Load file
read_files <- function(filename, loc){
  file <- read.csv(filename, header = FALSE)
  data <- file[loc,] %>%
    filter(V2 == "MeasB:Result" | (V2 == "MeasB:WavelengthEms" & V1 == "A01")) %>%
    select(-V2) %>%
    t %>%
    as.data.frame
  data[1,1] <- "emission"
  names(data) <- data[1,]
  
  data <- filter(data, emission != "emission")
  data <- mutate_all(data, function(x) as.numeric(as.character(x)))
}


############################# Rearrange data  ################################## 

# Read data
data_434 <- read_files("Emission_spectra_18h.csv", loc = 6:150)

# Dividing by reading at 477 nm
data_434_norm <- mutate_at(data_434, 2:49, function(x) x/x[23])

# Rearranging reads to a single column, adding replicate/treatment
data_434_norm <- data_434_norm %>%
  melt(id = "emission") %>%
  mutate(replicate = substr(variable, 3, 3)) %>%
  rename(well = variable) %>%
  mutate(treatment = ifelse(replicate %in% c(1, 3, 5), "Protease", "No protease"))

# Average data points for each replicate
replicate_data_points <- data_434_norm %>%
  group_by(treatment, replicate, emission) %>%
  summarise(average = mean(value))

# Average data across treatment groups
treatment_average <- data_434_norm %>%
  group_by(treatment, emission) %>%
  summarise(average = mean(value))

# Draw graph
ggplot(treatment_average, aes(emission, average, fill = treatment)) +
  geom_line() +
  geom_point(data = replicate_data_points, aes(emission, average, colour = treatment), alpha = 0.4) +
  labs(y = "Emission relative to 477 nm", x = "Wavelength for emission") +
  theme_classic() +
  scale_color_manual(values=c(complete_inhibition, no_inhibition)) +
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_blank())

ggsave("emission.svg", units = "in", width = 3, height = 3)

       