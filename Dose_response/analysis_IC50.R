################################################################################
#
# Author: Mark Calcott
# Date: 20200910
# Overview: Dose response curves for top hits
#
################################################################################

# Load relevant packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(broom)
library("readxl")
library(tidyr)
library(drc)

############################# Load data and rearrange ##########################
read_plates_EC50 <- function(treatment, replicate, plate){
  filename <- paste0("SARS-CoV-2_", treatment, "_", replicate, "_", plate, ".csv")
  data <- read.csv(filename, header = FALSE, skip = 9, nrows = 96)
  data$plate <- plate
  data$treatment <- treatment
  data$replicate <- replicate 
  data$well <- data$V1
  data$sample <- factor(as.numeric(substr(data$V1, 2,3)))
  data$dilution <- factor(substr(data$V1, 1,1), labels = c(20, 4, 0.8, 0.16, 0.032, 0.0064, 0.00128, 0))
  
  data$ratio <- data$V2
  data <- data %>%
    dplyr::select(plate, well, treatment, replicate, sample, dilution, ratio)

  data_summary <- data %>%
    group_by(dilution) %>%
    summarise(well, inhibition = (ratio-ratio[sample == 12])/
                (ratio[sample == 1]-ratio[sample == 12])*100)
  
  data$dilution <- as.numeric(as.character(data$dilution))
  data$inhibition <- data_summary$inhibition
  
  data
}


# Empty data frame to add each plate to
data <- data.frame()

# Go through each file and add to the dataframe
for(treatment in c("EC50", "TX50")){
  for(replicate in 1:3){
    for(plate in c(1:8)){
      data <- rbind(data, read_plates_EC50(treatment, replicate, plate))
    }
  }
}

############################ Add names of compounds ############################

naming_file <- read.csv("EC50_naming.csv")
lopacData <- read_excel("../LOPAC_screen_w_protease/AcTable.xls") 
lopacData$RackNo <- as.numeric(lopacData$RackNo)

sample_names <- merge(lopacData, naming_file) %>%
  dplyr::select(plate, sample, name)

sample_data <- merge(data, sample_names)



############################## Calculate EC50s ###############################

# Create a table with all the EC50 values
EC50_table <- data.frame()
for(name in unique(sample_data$name)){
  for(treatment in c("EC50", "TX50")){
    EC50_entry <- data.frame(name = rep(name, 3),
                             treatment = rep(treatment, 3),
                             replicate = 1:3)
    sample <- sample_data[sample_data$name == name & 
                            sample_data$treatment == treatment,]
    model <- drm(inhibition~dilution, replicate, data=sample, fct = LL.4())
    EC50_entry$EC50 <- model$coefficients[startsWith(names(model$coefficients), "e")]
    EC50_table <- rbind(EC50_table, EC50_entry)
  }
}

# Calculate the EC50 mean and standard deviation for each sample
data_summary <- EC50_table %>%
  group_by(name, treatment) %>%
  summarise(EC50.mean = mean(EC50), EC50.sd = sd(EC50))

# Sorted data to identify the best hits
sorted <-  data_summary %>%
  filter(treatment == "TX50") %>%
  arrange(EC50.mean)
sorted$rank = rownames(sorted)


############################## Create EC50 plots ###############################

sample_data$name <- factor(sample_data$name, levels = sorted$name)
levels(sample_data$name)
# Set the lowest dilution to 0.0001 for display purposes
#sample_data$dilution[sample_data$dilution < 0.00120] <- 0.00001
# Draw dose response curves
ggplot(sample_data, aes(x = dilution, y = inhibition, color = treatment)) +
  geom_point(alpha = 0.5) +
  stat_smooth(
    method = "drm",
    method.args = list(
      fct = LL.4(), logDose = 10
    ),
    se = FALSE
  ) +
  scale_x_log10() +
  labs(x = "Compound concentration (log (\u00b5M))", y = "Relative inhibition (%)") +
  facet_wrap(vars(as.numeric(name))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.position = "none")

################# Create EC50 plots of the top 20 hits #########################
# Select the hits
top_hits <- sample_data %>%
  filter(name %in% sorted$name[1:22]) %>%
  filter(!(name %in% c("Disulfiram", "SCH-202676")))

# Draw plot
ggplot(top_hits, aes(x = dilution, y = inhibition, color = treatment)) +
  geom_point(alpha = 0.5) +
  stat_smooth(
    method = "drm",
    method.args = list(
      fct = LL.4(), logDose = 10
    ),
    se = FALSE
  ) +
  scale_x_log10() +
  labs(x = "Compound concentration (log (\u00b5M))", y = "Relative inhibition (%)") +
  facet_wrap(vars(name), ncol = 4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.position = "none")

######################### Make table of EC50 values ############################
with_protease <- data_summary %>%
  filter(treatment == "TX50") %>%
  rename(with_mean = EC50.mean,
         with_sd = EC50.sd) %>%
  dplyr::select(-treatment)

without_protease <- data_summary %>%
  filter(treatment == "EC50") %>%
  rename(without_mean = EC50.mean,
         without_sd = EC50.sd) %>%
  dplyr::select(-treatment)

export_data <- merge(with_protease, without_protease) %>%
  mutate(ratio = with_mean/without_mean)

export_data <- export_data[order(export_data$with_mean),]
rownames(export_data) <- NULL
write.csv(export_data, "EC50_table.csv")
