################################################################################
#
# Author: Mark Calcott
# Date: 20200910
# Overview: Graph for testing background fluorescence of LOPAC compounds
#
################################################################################

# Load relevant packages
library(dplyr)
library(ggplot2)


############################# Load data and rearrange ###################### 

# Function for reading each plate and standardising data
read_plates <- function(replicate, plate){
  filename = paste0("DMSO_", replicate, "_", plate, ".csv")
  data <- read.csv(filename, header = FALSE, skip = 9, nrows = 96)
  data$RackPos <- data$V1
  data$RackNo <- plate 
  data$Replicate <- replicate 
  data$column <- substr(data$RackPos, 2,3)
  data$ratio <- data$V2
  data <- data %>%
    select(RackNo, RackPos, Replicate, column, ratio)
  #Get control values
  data_controls <- data %>% 
    filter(substr(RackPos, 2,3) == "12" | substr(RackPos, 2,3) == "01")
  
  mean_plate <- mean(data_controls$ratio)
  sd_plate <- sd(data_controls$ratio)
  
  #Add standard value
  data$standard <- (data$ratio-mean_plate)/sd_plate
  
  data
}

# Empty data frame to add each plate to
data <- data.frame()

# Go through each file and add to the dataframe
for(replicate in 1:2){
  for(plate in 1:16){
    data <- rbind(data, read_plates(replicate, plate))
  }
}


############################   Scatterplots  ############################
### Filter to remove controls
data <- data %>%
  filter(column != "01" & column != "12" )

### Rearrange data to correct format
sample_data <- data %>%
  filter(Replicate == 1) %>%
  select(-column, -ratio, -Replicate)
sample_data2 <- data %>%
  filter(Replicate == 2) %>%
  select(standard2 = standard)
sample_data <- cbind(sample_data, sample_data2)

### Rename columns, and add mean and standard deviation
names(sample_data)[names(sample_data) == "standard"] <- "Replicate_1"
names(sample_data)[names(sample_data) == "standard2"] <- "Replicate_2"
sample_data$mean <- rowMeans(sample_data[,c('Replicate_1', 'Replicate_2')])
sample_data$sd <- apply(sample_data[, c('Replicate_1', 'Replicate_2')],1,sd)

# Make a scatter graph for all data
ggplot(sample_data, aes(Replicate_1, Replicate_2)) +
  geom_point(alpha = 0.4) +
  labs(y = "Replicate 1 normalised", x = "Replicate 2 normalised") +
  geom_text(aes(label = "Dipyridamole"), x = -36.5, y = -43.7, 
            color = "black",
            hjust = 0) +
  geom_text(aes(label = "Anthrapyrazolone"), x = -12.1, y = -14.8, 
            color = "black",
            hjust = 0) +
  geom_text(aes(label = "NF 023"), x = -14, y = -5.5, 
            color = "black",
            hjust = 0) +
  theme_classic() 


