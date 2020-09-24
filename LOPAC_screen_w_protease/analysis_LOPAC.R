################################################################################
#
# Author: Mark Calcott
# Date: 20200910
# Overview: Graph for screening the LOPAC library in duplicate
#
################################################################################

# Load relevant packages
library(dplyr)
library(ggplot2)
library("readxl")
library(reshape2)

# Color scheme
no_inhibition <- "#0072B2"
complete_inhibition <- "#D55E00"


############################# Load data and rearrange ###################### 

# Function for reading each plate
read_plates <- function(replicate, plate){
  filename = paste0("Raw_data/SARS-CoV-2_", replicate, "_", plate, ".csv")
  data <- read.csv(filename, header = FALSE, skip = 9, nrows = 96)
  data$RackPos <- data$V1
  data$RackNo <- plate 
  data$Replicate <- replicate 
  data$column <- substr(data$RackPos, 2,3)
  data$ratio <- data$V2
  data <- data %>%
    select(RackNo, RackPos, Replicate, column, ratio)
  #Get neg value
  data_neg <- data %>% 
    filter(substr(RackPos, 2,3) == "12") %>%
    select(ratio) %>%
    colMeans()
  #Get pos value
  data_pos <- data %>% 
    filter(substr(RackPos, 2,3) == "01") %>%
    select(ratio) %>%
    colMeans()
  
  #Add relative inhibition value
  data$relative_inhibition <- (data$ratio-data_neg)/(data_pos-data_neg)*100
  
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


########### Calculate Z factor and graph positive/negative controls ############
# Select just the eight positive and eight negative controls from each plate
control_data <- data %>%
  filter(column == "01" | column == "12") %>%
  mutate(location = ceiling(seq(1, length(.$column))/2)) %>%
  mutate(Sample = factor(column, levels = c("01", "12"), labels = c("Positive", "Negative")))

# Calculate the mean and standard deviation for the controls
data_summary <- control_data %>%
  group_by(Sample) %>%
  summarise(average = mean(relative_inhibition), sd = sd(relative_inhibition))

# Calculate Z factor
sum_control_sd <- data_summary$sd[data_summary$Sample == "Positive"] + data_summary$sd[data_summary$Sample == "Negative"]
diff_control_mean <- data_summary$average[data_summary$Sample == "Positive"] - data_summary$average[data_summary$Sample == "Negative"]
z_factor <- 1-((3*sum_control_sd)/diff_control_mean)
z_factor_label <- paste("Z' = ", formatC(z_factor, digits = 2))

# Graph positives and negatives
ggplot(control_data, aes(location, relative_inhibition, color = Sample)) +
  geom_point(alpha = 0.6) + 
  geom_segment(
    x = 160, y = 10,
    xend = 160, yend = 85,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    size = 1, 
    arrow = arrow(length = unit(0.1, "inches"), ends = "both"),
    colour = "black" # Also accepts "red", "blue' etc
  ) +
  geom_text(aes(label = z_factor_label), x = 190, y = 50, color = "black") +
  labs(y = "Relative inhibition (%)", x = "Sample number") +
  theme_classic() +
  scale_color_manual(values=c(complete_inhibition, no_inhibition)) + 
  theme(legend.position = "none") +
  geom_text(aes(label = "No protease"), x = 80, y = 90, color = complete_inhibition) +
  geom_text(aes(label = "Protease"), x = 80, y = 10, color = no_inhibition)
  

############################   Calculate R2 values  ############################
### Filter to remove controls
data <- data %>%
  filter(column != "01" & column != "12" )

### Rearrange data to correct format
sample_data <- data %>%
  filter(Replicate == 1) %>%
  select(-column, -ratio, -Replicate)
sample_data2 <- data %>%
  filter(Replicate == 2) %>%
  select(relative_inhibition2 = relative_inhibition)
sample_data <- cbind(sample_data, sample_data2)

### Rename columns, and add mean and standard deviation
names(sample_data)[names(sample_data) == "relative_inhibition"] <- "Replicate_1"
names(sample_data)[names(sample_data) == "relative_inhibition2"] <- "Replicate_2"
sample_data$mean <- rowMeans(sample_data[,c('Replicate_1', 'Replicate_2')])
sample_data$sd <- apply(sample_data[, c('Replicate_1', 'Replicate_2')],1,sd)


### Calculate Pearson's correlation
r2 <- with(sample_data, cor(Replicate_1, Replicate_2, method = "pearson")^2)
r2 <- signif(r2, 2)

# Make a scatter graph
ggplot(sample_data, aes(Replicate_1, Replicate_2)) +
  geom_smooth(method='lm', se = FALSE, color = "grey") +
  geom_point(alpha = 0.4) +
  geom_text(aes(label = paste("R^2 ==", r2)), x = 20, y = 80, 
            color = "black", 
            parse = TRUE) +
  labs(y = "Replicate 1 relative inhibition (%)", x = "Replicate 2 relative inhibition (%)") +
  theme_classic() 

##############################  Make a QQ-plot  ################################

# QQ plot
ggplot(sample_data, aes(sample=mean)) +
  stat_qq(alpha = 0.4) +
  geom_hline(aes(yintercept = 19), linetype = "dashed", alpha = 0.8) +
  labs(y = "Mean relative inhibtion (%)", x = "Normal theoretical quantiles") +
  geom_text(aes(label = paste("Hit threshold")), x = -2, y = 23, 
            color = "grey15",
            size = 4) +
  theme_classic()



##################  Export replicate data combined with LOPAC  #################

lopacData <- read_excel("AcTable.xls")
lopacData$RackNo <- as.numeric(lopacData$RackNo)
export_data <- merge(sample_data, lopacData)
export_data <- export_data[order(-export_data$mean),]
rownames(export_data) <- NULL
write.csv(export_data, "LOPAC_screen_w_protease.csv")




