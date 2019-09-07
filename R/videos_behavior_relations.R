# Load packages.
library("tidyverse")
library("ggplot2")
library("brms")
library("brmstools")
library("bayesplot")
library("reshape2")

# Load misc. functions.
source("./R/misc_functions.R")

# Load recorded behaviors data.
data <- read_csv("./data/videos/recorded_behaviors.csv")

# Add "Before BA" flag. 
data$behavior_before = 0
for (i in 1:(nrow(data)-1)){
  if (data$behavior[i+1]=="BA"){
    data$behavior_before[i]=1
  }
  
}

# Select only cases where "Before BA" flag is 1.
data_subset = data[data$behavior_before==1,]

# Case counts.
cases = as.data.frame.matrix(table(data_subset$group,data_subset$behavior))
cases$group = row.names(cases)

# Write table to disk.
write_csv(cases,"./images4/4_freq_behav_before_ba/1_4_freq_behav_before_ba_table.csv")

# Melted cases data.frame.
melted_cases = melt(cases)

# Plot.
ggplot(melted_cases, aes(variable,value,fill=variable)) + 
       geom_bar(stat = "identity") + 
       facet_wrap(~ group)

