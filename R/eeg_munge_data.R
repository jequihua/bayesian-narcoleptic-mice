
# Load packages.
library("tidyverse")
library("ggplot2")
library("brms")
library("bayesplot")
library("readr")

# Load recorded behaviors data.
data <- read_csv("./data/eeg/stacked_eeg_slps.csv")

# Set timestamp and time fields 
data$time = as.POSIXct(data$time,format="%H:%M:%S")
data$timestamp = as.POSIXct(data$timestamp,format="%m/%d/%Y - %H:%M:%S")

# Note on Behavior labels.
# 3 - Wake
# 2 - No-REM
# 1 - REM

# Mice.
mice = unique(data$mouse)

mice_list = list()

for (i in 1:length(mice)){

  mouse = data[data$mouse==mice[i],]
  counter=0
  
  # Filter out unwanted hours.
  mouse = mouse[mouse$timestamp>=paste0(min(as.Date(mouse$timestamp, "%m/%d/%Y",tz=""))," 18:00:00")&
                mouse$timestamp<=paste0(max(as.Date(mouse$timestamp, "%m/%d/%Y",tz=""))," 02:00:00"),]
  
  rep_df = data.frame(behavior = mouse$Behavior, counts= sequence(rle(as.character(mouse$Behavior))$lengths))

  rep_df$mouse=mouse$mouse[1]
  
  # 1 silvestre
  # 2 txcb
  # 3 txox
  # 4 NB
  # 5 txsham
  
  if(mouse$group[1]==1){
    rep_df$group="silvestre"
  }
  else if (mouse$group[1]==2) {
    rep_df$group="txcb"
  }
  else if (mouse$group[1]==3) {
    rep_df$group="txox"
  }
  else if (mouse$group[1]==4) {
    rep_df$group="NB"
  }
  else if (mouse$group[1]==5) {
    rep_df$group="txsham"
  }
  
  idx_list=list()
  for (j in 1:(nrow(rep_df)-1)){
    if(rep_df$behavior[j]!=1 | 
       (rep_df$behavior[j]==1&rep_df$behavior[j+1]==1)){
      counter=counter+1
      idx_list[[counter]]=j
    }
  }
  
  idx_list = -1*unlist(idx_list)
  
  rep_df = rep_df[idx_list,]
  
  rep_df$counts = rep_df$counts*12
  
  rep_df = rep_df[rep_df$behavior==1,]
  
  mice_list[[i]]=rep_df
  
}

binded_rep_dfs = bind_rows(mice_list)

View(binded_rep_dfs)

write_excel_csv(binded_rep_dfs,"D:/repositories/bayesian-narcoleptic-mice/data/eeg/eeg_rem_epochs.csv")



