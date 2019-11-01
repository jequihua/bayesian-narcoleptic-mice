
# Load packages.
library("tidyverse")
library("ggplot2")
library("brms")
library("bayesplot")

# Load recorded behaviors data.
rem <- read_csv("./data/eeg/eeg_rem_epochs.csv")
names(rem)[2]="seconds"
norem <- read_csv("./data/eeg/eeg_norem_epochs.csv")
names(norem)[2]="seconds"
wake <- read_csv("./data/eeg/eeg_wake_epochs.csv")
names(wake)[2]="seconds"

### Total time in REM model.

# Poisson regression offset log(time)

# Table available number of seconds-data per mice.
mice_totalhours <- group_by(rem,mouse) %>%
                   summarise(totaltime = sum(seconds))

# Join data and mice_hours.
dataa <- right_join(rem,mice_totalhours,by="mouse")

rem_counts <- dataa[!duplicated(dataa[,c("mouse","group")]),]
rem_counts$group[rem_counts$group=="silvestre"]="WT"
rem_counts$group[rem_counts$group=="txsham"]="sham"

# Change variable types.
rem_counts$group <- as.factor(rem_counts$group)
rem_counts$mouse <- as.factor(rem_counts$mouse)

rem_counts$group = factor(as.character(rem_counts$group), levels=c("WT","NB","sham","txcb","txox"))

# Which priors can we specify for our model?
#get_prior(totaltime ~ group + (1|mouse), data=rem_counts)

# Fit brms model for time spent in REM state per hour.            
rem_counts_fit <- brm(totaltime ~ group, 
                     data = rem_counts, family = lognormal, control = list(adapt_delta = 0.9999,max_treedepth=15))


View(rem_counts)
# Chart raw data boxplots.
ggplot(rem_counts, aes(x=group, y=totaltime)) + geom_boxplot(lwd=1.2)+
  labs(x = NULL, y = "Total seconds")+
  ggtitle("REM")+
  theme_bw(base_size = 40)+
  xlab("Groups")

ggsave("./images5/1_rem_Tseconds/1_rem_seconds_rawdata_boxplots.png",
       width = 10, height = 8, dpi = 300,device="png")

# Summary stats. 
rem_counts_grouped = group_by(rem_counts,group)
rem_counts_summary = summarise(rem_counts_grouped,
                               n=n(),
                               n_mouse=n_distinct(mouse),
                               mean=mean(totaltime),
                               median = median(totaltime),
                               min=min(totaltime),
                               max=max(totaltime),
                               sd=sd(totaltime))

write_csv(rem_counts_summary,"./images5/1_rem_Tseconds/1_rem_seconds_rawdata_summary.csv")
summary(rem_counts_fit)
# Coefficient plot.
coeff_plot_data = interval_data_eeg(rem_counts_fit)
plot_intervals(coeff_plot_data)
ggsave("./images5/1_rem_Tseconds/2_rem_seconds_coefficientplot.png",
       width = 10, height = 8, dpi = 300,device="png")


# Evaluate model fit.
summary(rem_counts_fit)

stanplot(rem_counts_fit,pars = "^b_")

# PP plot
pp <- brms::pp_check(rem_counts_fit,nsamples=10)
pp + theme_bw()

# PP hist plot
yrep_poisson <- posterior_predict(ba_counts_fit, draws = 500)
ppc_hist(ba_counts$BA_counts, yrep_poisson[1:5, ],binwidth = 3)

mean(ba_counts$BA_counts == 0)
mean(yrep_poisson == 0)

brms::marginal_effects(rem_counts_fit,cex=5)


########################################################################################################


# Poisson regression offset log(time)

# Table available number of hours-data per mice.
mice_totalhours <- group_by(wake,mouse) %>%
  summarise(totaltime = sum(seconds))

# Join data and mice_hours.
dataa <- right_join(wake,mice_totalhours,by="mouse")

rem_counts <- dataa[!duplicated(dataa[,c("mouse","group")]),]
rem_counts$group[rem_counts$group=="silvestre"]="WT"
rem_counts$group[rem_counts$group=="txsham"]="sham"

# Change variable types.
rem_counts$group <- as.factor(rem_counts$group)
rem_counts$mouse <- as.factor(rem_counts$mouse)
rem_counts$group = factor(as.character(rem_counts$group), levels=c("WT","NB","sham","txcb","txox"))


# Fit brms model for absolute number of BA attacks per hour.            
rem_counts_fit <- brm(totaltime ~ group , 
                      data = rem_counts, family = lognormal,
                      control = list(adapt_delta = 0.9999,max_treedepth=15))

View(rem_counts)

# Chart raw data boxplots.
ggplot(rem_counts, aes(x=group, y=totaltime)) + geom_boxplot(lwd=1.2)+
  labs(x = NULL, y = "Total seconds")+
  ggtitle("Wake")+
  theme_bw(base_size = 40)+
  xlab("Groups")

ggsave("./images5/2_wake_Tseconds/1_wake_seconds_rawdata_boxplots.png",
       width = 10, height = 8, dpi = 300,device="png")

# Summary stats. 
rem_counts_grouped = group_by(rem_counts,group)
rem_counts_summary = summarise(rem_counts_grouped,
                               n=n(),
                               n_mouse=n_distinct(mouse),
                               mean=mean(totaltime),
                               median = median(totaltime),
                               min=min(totaltime),
                               max=max(totaltime),
                               sd=sd(totaltime))

write_csv(rem_counts_summary,"./images5/2_wake_Tseconds/1_wake_seconds_rawdata_summary.csv")
summary(rem_counts_fit)
# Coefficient plot.
coeff_plot_data = interval_data_eeg(rem_counts_fit)
plot_intervals(coeff_plot_data)
ggsave("./images5/2_wake_Tseconds/2_wake_seconds_coefficientplot.png",
       width = 10, height = 8, dpi = 300,device="png")


stanplot(rem_counts_fit,pars = "^b_")

# PP plot
pp <- brms::pp_check(rem_counts_fit,nsamples=10)
pp + theme_bw()

# PP hist plot
yrep_poisson <- posterior_predict(ba_counts_fit, draws = 500)
ppc_hist(ba_counts$BA_counts, yrep_poisson[1:5, ],binwidth = 3)

mean(ba_counts$BA_counts == 0)
mean(yrep_poisson == 0)

brms::marginal_effects(rem_counts_fit,cex=5)

dd[[1]]

########################################################################################################



########################################################################################################


# Poisson regression offset log(time)

# Table available number of hours-data per mice.
mice_totalhours <- group_by(norem,mouse) %>%
  summarise(totaltime = sum(seconds))

# Join data and mice_hours.
dataa <- right_join(norem,mice_totalhours,by="mouse")

dataa

rem_counts <- dataa[!duplicated(dataa[,c("mouse","group")]),]
rem_counts$group[rem_counts$group=="silvestre"]="WT"
rem_counts$group[rem_counts$group=="txsham"]="sham"
head(rem_counts)
# Change variable types.
rem_counts$group <- as.factor(rem_counts$group)
rem_counts$mouse <- as.factor(rem_counts$mouse)

rem_counts$group = factor(as.character(rem_counts$group), levels=c("WT","NB","sham","txcb","txox"))

# Fit brms model for absolute number of BA attacks per hour.            
rem_counts_fit <- brm(totaltime ~ group , 
                      data = rem_counts, family = lognormal, control = list(adapt_delta = 0.9999,max_treedepth=15))

# Evaluate model fit.
summary(rem_counts_fit)

View(rem_counts)
# Chart raw data boxplots.
ggplot(rem_counts, aes(x=group, y=totaltime)) + geom_boxplot(lwd=1.2)+
  labs(x = NULL, y = "Total seconds")+
  ggtitle("NOREM")+
  theme_bw(base_size = 40)+
  xlab("Groups")

ggsave("./images5/3_norem_Tseconds/1_norem_seconds_rawdata_boxplots.png",
       width = 10, height = 8, dpi = 300,device="png")

# Summary stats. 
rem_counts_grouped = group_by(rem_counts,group)
rem_counts_summary = summarise(rem_counts_grouped,
                               n=n(),
                               n_mouse=n_distinct(mouse),
                               mean=mean(totaltime),
                               median = median(totaltime),
                               min=min(totaltime),
                               max=max(totaltime),
                               sd=sd(totaltime))

write_csv(rem_counts_summary,"./images5/3_norem_Tseconds/1_norem_seconds_rawdata_summary.csv")
summary(rem_counts_fit)
# Coefficient plot.
coeff_plot_data = interval_data_eeg(rem_counts_fit)
plot_intervals(coeff_plot_data)
ggsave("./images5/3_norem_Tseconds/2_norem_Tseconds_coefficientplot.png",
       width = 10, height = 8, dpi = 300,device="png")



stanplot(rem_counts_fit,pars = "^b_")
?posterior_samples
# PP plot
pp <- brms::pp_check(rem_counts_fit,nsamples=10)
pp + theme_bw()

# PP hist plot
yrep_poisson <- posterior_predict(ba_counts_fit, draws = 500)
ppc_hist(ba_counts$BA_counts, yrep_poisson[1:5, ],binwidth = 3)

mean(ba_counts$BA_counts == 0)
mean(yrep_poisson == 0)

brms::marginal_effects(rem_counts_fit,cex=5)

dd[[1]]

########################################################################################################
