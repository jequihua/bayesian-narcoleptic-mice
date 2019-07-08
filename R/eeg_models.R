
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

data = bind_rows(list(rem,norem,wake))

names(data)[2]="seconds"

### Total time in REM model.

# Poisson regression offset log(time)

# Table available number of hours-data per mice.
mice_totalhours <- group_by(rem,mouse) %>%
                   summarise(totaltime = sum(seconds))

# Join data and mice_hours.
dataa <- right_join(rem,mice_totalhours,by="mouse")

rem_counts <- dataa[!duplicated(dataa[,c("mouse","group")]),]

# Change variable types.
rem_counts$group <- as.factor(rem_counts$group)
rem_counts$mouse <- as.factor(rem_counts$mouse)

rem_counts$group = factor(as.character(rem_counts$group), levels=c("silvestre","NB","txsham","txcb","txox"))

rem_counts$totaltime = rem_counts$totaltime/60

# Which priors can we specify for our model?
get_prior(totaltime ~ group + (1|mouse), data=rem_counts)

prior <- c(set_prior("student_t(3, 5, 10)", class = "b"),
           set_prior("normal(0,1)", class = "b", coef = "grouptxsham"),
           set_prior("normal(0,1)", class = "b", coef = "grouptxcb"),
           set_prior("normal(0,1)", class = "b", coef = "grouptxox"),
           set_prior("student_t(3,0,  1)", class = "sd", 
                     group = "mouse", coef = "Intercept"),
           set_prior("student_t(3,0,  1)", class = "sd"))

# Fit brms model for absolute number of BA attacks per hour.            
rem_counts_fit <- brm(totaltime ~ group + (1|mouse), 
                     data = rem_counts, family = Gamma(), control = list(adapt_delta = 0.9999,max_treedepth=15),
                     prior=prior)


View(rem_counts)

library("ggpubr")
ggboxplot(rem_counts, x = "group", y = "totaltime", 
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = unique(rem_counts$group),
          ylab = "totaltime", xlab = "Treatment")

kruskal.test(totaltime ~ group, data = rem_counts)

rem_counts_aux = rem_counts[rem_counts$group=="silvestre"|rem_counts$group=="NB",]
res.aov <- aov(totaltime ~ group, data = rem_counts_aux)

View(rem_counts_aux)

summary(res.aov)

TukeyHSD(res.aov)

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

# Change variable types.
rem_counts$group <- as.factor(rem_counts$group)
rem_counts$mouse <- as.factor(rem_counts$mouse)

rem_counts$group = factor(as.character(rem_counts$group), levels=c("silvestre","NB","txsham","txcb","txox"))

# Fit brms model for absolute number of BA attacks per hour.            
rem_counts_fit <- brm(seconds ~ group + (1|mouse), 
                      data = rem_counts, family = exponential(), control = list(adapt_delta = 0.9999,max_treedepth=15))

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

dd = brms::marginal_effects(rem_counts_fit,cex=5)

dd[[1]]

########################################################################################################

### Percentage of total sampled time spent in BA state.

# Total hours of sampling by mouse.
mice_totalseconds <- group_by(data,mouse) %>%
                     summarise(totaltime = sum(seconds))

View(mice_totalseconds)

# Join data and mice_hours.
dataa <- right_join(data,mice_totalseconds,by="mouse")

ba_percentage <- filter(dataa,behavior==1) %>%
                 mutate(mouse=paste0(group," ",mouse)) %>%
                 group_by(mouse) %>%
                 summarise(BA_totaltime = sum(seconds),hours=min(totaltime)) %>%
                 mutate(percentage_time_BA=(BA_totaltime/hours)) %>%
                 separate(mouse,c("group","mouse"))

# Change variable types.
ba_percentage$group <- as.factor(ba_percentage$group)
ba_percentage$mouse <- as.factor(ba_percentage$mouse)

# Which priors can we specify for our model?
get_prior(percentage_time_BA ~ group + (1|mouse), data=ba_percentage)

# Fit brms model for total time spent in BA state.
ba_totaltime_fit <- brm(percentage_time_BA ~ group + (1|mouse), 
                        data = ba_percentage, family = Beta, control = list(adapt_delta = 0.999,max_treedepth = 15))

# Evaluate model fit.
summary(ba_totaltime_fit)

plot(ba_totaltime_fit)

stanplot(ba_totaltime_fit,pars = "^b_")

?stanplot

pp = brms::pp_check(ba_totaltime_fit,nsamples=10)
pp + theme_bw()

brms::marginal_effects(ba_totaltime_fit)
