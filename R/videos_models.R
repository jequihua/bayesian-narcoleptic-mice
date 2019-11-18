
# Load packages.
library("tidyverse")
library("ggplot2")
library("brms")
library("brmstools")
library("bayesplot")
library("emmeans")
library("tidybayes")
library("ggmcmc")
library("coda")

# Load misc. functions.
source("./R/misc_functions.R")

# Load recorded behaviors data.
data <- read_csv("./data/videos/recorded_behaviors.csv")

### Duration of each BA (narcoleptic attack) model.
ba_seconds <- filter(data,behavior=="BA")

# Change variable types.
ba_seconds$group <- as.factor(ba_seconds$group)
ba_seconds$mouse <- as.factor(ba_seconds$mouse)

# Which priors can we specify for our model?
get_prior(seconds ~ group + (1|mouse), data=ba_seconds)

# Fit brms model for duration of each BA attack.
ba_seconds_fit <- brm(seconds ~ group + (1|mouse), 
                      data = ba_seconds,
                      family = lognormal,
                      control = list(adapt_delta = 0.999))

summary(ba_seconds_fit)

mcmc_df = ggs(as.mcmc(ba_seconds_fit))
ggs_geweke(mcmc_df,frac1 = 0.2, frac2 = 0.6)

?ggs_geweke

plot(ba_seconds_fit)

modelposterior <- as.mcmc(ba_seconds_fit)

geweke.plot(modelposterior)


# Chart raw data boxplots.
ggplot(ba_seconds, aes(x=group, y=seconds)) + geom_boxplot(lwd=1.2)+
labs(x = NULL, y = "Seconds")+
ggtitle("Duration of each BA")+
theme_bw(base_size = 40)+
xlab("Groups")
 
ggsave("./images4/1_ba_duration/1_ba_duration_rawdata_boxplots.png",
       width = 10, height = 8, dpi = 500,device="png")

# Summary stats. 
ba_seconds_grouped = group_by(ba_seconds,group)
ba_seconds_summary = summarise(ba_seconds_grouped,
                               n=n(),
                               n_mouse=n_distinct(mouse),
                               mean=mean(seconds),
                               median = median(seconds),
                               min=min(seconds),
                               max=max(seconds),
                               sd=sd(seconds),
                               Q0.25=quantile(seconds,0.25),
                               Q0.75=quantile(seconds,0.75),
                               IQR=IQR(seconds))

write_csv(ba_seconds_summary,"./images4/1_ba_duration/1_ba_duration_rawdata_summary.csv")

# Coefficient plot.
coeff_plot_data = interval_data(ba_seconds_fit)
plot_intervals(coeff_plot_data)
ggsave("./images4/1_ba_duration/2_ba_duration_modelfit_coefficientplot.png",
       width = 10, height = 8, dpi = 500,device="png")

# Summary of model fit.
fit_summary = lazerhawk::brms_SummaryTable(ba_seconds_fit,astrology=TRUE)
write_csv(fit_summary,"./images4/1_ba_duration/2_ba_duration_modelfit_fitsummary.csv")

# Contrasts 
credible_intervals_dat = credible_intervals_data(ba_seconds,ba_seconds_fit)

write_csv(credible_intervals_dat,"./images4/1_ba_duration/3_ba_duration_modelinference_pairwisecomp_crudetable.csv")

plot_contrasts(credible_intervals_dat,xlab="Effect size")
ggsave("./images4/1_ba_duration/3_ba_duration_modelinference_pairwisecomp_crude.png",
       width = 10, height = 8, dpi = 500,device="png")


##################################


# Evaluate model fit.
summary(ba_seconds_fit)

plot(ba_seconds_fit)

# Predictive posterior checks.
pp <- brms::pp_check(ba_seconds_fit,nsamples=20)
pp + theme_bw()

# PP hist plot.
y_simulated <- posterior_predict(ba_seconds_fit, draws = 500)
ppc_hist(ba_seconds$seconds, y_simulated[1:11,], binwidth = 4)

ppc_stat_freqpoly_grouped(ba_seconds$seconds, y_simulated, 
                          group=ba_seconds$group, stat = "median",
                          binwidth = 4,freq=FALSE)+
                          coord_cartesian(xlim = c(10, 70))

head(ba_seconds)

ba_seconds%>%
  add_residual_draws(ba_seconds_fit) %>%
  ggplot(aes(x = .row, y = .residual)) +
  stat_pointinterval()

# Conditional Predictions.
brms::marginal_effects(ba_seconds_fit)

### Total number of BA (narcoleptic attacks) cases model per hour.

# Poisson regression offset log(time)

# Table available number of hours-data per mice.
mice_totalhours <- group_by(data,mouse) %>%
                   summarise(totaltime = sum(seconds)/3600)

# Join data and mice_hours.
dataa <- right_join(data,mice_totalhours,by="mouse")


ba_counts <- dataa %>% 
             mutate(mouse=paste0(group," ",mouse)) %>%
             group_by(mouse) %>%
             summarise(BA_counts = sum(behavior == "BA"),hours=min(totaltime)) %>%
             separate(mouse,c("group","mouse"))

# Change variable types.
ba_counts$group <- as.factor(ba_counts$group)
ba_counts$mouse <- as.factor(ba_counts$mouse)

# Which priors can we specify for our model?
get_prior(BA_counts ~ group + (1|mouse) + offset(log(hours)), data=ba_counts)

# Fit brms model for number of BA attacks per hour.            
ba_counts_fit <- brm(BA_counts ~ group + offset(log(hours)), 
                     data = ba_counts, family = poisson(), 
                     control = list(adapt_delta = 0.9999,max_treedepth=15))


head(ba_counts)

posterior_data = contrasts_data(ba_counts,ba_counts_fit)

head(posterior_data)

# Chart raw data boxplots.
ggplot(ba_counts, aes(x=group, y=BA_counts)) + geom_boxplot(lwd=1.2)+
  labs(x = NULL, y = "counts")+
  ggtitle("BA counts")+
  theme_bw(base_size = 40)+
  xlab("Groups")

ggsave("./images4/2_ba_counts/1_ba_counts_rawdata_boxplots.png",
       width = 10, height = 8, dpi = 500,device="png")

# Summary stats. 
ba_counts_grouped = group_by(ba_counts,group)
ba_counts_summary = summarise(ba_counts_grouped,
                               n=n(),
                               n_mouse=n_distinct(mouse),
                               mean=mean(BA_counts),
                               median = median(BA_counts),
                               min=min(BA_counts),
                               max=max(BA_counts),
                               sd=sd(BA_counts),
                               Q0.25=quantile(BA_counts,0.25),
                               Q0.75=quantile(BA_counts,0.75),
                               IQR=IQR(BA_counts))

write_csv(ba_counts_summary,"./images4/2_ba_counts/1_ba_counts_rawdata_summary.csv")

summary(ba_counts_fit)

# Coefficient plot.
coeff_plot_data = interval_data(ba_counts_fit)
plot_intervals(coeff_plot_data)


ggsave("./images4/2_ba_counts/2_ba_counts_modelfit_coefficientplot.png",
       width = 10, height = 8, dpi = 500,device="png")

# Summary of model fit.
fit_summary = lazerhawk::brms_SummaryTable(ba_counts_fit,astrology=TRUE)
write_csv(fit_summary,"./images4/2_ba_counts/2_ba_counts_modelfit_fitsummary.csv")

head(ba_counts)

credible_intervals_dat = credible_intervals_data_hours(ba_counts,ba_counts_fit)
write_csv(credible_intervals_dat,"./images4/2_ba_counts/3_ba_counts_modelinference_pairwisecomp_crudetable.csv")

plot_contrasts(credible_intervals_dat,xlab="Effect size")
ggsave("./images4/2_ba_counts/3_ba_duration_modelinference_pairwisecomp_crude.png",
       width = 10, height = 8, dpi = 500,device="png")

# Evaluate model fit.
summary(ba_counts_fit)

plot(ba_counts_fit)

stanplot(ba_counts_fit,pars = "^b_")

# PP plot
pp <- brms::pp_check(ba_counts_fit,nsamples=15)
pp + theme_bw()

# PP hist plot
yrep_poisson <- posterior_predict(ba_counts_fit, draws = 1000)
ppc_hist(ba_counts$BA_counts, yrep_poisson[1:11, ],binwidth = 4)

# PP histogram plot (observed vs simulated values).
y_simulated <- posterior_predict(ba_counts_fit, draws = 1000)
ppc_hist(ba_counts$BA_counts, y_simulated[1:11,], binwidth = 4)

# PP frecuency polygons vs the observed median of each group.
ppc_stat_freqpoly_grouped(ba_counts$BA_counts, y_simulated, 
                          group=ba_counts$group, stat = "median", binwidth = 4)+
  coord_cartesian(xlim = c(0, 30))





mean(ba_counts$BA_counts == 0)
mean(yrep_poisson == 0)

brms::marginal_effects(ba_counts_fit,cex=2)

### Percentage of total sampled time spent in BA state.

# Total hours of sampling by mouse.
mice_totalseconds <- group_by(data,mouse) %>%
  summarise(totaltime = sum(seconds))

# Join data and mice_hours.
dataa <- right_join(data,mice_totalseconds,by="mouse")

ba_percentage <- filter(dataa,behavior=="BA") %>%
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
ba_totaltime_fit <- brm(percentage_time_BA ~ group , 
                        data = ba_percentage, 
                        family = Beta, 
                        control = list(adapt_delta = 0.999,max_treedepth = 15))

summary(ba_totaltime_fit)

# Chart raw data boxplots.
ggplot(ba_percentage, aes(x=group, y=percentage_time_BA)) + geom_boxplot(lwd=1.2)+
  labs(x = NULL, y = "% of time")+
  ggtitle("% of time spent in BA")+
  theme_bw(base_size = 40)+
  xlab("Groups")

ggsave("./images4/3_ba_percentage/1_ba_counts_rawdata_boxplots.png",
       width = 10, height = 8, dpi = 500,device="png")

# Summary stats. 
ba_percentage_grouped = group_by(ba_percentage,group)
ba_percentage_summary = summarise(ba_percentage_grouped,
                              n=n(),
                              n_mouse=n_distinct(mouse),
                              mean=mean(percentage_time_BA),
                              median = median(percentage_time_BA),
                              min=min(percentage_time_BA),
                              max=max(percentage_time_BA),
                              sd=sd(percentage_time_BA),
                              Q0.25=quantile(percentage_time_BA,0.25),
                              Q0.75=quantile(percentage_time_BA,0.75),
                              IQR=IQR(percentage_time_BA))

write_csv(ba_percentage_summary,"./images4/3_ba_percentage/1_ba_percentage_rawdata_summary.csv")

ba_percentage_fit = ba_totaltime_fit

# Coefficient plot.
coeff_plot_data = interval_data(ba_percentage_fit)
plot_intervals(coeff_plot_data)
ggsave("./images4/3_ba_percentage/2_ba_percentage_modelfit_coefficientplot.png",
       width = 10, height = 8, dpi = 500,device="png")

# Summary of model fit.
fit_summary = lazerhawk::brms_SummaryTable(ba_percentage_fit,astrology=TRUE)
write_csv(fit_summary,"./images4/3_ba_percentage/2_ba_percentage_modelfit_fitsummary.csv")

# Contrasts 
#copms = tukeys_mat(ba_seconds,ba_seconds_fit)
#write_csv(comps,"./images4/3_ba_percentage/3_ba_duration_modelinference_pairwisecomp_crudetable.csv")
#plot_contrasts(comps)

credible_intervals_dat = credible_intervals_data(ba_percentage,ba_percentage_fit)

write_csv(credible_intervals_dat,"./images4/3_ba_percentage/3_ba_percentage_modelinference_pairwisecomp_crudetable.csv")

plot_contrasts(credible_intervals_dat,xlab="Effect size")
ggsave("./images4/3_ba_percentage/3_ba_percentage_modelinference_pairwisecomp_crude.png",
       width = 10, height = 8, dpi = 500,device="png")


###############################################################################################################3

comps = tukeys_mat(ba_percentage,ba_totaltime_fit)
comps[,2:5]=-1*comps[,2:5]

write_csv(comps,"./images4/3_ba_percentage/3_ba_percentage_modelinference_pairwisecomp_table.csv")
plot_contrasts(comps)
ggsave("./images4/3_ba_percentage/3_ba_percentage_modelinference_pairwisecomp.png",
       width = 10, height = 8, dpi = 500,device="png")

# Evaluate model fit.
summary(ba_totaltime_fit)

plot(ba_totaltime_fit)

stanplot(ba_totaltime_fit,pars = "^b_")


pp = brms::pp_check(ba_totaltime_fit,nsamples=15)
pp + theme_bw()

# PP histogram plot (observed vs simulated values).
y_simulated <- posterior_predict(ba_totaltime_fit, draws = 500)
ppc_hist(ba_percentage$percentage_time_BA, y_simulated[1:11,], binwidth = 0.02)

# PP frecuency polygons vs the observed median of each group.
ppc_stat_freqpoly_grouped(ba_percentage$percentage_time_BA, y_simulated, 
                          group=ba_percentage$group, stat = "median", binwidth = 0.02)+
                          coord_cartesian(xlim = c(-0.01, 0.1))


brms::marginal_effects(ba_totaltime_fit)

