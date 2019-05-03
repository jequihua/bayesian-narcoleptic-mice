
# Load packages.
library("tidyverse")
library("ggplot2")
library("brms")

# Load recorded behaviors data.
data <- read_csv("./data/videos/recorded_behaviors.csv")

### Duration of each BA (narcoleptic attack) model.
ba_seconds <- filter(data,behavior=="BA")

# Change variable types.
ba_seconds$group <- as.factor(ba_seconds$group)
ba_seconds$mouse <- as.factor(ba_seconds$mouse)

# Fit brms model for duration of each BA attack.
ba_seconds_fit <- brm(seconds ~ group + (1|mouse), 
                      data = ba_seconds, family = lognormal(), control = list(adapt_delta = 0.999))

# Evaluate model fit.
summary(ba_seconds_fit)

pp <- brms::pp_check(ba_seconds_fit)
pp + theme_bw()

brms::marginal_effects(ba_seconds_fit)

# Total time spent in BA state.
ba_totaltime <- filter(data,behavior=="BA") %>%
                mutate(mouse=paste0(group," ",mouse)) %>%
                group_by(mouse) %>%
                summarise(BA_totaltime = sum(seconds)) %>%
                separate(mouse,c("group","mouse"))

# Change variable types.
ba_totaltime$group <- as.factor(ba_totaltime$group)
ba_totaltime$mouse <- as.factor(ba_totaltime$mouse)

# Fit brms model for total time spent in BA state.
ba_totaltime_fit <- brm(BA_totaltime ~ group + (1|mouse), 
                        data = ba_totaltime, family = lognormal(), control = list(adapt_delta = 0.999,max_treedepth = 13))

# Evaluate model fit.
summary(ba_totaltime_fit)

plot(ba_totaltime_fit)

pp <- brms::pp_check(ba_totaltime_fit)
pp + theme_bw()

brms::marginal_effects(ba_totaltime_fit)

### Total number of BA (narcoleptic attacks) cases model per hour.

# Table available number of hours-data per mice.
mice_totalhours <- group_by(data,mouse) %>%
  summarise(totaltime = sum(seconds)/3600)

# Join data and mice_hours.
dataa <- right_join(data,mice_totalhours,by="mouse")

ba_counts <- dataa %>% 
             mutate(mouse=paste0(group," ",mouse)) %>%
             group_by(mouse) %>%
             summarise(BA_counts = sum(behavior == "BA"),hours=min(totaltime)) %>%
             mutate(BA_counts_hour=BA_counts/hours) %>%
             separate(mouse,c("group","mouse"))

# Change variable types.
ba_counts$group <- as.factor(ba_counts$group)
ba_counts$mouse <- as.factor(ba_counts$mouse)

# Fit brms model for absolute number of BA attacks per hour.            
ba_counts_fit <- brm(BA_counts_hour ~ group + (1|mouse), 
                     data = ba_counts, family = lognormal(), control = list(adapt_delta = 0.9999,max_treedepth=15))

# Evaluate model fit.
summary(ba_counts_fit)

plot(ba_counts_fit)

pairs(ba_counts_fit)

pp <- brms::pp_check(ba_counts_fit)
pp + theme_bw()

brms::marginal_effects(ba_counts_fit)

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
ba_totaltime$group <- as.factor(ba_totaltime$group)
ba_totaltime$mouse <- as.factor(ba_totaltime$mouse)

# Which priors can we specify for our model?
get_prior(percentage_time_BA ~ group + (1|mouse), data=ba_percentage)

prior <- c(set_prior("student_t(3, 5, 10)", class = "b"),
           set_prior("normal(0,1)", class = "b", coef = "groupsham"),
           set_prior("normal(0,1)", class = "b", coef = "grouptxcb"),
           set_prior("normal(0,1)", class = "b", coef = "grouptxox"),
           set_prior("student_t(3,0,  1)", class = "sd", 
                     group = "mouse", coef = "Intercept"),
           set_prior("student_t(3,0,  1)", class = "sd"))

# Fit brms model for total time spent in BA state.
ba_totaltime_fit <- brm(percentage_time_BA ~ group + (1|mouse), 
                        data = ba_percentage, family = Beta, control = list(adapt_delta = 0.999,max_treedepth = 15))

# Evaluate model fit.
summary(ba_totaltime_fit)

plot(ba_totaltime_fit)

pp = brms::pp_check(ba_totaltime_fit)
pp + theme_bw()

brms::marginal_effects(ba_totaltime_fit)
