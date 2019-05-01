
# Load packages.
library("tidyverse")
library("ggplot2")
library("brms")

# Load recorded behaviors data.
data <- read_csv("./data/videos/recorded_behaviors.csv")

### Total number of BA (narcoleptic attacks) cases model.
ba_counts <- data %>% 
             mutate(mouse=paste0(group," ",mouse)) %>%
             group_by(mouse) %>%
             summarise(BA_counts = sum(behavior == "BA")) %>%
             separate(mouse,c("group","mouse"))
             
             
# Change variable types.
ba_counts$group <- as.factor(ba_counts$group)
ba_counts$mouse <- as.factor(ba_counts$mouse)

# Fit brms model for absolute number of BA attacks.            
ba_counts_fit <- brm(BA_counts ~ group + (1|mouse), 
                     data = ba_counts, family = poisson(), control = list(adapt_delta = 0.999))

# Evaluate model fit.
summary(ba_counts_fit)

plot(ba_counts_fit)

pp = brms::pp_check(ba_counts_fit)
pp + theme_bw()

brms::marginal_effects(ba_counts_fit)

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

pp = brms::pp_check(ba_seconds_fit)
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
ba_totaltime_fit <- brm(BA_totaltime ~ group + (group | mouse), 
                        data = ba_totaltime, family = lognormal(), control = list(adapt_delta = 0.999,max_treedepth = 13))

# Evaluate model fit.
summary(ba_totaltime_fit)

plot(ba_totaltime_fit)

pairs(ba_totaltime_fit)

pp = brms::pp_check(ba_totaltime_fit)
pp + theme_bw()

brms::marginal_effects(ba_totaltime_fit)
