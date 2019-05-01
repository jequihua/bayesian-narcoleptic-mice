
# Load packages.
library("tidyverse")
library("ggplot2")
library("brms")

# Load recorded behaviors data.
data <- read_csv("./data/videos/recorded_behaviors.csv")

# Available number of hours-data per mice table.
mice_totalhours <- group_by(data,mouse) %>%
                   summarise(totaltime = sum(seconds)/3600)

# Join data and mice_hours.
data <- right_join(data,mice_totalhours,by="mouse")

### Total number of BA (narcoleptic attacks) cases model per hour.
ba_counts <- data %>% 
             mutate(mouse=paste0(group," ",mouse)) %>%
             group_by(mouse) %>%
             summarise(BA_counts = sum(behavior == "BA"),hours=min(totaltime)) %>%
             mutate(BA_counts_hour=BA_counts/hours) %>%
             separate(mouse,c("group","mouse"))

head(ba_counts)

# Change variable types.
ba_counts$group <- as.factor(ba_counts$group)
ba_counts$mouse <- as.factor(ba_counts$mouse)

# Fit brms model for absolute number of BA attacks per hour.            
ba_counts_fit <- brm(BA_counts_hour ~ group + (1|mouse), 
                     data = ba_counts, family = lognormal(), control = list(adapt_delta = 0.9999,max_treedepth=15))

# Evaluate model fit.
summary(ba_counts_fit)

plot(ba_counts_fit)

pp = brms::pp_check(ba_counts_fit)
pp + theme_bw()

brms::marginal_effects(ba_counts_fit)

### Percentage of time spent in BA state.

# Total hours of sampling by mouse.
mice_totalseconds <- group_by(data,mouse) %>%
                     summarise(totaltime = sum(seconds))

# Join data and mice_hours.
data <- right_join(data,mice_totalseconds,by="mouse")
head(data)
ba_percentage <- filter(data,behavior=="BA") %>%
                 mutate(mouse=paste0(group," ",mouse)) %>%
                 group_by(mouse) %>%
                 summarise(BA_totaltime = sum(seconds),hours=min(totaltime)) %>%
                 mutate(percentage_time_BA=(BA_totaltime/hours)*100) %>%
                 separate(mouse,c("group","mouse"))
View(ba_percentage)
# Change variable types.
ba_totaltime$group <- as.factor(ba_totaltime$group)
ba_totaltime$mouse <- as.factor(ba_totaltime$mouse)

# Fit brms model for total time spent in BA state.
ba_totaltime_fit <- brm(percentage_time_BA ~ group + (1|mouse), 
                        data = ba_percentage, family = lognormal(), control = list(adapt_delta = 0.9999,max_treedepth = 15))

# Evaluate model fit.
summary(ba_totaltime_fit)

plot(ba_totaltime_fit)

pairs(ba_totaltime_fit)

pp = brms::pp_check(ba_totaltime_fit)
pp + theme_bw()

brms::marginal_effects(ba_totaltime_fit)

