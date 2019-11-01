
# Create a Tukeys contrast matrix.
tukeys_mat = function(model_input_data,brms_model_fit,columns=1:4,estimate.method="median"){
  mcmc = brms_model_fit
  coefs <- as.matrix(mcmc)[, columns]
  newdata <- data.frame(x = levels(model_input_data$group))
  tuk.mat <- contrMat(n = table(newdata$x), type = "Tukey")
  Xmat <- model.matrix(~x, data = newdata)
  pairwise.mat <- tuk.mat %*% Xmat
  mcmc_areas(coefs %*% t(pairwise.mat))
  (comps = broom::tidyMCMC(coefs %*% t(pairwise.mat), conf.int = TRUE, conf.method = "HPDinterval",estimate.method=estimate.method))
  comps$term = factor(comps$term, levels = rev(unique(comps$term)))
  return(comps)
}

mcmc <- ba_seconds_fit
head(as.matrix(mcmc))
coefs <- as.matrix(mcmc)[, columns]



ba_seconds_fit %>%
  spread_draws(b_group) %>%
  compare_levels(r_condition, by = condition) %>%
  ggplot(aes(y = condition, x = r_condition)) +
  stat_halfeyeh()

get_variables(ba_seconds_fit)

summary(ba_seconds_fit)


# Conditional Predictions.
brms::marginal_effects(ba_seconds_fit)

?marginal_effects



(hyp2 <- hypothesis(ba_seconds_fit, c("group + (1|mouse) = 0", "Intercept = 0"),
                    alpha = 0.05))


fixed_form <- brms:::fix_factor_contrasts(ba_seconds_fit$formula)$fixed


mm <- brms:::get_model_matrix(fixed_form, data = fit$data) 

population_level_effects <- fixef(ba_seconds_fit)

mm <- brms:::get_model_matrix(population_level_effects, data = population_level_effects$data)

?get_model

?posterior_predict


ba_seconds$group

data_A1 <- ba_seconds[ba_seconds$group == levels(ba_seconds$group)[1], ]
PPD_A1 <- posterior_predict(ba_seconds_fit, newdata = data_A1)

data_A2 <- ba_seconds[ba_seconds$group == levels(ba_seconds$group)[4], ]
PPD_A2 <- posterior_predict(ba_seconds_fit, newdata = data_A2)
PPD_diff <- PPD_A2 - PPD_A1

PPD_A2






contrasts_dat = contrasts_data(ba_seconds,ba_seconds_fit)
credible_intervals_dat = credible_intervals_data(ba_seconds,ba_seconds_fit)
credible_intervals_dat
plot_contrasts(credible_intervals_dat)


treatments = unique(ba_seconds$group)

NB_vs_txox = contrasts_dat$txox

quants = quantile(NB_vs_txox, probs = c(.5, .025, .975))

quants[1]

hist(NB_vs_txox)

broom::tidyMCMC(NB_vs_txox, conf.int = TRUE, conf.method = "HPDinterval")

ss = marginal_effects(ba_seconds_fit)

pred = predict(ba_seconds_fit,ba_seconds)


pred

nrow(ba_seconds)

h.Trt01 <- hypothesis(ba_seconds_fit, hypothesis=c('Trt1vs0'=c("ba_seconds_fit = 0")))
