# Load packages.
library("ggplot2")
library("multcomp")
library("broom")

###########################################################################################
# Function to create interval data for plotting.
interval_data = function(brms_model_fit,columns=1:4,
                                        param_names=c("Intercept","sham","txcb","txox"),
                                        prob=0.5,
                                        prob_outer=0.95,
                                        point_est="median"){
  post = posterior_samples(brms_model_fit)[,columns]
  colnames(post)=param_names
  interval_df = mcmc_intervals_data(x=post,prob=prob,
                                             prob_outer=prob_outer,
                                             point_est=point_est)
  return(interval_df)
}

###########################################################################################
# Function to produce a coefficient plot (credible intervals).
plot_intervals <- function(data, draw_points = TRUE, draw_ref_line = TRUE,
                           line_position = 0) {
  
  maybe_points <- if (draw_points) data else data[0, ]
  
  geom_maybe_vline <- if (draw_ref_line) geom_vline else geom_ignore
  
  ggplot(data) + 
    aes_(y = ~ parameter, yend = ~ parameter) +
    geom_maybe_vline(xintercept = line_position, size = 2, color = "gray") +
    geom_segment(aes_(x = ~ ll, xend = ~ hh), size = 2) + 
    geom_segment(aes_(x = ~ l, xend = ~ h), size = 4) + 
    geom_point(aes_(x = ~ m), data = maybe_points, size = 7) + 
    scale_y_discrete(limits = rev(data$parameter)) + 
    labs(x = "Coefficient values", y = NULL)+
    theme_bw(base_size = 20)
}

###########################################################################################  
# Function to create data for contrasts between treatments (using "posterior_predict").
contrasts_data = function(model_input_data,brms_model_fit){
      posterior_data = as.data.frame(posterior_predict(
      brms_model_fit,
      newdata = data.frame(group = levels(model_input_data$group)),
      re_formula = NA,
      summary = FALSE # extract the full MCMC)
    ))
  
  colnames(posterior_data) = unique(model_input_data$group)
  return(posterior_data)
}

########################################################################################### 
# Create a Tukeys percentage contrast matrix.
tukeys_Pmat = function(model_input_data,brms_model_fit,columns=1:4,estimate.method="median"){
  mcmc = brms_model_fit
  coefs <- as.matrix(mcmc)[, columns]
  newdata <- data.frame(x = levels(model_input_data$group))
  tuk.mat <- contrMat(n = table(newdata$x), type = "Tukey")
  Xmat <- model.matrix(~x, data = newdata)
  pairwise.mat <- tuk.mat %*% Xmat
  tuk.mat[tuk.mat == -1] = 0
  comp.mat <- tuk.mat %*% Xmat
  comp.mcmc = 100 * (coefs %*% t(pairwise.mat))/coefs %*% t(comp.mat)
  (comps = broom::tidyMCMC(comp.mcmc, conf.int = TRUE, conf.method = "HPDinterval"))
  comps$term = factor(comps$term, levels = rev(unique(comps$term)))
  return(comps)
}

########################################################################################### 
# Create a Tukeys contrast matrix.
tukeys_mat = function(model_input_data,brms_model_fit,columns=1:4,estimate.method="median"){
  mcmc <- brms_model_fit
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

########################################################################################### 
# Plot pairwise comparisons (contrasts).
plot_contrasts = function(comps,xlab="Effect size"){
  ggplot(comps, aes(y = estimate, x = term)) + 
    geom_pointrange(aes(ymin = conf.low,ymax = conf.high),size=2) +
    geom_hline(yintercept = 0, size = 2, color = "gray") +
    scale_y_continuous(xlab) + scale_x_discrete("") + coord_flip() +
    theme_bw(base_size = 25)
}






