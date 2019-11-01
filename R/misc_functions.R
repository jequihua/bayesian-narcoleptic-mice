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
# Function to create interval data for plotting.
interval_data_eeg = function(brms_model_fit,columns=1:5,
                         param_names=c("Intercept","NB","sham","txcb","txox"),
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
    geom_segment(aes_(x = ~ ll, xend = ~ hh), size = 3) + 
    #geom_segment(aes_(x = ~ l, xend = ~ h), size = 4) + 
    geom_point(aes_(x = ~ m), data = maybe_points, size = 8) + 
    scale_y_discrete(limits = rev(data$parameter)) + 
    labs(x = "Coefficient values", y = NULL)+
    theme_bw(base_size = 40)
}

###########################################################################################  
# Function to create data for contrasts between treatments (using "posterior_predict" or "fitted").
contrasts_data = function(model_input_data,brms_model_fit){
      posterior_data = as.data.frame(fitted(
      brms_model_fit,
      newdata = data.frame(group = levels(model_input_data$group)),
      re_formula = NA,
      summary = FALSE # extract the full MCMC)
    ))
  
  colnames(posterior_data) = unique(model_input_data$group)
  return(posterior_data)
}

###########################################################################################  
# Function to create data for contrasts between treatments (using "posterior_predict" or "fitted").
contrasts_data_hours = function(model_input_data,brms_model_fit){
  posterior_data = as.data.frame(fitted(
    brms_model_fit,
    newdata = data.frame(group = levels(model_input_data$group),hours=rep(4,4)),
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
# Create data for interval plots.
credible_intervals_data = function(data,fit,untreated="NB",columns=1:3){
  contrasts_dat = contrasts_data(data[,columns],fit)
  
  treatments = as.character(unique(data$group))
  treatments = treatments[treatments!=untreated]
  
  output = data.frame(matrix(0,length(treatments),8))
  names(output)=c("term", "estimate", "conf.low", "conf.high","mean",
                  "ctrl_estimate","treatment_estimate","reduction_percentage")
  
  for (i in 1:nrow(output)){
    
    treatment = treatments[i]
    output[i,"term"]=paste0(treatment," - ",untreated)
    
    NB_vs_treatment = contrasts_dat[,treatment]-contrasts_dat[,untreated]
    
    NB = contrasts_dat[,untreated]
    treat =  contrasts_dat[,treatment]
    
    quants = quantile(NB_vs_treatment, probs = c(.5, .025, .975))
    
    output[i,"estimate"]=quants[1]
    output[i,"conf.low"]=quants[2]
    output[i,"conf.high"]=quants[3]
    output[i,"mean"]=mean(NB_vs_treatment)
    output[i,"ctrl_estimate"]=mean(NB)
    output[i,"treatment_estimate"]=mean(treat)
    output[i,"reduction_percentage"]=output[i,"mean"]/output[i,"ctrl_estimate"]
    
  }
  
  output$term = factor(output$term,levels=rev(output$term))
  
  return(output)
}

########################################################################################### 
# Create data for interval plots (hours variable, needs 1+).
credible_intervals_data_hours = function(data,fit,untreated="NB",columns=1:3){
  contrasts_dat = contrasts_data_hours(data[,columns],fit)
  
  treatments = as.character(unique(data$group))
  treatments = treatments[treatments!=untreated]
  
  output = data.frame(matrix(0,length(treatments),8))
  names(output)=c("term", "estimate", "conf.low", "conf.high","mean",
                  "ctrl_estimate","treatment_estimate","reduction_percentage")
  
  for (i in 1:nrow(output)){
    
    treatment = treatments[i]
    output[i,"term"]=paste0(treatment," - ",untreated)
    
    NB_vs_treatment = contrasts_dat[,treatment]-contrasts_dat[,untreated]
    
    NB = contrasts_dat[,untreated]
    treat =  contrasts_dat[,treatment]
    
    quants = quantile(NB_vs_treatment, probs = c(.5, .025, .975))
    
    output[i,"estimate"]=quants[1]
    output[i,"conf.low"]=quants[2]
    output[i,"conf.high"]=quants[3]
    output[i,"mean"]=mean(NB_vs_treatment)
    output[i,"ctrl_estimate"]=mean(NB)
    output[i,"treatment_estimate"]=mean(treat)
    output[i,"reduction_percentage"]=output[i,"mean"]/output[i,"ctrl_estimate"]
    
  }
  
  output$term = factor(output$term,levels=rev(output$term))
  
  return(output)
}

########################################################################################### 
# Plot pairwise comparisons (contrasts).
plot_contrasts = function(comps,xlab="Effect size"){
  ggplot(comps, aes(y = estimate, x = term)) + 
    geom_pointrange(aes(ymin = conf.low,ymax = conf.high),size=2) +
    geom_hline(yintercept = 0, size = 3, color = "gray") +
    scale_y_continuous(xlab) + scale_x_discrete("") + coord_flip() +
    theme_bw(base_size = 40)
}








