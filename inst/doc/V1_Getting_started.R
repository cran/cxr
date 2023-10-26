## ----setup,echo=FALSE---------------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## -----------------------------------------------------------------------------
library(cxr)
data("neigh_list") 

## -----------------------------------------------------------------------------
my.sp <- "HOMA" 
# get data from the list
obs_homa <- neigh_list[[my.sp]]
# no need for ID column
obs_homa <- subset(obs_homa,select = -c(obs_ID))
# For each observation, we need the individual plant fitness 
# and the number of neighbours per species (in columns).
head(obs_homa)

#There are a couple of outliers, eliminate from analyses
obs_homa<-subset(obs_homa, fitness<500)

## -----------------------------------------------------------------------------
#?cxr_pm_fit #check the help file for a description of the arguments
fit_homa <- cxr_pm_fit(data = obs_homa,
                       focal_column = my.sp,
                       model_family = "BH",
                       covariates = NULL,
                       optimization_method = "Nelder-Mead",
                       alpha_form = "pairwise",
                       lambda_cov_form = "none",
                       alpha_cov_form = "none",
                       initial_values = list(lambda = 45,
                                             alpha_intra = .1,
                                             alpha_inter = .1),
                       #not aplicable to this optimazation method
                       # lower_bounds = list(lambda = 0, 
                       #                     alpha_intra = 0,
                       #                     alpha_inter = 0),
                       # upper_bounds = list(lambda = 10,
                       #                     alpha_intra = 1,
                       #                     alpha_inter = 1),
                       fixed_terms = NULL,
                       # a low number of bootstrap samples
                       # for demonstration purposes, 
                       # increase it for robust results.
                       bootstrap_samples = 3) 

## -----------------------------------------------------------------------------
summary(fit_homa)

## -----------------------------------------------------------------------------
names(fit_homa) #list of all available elements.

#reproduction success in the absence of neighbors
fit_homa$lambda
# intraspecific interaction
fit_homa$alpha_intra
# interspecific interactions
fit_homa$alpha_inter

## -----------------------------------------------------------------------------
require(ggplot2)
ggplot(obs_homa, aes(HOMA , fitness)) + 
  geom_point() +
  stat_function(fun = function(x) fit_homa$lambda/(1+fit_homa$alpha_intra*x), lwd = 1.5, colour = "blue")

## -----------------------------------------------------------------------------
my.sp <- c("BEMA","CETE","LEMA")
obs_3sp <- neigh_list[my.sp]
# discard ID column
for(i in 1:length(obs_3sp)){
  obs_3sp[[i]] <- obs_3sp[[i]][,2:length(obs_3sp[[i]])]
}
# load covariates: salinity
data("salinity_list")
salinity <- salinity_list[my.sp]
# keep only salinity column
for(i in 1:length(salinity)){
  salinity[[i]] <- as.matrix(salinity[[i]][,2:length(salinity[[i]])])
  colnames(salinity[[i]]) <- "salinity"
}

## -----------------------------------------------------------------------------
names(obs_3sp)

# observation data
head(obs_3sp[[1]])
# number of fitness observations
nrow(obs_3sp[[1]])

# salinity data
head(salinity[[1]])
# number of covariate observations
nrow(salinity[[1]])

## -----------------------------------------------------------------------------

fit_3sp <- list()
all_sp <- names(obs_3sp)

for (i in 1:length(obs_3sp)){
  obs <- obs_3sp[i]
  my.sp <- all_sp[i]
  salinity.sp <- salinity[i]
  lambda <- mean(obs[[1]]$fitness)
  upper_lambda <- as.numeric(max(obs[[1]]$fitness))


fit <- cxr_pm_multifit(data = obs,
                           focal_column = my.sp,
                           model_family = "BH",
                           # here we use a bounded method for demonstration purposes
                           optimization_method = "bobyqa", 
                           covariates = salinity.sp,
                           alpha_form = "pairwise",
                           lambda_cov_form = "global", # effect of covariates over lambda
                           alpha_cov_form = "global", # effect of covariates over alpha
                           initial_values = list(lambda = lambda,
                                                 alpha_intra = 1,
                                                 alpha_inter = 1,
                                                 lambda_cov = 0.1,
                                                 alpha_cov = 0.1),
                           lower_bounds = list(lambda = 0,
                                               alpha_intra = -1,
                                               alpha_inter = -1,
                                               lambda_cov = -2,
                                               alpha_cov = -2),
                           upper_bounds = list(lambda = upper_lambda,
                                               alpha_intra = 3,
                                               alpha_inter = 3,
                                               lambda_cov = 2,
                                               alpha_cov = 2),
                           # no standard errors
                           bootstrap_samples = 0)

fit_3sp[[i]] <- fit
}
names(fit_3sp) <- all_sp

## -----------------------------------------------------------------------------
fit_3sp$BEMA

## -----------------------------------------------------------------------------
fit_3sp$BEMA$log_likelihood

## -----------------------------------------------------------------------------
fit_3sp$BEMA$lambda_cov
fit_3sp$BEMA$alpha_cov

## -----------------------------------------------------------------------------
lower_bounds = list(lambda = 0,
                    alpha_intra = 0,
                    alpha_inter = -1,
                    lambda_cov = 0,
                    alpha_cov = 0)
upper_bounds = list(lambda = 100,
                    alpha_intra = 1,
                    alpha_inter = 1,
                    lambda_cov = 1,
                    alpha_cov = 1)

## -----------------------------------------------------------------------------
# fit three species, as in the previous example
data <- neigh_list[1:3]

# keep only fitness and neighbours columns
for(i in 1:length(data)){
  data[[i]] <- data[[i]][,2:length(data[[i]])]
}
focal_column <- names(data)

# covariates
# in this example, we use salinity
# and two more covariates, randomly generated

data("salinity_list")

# observations for the first three species
cov_list <- salinity_list[1:3]

for(i in 1:length(cov_list)){
  # keep only salinity column
  cov_list[[i]] <- as.matrix(cov_list[[i]][,2:length(cov_list[[i]])])
  
  # add two random covariates
  cov_list[[i]] <- cbind(cov_list[[i]],
                         runif(nrow(cov_list[[i]]),1,10),
                         runif(nrow(cov_list[[i]]),10,100))
  colnames(cov_list[[i]]) <- c("salinity","cov2","cov3")
}

# this is how each element of the covariates list looks like
head(cov_list[[1]])

# function parameters
model_family <- "BH"
covariates <- cov_list
# bobyqa is generally more robust than other bounded methods
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "global"
alpha_cov_form <- "pairwise"

# note how lambda_cov and alpha_cov
# have different initial values for each covariate effect
# the commented assignations are also possible, 
# giving equal initial values to all parameters
initial_values = list(lambda = 1,
                      alpha_intra = 0.1,
                      alpha_inter = 0.1,
                      lambda_cov = c(0.1,0.2,0.1),
                      alpha_cov = c(0.1,0.2,0.1))
# lambda_cov = c(0.1),
# alpha_cov = c(0.1))

# same with boundaries
lower_bounds = list(lambda = 0,
                    alpha_intra = 0,
                    alpha_inter = -1,
                    lambda_cov = c(-1,0,-1),
                    alpha_cov = c(-1,0,-1))
# lambda_cov = c(-1),
# alpha_cov = c(-1))

upper_bounds = list(lambda = 100,
                    alpha_intra = 1,
                    alpha_inter = 1,
                    lambda_cov = c(1,2,1),
                    alpha_cov = c(1,2,1))
# lambda_cov = c(1),
# alpha_cov = c(1))

fixed_terms <- NULL
bootstrap_samples <- 3

## ----eval=FALSE---------------------------------------------------------------
#  # this fit is fairly complex, it may take a while,
#  # it also raises warnings, suggesting that either the data,
#  # model, or initial values/boundaries can be improved.
#  # This is consistent with having observational data
#  # and, furthermore, random covariates
#  
#  fit_multi_cov <- cxr_pm_multifit(data = data,
#                                   focal_column = focal_column,
#                                   model_family = model_family,
#                                   covariates = covariates,
#                                   optimization_method = optimization_method,
#                                   alpha_form = alpha_form,
#                                   lambda_cov_form = lambda_cov_form,
#                                   alpha_cov_form = alpha_cov_form,
#                                   initial_values = initial_values,
#                                   lower_bounds = lower_bounds,
#                                   upper_bounds = upper_bounds,
#                                   fixed_terms = fixed_terms,
#                                   bootstrap_samples = bootstrap_samples)

## -----------------------------------------------------------------------------
fixed_terms <- list(lambda = 1)

# now lambda does not appear in 'initial_values'
initial_values <- list(alpha_intra = 0,
                       alpha_inter = 0,
                       lambda_cov = 0,
                       alpha_cov = 0)
# lower and upper bounds should, likewise, 
# not contain lambda


## -----------------------------------------------------------------------------
fixed_terms <- list(list(lambda = 1), # focal sp 1
                    list(lambda = 1.2), # focal sp 2
                    list(lambda = 1.3)) # focal sp 3

