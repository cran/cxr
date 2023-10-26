## ----setup,echo=FALSE,message=FALSE,warning=FALSE-----------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

# library(tidyverse,quietly = TRUE)
# library(Matrix) # needed for fill_demography_matrix

## ----echo=FALSE, out.width = "99%", fig.cap="Conceptual overview of metapopulation functions. We model stage-specific demographic transitions from t to t+1 (here, t = one year) using a matrix population model for populations of n species in h sites. Stages are juveniles, J, non-reproductive, N, and reproductive, R. Stage transitions depend on demographic rates of survival (S), transition probability to the reproductive stage (T), and recruitment (O). Non-reproductive individuals can also disperse (D~N~) and have a probabaility of S~DN~ to survive dispersal. The example logistic regression for species 1 shows that rates are modeled as functions of an environmental variable (E) and intra- (here N~S1~) and interspecific (here N~S2~, N~S3~, N~S4~) density.The different colors highlight that the mean demographic rates (the α in the models) differ among species and habitats."----
# All defaults
knitr::include_graphics("Figure1.png")

## -----------------------------------------------------------------------------
# load the package
library(cxr)

# define species names
sp <- c("s1","s2","s3")
# number of species
num.sp <- length(sp)

# define site names
sites <- c("sa","sb")
# number of sites
num.sites <- length(sites)

# number of demographic stages - this should be always fixed
num.stages <- 3

# vital rates to account for - these names should be fixed
rates <- c("Sj","Sn","Sr","Rn","Rr","D","Ds","O")

# years (or time steps) of simulations
years <- 100 

# simulate environment from normal distribution
set.seed(123)
env <- rnorm(years, mean=0, sd=1)


## -----------------------------------------------------------------------------
data("metapopulation_example_param", package = "cxr")

# coefficients for species "s1" in site "sa"
metapopulation_example_param[["s1"]][["sa"]]

## -----------------------------------------------------------------------------
# build the templates for the vec-permutation matrices
# returns a nested list: [[matrix.type]][[species number]]
# where matrix.type is "demography", "dispersal", or "permutation"
vpm <- vec_permutation_matrices(num.sp,num.sites,num.stages)
# for example
vpm[["demography"]][[1]]

## -----------------------------------------------------------------------------
# initial.densities: list: [[species]][sites*stages]
# this needs to be filled manually
initial.densities <- list()
# sp1
initial.densities[[1]] <- matrix(c(10,10,15,10,10,13),
                                 nrow = num.sites,
                                 ncol = num.stages,
                                 byrow = TRUE)
# sp2
initial.densities[[2]] <- matrix(c(15,5,3,15,5,3),
                                 nrow = num.sites,
                                 ncol = num.stages,
                                 byrow = TRUE)
# sp2
initial.densities[[3]] <- matrix(c(5,5,2,3,4,2),
                                 nrow = num.sites,
                                 ncol = num.stages,
                                 byrow = TRUE)

## -----------------------------------------------------------------------------
initial.densities[[1]]

## -----------------------------------------------------------------------------
# projected.densities list: [[species]][[years]][[sites*stages]
projected.densities <- list()
for(i.sp in 1:num.sp){
  projected.densities[[i.sp]] <- list()
  for(i.year in 1:years){
    projected.densities[[i.sp]][[i.year]] <- matrix(0,
                                                    nrow = num.sites,
                                                    ncol = num.stages)
  }
}

## -----------------------------------------------------------------------------
# create an auxiliary list, to keep track of the densities per timestep
current.densities <- initial.densities

for(i.year in 1:years){

  # Update projected densities at this timestep
  for(i.sp in 1:num.sp){
    projected.densities[[i.sp]][[i.year]] <- current.densities[[i.sp]]
  }

  # update transition matrices ----------------------------------------------

  # this is a list per species and site
  transition_matrices <- list()

  for(i.sp in 1:length(sp)){
    transition_matrices[[i.sp]] <- list()
    for(i.site in 1:length(sites)){

      # store the transition matrix for this sp and site
      transition_matrices[[i.sp]][[i.site]] <- fill_transition_matrix(focal.sp = i.sp,
                                                                      site = i.site,
                                                                      param = metapopulation_example_param,
                                                                      env = env[i.year],
                                                                      current.densities = current.densities)

    }# for each site
    names(transition_matrices[[i.sp]]) <- sites
  }# for each species
  names(transition_matrices) <- sp

  # update demography and dispersal matrices --------------------------------

  for(i.sp in 1:length(sp)){

    # demography
    vpm[["demography"]][[i.sp]] <- fill_demography_matrix(focal.sp = i.sp,
                                                          vpm = vpm,
                                                          transition_matrices = transition_matrices)
    # dispersal
    vpm[["dispersal"]][[i.sp]] <- fill_dispersal_matrix(focal.sp = i.sp,
                                                        num.sites = num.sites,
                                                        param = metapopulation_example_param,
                                                        vpm = vpm,
                                                        env = env[i.year],
                                                        current.densities = current.densities)
  }# for i.sp

  # update densities --------------------------------------------------------
  # all stages and sites for each species
  for(i.sp in 1:length(sp)){
    current.densities[[i.sp]] <- calculate_densities(focal.sp = i.sp,
                                                     vpm,
                                                     current.densities)
  }# for i.sp
}# for i.year

## -----------------------------------------------------------------------------
# transform list of projected densities to dataframe
df <- densities_to_df(projected.densities)

# tidy
df$species <- dplyr::recode(as.character(df$species), "1" = "S1", "2" = "S2", "3" = "S3")

# plot
dynamics.plot <- ggplot2::ggplot(df,ggplot2::aes(year,density,col=species))+
  ggplot2::geom_line()+
  ggplot2::facet_grid(stage~site,scales = "free")+
  ggplot2::scale_color_manual(name="",values=c("darkgreen","orange","darkred"))+
  ggplot2::xlab("Simulation year")+ggplot2::ylab("Total density")+
  ggplot2::theme_bw(base_size=20)+
  ggplot2::theme(panel.grid = ggplot2::element_blank())+
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "black"))

## ----echo=FALSE, fig.height=9, fig.width=7------------------------------------
dynamics.plot

## -----------------------------------------------------------------------------
# this builds an empty data structure
example_param <- build_param(sp = sp,
                     sites = sites,
                     rates = rates,
                     env = env,
                     num.params = 8)

## -----------------------------------------------------------------------------
example_param[["s1"]][["sa"]]

## -----------------------------------------------------------------------------
data("glm_example_coefs",package = "cxr")

glm_example_coefs

## -----------------------------------------------------------------------------
glm.coef.equivalence <- list("alpha" = "(Intercept)",
                             "beta1" = "env",
                             "beta2" = "dens", 
                             "beta3" = "densS2", 
                             "beta4" = "densS3",
                             "beta5" = "dens:env",
                             "beta6" = "env:densS2",
                             "beta7" = "env:densS3")

## -----------------------------------------------------------------------------
example_param <- generate_vital_rate_coefs(example_param,
                                           sp = "s1",
                                           sites = "sa",
                                           vital.rate = "Sj",
                                           glm.object = glm_example_coefs,
                                           glm.coef.equivalence = glm.coef.equivalence)

# the resulting table. Note that we only entered the coefficients for survival of juveniles (Sj)
example_param[["s1"]][["sa"]]

## -----------------------------------------------------------------------------
example_param$s2$sa["Sr",5] <- 0.001

## ----echo=FALSE, out.width="20%", out.height="20%", fig.cap="The prey (S1)"----
knitr::include_graphics("rabbit.png")

## ----echo=FALSE, out.width="20%", out.height="20%", fig.cap="The meso-predator (S2)"----
knitr::include_graphics("fox.png")

## ----echo=FALSE, out.width="20%", out.height="20%", fig.cap="The top-predator (S2)"----
knitr::include_graphics("lynx.png")

## -----------------------------------------------------------------------------
example_param <- generate_vital_rate_coefs(param = example_param,
                                   sp = "s1",
                                   # sites = "sa",
                                   vital.rate = c("Sj"),#,"Sn","Sr","Rn","Rr","D","Ds","O"),
                                   vr.coef = "alpha",
                                   mean.coef = .1,sd.coef = 0)

