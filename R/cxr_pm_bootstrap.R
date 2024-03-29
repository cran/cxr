#' Standard error estimates for model parameters
#' 
#' Computes bootstrap standard errors for a given population dynamics model.
#' This function is provided for completeness, but error calculation is
#' integrated in the function \code{cxr_pm_fit}.
#'
#' @param fitness_model function returning a single value to minimize, given a set of parameters and a fitness metric
#' @param optimization_method numerical optimization method
#' @param data dataframe with observations in rows and two sets of columns:
#' * fitness: fitness metric for the focal individual
#' * neighbours: columns with user-defined names with number of neighbours for each group
#' @param focal_column optional integer value giving the position, or name, of the column
#' with neighbours from the same species as the focal one. This is necessary if "alpha_intra" is specified. 
#' @param covariates optional matrix with observations in rows and covariates in columns. Each cell is the value of a covariate
#' in a given observation.
#' @param init_par 1d vector of initial parameters
#' @param lower_bounds 1d vector of lower bounds
#' @param upper_bounds 1d vector of upper bounds
#' @param fixed_parameters optional list specifying values of fixed parameters, 
#' with components "lambda","alpha_intra","alpha_inter","lambda_cov", and "alpha_cov".
#' @param bootstrap_samples how many bootstrap samples to compute.
#'
#' @return 1d vector, the standard error of each parameter in init_par
#' @import stats 
#' @md
#' @export
cxr_pm_bootstrap <- function(fitness_model,
                             optimization_method,
                             data,
                             focal_column,
                             covariates,
                             init_par,
                             lower_bounds,
                             upper_bounds,
                             fixed_parameters,
                             bootstrap_samples){ 
  if(bootstrap_samples<2){
    print("cxr_pm_bootstrap: number of bootstrap samples cannot be < 2. Setting bootstrap samples to 2.")
    bootstrap_samples <- 2
  }
  
  bootres <- matrix(nrow = bootstrap_samples, ncol = length(init_par))
  
  for(i.sample in 1:bootstrap_samples){
    
    bsample <- sample(nrow(data),nrow(data),replace = TRUE)
    # bsample <- sample(length(log.fitness),length(log.fitness),replace = T)
    
    # sample data
    bdata <- data[bsample,]
    
    # same as in cxr_pm_fit
    # just to avoid a note in R CMD CHECK
    dropname <- "fitness"
    bneigh_matrix <- as.matrix(bdata[ , !(names(bdata) %in% dropname)])

    
    # initial check to remove any neighbour with no presences
    # also ensure there is at least one neighbour with presences
    if(any(colSums(bneigh_matrix) == 0)){
      bempty_neigh <- colnames(bneigh_matrix)[which(colSums(bneigh_matrix) == 0)]
      bpresent_neigh <- colnames(bneigh_matrix)[which(colSums(bneigh_matrix) != 0)]
      
      if(length(bpresent_neigh) == 0){
        message("cxr_pm_bootstrap ERROR: No neighbours with densities > 0 in any observation.")      
        return(NULL)
      }else{
        bneigh_matrix <- bneigh_matrix[,bpresent_neigh]
      }
      
      # neigh_matrix <- neigh_matrix[,present_neigh]
    }else{
      bempty_neigh <- NULL
      bpresent_neigh <- colnames(bneigh_matrix)
    }# if-else any neighbour with no presences
    
    
    if(is.null(focal_column)){
      # no alpha_intra
      bneigh_inter_matrix <- bneigh_matrix
      bneigh_intra_matrix <- NULL
      # set also names
      bneigh_inter <- colnames(bneigh_inter_matrix)
      bneigh_intra <- NULL
    }else{
      
      # which column number
      if(inherits(focal_column,"character")){
        bfocal_column_num <- which(colnames(bneigh_matrix) == focal_column)
      }else{
        bfocal_column_num <- focal_column -1
      }
      # intra and inter observations in different matrices
      bneigh_inter_matrix <- bneigh_matrix[,-bfocal_column_num]
      bneigh_intra_matrix <- as.matrix(bneigh_matrix[,bfocal_column_num])
      colnames(bneigh_intra_matrix) <- colnames(bneigh_matrix)[bfocal_column_num]
      # set also names
      bneigh_inter <- colnames(bneigh_inter_matrix)
      bneigh_intra <- colnames(bneigh_matrix)[bfocal_column_num]
    }
    
    if(is.data.frame(covariates)){
      bcov <- as.data.frame(covariates[bsample,])
    }else if(is.matrix(covariates)){
      bcov <- as.data.frame(covariates[bsample,])
    }else{
      bcov <- 0
    }
    
    bpar <- NULL
    
    ############
    if(optimization_method %in% c("BFGS", "CG", "Nelder-Mead", "ucminf")){
      tryCatch({
        bpar <- optimx::optimx(par = init_par, 
                               fn = fitness_model, 
                               gr = NULL, 
                               method = optimization_method,
                               # lower = lower_bounds,
                               # upper = upper_bounds,
                               control = list(), 
                               hessian = F,
                               fitness = log(bdata$fitness), 
                               neigh_intra_matrix = bneigh_intra_matrix,
                               neigh_inter_matrix = bneigh_inter_matrix,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        
        par.pos <- which(!names(bpar) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
        bpar <- as.numeric(bpar[,par.pos])
        row.names(bpar) <- NULL
        
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method %in% c("L-BFGS-B", "nlm", "nlminb", 
                                        "Rcgmin", "Rvmmin", "spg", 
                                        "bobyqa", "nmkb", "hjkb")){
      tryCatch({
        bpar <- optimx::optimx(par = init_par, 
                               fn = fitness_model, 
                               gr = NULL, 
                               method = optimization_method,
                               lower = lower_bounds,
                               upper = upper_bounds,
                               control = list(), 
                               hessian = F,
                               fitness = log(bdata$fitness), 
                               neigh_intra_matrix = bneigh_intra_matrix,
                               neigh_inter_matrix = bneigh_inter_matrix,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        
        par.pos <- which(!names(bpar) %in% c("value","fevals","gevals","niter","convcode","kkt1","kkt2","xtime"))
        bpar <- as.numeric(bpar[,par.pos])
        row.names(bpar) <- NULL
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_CRS2_LM"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                               eval_f = fitness_model,
                               opts = list("algorithm"="NLOPT_GN_CRS2_LM", "maxeval"=1e4),
                               lb = lower_bounds,
                               ub = upper_bounds,
                               fitness = log(bdata$fitness), 
                               neigh_intra_matrix = bneigh_intra_matrix,
                               neigh_inter_matrix = bneigh_inter_matrix,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_ISRES"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                               eval_f = fitness_model,
                               opts = list("algorithm"="NLOPT_GN_ISRES", "maxeval"=1e4),
                               lb = lower_bounds,
                               ub = upper_bounds,
                               fitness = log(bdata$fitness), 
                               neigh_intra_matrix = bneigh_intra_matrix,
                               neigh_inter_matrix = bneigh_inter_matrix,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "nloptr_DIRECT_L_RAND"){
      tryCatch({
        bpar <- nloptr::nloptr(x0 = init_par,
                               eval_f = fitness_model,
                               opts = list("algorithm"="NLOPT_GN_DIRECT_L_RAND", "maxeval"=1e4),
                               lb = lower_bounds,
                               ub = upper_bounds,
                               fitness = log(bdata$fitness), 
                               neigh_intra_matrix = bneigh_intra_matrix,
                               neigh_inter_matrix = bneigh_inter_matrix,
                               covariates = covariates, 
                               fixed_parameters = fixed_parameters)
        bpar <- bpar$solution
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "GenSA"){
      tryCatch({
        bpar <- GenSA::GenSA(par = init_par,
                             fn = fitness_model,
                             lower = lower_bounds,
                             upper = upper_bounds, 
                             control = list(maxit = 1e3), 
                             fitness = log(bdata$fitness), 
                             neigh_intra_matrix = bneigh_intra_matrix,
                             neigh_inter_matrix = bneigh_inter_matrix,
                             covariates = covariates, 
                             fixed_parameters = fixed_parameters)
        bpar <- bpar$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    # }else if(optimization_method == "hydroPSO"){
    #   tryCatch({
    #     # suppress annoying output
    #     bpar <- hydroPSO::hydroPSO(par = init_par,
    #                                fn = fitness_model,
    #                                lower = lower_bounds,
    #                                upper = upper_bounds, 
    #                                control=list(write2disk=FALSE, maxit = 1e3, MinMax = "min", verbose = F),
    #                                fitness = log(bdata$fitness), 
    #                                neigh_intra_matrix = bneigh_intra_matrix,
    #                                neigh_inter_matrix = bneigh_inter_matrix,
    #                                covariates = covariates, 
    #                                fixed_parameters = fixed_parameters)
    #     
    #     bpar <- bpar$par
    #   }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }else if(optimization_method == "DEoptimR"){
      tryCatch({
        bpar <- DEoptimR::JDEoptim(lower = lower_bounds,
                                   upper = upper_bounds,
                                   fn = fitness_model,
                                   fitness = log(bdata$fitness), 
                                   neigh_intra_matrix = bneigh_intra_matrix,
                                   neigh_inter_matrix = bneigh_inter_matrix,
                                   covariates = covariates, 
                                   fixed_parameters = fixed_parameters)
        bpar <- bpar$par
      }, error=function(e){cat("cxr_pm_bootstrap ERROR :",conditionMessage(e), "\n")})
    }
    
    if(!is.null(bpar)){
      if(sum(is.na(bpar)) == 0){
        bootres[i.sample,] <- bpar
      }
    }else{
      return(NULL)
    }
    
  }# for i.sample  
  
  bootres <- bootres[which(!is.na(rowSums(bootres))),]
  valid.samples <- nrow(bootres)
  if(!is.null(valid.samples)){
    if(valid.samples == 1){
      boot.se <- bootres
    }else{
      boot.se <- apply(bootres,2,sd)
    }
  }else{
    boot.se <- rep(NA_integer_,length(init_par))
  }
  # 
  # if(nrow(bootres)>2){
  #   boot.se <- apply(bootres,2,sd)
  # }else if(nrow(bootres)!= 0){
  #   boot.se <- bootres[1,]
  # }else{
  #   boot.se <- rep(NA_integer_,length(init_par))
  # }
  
  if(!is.null(names(init_par))){
    names(boot.se) <- paste(names(init_par),"_se",sep="")
  }
  boot.se

}


