# Function to create calls for the regression transition model
#
# Inputs:
# response - response variables
# model_type - the type of model to be fitted to each response variable, e.g. lm, binomial, gls
# df - Optional argument for extracting class of response variables and number of levels for factor variables 
# covs - additional non-static covariates, e.g. age
# static - static covariates, e.g. ethnicity
mvcall <- function(response,model_type = NULL,df = NULL, covs = NULL,static = NULL,cov_order = "current",prev_measure = "p",mode = "long",var_select_val=NULL,sigma_vars = NULL){
  # Create an empty list
  mv_list <- vector(mode = "list", length = length(response))
  response_init <- response
  if(mode == "long"){
    # If covs is non-empty, append two to its elements
    if(!is.null(covs)){
      covs <- paste(covs,2,sep = "")
    }
    if(!(var_select_val %in% c("determ","determ_step"))){
      for(i in 1:length(response)){
        if(model_type[[i]] == "binomial" & prev_measure == "i"){
          if(cov_order == "both"){ 
            # Create model formula
            mv_list[[i]]$call <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(unique(c(response[-i],response_init[-i])),covs,static), collapse = " + ")),collapse = " ~ ")
            mv_list[[i]]$varlist <- c(paste(c(response[i],2), collapse = ""),unique(c(response[-i],response_init[-i])),covs,static)
          } else {
            # Create model formula
            mv_list[[i]]$call <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(response[-i],covs,static), collapse = " + ")),collapse = " ~ ")
            mv_list[[i]]$varlist <- c(paste(c(response[i],2), collapse = ""),response[-i],covs,static)
          }
        } else {
          if(cov_order == "both"){ 
            # Create model formula
            mv_list[[i]]$call <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(unique(c(response,response_init)),covs,static), collapse = " + ")),collapse = " ~ ")
            mv_list[[i]]$varlist <- c(paste(c(response[i],2), collapse = ""),unique(c(response,response_init)),covs,static)
          } else {
            # Create model formula
            mv_list[[i]]$call <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(response,covs,static), collapse = " + ")),collapse = " ~ ")
            mv_list[[i]]$varlist <- c(paste(c(response[i],2), collapse = ""),response,covs,static)
          }
        }
        if(model_type[[i]] == "binomial"){
          mv_list[[i]]$prev_measure <- prev_measure
        }
        # If using lmvar model remove intercept from model formula
        if(model_type[[i]] %in% c("gls","lmvar")){
          mv_list[[i]]$call <- paste(c(mv_list[[i]]$call,"- 1"), collapse = " ")
        }
        # Add model type
        if(length(model_type) >= 1){
          mv_list[[i]]$modeltype <- model_type[i]
        }
        # If a data frame is supplied, add information about class of response variables and number of levels for factor variables
        if(length(df) >= 1){
          mv_list[[i]]$class <- class(df[,response[i]])
          if(class(df[,response[i]]) == "factor"){
            mv_list[[i]]$levels <- nlevels(df[,response[i]])
          }
        }
        if(cov_order %in% c("current","both")){
          # Append two to the name of the ith response variable to create multivariate structure
          response[i] <- paste(c(response[i],2),collapse = "")
        }
      }
    } else {
      for(i in 1:length(response)){
        if(model_type[[i]] == "binomial" & prev_measure == "i"){
          if(cov_order == "both"){ 
            # Create model formula
            mv_list[[i]]$call <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(unique(c(response[-i],response_init[-i])),covs,static), collapse = " + ")),collapse = " ~ ")
            mv_list[[i]]$varlist <- c(paste(c(response[i],2), collapse = ""),unique(c(response[-i],response_init[-i])),covs,static)
            # Create model variance formula
            if(!is.null(sigma_vars)){
              mv_list[[i]]$varlist_sigma <- c(paste(c(response[i],2), collapse = ""),sigma_vars)
              mv_list[[i]]$call_sigma <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(sigma_vars), collapse = " + ")),collapse = " ~ ")
            } else{
              mv_list[[i]]$varlist_sigma <- c(paste(c(response[i],2), collapse = ""),covs,static)
              mv_list[[i]]$call_sigma <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(covs,static), collapse = " + ")),collapse = " ~ ")
            }
          } else {
            # Create model formula
            mv_list[[i]]$call <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(response[-i],covs,static), collapse = " + ")),collapse = " ~ ")
            mv_list[[i]]$varlist <- c(paste(c(response[i],2), collapse = ""),response[-i],covs,static)
            # Create variance formula
            if(!is.null(sigma_vars)){
              mv_list[[i]]$varlist_sigma <- c(paste(c(response[i],2), collapse = ""),sigma_vars)
              mv_list[[i]]$call_sigma <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(sigma_vars), collapse = " + ")),collapse = " ~ ")
            } else{
              mv_list[[i]]$varlist_sigma <- c(paste(c(response[i],2), collapse = ""),covs,static)
              mv_list[[i]]$call_sigma <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(covs,static), collapse = " + ")),collapse = " ~ ")
            }
          }
        } else {
          if(cov_order == "both"){ 
            # Create model formula
            mv_list[[i]]$call <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(unique(c(response,response_init)),covs,static), collapse = " + ")),collapse = " ~ ")
            mv_list[[i]]$varlist <- c(paste(c(response[i],2), collapse = ""),unique(c(response,response_init)),covs,static)
            if(!is.null(sigma_vars)){
              mv_list[[i]]$varlist_sigma <- c(paste(c(response[i],2), collapse = ""),sigma_vars)
              mv_list[[i]]$call_sigma <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(sigma_vars), collapse = " + ")),collapse = " ~ ")
            } else{
              mv_list[[i]]$varlist_sigma <- c(paste(c(response[i],2), collapse = ""),covs,static)
              mv_list[[i]]$call_sigma <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(covs,static), collapse = " + ")),collapse = " ~ ")
            }
          } else {
            # Create model formula
            mv_list[[i]]$call <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(response,covs,static), collapse = " + ")),collapse = " ~ ")
            mv_list[[i]]$varlist <- c(paste(c(response[i],2), collapse = ""),response,covs,static)
            if(!is.null(sigma_vars)){
              mv_list[[i]]$varlist_sigma <- c(paste(c(response[i],2), collapse = ""),sigma_vars)
              mv_list[[i]]$call_sigma <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(sigma_vars), collapse = " + ")),collapse = " ~ ")
            } else{
              mv_list[[i]]$varlist_sigma <- c(paste(c(response[i],2), collapse = ""),covs,static)
              mv_list[[i]]$call_sigma <- paste(c(paste(c(response[i],2), collapse = ""),paste(c(covs,static), collapse = " + ")),collapse = " ~ ")
            }
          }
        }
        if(model_type[[i]] == "binomial"){
          mv_list[[i]]$prev_measure <- prev_measure
        }
        # If using lmvar model remove intercept from model formula
        if(model_type[[i]] %in% c("gls","lmvar")){
          mv_list[[i]]$call <- paste(c(mv_list[[i]]$call,"- 1"), collapse = " ")
          mv_list[[i]]$call_sigma <- paste(c(mv_list[[i]]$call_sigma,"- 1"), collapse = " ")
        }
        # Add model type
        if(length(model_type) >= 1){
          mv_list[[i]]$modeltype <- model_type[i]
        }
        # If a data frame is supplied, add information about class of response variables and number of levels for factor variables
        if(length(df) >= 1){
          mv_list[[i]]$class <- class(df[,response[i]])
          if(class(df[,response[i]]) == "factor"){
            mv_list[[i]]$levels <- nlevels(df[,response[i]])
          }
        }
        if(cov_order %in% c("current","both")){
          # Append two to the name of the ith response variable to create multivariate structure
          response[i] <- paste(c(response[i],2),collapse = "")
        }
      }
    }
  } else {
    for(i in 1:length(response)){
      if(model_type[[i]] == "binomial" & prev_measure == "i"){
        # Create model formula
        mv_list[[i]]$call <- paste(c(response[i],paste(c(response[0:(i-1)],covs,static), collapse = " + ")),collapse = " ~ ")
        mv_list[[i]]$varlist <- c(response[i],response[0:(i-1)],covs,static)
      } else {
        # Create model formula
        mv_list[[i]]$call <- paste(c(response[i],paste(c(response[0:(i-1)],covs,static), collapse = " + ")),collapse = " ~ ")
        mv_list[[i]]$varlist <- c(response[i],response[0:(i-1)],covs,static)
      }
      if(model_type[[i]] == "binomial"){
        mv_list[[i]]$prev_measure <- prev_measure
      }
      # If using lmvar model remove intercept from model formula
      if(model_type[[i]] %in% c("gls","lmvar")){
        mv_list[[i]]$call <- paste(c(mv_list[[i]]$call,"- 1"), collapse = " ")
        mv_list[[i]]$call_sigma <- paste(c(mv_list[[i]]$call_sigma,"- 1"), collapse = " ")
      }
      # Add model type
      if(length(model_type) >= 1){
        mv_list[[i]]$modeltype <- model_type[i]
      }
      # If a data frame is supplied, add information about class of response variables and number of levels for factor variables
      if(length(df) >= 1){
        mv_list[[i]]$class <- class(df[,response[i]])
        if(class(df[,response[i]]) == "factor"){
          mv_list[[i]]$levels <- nlevels(df[,response[i]])
        }
      }
    }
  }
  names(mv_list) <- response
  return(mv_list)
}
