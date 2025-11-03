# Function which takes the input of the mvcall function to fit multivariate model for use in rtmodel
#
# Inputs: 
# mv_calls - output from mvcall function
# df - data frame to fit models on 
# var_select - variable selection method
library(lmvar,MASS)
mvfit <- function(mv_calls,df,var_select = "full", sigma_covs_rand = FALSE){
  # Create empty list
  mv_fits <- vector(mode = "list", length = length(mv_calls))
  # If full model
  if(var_select == "full"){
    # Fit models
    for(i in 1:length(mv_calls)){
      # If model type is binomial glm
      if(mv_calls[[i]]$modeltype == "binomial"){
        if(mv_calls[[i]]$prev_measure == "i"){
          ind_var <- substr(mv_calls[[i]]$varlist[1], 1, nchar(mv_calls[[i]]$varlist[1])-1)
          mv_fits[[i]] <- glm(mv_calls[[i]]$call,family = binomial,data = df[df[,ind_var] == 0,])
          mv_fits[[i]]$prev_measure <- "i"
        }      
        else {
          mv_fits[[i]] <- glm(mv_calls[[i]]$call,family = binomial,data = df)
          mv_fits[[i]]$prev_measure <- "p"
        }
      } # If model type is ordered logistic
      else if (mv_calls[[i]]$modeltype == "polr"){
        
        mv_fits[[i]] <- polr(mv_calls[[i]]$call,data = df)
        
      } # If model type is a linear regression with non-constant variances
      else if (mv_calls[[i]]$modeltype %in% c("gls","lmvar")){
        model_call <- paste("~",gsub("~","+",mv_calls[[i]]$call))
        df_mat <- model.matrix(formula(model_call), data = df)
        mv_fits[[i]] <- lmvar(y = df_mat[,1],X_mu = df_mat[,-1,drop=FALSE], X_sigma = df_mat[,-1,drop=FALSE],intercept_mu = TRUE,intercept_sigma = TRUE)
        mv_fits[[i]]$model_call <- mv_calls[[i]]$call
        
      } # If model type is a linear model
      else {
        
        mv_fits[[i]] <- lm(mv_calls[[i]]$call,data = df)
        
      }
    }
  }
  # If stepwise model
  else if (var_select == "stepwise"){
    # Fit models
    for(i in 1:length(mv_calls)){
      if(mv_calls[[i]]$modeltype == "binomial"){
        if(mv_calls[[i]]$prev_measure == "i"){
          ind_var <- substr(mv_calls[[i]]$varlist[1], 1, nchar(mv_calls[[i]]$varlist[1])-1)
          stepcall <- step(glm(mv_calls[[i]]$call,family = binomial,data = na.omit(df[df[,ind_var] == 0,mv_calls[[i]]$varlist])),trace = FALSE)$call
          mv_fits[[i]] <- glm(stepcall,family = binomial,data = df[df[ind_var,] == 0,])
          mv_fits[[i]]$prev_measure <- "i"
        } else{
          stepcall <- step(glm(mv_calls[[i]]$call,family = binomial,data = na.omit(df[,mv_calls[[i]]$varlist])),trace = FALSE)$call
          mv_fits[[i]] <- glm(stepcall,family = binomial,data = df)
          mv_fits[[i]]$prev_measure <- "p"
        }
        
        
      } else if (mv_calls[[i]]$modeltype == "polr"){
        stepcall <- step(polr(formula(mv_calls[[i]]$call),data = na.omit(df[,mv_calls[[i]]$varlist]),Hess = TRUE),trace = FALSE)$call 
        mv_fits[[i]] <- polr(stepcall,data = df)
        
      } else if (mv_calls[[i]]$modeltype %in% c("gls","lmvar")){
        model_call <- paste("~",gsub("~","+",mv_calls[[i]]$call))
        df_mat <- na.omit(model.matrix(formula(model_call), data = df))
        stepcall <- fwbw(lmvar(y = df_mat[,1],X_mu = df_mat[,-1,drop=FALSE], X_sigma = df_mat[,-1,drop=FALSE],intercept_mu = TRUE,intercept_sigma = TRUE),fun = AIC)$object
        df_mat <- model.matrix(formula(model_call), data = df)
        # Remove intercepts from model matrices
        mv_fits[[i]] <- lmvar(y = df_mat[,1],X_mu = df_mat[,names(stepcall$coefficients_mu)[-1]], X_sigma = df_mat[,names(stepcall$coefficients_sigma)[-1]],intercept_mu = FALSE,,intercept_sigma = FALSE)
        mv_fits[[i]]$model_call <- mv_calls[[i]]$call
        
      } else {
        
        stepcall <- step(lm(mv_calls[[i]]$call,data = na.omit(df[,mv_calls[[i]]$varlist])),trace = FALSE)$call 
        mv_fits[[i]] <- lm(stepcall,data = df)
      }
      
    }
  } else if(var_select == "determ"){
    # Fit models
    for(i in 1:length(mv_calls)){
      # If model type is binomial glm
      if(mv_calls[[i]]$modeltype == "binomial"){
        if(mv_calls[[i]]$prev_measure == "i"){
          ind_var <- substr(mv_calls[[i]]$varlist[1], 1, nchar(mv_calls[[i]]$varlist[1])-1)
          mv_fits[[i]] <- glm(mv_calls[[i]]$call,family = binomial,data = df[df[,ind_var] == 0,])
          mv_fits[[i]]$prev_measure <- "i"
        }      
        else {
          mv_fits[[i]] <- glm(mv_calls[[i]]$call,family = binomial,data = df)
          mv_fits[[i]]$prev_measure <- "p"
        }
      } # If model type is ordered logistic
      else if (mv_calls[[i]]$modeltype == "polr"){
        
        mv_fits[[i]] <- polr(mv_calls[[i]]$call,data = df)
        
      } # If model type is a linear regression with non-constant variances
      else if (mv_calls[[i]]$modeltype %in% c("gls","lmvar")){
        model_call <- paste("~",gsub("~","+",mv_calls[[i]]$call))
        df_mat <- model.matrix(formula(model_call), data = df)
        model_call_sigma <- paste("~",gsub("~","+",mv_calls[[i]]$call_sigma))
        df_mat_sigma <- model.matrix(formula(model_call_sigma), data = df)
        mv_fits[[i]] <- lmvar(y = df_mat[,1],X_mu = df_mat[,-1,drop=FALSE], X_sigma = df_mat_sigma[,-1,drop=FALSE],intercept_mu = TRUE,intercept_sigma = TRUE)
        mv_fits[[i]]$model_call <- mv_calls[[i]]$call
        mv_fits[[i]]$model_call_sigma <- mv_calls[[i]]$call_sigma
        
      } # If model type is a linear model
      else {
        
        mv_fits[[i]] <- lm(mv_calls[[i]]$call,data = df)
        
      }
    }
  } else if(var_select == "determ_step"){
    # Fit models
    for(i in 1:length(mv_calls)){
      # If model type is binomial glm
      if(mv_calls[[i]]$modeltype == "binomial"){
        if(mv_calls[[i]]$prev_measure == "i"){
          ind_var <- substr(mv_calls[[i]]$varlist[1], 1, nchar(mv_calls[[i]]$varlist[1])-1)
          mv_fits[[i]] <- glm(mv_calls[[i]]$call,family = binomial,data = df[df[,ind_var] == 0,])
          mv_fits[[i]]$prev_measure <- "i"
        }      
        else {
          mv_fits[[i]] <- glm(mv_calls[[i]]$call,family = binomial,data = df)
          mv_fits[[i]]$prev_measure <- "p"
        }
      } # If model type is ordered logistic
      else if (mv_calls[[i]]$modeltype == "polr"){
        
        mv_fits[[i]] <- polr(mv_calls[[i]]$call,data = df)
        
      } # If model type is a linear regression with non-constant variances
      else if (mv_calls[[i]]$modeltype %in% c("gls","lmvar")){
        model_call <- paste("~",gsub("~","+",mv_calls[[i]]$call))
        df_mat <- model.matrix(formula(model_call), data = df)
        model_call_sigma <- paste("~",gsub("~","+",mv_calls[[i]]$call_sigma))
        df_mat_sigma <- model.matrix(formula(model_call_sigma), data = df)
        mv_fits[[i]]$model_call_sigma <- mv_calls[[i]]$call_sigma
        stepcall <- fwbw(lmvar(y = df_mat[,1],X_mu = df_mat[,-1,drop=FALSE], X_sigma = df_mat_sigma[,-1,drop=FALSE],intercept_mu = TRUE,intercept_sigma = TRUE),fun = AIC)$object
        df_mat <- model.matrix(formula(model_call), data = df)
        # Remove intercepts from model matrices
        mv_fits[[i]] <- lmvar(y = df_mat[,1],X_mu = df_mat[,names(stepcall$coefficients_mu)[-1]], X_sigma = df_mat_sigma[,-1,drop=FALSE],intercept_mu = FALSE,,intercept_sigma = FALSE)
        mv_fits[[i]]$model_call <- mv_calls[[i]]$call
      } # If model type is a linear model
      else {
        
        mv_fits[[i]] <- lm(mv_calls[[i]]$call,data = df)
        
      }
    }
  }
  for(i in 1:length(mv_calls)){
    names(mv_fits)[i] <- mv_calls[[i]]$varlist[1]
  }
  return(mv_fits)
}
