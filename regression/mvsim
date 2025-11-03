# Function which takes the input of the mvfits function to simulate synthetic individuals in a regression transition model
#
# Inputs: 
# mvfits - output from mvfits function
# simpop - Baseline simulated population
# stochastic 
# bincutoff = cutoff in binomial predictions

require(stringr)
mvsim <- function(mvfits, simpop, stochastic = FALSE,bincutoff = 0.5,preserve_old = FALSE, preserve_new = TRUE){
  # If deterministic
  if(stochastic == FALSE){
    for(i in 1:length(mvfits)){
      if(class(mvfits[[i]])[1] == "glm"){
        if (mvfits[[i]]$family$family == "binomial"){ 
          if(mvfits[[i]]$prev_measure == "i"){
            ind_var <- substr(names(mvfits)[i], 1, nchar(names(mvfits)[i])-1) 
            simpop[,names(mvfits)[i]] <- simpop[,ind_var]
            simpop[,names(mvfits)[i]][simpop[,names(mvfits)[i]] == 0] <- as.numeric(predict(mvfits[[i]],simpop[simpop[,names(mvfits)[i]] == 0,] , type = "response") > bincutoff)
            simpop[,names(mvfits)[i]] <- factor(simpop[,names(mvfits)[i]])
          } else {
          simpop[,names(mvfits)[i]] <- factor(as.numeric(predict(mvfits[[i]],simpop, type = "response") > bincutoff))
          }
        } else {
          simpop[,names(mvfits)[i]] <- predict(mvfits[[i]],simpop) 
        }
      } 
      else if (class(mvfits[[i]])[1] == "lmvar"){
        mu_cols <- names(mvfits[[i]]$coefficients_mu)[-1]
        sigma_cols <- names(mvfits[[i]]$coefficients_sigma)[-1]
        model_call <- as.formula(gsub(".*~","~", mvfits[[i]]$model_call))
        if("model_call_sigma" %in% names(mvfits[[i]])){
          model_call_sigma <- as.formula(gsub(".*~","~", mvfits[[i]]$model_call_sigma))
          simpop[,names(mvfits)[i]] <- predict(mvfits[[i]], X_mu = model.matrix(model_call,simpop)[,mu_cols],X_sigma  =  model.matrix(model_call_sigma,simpop)[,sigma_cols])[,1]
        } else{
          simpop[,names(mvfits)[i]] <- predict(mvfits[[i]], X_mu = model.matrix(model_call,simpop)[,mu_cols],X_sigma  =  model.matrix(model_call,simpop)[,sigma_cols])[,1]
        }
      } 
      else {
        simpop[,names(mvfits)[i]] <- predict(mvfits[[i]],simpop)
      }
    }
  }
  else { 
    for(i in 1:length(mvfits)){
      if(class(mvfits[[i]])[1] == "glm"){
        if (mvfits[[i]]$family$family == "binomial"){ 
          if(mvfits[[i]]$prev_measure == "i"){
            ind_var <- substr(names(mvfits)[i], 1, nchar(names(mvfits)[i])-1) 
            simpop[,names(mvfits)[i]] <- simpop[,ind_var]
            simpop[,names(mvfits)[i]][simpop[,names(mvfits)[i]] == 0] <- rbinom(n = length(simpop[,names(mvfits)[i]][simpop[,names(mvfits)[i]] == 0]),size = 1, prob = predict(mvfits[[i]],simpop[simpop[,names(mvfits)[i]] == 0,], type = "response"))
            simpop[,names(mvfits)[i]] <- factor(simpop[,names(mvfits)[i]])
          } else{
            simpop[,names(mvfits)[i]] <- factor(rbinom(n = nrow(simpop),size = 1, prob = predict(mvfits[[i]],simpop, type = "response")))
          }
        } else {
          simpop[,names(mvfits)[i]] <- predict(mvfits[[i]],simpop) 
        }
      } 
      else if (class(mfits[[i]])[1] == "polr"){
        pred_vec <- predict(mvfits[[i]],simpop,type = "p")
        simpop[,names(mvfits)[i]] <- apply(pred_vec,1, function(x) sample(colnames(pred_vec),1, prob = x))
        simpop[,names(mvfits)[i]] <- factor(simpop[,names(mvfits)[i]])
      } 
      else if ("lmvar" %in% class(mvfits[[i]])){
        mu_cols <- names(mvfits[[i]]$coefficients_mu)
        sigma_cols <- names(mvfits[[i]]$coefficients_sigma)
        model_call <- as.formula(gsub(".*~","~", mvfits[[i]]$model_call))
        if("model_call_sigma" %in% names(mvfits[[i]])){
          model_call_sigma <- as.formula(gsub(".*~","~", mvfits[[i]]$model_call_sigma))
        } else{
          model_call_sigma <- as.formula(gsub(".*~","~", mvfits[[i]]$model_call))
        }
        pred_vec <- predict(mvfits[[i]],X_mu = model.matrix(model_call,simpop)[,mu_cols,drop=FALSE],X_sigma  = model.matrix(model_call_sigma,simpop)[,sigma_cols,drop=FALSE])[,1] 
        rand_vec <- rnorm(nrow(simpop), mean = 0, sd = predict(mvfits[[i]],
                                                     X_mu = model.matrix(model_call,simpop)[,mu_cols,drop=FALSE],
                                                     X_sigma  = model.matrix(model_call_sigma,simpop)[,sigma_cols,drop=FALSE])[,2])
        simpop[,names(mvfits)[i]] <- pred_vec + rand_vec
      } 
      else {
        simpop[,names(mvfits)[i]] <- predict(mvfits[[i]],simpop) + rnorm(nrow(simpop), mean = 0, sd = sigma(mvfits[[i]]))
      }
    }
  }
  if(preserve_old == FALSE){
    simpop[,str_sub(names(mvfits),end = -2)] <- simpop[,names(mvfits)]
    if(preserve_new == FALSE){
      simpop[,names(mvfits)] <- NULL
    }
  }
  return(simpop)
}
