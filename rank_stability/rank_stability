# Function for applying the rank stability method to the data
#
# Inputs:
# df - Data with previous and next risk factor values
# prev_cols - A vector with the names of the columns containing the previous risk factor values
# next_cols - A vector with the names of the columns containing the next risk factor values
# strata_cols - A vector with the names of the columns with the ranking strata
# replace_cols - Replace values in next cols or create new columns
# rank_cols <- a pair of columns used to match an entire row
# inf_cols - Strata

rank_stab<- function(df,next_df,id_col, vars = NULL, strata_cols = NULL,
                     prev_cols = NULL,next_cols = NULL,replace_cols = FALSE, inf_cols = NULL,rank_cols = NULL){
  #### If ranking independently  ####
  if(is.null(rank_cols)){
    if(is.null(prev_cols)){
      if(class(prev_dfs) == "list"){
        prev_df <- df[[length(prev_dfs)]]
      }
      if(length(strata_cols) == 0 ){
        for(i in 1:length(vars)){
          prev_rank <- order(prev_df[,vars[i]])
          next_rank <-  order(next_df[,vars[i]])
          next_df[vars[i]] <- next_df[,vars[i]][next_rank][prev_rank]
          if(i == 1){
            next_df$id <- 1:nrow(next_df)
            next_df$id <- next_df$id[next_rank][prev_rank]
          }
        }
      }
      else {
        for(i in 1:length(vars)){
          prev_rank <- do.call("order",prev_df[,c(strata_cols,vars[i])])
          next_rank <-  do.call("order",df[,c(strata_cols,vars[i])])
          next_df[vars[i]] <- next_df[,vars[i]][next_rank][prev_rank]
          if(i == 1){
            next_df$id <- 1:nrow(next_df)
            next_df$id <- next_df$id[next_rank][prev_rank]
          }
        }
      }
      next_df <- next_df[order(next_df$id),]
      df <- list(df,next_df)
    } # If previous columns specified 
    else {
      if(replace_cols == FALSE){
        new_cols <- paste0("pred_",next_cols)
      } else {
        new_cols <- next_cols
      }
      if(length(strata_cols) == 0 ){
        for(i in 1:length(prev_cols)){
          prev_rank <- order(df[,prev_cols[i]])
          next_rank <-  order(df[,next_cols[i]])
          df[prev_rank,new_cols[i]] <- df[next_rank,next_cols[i]]  
        }
      }
      else {
        for(i in 1:length(prev_cols)){
          prev_rank <- do.call("order",df[,c(strata_cols,prev_cols[i])])
          next_rank <-  do.call("order",df[,c(strata_cols,next_cols[i])])
          df[prev_rank,new_cols[i]] <- df[next_rank,next_cols[i]]
        }
      }
      if(length(inf_cols) == 0){
      } else{
        inf_combs <- do.call("c", lapply(seq(inf_cols), function(i) combn(inf_cols, i, FUN = list)))
        inf_list <- vector(mode = "list", length = (length(inf_combs)+1))
        rank_diff <- abs(df[,next_cols] - df[,new_cols])
        for(j in 1:(length(inf_combs)+1))
          if(j == 1){
            inf_list[[j]] <- colMeans(rank_diff)
          } else {
            by_list <- list()
            for(i in 1:length(inf_combs[[j-1]])){
              by_list[[i]] <- df[,inf_combs[[j-1]][i]]
            }
            inf_list[[j]] <- aggregate(rank_diff, by = by_list,FUN = mean)
          }
        df <- list(df,inf_list)
      }
    }
  } 
  #### If sorting all columns one variable ####
  else {
    if(is.null(prev_cols)){
      if(class(prev_dfs) == "list"){
        prev_df <- df[[length(prev_dfs)]]
      }
      if(length(strata_cols) == 0 ){
        for(i in 1:length(vars)){
          prev_rank <- order(prev_df[,vars[i]])
          next_rank <-  order(next_df[,vars[i]])
          next_df[vars[i]] <- next_df[,vars[i]][next_rank][prev_rank]
          if(i == 1){
            next_df$id <- 1:nrow(next_df)
            next_df$id <- next_df$id[next_rank][prev_rank]
          }
        }
      }
      else {
        for(i in 1:length(vars)){
          prev_rank <- do.call("order",prev_df[,c(strata_cols,vars[i])])
          next_rank <-  do.call("order",df[,c(strata_cols,vars[i])])
          next_df[vars[i]] <- next_df[,vars[i]][next_rank][prev_rank]
          if(i == 1){
            next_df$id <- 1:nrow(next_df)
            next_df$id <- next_df$id[next_rank][prev_rank]
          }
        }
      }
      next_df <- next_df[order(next_df$id),]
      df <- list(df,next_df)
    } # If previous columns specified 
    else {
      if(replace_cols == FALSE){
        new_cols <- paste0("pred_",next_cols)
      } else {
        new_cols <- next_cols
      }
      if(length(strata_cols) == 0 ){
          prev_rank <- order(df[,rank_cols[1]])
          next_rank <-  order(df[,rank_cols[2]])
          df[prev_rank,new_cols] <- df[next_rank,next_cols]  
      }
      else {
          prev_rank <- do.call("order",df[,c(strata_cols,rank_cols[1])])
          next_rank <-  do.call("order",df[,c(strata_cols,rank_cols[2])])
          df[prev_rank,new_cols] <- df[next_rank,next_cols]
      }
      if(length(inf_cols) == 0){
      } else{
        inf_combs <- do.call("c", lapply(seq(inf_cols), function(i) combn(inf_cols, i, FUN = list)))
        inf_list <- vector(mode = "list", length = (length(inf_combs)+1))
        rank_diff <- abs(df[,next_cols] - df[,new_cols])
        for(j in 1:(length(inf_combs)+1))
          if(j == 1){
            inf_list[[j]] <- colMeans(rank_diff)
          } else {
            by_list <- list()
            for(i in 1:length(inf_combs[[j-1]])){
              by_list[[i]] <- df[,inf_combs[[j-1]][i]]
            }
            inf_list[[j]] <- aggregate(rank_diff, by = by_list,FUN = mean)
          }
        df <- list(df,inf_list)
      }
    }
  }
  return(df)
}
