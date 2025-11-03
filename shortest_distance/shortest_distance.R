# Function for applying shortest distance
#
# Inputs:
# df -
# prev_cols -
# next_cols - 
# rank_cols -
# replace_cols -
# inf_cols - 
# scale_within_cohorts -
sd_method <- function(df, prev_cols = NULL,next_cols = NULL, rank_cols = NULL,replace_cols = FALSE, inf_cols = NULL,scale_within_cohorts = FALSE){
  library("tidyr","fields","clue")
  if(replace_cols == FALSE){
    new_cols <- paste0("pred_",next_cols)
  } else {
    new_cols <- next_cols
    df[,new_cols] <- NA
  }
  if(length(rank_cols) == 0){
    distmat <- rdist(x1 = scale(df[,prev_cols]), x2 = scale(df[,next_cols]), compact = FALSE)
    distmat <- solve_LSAP(distmat, maximum = FALSE)
    df[distmat,new_cols] <- df[,next_cols]
  } 
  else if (length(rank_cols) == 1) {
    rank_combs <- unique(df[,rank_cols])
    for (i in 1:length(rank_combs)){
      row_inds <- which(df[,rank_cols] == rank_combs[i])
      df_prev <- df[row_inds,prev_cols]
      df_next <- df[row_inds,next_cols]
      if(scale_within_cohorts == TRUE){
        if(length(row_inds) >= 3){
          distmat <- rdist(x1 = scale(df_prev), x2 = scale(df_next), compact = FALSE)
        } else{
          distmat <- rdist(x1 = df_prev, x2 = df_next, compact = FALSE)
        }
      } else {
        distmat <- rdist(x1 = scale(df[,prev_cols])[row_inds,], x2 = scale(df[,next_cols])[row_inds,], compact = FALSE) 
      }
      distmat <- solve_LSAP(distmat, maximum = FALSE)
      df[row_inds[distmat],new_cols] <- df_next
    }
  }
  else {
    rank_combs <- as.data.frame(crossing(df[,rank_cols]))
    for (i in 1:nrow(rank_combs)){
      print(i)
      row_inds <- which(rowSums(df[,rank_cols] == c(rank_combs[i,])) == length(rank_cols))
      df_prev <- df[row_inds,prev_cols]
      df_next <- df[row_inds,next_cols]      
      if(scale_within_cohorts == TRUE){
        if(length(row_inds) >= 3){
          distmat <- rdist(x1 = scale(df_prev), x2 = scale(df_next), compact = FALSE)
        } else{
          distmat <- rdist(x1 = df_prev, x2 = df_next, compact = FALSE)
        }
      } else {
        distmat <- rdist(x1 = scale(df[,prev_cols])[row_inds,], x2 = scale(df[,next_cols])[row_inds,], compact = FALSE) 
      }
      distmat <- solve_LSAP(distmat, maximum = FALSE)
      df[row_inds[distmat],new_cols] <- df_next
    }
  }
  return(df)
}
