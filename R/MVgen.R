#' @title Amputation of a dataset
#' 
#' @description This function is designed to ampute datasets.
#' 
#' @param dataset dataset to be amputed
#' @param prop_NA desired proportion of missing values in the amputed dataset
#'
#' @return A dataset with missing values.
#' @export
#'
#' @examples
#' library(mi4p)
#' data(datasim)
#' datasim_amp <- MVgen(datasim, .2)
#' sum(is.na(datasim_amp))/prod(dim(datasim_amp))
MVgen <- function(dataset, prop_NA){
  data_NA <- dataset
  n_NA <- prop_NA * prod(dim(dataset))
  row_NA <- sample(1:nrow(dataset), n_NA, replace = TRUE)
  col_NA <- sample(1:ncol(dataset), n_NA, replace = TRUE)
  for (i in 1:n_NA) data_NA[row_NA[i], col_NA[i]] <- NA
  return(data_NA)
}

