
# Amputation of the initial (simulated) dataset
# dataset = dataset to be amputed
# prop_NA = desired proportion of missing values in the amputed dataset

#' Title
#'
#' @param dataset 
#' @param prop_NA 
#'
#' @return
#' @export
#'
#' @examples
#' 1+1
MVgen <- function(dataset, prop_NA){
  data_NA <- dataset
  n_NA <- prop_NA * prod(dim(dataset))
  row_NA <- sample(1:nrow(dataset), n_NA, replace = TRUE)
  col_NA <- sample(1:ncol(dataset), n_NA, replace = TRUE)
  for (i in 1:n_NA) data_NA[row_NA[i], col_NA[i]] <- NA
  return(data_NA)
}

