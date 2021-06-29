#' @export
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
partial_bf <- function(ztest, ztrain, nsamples = 20000, burnin = 10000,
                       Hmax = 8, alpha = 1/2, logpenalty = 2, renewal = NULL){
  if(is.null(renewal)) renewal <- numeric()
  m <- length(unique(ztest))
  outcpp <- ibf(ztest, ztrain, renewal, alpha,
                logpenalty, Hmax, m, burnin, nsamples)
  out <- as_tibble(outcpp$logQ) %>%
    left_join(as_tibble(outcpp$posterior), by = "tree")
  return(out)
}

# partial_bf(binchain[[1]], binchain[2:4])
