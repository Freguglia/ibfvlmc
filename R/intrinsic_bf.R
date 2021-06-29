#' @export
#' @importFrom furrr future_map_dbl furrr_options
intrinsic_bf <- function(z, renewal0, renewal1 = numeric(),
                         nsamples = 20000, burnin = 10000,
                         Hmax = 5,
                         alpha0 = 1/2, alpha1 = 1/2,
                         logpenalty0 = NULL, logpenalty1 = NULL,
                         type = "arithmetic", seed = NULL){
  I <- length(z)
  m <- length(unique(z[[1]]))
  if(is.null(seed)) seed <- floor(runif(1)*100000)
  if(is.null(logpenalty0)) logpenalty0 <- m - length(renewal0)
  if(is.null(logpenalty1)) logpenalty1 <- m - length(renewal1)
  
  partials <- future_map_dbl(1:I, function(i) {
    set.seed(seed + i)
    partial_bf(ztest = z[-i], ztrain = z[[i]], 
               nsamples = nsamples, burnin = burnin, 
               Hmax = Hmax, alpha0 = alpha0, alpha1 = alpha1, 
               logpenalty0 = logpenalty0, logpenalty1 = logpenalty1, 
               renewal0 = renewal0, renewal1 = renewal1)
  }, .options = furrr_options(seed = NULL))
  if(type == "arithmetic"){
    return(list(mean(partials), partials))
  } else if (type == "geometric"){
    return(list(prod(partials)^{1/I}, partials))
  }
}
