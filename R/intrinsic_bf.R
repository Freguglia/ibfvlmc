#' @export
#' @importFrom furrr future_map furrr_options
#' @importFrom purrr map_dbl
intrinsic_bf <- function(z, renewal0, renewal1 = numeric(),
                         nsamples = 20000, burnin = 10000,
                         Hmax = 5,
                         alpha0 = 1/2, alpha1 = 1/2,
                         logpenalty0 = NULL, logpenalty1 = NULL,
                         seed = NULL){
  init_time <- Sys.time()
  I <- length(z)
  m <- length(unique(z[[1]]))
  if(is.null(seed)) seed <- floor(runif(1)*100000)
  if(is.null(logpenalty0)) logpenalty0 <- m - length(renewal0)
  if(is.null(logpenalty1)) logpenalty1 <- m - length(renewal1)
  
  partials <- future_map(1:I, function(i) {
    set.seed(seed + i)
    partial_bf(ztest = z[-i], ztrain = z[[i]],
               nsamples = nsamples, burnin = burnin,
               Hmax = Hmax, alpha0 = alpha0, alpha1 = alpha1,
               logpenalty0 = logpenalty0, logpenalty1 = logpenalty1,
               renewal0 = renewal0, renewal1 = renewal1)
  }, .options = furrr_options(seed = NULL))
  pbfs <- map_dbl(partials, "pbf")
  ibf_arithmetic <- mean(pbfs)
  ibf_geometric <- prod(pbfs^(1/I))
  total_time <- Sys.time() - init_time
  return(list(ibf_arithmetic = ibf_arithmetic,
              ibf_geometric = ibf_geometric,
              partials = partials,
              total_time = total_time))}
