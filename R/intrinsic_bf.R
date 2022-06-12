#' @importFrom furrr future_map furrr_options
#' @importFrom purrr map_dbl
intrinsic_bf <- function(z, renewal0, renewal1 = numeric(),
                         nsamples = 20000, burnin = 10000,
                         Hmax = 5,
                         allowedMatrix = NULL,
                         alpha0 = 1/2, alpha1 = 1/2,
                         logpenalty0 = 0, logpenalty1 = 0,
                         seed = NULL, subset_size = 2){
  init_time <- Sys.time()
  I <- length(z)
  m <- length(unique(z[[1]]))
  if(is.null(seed)) seed <- floor(runif(1)*100000)
  if(is.null(logpenalty0)) logpenalty0 <- m - length(renewal0)
  if(is.null(logpenalty1)) logpenalty1 <- m - length(renewal1)
  
  cbn <- combn(1:I, (I-subset_size))
  cbn <- unname(split(cbn, rep(1:ncol(cbn), each = nrow(cbn))))
  progressr::with_progress({
    p <- progressr::progressor(steps = length(cbn))
    partials <- future_map(1:length(cbn), function(i) {
      p()
      set.seed(seed + i)
      partial_bf(ztest = z[-cbn[[i]]], ztrain = z[cbn[[i]]],
                 nsamples = nsamples, burnin = burnin,
                 Hmax = Hmax, alpha0 = alpha0, alpha1 = alpha1,
                 logpenalty0 = logpenalty0, logpenalty1 = logpenalty1,
                 renewal0 = renewal0, renewal1 = renewal1,
                 allowedMatrix = allowedMatrix)
    }, .options = furrr_options(seed = NULL)) 
  })
  pbfs <- map_dbl(partials, "pbf")
  ibf_arithmetic <- mean(pbfs)
  ibf_geometric <- prod(pbfs^(1/length(cbn)))
  total_time <- Sys.time() - init_time
  return(list(ibf_arithmetic = ibf_arithmetic,
              ibf_geometric = ibf_geometric,
              partials = partials,
              total_time = total_time))}

#' @title Intrinsic Bayes Factor for Renewal State Testing
#' 
#' @param z A `list` object containing independent VLMC sequences in each element.
#' @param nsamples Size of context-tree chains in each internal Metropolis-Hastings algorithm run.
#' @param burnin Number of Metropolis-Hastings algorithm steps executed before starting to record states.
#' @param Hmax Maximum height of the context tree considered.
#' @param allowedMatrix A logical `matrix` specifying which transitions are allowed (i.e., have positive probability). 
#' The row `i` and column `j` indicates whether a transition from (i-1) to (j-1) is allowed. 
#' `NULL` can be used to indicate that every transition is allowed. 
#' @param renewal Which states are considered renewal states under the H0 hypothesis.
#' @param alpha0,alpha1 The Dirichlet parameter to be considered in the prior distribution of probabilities, under
#' the H0 and H1 hypothesis respectively. 
#' @param logpenalty0,logpenalty1 Penalty considered in the context tree prior distribution for growing a new branch for the current tree. `0` represents a uniform context tree prior distribution.
#' @param seed Random seed to be used as a reference for the Partial Bayes Factors computations. Passing this argument instead of using the `set.seed()` function is required to ensure reproducibility because of parallel computations carried in the function.
#' @param subset_size Number of samples used in each training set. Increasing this value may drastically increases the number of combinations of sequences used as training samples, it is not recommended to use values larger than `2`.
#' 
#' @return A `list` object containing:
#'   * `ibf_arithmetic`: Arithmetic Intrinsic Bayes Factor computed.
#'   * `ibf_geometric`: Geometric Intrinsic Bayes Factor computed.
#'   * `partials`: A list with the output of each Partial Bayes Factor computed.
#'   * `total_time`: Total time elapsed. 
#' 
#' @details 
#' 
#' In order to fully take advantage of parallel computing in this function, the user should set a `plan` from the `future` package.
#' 
#' ```
#' library(future)
#' plan(multisession)
#' ```
#' 
#' @author Victor Freguglia
#' @export
intrinsic_bf_cmp <- function(z, renewal,
                         nsamples = 20000, burnin = 10000,
                         Hmax = 5,
                         allowedMatrix = NULL,
                         alpha0 = 1/2, alpha1 = 1/2,
                         logpenalty0 = 0, logpenalty1 = 0,
                         seed = NULL, subset_size = 2,
                         max_subsets_evaluated = NULL){
  init_time <- Sys.time()
  I <- length(z)
  m <- length(unique(z[[1]]))
  if(is.null(seed)) seed <- floor(runif(1)*100000)
  if(is.null(logpenalty0)) logpenalty0 <- m - length(renewal)
  if(is.null(logpenalty1)) logpenalty1 <- m
  
  subset_number <- choose(I, subset_size)
  if(subset_number > 10000 && is.null(max_subsets_evaluated)){
    stop("Too many combinations and no max_subsets_evaluated set.")
  } else if(subset_number > 10000){
    cbn <- lapply(seq_len(max_subsets_evaluated),
                  function(x) sample(seq_len(I), size = min(I,subset_size)))
  } else {
    cbn <- combn(1:I, (I-subset_size))
    cbn <- unname(split(cbn, rep(1:ncol(cbn), each = nrow(cbn))))
    if(!is.null(max_subsets_evaluated)){
      cbn <- cbn[sample(seq_along(cbn), size = min(length(cbn), max_subsets_evaluated))]
    }
  }
  progressr::with_progress({
    p <- progressr::progressor(steps = length(cbn))
    partials <- future_map(1:length(cbn), function(i) {
      set.seed(seed + i)
      out <- 
        partial_bf_cmp(ztest = z[-cbn[[i]]], ztrain = z[cbn[[i]]],
                       nsamples = nsamples, burnin = burnin,
                       Hmax = Hmax, alpha0 = alpha0, alpha1 = alpha1,
                       logpenalty0 = logpenalty0, logpenalty1 = logpenalty1,
                       renewal = renewal, allowedMatrix = allowedMatrix)
      p()
      return(out)
    }, .options = furrr_options(seed = NULL))})
  pbfs <- map_dbl(partials, "pbf")
  ibf_arithmetic <- mean(pbfs)
  ibf_geometric <- prod(pbfs^(1/length(cbn)))
  total_time <- Sys.time() - init_time
  return(list(ibf_arithmetic = ibf_arithmetic,
              ibf_geometric = ibf_geometric,
              partials = partials))}
