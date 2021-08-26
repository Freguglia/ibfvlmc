#' @export
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate arrange desc
#' @importFrom Brobdingnag as.brob
expected_Q <- function(ztest, ztrain, nsamples = 20000, burnin = 10000,
                       Hmax = 8, alpha = 1/2, logpenalty = 2, renewal = NULL,
                       prohibited = NULL){
  if(is.null(renewal)) renewal <- numeric()
  m <- length(unique(ztrain))
  if(is.null(prohibited)){
    prohibited = lapply(1:m, function(x) integer(0))}
  outcpp <- ibf(ztest, ztrain, renewal, prohibited, alpha,
                logpenalty, Hmax, m, burnin, nsamples)
  out <- as_tibble(outcpp$logQ) %>%
    left_join(as_tibble(outcpp$posterior), by = "tree") %>%
    mutate(prob = count/sum(count)) %>%
    arrange(desc(prob))
  Q <- as.brob(out$logq)
  Q <- exp(Q)
  Q <- Q*out$prob
  return(list(EQ = Brobdingnag::sum(Q), posterior = out))
}

#' @export
expected_Q_cmp <- function(ztest, ztrain, nsamples = 20000, burnin = 10000,
                           Hmax = 8, alpha = 1/2, logpenalty = 2, renewal = NULL,
                           prohibited = NULL){
  if(is.null(renewal)) renewal <- numeric()
  m <- length(unique(ztrain))
  if(is.null(prohibited)){
    prohibited = lapply(1:m, function(x) integer(0))}
  outcpp <- ibf_comp(ztest, ztrain, renewal, prohibited, alpha,
                logpenalty, Hmax, m, burnin, nsamples)
  out <- as_tibble(outcpp$logQ) %>%
    left_join(as_tibble(outcpp$posterior), by = "tree") %>%
    mutate(prob = count/sum(count)) %>%
    arrange(desc(prob))
  Q <- as.brob(out$logq)
  Q <- exp(Q)
  Q <- Q*out$prob
  return(list(EQ = Brobdingnag::sum(Q), posterior = out))
}

#' @export
partial_bf <- function(ztest, ztrain,
                          nsamples = 20000, burnin = 10000,
                          Hmax = 8,
                          alpha0 = 1/2, alpha1 = 1/2,
                          logpenalty0 = 0, logpenalty1 = 0,
                          renewal0 = 0, renewal1 = numeric(),
                          prohibited = NULL){
  Q0 <- expected_Q(ztest, ztrain, nsamples, burnin, Hmax, alpha0,
                    logpenalty0, renewal = renewal0, prohibited)
  Q1 <- expected_Q(ztest, ztrain, nsamples, burnin, Hmax, alpha1,
                    logpenalty1, renewal = renewal1, prohibited)
  return(list(pbf = as.numeric(Q0$EQ/Q1$EQ),
         Q0 = Q0,
         Q1 = Q1))
}

#' @export
partial_bf_cmp<- function(ztest, ztrain,
                          nsamples = 20000, burnin = 10000,
                          Hmax = 8,
                          alpha0 = 1/2, alpha1 = 1/2,
                          logpenalty0 = 0, logpenalty1 = 0,
                          renewal = 0,
                          prohibited = NULL){
  Q0 <- expected_Q(ztest, ztrain, nsamples, burnin, Hmax, alpha0,
                    logpenalty0, renewal = renewal, prohibited = prohibited)
  Q1 <- expected_Q_cmp(ztest, ztrain, nsamples, burnin, Hmax, alpha1,
                    logpenalty1, renewal = renewal, prohibited = prohibited)
  return(list(pbf = as.numeric(Q0$EQ/Q1$EQ),
         Q0 = Q0,
         Q1 = Q1))
}

# partial_bf_cmp(binchain[-1], binchain[[1]], Hmax = 8)
