#' @export
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate arrange desc
#' @importFrom Brobdingnag as.brob
expected_Q <- function(ztest, ztrain, nsamples = 20000, burnin = 10000,
                       Hmax = 8, alpha = 1/2, logpenalty = 2, renewal = NULL,
                       allowedMatrix = NULL){
  if(is.null(renewal)) renewal <- integer()
  m <- length(unique(ztrain))
  if(is.null(allowedMatrix)) allowedMatrix <- matrix(TRUE, nrow = m, ncol = m)
  outcpp <- ibf(ztest, ztrain, renewal, allowedMatrix, alpha,
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
                           allowedMatrix = NULL){
  if(is.null(renewal)) renewal <- integer()
  m <- length(unique(ztrain))
  if(is.null(allowedMatrix)) allowedMatrix <- matrix(TRUE, nrow = m, ncol = m)
  outcpp <- ibf_comp(ztest, ztrain, renewal, allowedMatrix, alpha,
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
                          allowedMatrix = NULL){
  Q0 <- expected_Q(ztest, ztrain, nsamples, burnin, Hmax, alpha0,
                    logpenalty0, renewal = renewal0, allowedMatrix)
  Q1 <- expected_Q(ztest, ztrain, nsamples, burnin, Hmax, alpha1,
                    logpenalty1, renewal = renewal1, allowedMatrix)
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
                          allowedMatrix = NULL){
  Q0 <- expected_Q(ztest, ztrain, nsamples, burnin, Hmax, alpha0,
                    logpenalty0, renewal = renewal, allowedMatrix = allowedMatrix)
  Q1 <- expected_Q_cmp(ztest, ztrain, nsamples, burnin, Hmax, alpha1,
                    logpenalty1, renewal = renewal, allowedMatrix = allowedMatrix)
  return(list(pbf = as.numeric(Q0$EQ/Q1$EQ),
         Q0 = Q0,
         Q1 = Q1))
}

# partial_bf_cmp(binchain[-1], binchain[[1]], Hmax = 8)
