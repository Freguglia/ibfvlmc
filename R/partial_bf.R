#' @export
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join mutate arrange desc
#' @importFrom Brobdingnag as.brob
expected_Q <- function(ztest, ztrain, nsamples = 20000, burnin = 10000,
                       Hmax = 8, alpha = 1/2, logpenalty = 2, renewal = NULL){
  if(is.null(renewal)) renewal <- numeric()
  m <- length(unique(ztrain))
  outcpp <- ibf(ztest, ztrain, renewal, alpha,
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
                          logpenalty0 = 2, logpenalty1 = 2,
                          renewal0 = 0, renewal1 = numeric()){
  Q0 <- expected_Q(ztest, ztrain, nsamples, burnin, Hmax, alpha0,
                    logpenalty0, renewal = renewal0)$EQ
  Q1 <- expected_Q(ztest, ztrain, nsamples, burnin, Hmax, alpha1,
                    logpenalty1, renewal = renewal1)$EQ
  print(Q0)
  print(Q1)
  return(as.numeric(Q0/Q1))
}

# partial_bf(binchain[-1], binchain[[1]], Hmax = 5)
