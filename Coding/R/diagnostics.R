## diagnostics.R -- MCMC diagnostics for RQ6 (T04).

#' Effective sample size of a single chain, by the initial positive sequence
#' estimator used in coda. Falls back to a simple autocorrelation sum when coda
#' is unavailable.
ess <- function(chain) {
  chain <- as.numeric(chain)
  if (requireNamespace("coda", quietly = TRUE))
    return(as.numeric(coda::effectiveSize(coda::as.mcmc(chain))))
  n <- length(chain)
  a <- acf(chain, plot = FALSE, lag.max = min(n - 1, 200))$acf[-1]
  keep <- which(a < 0)[1]
  if (is.na(keep)) keep <- length(a)
  rho <- if (keep > 1) sum(a[seq_len(keep - 1)]) else 0
  n / (1 + 2 * rho)
}

#' Geweke z statistic comparing the first 10 per cent and last 50 per cent of a
#' chain. Values far from zero suggest the chain has not settled.
geweke_z <- function(chain, frac1 = 0.1, frac2 = 0.5) {
  chain <- as.numeric(chain)
  if (requireNamespace("coda", quietly = TRUE))
    return(as.numeric(coda::geweke.diag(coda::as.mcmc(chain), frac1, frac2)$z))
  n <- length(chain)
  a <- chain[seq_len(floor(frac1 * n))]
  b <- chain[(n - floor(frac2 * n) + 1):n]
  (mean(a) - mean(b)) / sqrt(var(a) / length(a) + var(b) / length(b))
}

#' Summary of a sigma chain: ESS, ESS per second, Geweke z, posterior mean.
sigma_diagnostics <- function(sigma_draws, elapsed_sec) {
  e <- ess(sigma_draws)
  list(ess = e,
       ess_per_sec = if (is.finite(elapsed_sec) && elapsed_sec > 0) e / elapsed_sec else NA_real_,
       geweke_z = geweke_z(sigma_draws),
       sigma_mean = mean(sigma_draws))
}

#' Agreement between chains: the spread of per-chain posterior means of f at the
#' test points, relative to the overall posterior sd. Values near zero indicate
#' the chains agree.
chain_disagreement <- function(f_draws, nchain) {
  M <- nrow(f_draws)
  if (nchain < 2 || M %% nchain != 0) return(NA_real_)
  idx <- split(seq_len(M), rep(seq_len(nchain), each = M / nchain))
  means <- vapply(idx, function(i) colMeans(f_draws[i, , drop = FALSE]), numeric(ncol(f_draws)))
  between <- apply(means, 1, sd)
  within <- apply(f_draws, 2, sd)
  mean(between / pmax(within, .Machine$double.eps))
}
