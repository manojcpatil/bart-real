## metrics.R -- performance measures (T04).
## Conventions: `draws` is an M x Q matrix of posterior predictive draws
## (M draws, Q test points); `y` is the observed test response; `f` the true
## mean at the test points. Lower is better for every score except coverage,
## which should sit at its nominal level.

rmse <- function(pred, truth) sqrt(mean((pred - truth)^2))

#' Equal-tailed intervals from draws.
draw_intervals <- function(draws, level = 0.95) {
  a <- (1 - level) / 2
  q <- apply(draws, 2, quantile, probs = c(a, 1 - a), names = FALSE)
  list(lower = q[1, ], upper = q[2, ])
}

coverage <- function(lower, upper, truth) mean(truth >= lower & truth <= upper)

mean_width <- function(lower, upper) mean(upper - lower)

#' Interval score (Gneiting and Raftery 2007, eq. 43), negatively oriented.
#' alpha = 1 - level, e.g. 0.05 for a 95% interval.
interval_score <- function(lower, upper, y, level = 0.95) {
  alpha <- 1 - level
  width <- upper - lower
  below <- (2 / alpha) * (lower - y) * (y < lower)
  above <- (2 / alpha) * (y - upper) * (y > upper)
  mean(width + below + above)
}

#' CRPS from posterior draws, using
#'   CRPS = E|X - y| - 0.5 E|X - X'|.
#' The second term is computed from the sorted draws in O(M log M) rather than
#' by forming all M^2 pairs.
crps_sample_one <- function(x, y) {
  M <- length(x)
  xs <- sort(x)
  term1 <- mean(abs(xs - y))
  i <- seq_len(M)
  term2 <- 2 * sum((2 * i - M - 1) * xs) / (M^2)
  term1 - 0.5 * term2
}

crps_draws <- function(draws, y) {
  stopifnot(ncol(draws) == length(y))
  mean(vapply(seq_along(y), function(j) crps_sample_one(draws[, j], y[j]), numeric(1)))
}

#' Variable-selection scores. `selected` and `truth` are vectors of indices.
selection_scores <- function(selected, truth, p) {
  sel <- seq_len(p) %in% selected
  act <- seq_len(p) %in% truth
  tp <- sum(sel & act); fp <- sum(sel & !act); fn <- sum(!sel & act)
  precision <- if (tp + fp == 0) NA_real_ else tp / (tp + fp)
  recall    <- if (tp + fn == 0) NA_real_ else tp / (tp + fn)
  f1 <- if (is.na(precision) || is.na(recall) || precision + recall == 0) 0
        else 2 * precision * recall / (precision + recall)
  list(tp = tp, fp = fp, fn = fn, precision = precision, recall = recall, f1 = f1)
}

#' Monte Carlo standard errors across replications (Morris et al. 2019).
mcse_mean <- function(v) {
  v <- v[is.finite(v)]
  sd(v) / sqrt(length(v))
}

#' Number of replications needed for a target MCSE of a mean, given a pilot sd.
nsim_for_mcse <- function(pilot_sd, target_mcse) ceiling((pilot_sd / target_mcse)^2)

#' All predictive scores for one fitted method on one test set.
#' `f_draws` (optional) are draws of the regression function, used for the
#' secondary credible-interval coverage of f.
score_prediction <- function(pred_draws, y, f = NULL, f_draws = NULL, level = 0.95) {
  pi <- draw_intervals(pred_draws, level)
  point <- colMeans(pred_draws)
  out <- list(
    rmse_y   = rmse(point, y),
    rmse_f   = if (is.null(f)) NA_real_ else rmse(point, f),
    coverage = coverage(pi$lower, pi$upper, y),
    width    = mean_width(pi$lower, pi$upper),
    int_score = interval_score(pi$lower, pi$upper, y, level),
    crps     = crps_draws(pred_draws, y)
  )
  if (!is.null(f_draws) && !is.null(f)) {
    ci <- draw_intervals(f_draws, level)
    out$coverage_f <- coverage(ci$lower, ci$upper, f)
    out$width_f <- mean_width(ci$lower, ci$upper)
  } else {
    out$coverage_f <- NA_real_
    out$width_f <- NA_real_
  }
  out
}
