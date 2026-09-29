## dgp.R -- data-generating mechanisms for the simulation studies (T04).
## Every generator returns a list with x (matrix), y (vector), f (true mean
## vector) and sigma (the error sd used, a vector when heteroscedastic).

## Friedman #1 (Friedman 1991). Only x1..x5 affect the response; the remaining
## p - 5 columns are irrelevant.
f_friedman <- function(x) {
  10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 + 10 * x[, 4] + 5 * x[, 5]
}

## Sparse linear function: three active predictors.
f_linear <- function(x) {
  3 * x[, 1] + 2 * x[, 2] - 4 * x[, 3]
}

## Step function with an interaction: piecewise constant, favours trees.
f_step <- function(x) {
  5 * (x[, 1] > 0.5) + 3 * (x[, 2] > 0.7) + 6 * (x[, 1] > 0.5) * (x[, 3] > 0.4)
}

mean_functions <- list(friedman = f_friedman, linear = f_linear, step = f_step)

## Error models.
##   gaussian : e ~ N(0, sigma^2)                      (correctly specified)
##   hetero   : e ~ N(0, (sigma * (0.5 + |x1 - 0.5| * 2))^2)  (variance depends on x)
##   t3       : scaled t with 3 df, same variance as gaussian (heavy tails)
draw_errors <- function(n, x, sigma, error = c("gaussian", "hetero", "t3")) {
  error <- match.arg(error)
  if (error == "gaussian") {
    s <- rep(sigma, n)
    return(list(e = rnorm(n, 0, s), s = s))
  }
  if (error == "hetero") {
    s <- sigma * (0.5 + 2 * abs(x[, 1] - 0.5))
    return(list(e = rnorm(n, 0, s), s = s))
  }
  ## t3 has variance 3; divide by sqrt(3) so the marginal variance is sigma^2.
  s <- rep(sigma, n)
  list(e = sigma * rt(n, df = 3) / sqrt(3), s = s)
}

#' Generate one data set.
#' @param n number of observations
#' @param p total number of predictors (>= 5 for friedman)
#' @param sigma noise scale
#' @param mean_fn one of names(mean_functions)
#' @param error one of "gaussian", "hetero", "t3"
gen_data <- function(n, p, sigma, mean_fn = "friedman", error = "gaussian") {
  stopifnot(mean_fn %in% names(mean_functions))
  if (mean_fn == "friedman" && p < 5) stop("friedman needs p >= 5")
  x <- matrix(runif(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("x", seq_len(p))
  f <- mean_functions[[mean_fn]](x)
  err <- draw_errors(n, x, sigma, error)
  list(x = x, y = f + err$e, f = f, sigma = err$s)
}

## Which predictors are truly active, for scoring variable selection (RQ5).
active_vars <- function(mean_fn = "friedman") {
  switch(mean_fn,
         friedman = 1:5,
         linear   = 1:3,
         step     = 1:3,
         stop("unknown mean function"))
}

#' Training and test data from the same mechanism, for one replication.
gen_train_test <- function(n, p, sigma, mean_fn = "friedman",
                           error = "gaussian", n_test = 1000) {
  list(train = gen_data(n, p, sigma, mean_fn, error),
       test  = gen_data(n_test, p, sigma, mean_fn, error))
}
