## methods.R -- fitting wrappers (T04).
## Every fit_* function returns a list with
##   pred_draws : M x Q matrix of draws from the predictive distribution of a
##                new y at the test points (for RQ2), or a 2-row matrix of
##                interval endpoints for methods that produce intervals
##                directly (flagged by interval_only = TRUE)
##   f_draws    : M x Q draws of the regression function, or NULL
##   importance : named numeric vector of variable importance, or NULL
##   selected   : indices selected as relevant, or NULL
##   time       : elapsed seconds
##   extra      : method-specific diagnostics (sigma draws, acceptance, etc.)

need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("package '%s' is not installed; see Coding/R/install_packages.R", pkg))
}

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, time = proc.time()[["elapsed"]] - t0)
}

## ---------------------------------------------------------------- BART ----
#' BART via the BART package (birth and death proposals only).
#' @param sparse TRUE gives DART, the Dirichlet splitting prior (Linero 2018).
fit_bart <- function(x, y, x_test, ndpost = 2000, nskip = 1000, ntree = 200,
                     sparse = FALSE, keepevery = 1, seed = NULL) {
  need("BART")
  if (!is.null(seed)) set.seed(seed)
  r <- timed(BART::wbart(x.train = x, y.train = y, x.test = x_test,
                         ndpost = ndpost, nskip = nskip, ntree = ntree,
                         sparse = sparse, keepevery = keepevery,
                         printevery = ndpost + nskip + 1))
  fit <- r$value
  f_draws <- fit$yhat.test                       # M x Q draws of f(x)
  ## wbart returns sigma for every sweep (nskip + ndpost * keepevery), while the
  ## tree draws are thinned. Take the sweeps that the kept draws correspond to,
  ## so sigma and f come from the same states of the chain.
  sigma_post <- if (length(fit$sigma) >= nskip + ndpost * keepevery)
    fit$sigma[nskip + seq_len(ndpost) * keepevery]
  else tail(fit$sigma, ndpost)
  ## Predictive draws for a new y: add noise using each draw's own sigma.
  pred_draws <- f_draws + matrix(rnorm(length(f_draws), 0, rep(sigma_post, ncol(f_draws))),
                                 nrow = nrow(f_draws))
  imp <- if (!is.null(fit$varcount.mean)) fit$varcount.mean / sum(fit$varcount.mean) else NULL
  sel <- NULL
  if (sparse && !is.null(fit$varprob)) {
    ## Median probability model: inclusion probability above 0.5.
    incl <- colMeans(fit$varcount > 0)
    sel <- which(incl > 0.5)
  }
  list(pred_draws = pred_draws, f_draws = f_draws, importance = imp,
       selected = sel, time = r$time,
       extra = list(sigma = sigma_post, sigma_all = fit$sigma,
                    varprob = fit$varprob.mean))
}

## --------------------------------------------------------------- dbarts ----
#' BART via dbarts, which also uses the change proposal and multiple chains.
fit_dbarts <- function(x, y, x_test, ndpost = 2000, nskip = 1000, ntree = 200,
                       nchain = 4, keepevery = 1, seed = NULL) {
  need("dbarts")
  if (!is.null(seed)) set.seed(seed)
  ## When p >= n the default sigma estimate from least squares is undefined, so
  ## fall back to sd(y), which is what the BART package does internally.
  args <- list(x.train = x, y.train = y, x.test = x_test,
               ndpost = ndpost, nskip = nskip, ntree = ntree, keepevery = keepevery,
               nchain = nchain, keeptrees = FALSE, verbose = FALSE)
  if (ncol(x) >= nrow(x) - 1) args$sigest <- sd(y)
  r <- timed(do.call(dbarts::bart, args))
  fit <- r$value
  f_draws <- fit$yhat.test
  if (length(dim(f_draws)) == 3)                   # chains x draws x points
    f_draws <- matrix(aperm(f_draws, c(2, 1, 3)), ncol = dim(f_draws)[3])
  sigma_post <- as.vector(fit$sigma)
  sigma_post <- tail(sigma_post, nrow(f_draws))
  pred_draws <- f_draws + matrix(rnorm(length(f_draws), 0, rep(sigma_post, ncol(f_draws))),
                                 nrow = nrow(f_draws))
  list(pred_draws = pred_draws, f_draws = f_draws,
       importance = if (!is.null(fit$varcount)) colMeans(as.matrix(fit$varcount)) else NULL,
       selected = NULL, time = r$time,
       extra = list(sigma = sigma_post, nchain = nchain))
}

## ------------------------------------------------- random forest / QRF ----
#' Random forest point predictions with quantile regression forest intervals.
fit_rf <- function(x, y, x_test, num.trees = 500, level = 0.95, seed = NULL) {
  need("ranger")
  if (!is.null(seed)) set.seed(seed)
  df <- data.frame(y = y, x)
  r <- timed(ranger::ranger(y ~ ., data = df, num.trees = num.trees,
                            quantreg = TRUE, importance = "permutation"))
  fit <- r$value
  a <- (1 - level) / 2
  qs <- predict(fit, data = data.frame(x_test), type = "quantiles",
                quantiles = c(a, 0.5, 1 - a))$predictions
  point <- predict(fit, data = data.frame(x_test))$predictions
  list(pred_draws = NULL, f_draws = NULL,
       interval_only = TRUE,
       point = point, lower = qs[, 1], upper = qs[, 3],
       importance = fit$variable.importance, selected = NULL,
       time = r$time, extra = list())
}

## ------------------------------------------- XGBoost with conformal PI ----
#' Tuned XGBoost with split-conformal prediction intervals.
#' A calibration fraction of the training data is held out to set the width.
fit_xgb <- function(x, y, x_test, level = 0.95, cal_frac = 0.3,
                    grid = NULL, nrounds_max = 500, seed = NULL) {
  need("xgboost")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(grid))
    grid <- expand.grid(eta = c(0.05, 0.1, 0.3), max_depth = c(3, 6),
                        subsample = 0.8, colsample_bytree = 0.8)
  n <- nrow(x)
  idx_cal <- sample.int(n, max(20, floor(cal_frac * n)))
  idx_fit <- setdiff(seq_len(n), idx_cal)
  r <- timed({
    ## xgboost >= 3.0 requires an xgb.DMatrix and xgb.train.
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x[idx_fit, , drop = FALSE]),
                                   label = y[idx_fit])
    best <- NULL
    for (i in seq_len(nrow(grid))) {
      par <- list(eta = grid$eta[i], max_depth = grid$max_depth[i],
                  subsample = grid$subsample[i],
                  colsample_bytree = grid$colsample_bytree[i],
                  objective = "reg:squarederror", nthread = 1)
      cv <- xgboost::xgb.cv(params = par, data = dtrain, nrounds = nrounds_max,
                            nfold = 5, early_stopping_rounds = 20, verbose = 0)
      log <- as.data.frame(cv$evaluation_log)
      col <- grep("^test_.*mean$", names(log), value = TRUE)[1]
      score <- min(log[[col]])
      if (is.null(best) || score < best$score)
        best <- list(score = score, par = par,
                     nrounds = if (!is.null(cv$best_iteration)) cv$best_iteration
                               else which.min(log[[col]]))
    }
    model <- xgboost::xgb.train(params = best$par, data = dtrain,
                                nrounds = best$nrounds, verbose = 0)
    ## Split-conformal half width from absolute residuals on the calibration set.
    res <- abs(y[idx_cal] -
               predict(model, xgboost::xgb.DMatrix(as.matrix(x[idx_cal, , drop = FALSE]))))
    k <- ceiling((length(res) + 1) * level)
    half <- sort(res)[min(k, length(res))]
    list(model = model, half = half, best = best)
  })
  point <- predict(r$value$model, xgboost::xgb.DMatrix(as.matrix(x_test)))
  imp <- xgboost::xgb.importance(model = r$value$model)
  importance <- setNames(rep(0, ncol(x)), colnames(x))
  if (!is.null(imp)) importance[imp$Feature] <- imp$Gain
  list(pred_draws = NULL, f_draws = NULL, interval_only = TRUE,
       point = point, lower = point - r$value$half, upper = point + r$value$half,
       importance = importance, selected = NULL, time = r$time,
       extra = list(params = r$value$best$par, nrounds = r$value$best$nrounds))
}

## ------------------------------------------------ variable selection -----
#' Permutation thresholds for BART inclusion proportions (Bleich et al. 2014).
#' Implemented through bartMachine, which provides the published procedure.
select_bart_permute <- function(x, y, num_permute = 50, alpha = 0.05,
                                num_trees = 20, type = c("global_se", "local", "global_max")) {
  need("bartMachine")
  type <- match.arg(type)
  bm <- bartMachine::bartMachine(data.frame(x), y, num_trees = num_trees, verbose = FALSE)
  vs <- bartMachine::var_selection_by_permute(bm, num_permute_samples = num_permute,
                                              alpha = alpha, plot = FALSE)
  nm <- switch(type,
               local = vs$important_vars_local_names,
               global_max = vs$important_vars_global_max_names,
               global_se = vs$important_vars_global_se_names)
  match(nm, colnames(x))
}

#' Top-k importance selection, used for random forests and XGBoost.
select_top_k <- function(importance, k) {
  if (is.null(importance)) return(integer(0))
  order(importance, decreasing = TRUE)[seq_len(min(k, length(importance)))]
}
