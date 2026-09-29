## pdp_residuals.R -- T11, second part. Partial dependence and residual checks.
##
## Partial dependence. BART is fitted once to each full dataset. For every
## predictor, the partial dependence function PD_j(v) = mean_i f(v, x_i,-j) is
## computed for every posterior draw of f over a grid of 20 quantiles of x_j (a
## random 300 rows stand in for the average over i). The band is the pointwise
## 95% credible interval of PD_j: it describes uncertainty about the average
## effect, not the spread of new observations. The figure shows the four
## predictors whose PD has the largest range. They are not chosen by share of
## splits: on Energy that share put orientation, whose PD is flat, in the top
## four, because split counts favour predictors with few distinct values that
## are correlated with others. Both measures are kept in pdp_<data>.csv.
##
## Residual checks. Five-fold out-of-fold predictions from BART with the Study A
## settings. In-sample residuals from BART are shrunk towards zero (the trees fit
## the training points), so every check uses out-of-fold residuals. Reported per
## dataset: the Spearman correlation of |residual| with the fitted value and the
## ratio of the residual sd in the top and bottom deciles of the fitted value
## (heteroscedasticity), skewness and excess kurtosis of the standardised
## residuals (tails), and the coverage of the 95% predictive interval within each
## decile of the fitted value. The last is the real-data counterpart of Study B:
## a model with a single sigma can cover 95% on average while over-covering
## where the noise is small and under-covering where it is large.
##
## Usage (from the BART_Research folder):
##   Rscript Application/analysis/pdp_residuals.R                    # three small datasets
##   Rscript Application/analysis/pdp_residuals.R --data california  # about 4 hours
##   Rscript Application/analysis/pdp_residuals.R --pdp-only         # skip the 5-fold fits
##   Rscript Application/analysis/pdp_residuals.R --smoke
## Out: Application/Thesis_LaTeX/figures/{pdp,resid}_<data>.pdf,
##      Application/results/pdp_<data>.csv, coverage_by_decile.csv,
##      residual_checks.csv (rows are replaced per dataset, so runs can be split)

if (!file.exists("Coding/R/methods.R")) stop("run this from the BART_Research folder")
source("Coding/R/methods.R"); source("Coding/R/metrics.R")
source("Coding/sims/design.R")
tools::psnice(value = 19)                        # keep the laptop usable
options(warn = 1)

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
smoke <- "--smoke" %in% args
data_names <- strsplit(getopt("--data", "concrete,energy,wine_red"), ",")[[1]]
MC <- DESIGN$mcmc
K_VARS <- 4; N_GRID <- 20; N_ROWS <- 300; N_FOLDS <- 5
RES <- "Application/results"; FIG <- "Application/Thesis_LaTeX/figures"
if (smoke) {
  MC <- modifyList(MC, list(nskip = 100, ndpost = 100, ntree = 50))
  N_ROWS <- 50; N_FOLDS <- 2
  RES <- FIG <- file.path(tempdir(), "pdp_smoke")
}
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

LABELS <- list(concrete = "Compressive strength (MPa)", energy = "Heating load",
               wine_red = "Quality score", california = "Median house value ($100k)")

upsert <- function(df, file, key) {
  if (file.exists(file)) {
    old <- read.csv(file, stringsAsFactors = FALSE)
    df <- rbind(old[!old[[key]] %in% unique(df[[key]]), names(df)], df)
  }
  write.csv(df, file, row.names = FALSE)
}
skew <- function(z) mean((z - mean(z))^3) / sd(z)^3
exkurt <- function(z) mean((z - mean(z))^4) / sd(z)^4 - 3

for (nm in data_names) {
  t0 <- proc.time()[["elapsed"]]
  d <- read.csv(file.path("Application/data/clean", paste0(nm, ".csv")))
  x <- as.matrix(d[setdiff(names(d), "y")]); y <- d$y; n <- nrow(x)
  seed <- 7e6 + match(nm, c("concrete", "energy", "wine_red", "california"))

  ## ---- partial dependence from one full fit ---------------------------------
  set.seed(seed)
  fit <- BART::wbart(x.train = x, y.train = y, ndpost = MC$ndpost, nskip = MC$nskip,
                     ntree = MC$ntree, nkeeptrain = 0,
                     printevery = MC$ndpost + MC$nskip + 1)
  share <- fit$varcount.mean / sum(fit$varcount.mean)
  rows <- x[sample.int(n, min(N_ROWS, n)), , drop = FALSE]
  pd <- do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
    v <- unique(quantile(x[, j], seq(0.025, 0.975, length.out = N_GRID), type = 1))
    nd <- rows[rep(seq_len(nrow(rows)), length(v)), , drop = FALSE]
    nd[, j] <- rep(v, each = nrow(rows))
    f <- predict(fit, nd, mc.cores = 1)          # draws x (rows * grid)
    g <- rep(seq_along(v), each = nrow(rows))
    pdd <- vapply(seq_along(v), function(k) rowMeans(f[, g == k, drop = FALSE]),
                  numeric(nrow(f)))              # draws x grid
    data.frame(dataset = nm, variable = colnames(x)[j], split_share = unname(share[j]),
               value = v, pd = colMeans(pdd),
               lower = apply(pdd, 2, quantile, 0.025),
               upper = apply(pdd, 2, quantile, 0.975))
  }))
  rng <- tapply(pd$pd, pd$variable, function(v) diff(range(v)))
  pd$pd_range <- rng[pd$variable]
  write.csv(pd, file.path(RES, sprintf("pdp_%s.csv", nm)), row.names = FALSE)
  rm(fit); gc()
  top_vars <- names(sort(rng, decreasing = TRUE))[seq_len(min(K_VARS, length(rng)))]

  pdf(file.path(FIG, sprintf("pdp_%s.pdf", nm)), width = 7, height = 5.5)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), mgp = c(2.3, 0.7, 0))
  yl <- range(pd$lower[pd$variable %in% top_vars], pd$upper[pd$variable %in% top_vars])
  for (v in top_vars) {
    s <- pd[pd$variable == v, ]
    plot(s$value, s$pd, type = "n", ylim = yl, xlab = v, ylab = LABELS[[nm]],
         main = sprintf("%s (%.0f%% of splits)", v, 100 * s$split_share[1]),
         cex.main = 0.95, font.main = 1)
    polygon(c(s$value, rev(s$value)), c(s$lower, rev(s$upper)),
            col = "grey85", border = NA)
    lines(s$value, s$pd, lwd = 2)
    rug(quantile(x[, v], seq(0.05, 0.95, 0.1)), ticksize = 0.03)
  }
  par(op); dev.off()

  if ("--pdp-only" %in% args) {
    cat(sprintf("%s: pdp only, %.1f min\n", nm, (proc.time()[["elapsed"]] - t0) / 60))
    next
  }

  ## ---- out-of-fold predictions ------------------------------------------------
  set.seed(seed + 1)
  fold <- sample(rep(seq_len(N_FOLDS), length.out = n))
  oof <- data.frame(fitted = NA_real_, lower = NA_real_, upper = NA_real_,
                    sigma = NA_real_)[rep(1, n), ]
  ## The folds run in parallel where forking is available (Linux, e.g. GitHub
  ## Actions, where California would otherwise approach the 6-hour job limit)
  ## and one after another on Windows. Each fold has its own seed, so the
  ## result does not depend on which.
  cores <- if (.Platform$OS.type == "windows") 1L else
    min(N_FOLDS, parallel::detectCores())
  fits <- parallel::mclapply(seq_len(N_FOLDS), function(k) {
    te <- which(fold == k)
    b <- fit_bart(x[-te, ], y[-te], x[te, , drop = FALSE], ndpost = MC$ndpost,
                  nskip = MC$nskip, ntree = MC$ntree, keepevery = MC$keepevery,
                  seed = seed + 10 + k)
    pi <- draw_intervals(b$pred_draws, DESIGN$level)
    data.frame(row = te, fitted = colMeans(b$f_draws), lower = pi$lower,
               upper = pi$upper, sigma = mean(b$extra$sigma))
  }, mc.cores = cores, mc.set.seed = FALSE)
  for (f in fits) {
    if (!is.data.frame(f)) stop("a fold failed: ", paste(f, collapse = " "))
    oof[f$row, ] <- f[c("fitted", "lower", "upper", "sigma")]
  }
  rm(fits); gc()
  res <- y - oof$fitted
  z <- res / oof$sigma
  dec <- cut(oof$fitted, unique(quantile(oof$fitted, 0:10 / 10)), include.lowest = TRUE,
             labels = FALSE)
  cov_in <- y >= oof$lower & y <= oof$upper
  bydec <- data.frame(dataset = nm, decile = sort(unique(dec)),
                      fitted_mid = tapply(oof$fitted, dec, median),
                      n = as.vector(table(dec)),
                      coverage = tapply(cov_in, dec, mean),
                      width = tapply(oof$upper - oof$lower, dec, mean),
                      resid_sd = tapply(res, dec, sd),
                      sigma_post = tapply(oof$sigma, dec, mean))
  upsert(bydec, file.path(RES, "coverage_by_decile.csv"), "dataset")
  sp <- suppressWarnings(cor.test(abs(res), oof$fitted, method = "spearman"))
  chk <- data.frame(dataset = nm, n = n, oof_rmse = sqrt(mean(res^2)),
                    sigma_post = mean(oof$sigma), oof_coverage = mean(cov_in),
                    coverage_min_decile = min(bydec$coverage),
                    coverage_max_decile = max(bydec$coverage),
                    spearman_absres_fitted = unname(sp$estimate),
                    spearman_p = sp$p.value,
                    sd_ratio_top_bottom = tail(bydec$resid_sd, 1) / bydec$resid_sd[1],
                    skewness = skew(z), excess_kurtosis = exkurt(z),
                    top_variables = paste(top_vars, collapse = ";"),
                    minutes = (proc.time()[["elapsed"]] - t0) / 60)
  upsert(chk, file.path(RES, "residual_checks.csv"), "dataset")

  pdf(file.path(FIG, sprintf("resid_%s.pdf", nm)), width = 7, height = 3.4)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 1.5, 1), mgp = c(2.3, 0.7, 0))
  plot(oof$fitted, res, pch = 16, cex = 0.4, col = adjustcolor("black", 0.35),
       xlab = "Out-of-fold fitted value", ylab = "Residual")
  abline(h = 0, col = "grey50")
  lines(bydec$fitted_mid, 1.96 * bydec$resid_sd, lty = 2)
  lines(bydec$fitted_mid, -1.96 * bydec$resid_sd, lty = 2)
  abline(h = c(-1.96, 1.96) * mean(oof$sigma), col = "grey40", lty = 3)
  plot(bydec$decile, bydec$coverage, type = "b", pch = 16,
       ylim = range(0.8, 1, bydec$coverage), xlab = "Decile of fitted value",
       ylab = "Coverage of 95% interval")
  abline(h = DESIGN$level, lty = 2, col = "grey40")
  par(op); dev.off()

  cat(sprintf("%s: %.1f min; oof coverage %.3f (deciles %.2f-%.2f); spearman %.2f\n",
              nm, chk$minutes, chk$oof_coverage, chk$coverage_min_decile,
              chk$coverage_max_decile, chk$spearman_absres_fitted))
}
