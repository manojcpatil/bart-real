## real_data.R -- T11. The four real datasets under repeated random splits.
##
## For each dataset, 20 random 80/20 train/test splits. On every split the five
## methods of the simulation (BART, DART, dbarts, quantile random forest, XGBoost
## with split-conformal intervals) are fitted with the same MCMC settings as
## Study A and scored on the held-out 20 per cent with the same metrics. There is
## no true f here, so coverage and every other score refer to the observed y.
##
## Noise columns (RQ3 and RQ5 on real data). Each split is also run with 50 and
## 200 extra columns that are known to be irrelevant. Each extra column is a copy
## of a randomly chosen real predictor with its rows randomly permuted: it has
## exactly the marginal distribution of a real predictor, so a method cannot tell
## it apart by its values, but it is independent of y and of every real predictor.
## The same split is used at every noise level, so the comparison is paired.
## Reported: the share of importance given to noise columns, and for DART the
## number of noise and real columns selected by the median probability model.
##
## California Housing (n = 20640) is the runtime stress test. It runs without
## noise columns only (CFG$noise_california), and its coverage is reported with
## and without the rows top-coded at $500,001, whose true values are unknown.
##
## Usage (from the BART_Research folder):
##   Rscript Application/analysis/real_data.R                  # everything missing
##   Rscript Application/analysis/real_data.R --data concrete  # one dataset
##   Rscript Application/analysis/real_data.R --splits 1:5
##   Rscript Application/analysis/real_data.R --dry
##   Rscript Application/analysis/real_data.R --smoke          # tiny settings
##   Rscript Application/analysis/real_data.R --timing         # seconds per sweep
##   Rscript Application/analysis/real_data.R --gather
##
## Out: Application/results/raw/<data>_<noise>_<split>.rds, then
##      Application/results/real_data.csv and real_summary.csv

if (!file.exists("Coding/R/methods.R")) stop("run this from the BART_Research folder")
source("Coding/R/methods.R"); source("Coding/R/metrics.R")
source("Coding/R/diagnostics.R"); source("Coding/sims/design.R")

args <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
smoke <- "--smoke" %in% args

CFG <- list(datasets = c("concrete", "energy", "wine_red", "california"),
            splits = 20, test_frac = 0.2, noise = c(0, 50, 200),
            noise_california = 0, mcmc = DESIGN$mcmc, level = DESIGN$level,
            rf_trees = DESIGN$rf_trees, xgb_grid = DESIGN$xgb_grid,
            workers = DESIGN$workers)
RAW <- "Application/results/raw"
if (smoke) {
  CFG <- modifyList(CFG, list(splits = 1, noise = c(0, 50)))
  CFG$mcmc <- modifyList(CFG$mcmc, list(nskip = 100, ndpost = 100, ntree = 50,
                                        nchain_dbarts = 2))
  CFG$rf_trees <- 100
  CFG$xgb_grid <- CFG$xgb_grid[1, ]
  CFG$workers <- 1
  RAW <- file.path(tempdir(), "real_smoke")
}
## --workers overrides the machine setting in design.R (the GitHub Actions
## workflow passes the runner's core count).
if (!is.na(getopt("--workers", NA))) CFG$workers <- as.integer(getopt("--workers"))
dir.create(RAW, showWarnings = FALSE, recursive = TRUE)

load_data <- function(name) {
  d <- read.csv(file.path("Application/data/clean", paste0(name, ".csv")))
  list(x = as.matrix(d[setdiff(names(d), "y")]), y = d$y,
       censored = if (name == "california")
         read.csv("Application/data/clean/california_censored.csv")$row else integer(0))
}

## k irrelevant columns, each a row-permuted copy of a random real predictor.
add_noise <- function(x, k) {
  src <- sample.int(ncol(x), k, replace = TRUE)
  z <- vapply(src, function(j) x[sample.int(nrow(x)), j], numeric(nrow(x)))
  z <- matrix(z, nrow = nrow(x))
  colnames(z) <- sprintf("noise%03d_%s", seq_len(k), colnames(x)[src])
  z
}

## ---- gather ---------------------------------------------------------------
if ("--gather" %in% args) {
  files <- list.files(RAW, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("no results yet in ", RAW)
  d <- do.call(rbind, lapply(files, readRDS))
  write.csv(d, "Application/results/real_data.csv", row.names = FALSE)
  scores <- c("rmse", "coverage", "width", "int_score", "crps", "coverage_uncens",
              "seconds", "ess_sigma", "noise_importance", "dart_noise_sel",
              "dart_real_sel")
  agg <- do.call(rbind, lapply(scores, function(v) {
    if (all(is.na(d[[v]]))) return(NULL)
    a <- aggregate(d[[v]], by = list(dataset = d$dataset, noise = d$noise,
                                     method = d$method),
                   FUN = function(x) c(mean(x, na.rm = TRUE), mcse_mean(x),
                                       sum(is.finite(x))))
    data.frame(score = v, dataset = a$dataset, noise = a$noise, method = a$method,
               mean = a$x[, 1], se = a$x[, 2], splits = a$x[, 3])
  }))
  write.csv(agg, "Application/results/real_summary.csv", row.names = FALSE)
  wide <- reshape(agg[agg$score %in% c("rmse", "coverage", "width", "crps"),
                      c("score", "dataset", "noise", "method", "mean")],
                  idvar = c("dataset", "noise", "method"), timevar = "score",
                  direction = "wide")
  names(wide) <- sub("^mean\\.", "", names(wide))
  wide <- wide[order(wide$dataset, wide$noise, wide$method), ]
  num <- vapply(wide, is.numeric, logical(1))
  wide[num] <- lapply(wide[num], round, 3)
  print(wide, row.names = FALSE)
  quit(save = "no")
}

## ---- timing: seconds per sweep with the full settings ------------------------
## A short chain of each sampler on each dataset, to budget the full run.
if ("--timing" %in% args) {
  tools::psnice(value = 19)
  sweeps <- 200
  rows <- list()
  for (nm in CFG$datasets) for (k in unique(c(0, max(CFG$noise)))) {
    if (nm == "california" && k > max(CFG$noise_california)) next
    D <- load_data(nm); set.seed(1)
    x <- D$x
    if (k > 0) x <- cbind(x, add_noise(x, k))
    i <- sample.int(nrow(x), floor(0.8 * nrow(x)))
    b <- fit_bart(x[i, ], D$y[i], x[-i, ], ndpost = sweeps / 2, nskip = sweeps / 2,
                  ntree = CFG$mcmc$ntree)
    db <- fit_dbarts(x[i, ], D$y[i], x[-i, ], ndpost = sweeps / 2, nskip = sweeps / 2,
                     ntree = CFG$mcmc$ntree, nchain = 1)
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, noise = k, n_train = length(i),
      p = ncol(x), bart_ms = 1000 * b$time / sweeps, dbarts_ms = 1000 * db$time / sweeps)
    print(tail(rows, 1)[[1]], row.names = FALSE)
  }
  tm <- do.call(rbind, rows)
  write.csv(tm, "Application/results/real_timing.csv", row.names = FALSE)
  quit(save = "no")
}

## ---- job list -------------------------------------------------------------
jobs <- do.call(rbind, lapply(CFG$datasets, function(nm) {
  nz <- if (nm == "california") CFG$noise_california else CFG$noise
  expand.grid(dataset = nm, noise = nz, split = seq_len(CFG$splits),
              stringsAsFactors = FALSE)
}))
opt_data <- getopt("--data", NA)
opt_splits <- getopt("--splits", NA)
if (!is.na(opt_data)) jobs <- jobs[jobs$dataset %in% strsplit(opt_data, ",")[[1]], ]
if (!is.na(opt_splits)) jobs <- jobs[jobs$split %in% eval(parse(text = opt_splits)), ]
opt_noise <- getopt("--noise", NA)
if (!is.na(opt_noise)) jobs <- jobs[jobs$noise %in% as.integer(strsplit(opt_noise, ",")[[1]]), ]
jobs$file <- file.path(RAW, sprintf("%s_%03d_%02d.rds", jobs$dataset, jobs$noise, jobs$split))
done <- file.exists(jobs$file)
cat(sprintf("%d jobs requested, %d already done, %d to run\n",
            nrow(jobs), sum(done), sum(!done)))
jobs <- jobs[!done, , drop = FALSE]
if ("--dry" %in% args || nrow(jobs) == 0) {
  if (nrow(jobs)) print(table(jobs$dataset, jobs$noise))
  quit(save = "no")
}

## ---- one split ----------------------------------------------------------------
run_one <- function(j, cfg) {
  D <- load_data(j$dataset)
  ds <- match(j$dataset, c("concrete", "energy", "wine_red", "california"))
  ## The split depends on the dataset and split number only, so it is the same
  ## at every noise level; the noise columns get their own stream.
  seed <- as.integer(5e6 + 1e4 * ds + j$split)
  set.seed(seed)
  n <- nrow(D$x)
  te <- sort(sample.int(n, round(cfg$test_frac * n)))
  tr <- setdiff(seq_len(n), te)
  x <- D$x
  if (j$noise > 0) {
    set.seed(seed + 1e3 * j$noise)
    x <- cbind(x, add_noise(x, j$noise))
  }
  is_noise <- grepl("^noise", colnames(x))
  uncens <- !(te %in% D$censored)
  y_te <- D$y[te]
  mc <- cfg$mcmc
  rows <- list()
  add <- function(method, pt, lo, hi, draws, seconds, extra = list()) {
    `%||%` <- function(a, b) if (is.null(a)) b else a
    rows[[length(rows) + 1]] <<- data.frame(
      dataset = j$dataset, noise = j$noise, split = j$split, seed = seed,
      n_train = length(tr), n_test = length(te), p = ncol(x), method = method,
      seconds = seconds, rmse = rmse(pt, y_te),
      coverage = coverage(lo, hi, y_te), width = mean_width(lo, hi),
      int_score = interval_score(lo, hi, y_te, cfg$level),
      crps = if (is.null(draws)) NA_real_ else crps_draws(draws, y_te),
      coverage_uncens = if (all(uncens)) NA_real_
                        else coverage(lo[uncens], hi[uncens], y_te[uncens]),
      ess_sigma = extra$ess %||% NA_real_, geweke_z = extra$geweke %||% NA_real_,
      chain_gap = extra$chain_gap %||% NA_real_,
      noise_importance = extra$noise_imp %||% NA_real_,
      dart_noise_sel = extra$noise_sel %||% NA_real_,
      dart_real_sel = extra$real_sel %||% NA_real_,
      stringsAsFactors = FALSE)
  }
  noise_share <- function(imp) {
    if (is.null(imp) || !any(is_noise)) return(NA_real_)
    imp <- pmax(imp, 0)
    if (sum(imp) == 0) NA_real_ else sum(imp[is_noise]) / sum(imp)
  }
  bayes <- function(method, fit) {
    pi <- draw_intervals(fit$pred_draws, cfg$level)
    dg <- sigma_diagnostics(fit$extra$sigma, fit$time)
    ex <- list(ess = dg$ess, geweke = dg$geweke_z, noise_imp = noise_share(fit$importance))
    if (method == "dart" && any(is_noise))
      ex <- c(ex, list(noise_sel = sum(is_noise[fit$selected]),
                       real_sel = sum(!is_noise[fit$selected])))
    if (method == "dbarts")
      ex$chain_gap <- chain_disagreement(fit$f_draws, mc$nchain_dbarts)
    add(method, colMeans(fit$pred_draws), pi$lower, pi$upper, fit$pred_draws,
        fit$time, ex)
  }

  bayes("bart", fit_bart(x[tr, ], D$y[tr], x[te, ], ndpost = mc$ndpost, nskip = mc$nskip,
                         ntree = mc$ntree, keepevery = mc$keepevery, seed = seed + 1))
  bayes("dart", fit_bart(x[tr, ], D$y[tr], x[te, ], ndpost = mc$ndpost, nskip = mc$nskip,
                         ntree = mc$ntree, keepevery = mc$keepevery, sparse = TRUE,
                         seed = seed + 2))
  bayes("dbarts", fit_dbarts(x[tr, ], D$y[tr], x[te, ], ndpost = mc$ndpost,
                             nskip = mc$nskip, ntree = mc$ntree,
                             nchain = mc$nchain_dbarts, keepevery = mc$keepevery,
                             seed = seed + 3))
  rf <- fit_rf(x[tr, ], D$y[tr], x[te, ], num.trees = cfg$rf_trees, level = cfg$level,
               seed = seed + 4)
  add("rf", rf$point, rf$lower, rf$upper, NULL, rf$time,
      list(noise_imp = noise_share(rf$importance)))
  xg <- fit_xgb(x[tr, ], D$y[tr], x[te, ], level = cfg$level, grid = cfg$xgb_grid,
                nrounds_max = 500, seed = seed + 5)
  add("xgb", xg$point, xg$lower, xg$upper, NULL, xg$time,
      list(noise_imp = noise_share(xg$importance)))

  out <- do.call(rbind, rows)
  saveRDS(out, j$file)
  sprintf("%s noise=%d split %d done", j$dataset, j$noise, j$split)
}

## Largest first, so no worker is left with California at the end.
size <- c(concrete = 1, energy = 1, wine_red = 1.5, california = 20)
jobs <- jobs[order(-size[jobs$dataset] * (8 + jobs$noise), jobs$split), , drop = FALSE]
tools::psnice(value = 19)                       # keep the laptop usable
t0 <- proc.time()[["elapsed"]]
safe <- function(j) tryCatch(run_one(j, CFG), error = function(e)
  paste("FAILED:", j$file, conditionMessage(e)))
job_list <- split(jobs, seq_len(nrow(jobs)))
if (CFG$workers == 1) {
  res <- lapply(job_list, safe)
} else {
  cl <- parallel::makeCluster(CFG$workers, outfile = "Application/results/real_data.log")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  ROOT <- getwd()
  parallel::clusterExport(cl, c("ROOT", "CFG", "load_data", "add_noise", "run_one", "safe"))
  parallel::clusterEvalQ(cl, {
    setwd(ROOT); tools::psnice(value = 19)
    source("Coding/R/methods.R"); source("Coding/R/metrics.R")
    source("Coding/R/diagnostics.R")
  })
  res <- parallel::parLapplyLB(cl, job_list, safe, chunk.size = 1)
}
cat(sprintf("%d jobs in %.2f hours\n", length(res), (proc.time()[["elapsed"]] - t0) / 3600))
bad <- grep("^FAILED", unlist(res), value = TRUE)
if (length(bad)) cat(bad, sep = "\n") else cat("no failures\n")

if (smoke) {
  d <- do.call(rbind, lapply(list.files(RAW, full.names = TRUE), readRDS))
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], round, 3)
  print(d[, c("dataset", "noise", "p", "method", "seconds", "rmse", "coverage",
              "width", "crps", "coverage_uncens", "noise_importance",
              "dart_noise_sel", "dart_real_sel")], row.names = FALSE)
}
