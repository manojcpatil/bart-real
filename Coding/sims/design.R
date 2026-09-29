## design.R -- the agreed simulation design. Every later script sources this so
## that the grid, the replication count and the MCMC settings are defined once.
##
## Decision of 2026-09-23 (budget option B in Coding/results/budget.csv):
##   keep the full 2000 kept draws and 200 replications, and pay for that by
##   dropping the single most expensive cell, n = 5000 with p = 500, which alone
##   accounted for about 40 per cent of the projected 3.3 days. That cell is run
##   once, as a single-replication check reported in the discussion, rather than
##   200 times.
## Sampler settings, decided 2026-09-23 by three experiments in this folder:
##   thinning_check.R  thinning does not raise the effective sample size, so the
##                     poor mixing seen in the timing pilot is drift, not
##                     autocorrelation. keepevery stays at 1.
##   burnin_check.R    at equal cost, sweeps spent on burn-in give a converged
##                     chain (Geweke z near 0, ESS above 220 of 500, sigma within
##                     one per cent of the truth) where thinning gives neither.
##   trace_check.R     20000-sweep chains show sigma settling by roughly 4000 to
##                     6000 sweeps at n = 5000, and still drifting downwards at
##                     20000 sweeps when n = 1000 with p >= 100.
## Hence a long burn-in and no thinning. The replication count is halved from the
## 200 of option B to pay for it, which leaves the total cost close to the option
## B budget; the pilot confirms that 100 replications still meet the Monte Carlo
## error target.

DESIGN <- list(

  ## ---- the main grid (Study A) -------------------------------------------
  grid = subset(expand.grid(n = c(200, 1000, 5000), p = c(10, 100, 500)),
                !(n == 5000 & p == 500)),
  dropped_cell = data.frame(n = 5000, p = 500),   # run once, not replicated
  reps = 100,                                     # provisional, confirmed by the pilot
  sigma = 1,
  mean_fn = "friedman",
  errors = c("gaussian", "hetero", "t3"),         # gaussian for A, all three for B
  n_test = 1000,
  level = 0.95,

  ## ---- MCMC ---------------------------------------------------------------
  mcmc = list(nskip = 5000, ndpost = 2000, ntree = 200,
              keepevery = 1,                      # thinning does not help: see above
              nchain_dbarts = 4),

  ## Convergence is itself an experimental factor, because the evidence above
  ## says that a chain which has not settled produces wider intervals and so
  ## looks better calibrated than it is. Study C repeats a middle cell at each
  ## of these burn-in lengths to show the effect directly.
  burnin_ladder = c(500, 2000, 5000, 20000),

  ## Study C must separate two explanations of the pilot result, in which the
  ## BART package covered 0.89 and dbarts covered 0.94 on identical data. The
  ## packages differ in the sampler (dbarts also proposes CHANGE moves) but the
  ## pilot also ran dbarts with two chains, and pooling draws across chains that
  ## have not mixed widens the intervals on its own, which is the same mechanism
  ## that made thinning look helpful. Crossing the two factors settles it.
  study_c = list(cell = data.frame(n = 1000, p = 100),
                 sampler = c("bart", "dbarts"),
                 nchain = c(1, 4),
                 reps = 30),

  ## ---- competitors --------------------------------------------------------
  rf_trees = 500,
  xgb_grid = expand.grid(eta = c(0.05, 0.1, 0.3), max_depth = c(3, 6),
                         subsample = 0.8, colsample_bytree = 0.8),

  ## ---- machine ------------------------------------------------------------
  workers = 2                                     # 2 physical / 4 logical cores; 3 left
                                                  # the laptop unusable (2026-09-24)
)

#' Pilot settings: shorter chains and a subset of cells, used only to estimate
#' the between-replication spread of each score so that the replication count
#' can be justified rather than assumed.
PILOT <- list(
  grid = data.frame(n = c(200, 1000, 1000), p = c(10, 100, 500)),
  reps = 50,
  mcmc = list(nskip = 2000, ndpost = 1000, ntree = 200, keepevery = 1,
              nchain_dbarts = 2),
  xgb_grid = expand.grid(eta = 0.1, max_depth = 3,
                         subsample = 0.8, colsample_bytree = 0.8),
  target_mcse = list(coverage = 0.005, rmse_f = 0.01, crps = 0.01)
)
