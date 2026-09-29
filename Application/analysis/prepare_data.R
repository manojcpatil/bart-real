## prepare_data.R -- T10. Obtain and clean the four real datasets.
##
## Every dataset is downloaded from its original public source if the raw file is
## not already in data/raw, checked (dimensions, missing values, constant columns,
## duplicate rows, ranges), given short column names, and written to data/clean
## as a CSV whose last column is the response y. A registry, data/clean/datasets.csv,
## records the source, checksum and every cleaning decision, so Chapter 5 can cite
## exactly what was done.
##
## Cleaning decisions (each is reversible by editing the list DECISIONS below):
##  * Exact duplicate rows are removed. The analysis draws random 80/20 splits;
##    a row with an identical twin can land in both halves, and the model is then
##    scored on a point it has already seen, which flatters both accuracy and
##    coverage. Concrete has 25 such rows and red wine 240.
##  * Energy: the response is heating load (Y1). Cooling load (Y2) is dropped; it
##    is a second response, not a predictor.
##  * California Housing: y is median house value in units of $100,000. Values
##    were top-coded by the Census at $500,001; those rows are kept (they are part
##    of the data a user would face) but flagged in data/clean/california_censored.csv
##    so the analysis can report coverage with and without them.
##  * No transformation of y or of the predictors: BART is invariant to monotone
##    transformations of the predictors, and y is analysed on its natural scale
##    because calibration on that scale is the question.
##
## Usage (from the BART_Research folder):
##   Rscript Application/analysis/prepare_data.R
## Out: Application/data/clean/{concrete,energy,wine_red,california}.csv,
##      Application/data/clean/california_censored.csv, datasets.csv

if (!dir.exists("Application/data")) stop("run this from the BART_Research folder")
RAW <- "Application/data/raw"; CLEAN <- "Application/data/clean"
dir.create(RAW, showWarnings = FALSE, recursive = TRUE)
dir.create(CLEAN, showWarnings = FALSE, recursive = TRUE)

DECISIONS <- list(drop_duplicates = TRUE)

SOURCES <- list(
  concrete = list(file = "Concrete_Data.xls",
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/concrete/compressive/Concrete_Data.xls",
    cite = "Yeh (1998), UCI Machine Learning Repository"),
  energy = list(file = "ENB2012_data.xlsx",
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00242/ENB2012_data.xlsx",
    cite = "Tsanas and Xifara (2012), UCI Machine Learning Repository"),
  wine_red = list(file = "winequality-red.csv",
    url = "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv",
    cite = "Cortez et al. (2009), UCI Machine Learning Repository"),
  california = list(file = "houses.zip",
    url = "http://lib.stat.cmu.edu/datasets/houses.zip",
    cite = "Pace and Barry (1997), StatLib")
)

for (s in SOURCES) {
  f <- file.path(RAW, s$file)
  if (!file.exists(f)) download.file(s$url, f, mode = "wb", quiet = TRUE)
}

## ---- readers ----------------------------------------------------------------
read_concrete <- function() {
  d <- as.data.frame(readxl::read_excel(file.path(RAW, "Concrete_Data.xls")))
  names(d) <- c("cement", "slag", "fly_ash", "water", "superplasticizer",
                "coarse_agg", "fine_agg", "age", "y")
  d
}
read_energy <- function() {
  d <- as.data.frame(readxl::read_excel(file.path(RAW, "ENB2012_data.xlsx")))
  d <- d[rowSums(!is.na(d)) > 0, colSums(!is.na(d)) > 0]   # trailing blank cells
  names(d) <- c("rel_compactness", "surface_area", "wall_area", "roof_area",
                "height", "orientation", "glazing_area", "glazing_dist",
                "y", "cooling_load")
  d$cooling_load <- NULL
  d
}
read_wine <- function() {
  d <- read.csv(file.path(RAW, "winequality-red.csv"), sep = ";", check.names = FALSE)
  names(d) <- c("fixed_acidity", "volatile_acidity", "citric_acid", "residual_sugar",
                "chlorides", "free_so2", "total_so2", "density", "ph", "sulphates",
                "alcohol", "y")
  d
}
read_california <- function() {
  txt <- readLines(unz(file.path(RAW, "houses.zip"), "cadata.txt", encoding = "latin1"),
                   warn = FALSE)
  txt <- iconv(txt, "latin1", "UTF-8")          # the notes use Windows quote marks
  txt <- txt[vapply(strsplit(trimws(txt), "\\s+"), length, 0L) == 9]  # skip the notes
  d <- read.table(text = txt)
  names(d) <- c("value", "median_income", "housing_age", "total_rooms",
                "total_bedrooms", "population", "households", "latitude", "longitude")
  d$y <- d$value / 1e5
  d$value <- NULL
  d
}

## ---- clean and check ----------------------------------------------------------
clean_one <- function(name, d) {
  n_raw <- nrow(d)
  stopifnot(names(d)[ncol(d)] == "y", all(vapply(d, is.numeric, TRUE)))
  n_na <- sum(!complete.cases(d))
  d <- d[complete.cases(d), ]
  const <- names(d)[vapply(d, function(v) length(unique(v)) == 1, TRUE)]
  if (length(const)) d <- d[setdiff(names(d), const)]
  dup <- duplicated(d)
  n_dup <- sum(dup)
  if (DECISIONS$drop_duplicates) d <- d[!dup, ]
  ## Rows whose predictors repeat but whose response differs are kept: they are
  ## genuine replicate measurements and carry information about the noise.
  n_xrep <- sum(duplicated(d[setdiff(names(d), "y")]))
  rownames(d) <- NULL
  write.csv(d, file.path(CLEAN, paste0(name, ".csv")), row.names = FALSE)
  list(d = d, row = data.frame(
    dataset = name, n_raw = n_raw, n = nrow(d), p = ncol(d) - 1,
    rows_missing = n_na, duplicate_rows = n_dup,
    duplicates_removed = DECISIONS$drop_duplicates,
    repeated_x_different_y = n_xrep, constant_cols = length(const),
    y_min = min(d$y), y_median = median(d$y), y_max = max(d$y), y_sd = sd(d$y),
    y_distinct = length(unique(d$y)),
    source = SOURCES[[name]]$url, cite = SOURCES[[name]]$cite,
    raw_md5 = unname(tools::md5sum(file.path(RAW, SOURCES[[name]]$file))),
    prepared = format(Sys.Date()), stringsAsFactors = FALSE))
}

out <- list(concrete = clean_one("concrete", read_concrete()),
            energy = clean_one("energy", read_energy()),
            wine_red = clean_one("wine_red", read_wine()),
            california = clean_one("california", read_california()))

## Top-coded California rows, matched on the cleaned file's row order.
ca <- out$california$d
write.csv(data.frame(row = which(ca$y >= 5.00001 - 1e-9)),
          file.path(CLEAN, "california_censored.csv"), row.names = FALSE)

reg <- do.call(rbind, lapply(out, `[[`, "row"))
reg$note <- c(
  "response: 28-day-equivalent compressive strength (MPa); age is in days and takes 14 values",
  "response: heating load; 12 building shapes x 4 orientations x 16 glazing settings; cooling load dropped",
  "response: sensory quality score, integer 3-8, so intervals are for a discrete y",
  sprintf("response: median house value ($100k); %d rows top-coded at 5.00001 kept and flagged; runtime stress test",
          sum(ca$y >= 5.00001 - 1e-9)))
write.csv(reg, file.path(CLEAN, "datasets.csv"), row.names = FALSE)
rownames(reg) <- NULL
print(reg[, c("dataset", "n_raw", "n", "p", "rows_missing", "duplicate_rows",
              "repeated_x_different_y", "constant_cols", "y_min", "y_median",
              "y_max", "y_distinct")])
