# Measuring calibration, random forest
# Using real-world data
# Here: measuring calibration of both regression and classification RF
# Note: results presented in `rf-simul-results.R`

library(tidyverse)
library(locfit)
library(randomForest)
library(future)
library(binom)


# The functions used to compute calibration for given scores are defined in the
# calibrationBinary package.

# Some additional functions, specific to the simulation exercise:
source("functions/rf-functions.R")


# 1. Load Results----

# See `rf-simuls.R`
load("output/metrics_rf_tuned_reg.rda")
load("output/metrics_rf_tuned_class.rda")


# 2. Extract Metrics----

summary_metrics_tuned_reg <-
  map(metrics_rf_tuned_reg, "summary_metrics") |>
  list_rbind()
summary_metrics_tuned_class <- map(metrics_rf_tuned_class, "summary_metrics") |>
  list_rbind()

# Bind in a single tibble
metrics_rf <- summary_metrics_tuned_reg |>
  mutate(model = "Regression") |>
  bind_rows(
    summary_metrics_tuned_class |>
      mutate(model = "Classification")
  ) |>
  filter(sample %in% c("calibration", "test")) |>
  mutate(
    AUC = 1-AUC,
    model = factor(
      model,
      levels = c("Regression", "Classification")
    ),
    method = factor(
      method,
      levels = c("No Calibration", "platt", "isotonic", "beta",
                 "locfit_0", "locfit_1", "locfit_2") |> rev(),
      labels = c("No Calibration", "Platt Scaling", "Isotonic", "Beta",
                 "Locfit (deg=0)", "Locfit (deg=1)", "Locfit (deg=2)") |> rev()
    ),
    sample = factor(
      sample,
      levels = c("calibration", "test") |> rev(),
      labels = c("Calibration", "Test") |> rev())
  )

# 3. Visualize Metrics----

colours <- c("Train" = "#0072B2", "Calibration" = "#D55E00", "Test" = "#009E73")
library(tidyverse)
colours_samples <- c("Train" = "#0072B2", "Calibration" = "#D55E00", "Test" = "#009E73")
colours_calib <- c(
  # "#332288",
  "#117733",
  "#44AA99", "#88CCEE",
  "#DDCC77", "#CC6677", "#AA4499", "#882255") |> rev()

colours_test <- adjustcolor(colours_calib, alpha.f = .5)

colours <- NULL
for (k in 1:length(colours_calib))
  colours <- c(colours, colours_calib[k], colours_test[k])

boxplot_metrics_gof(metrics_rf)
boxplot_metrics_calib(metrics_rf)

# 4. Calibration Curves----

curves_reg <-
  map(metrics_rf_tuned_reg, "recalibration_curves") |>
  list_rbind()
curves_class <-
  map(metrics_rf_tuned_class, "recalibration_curves") |>
  list_rbind()

calib_curves <-
  curves_reg |>
  bind_rows(curves_class) |>
  group_by(xlim, sample, type, method) |>
  summarise(
    mean = mean(locfit_pred),
    lower = quantile(locfit_pred, probs = .025),
    upper = quantile(locfit_pred, probs = 0.975),
    .groups = "drop"
  )

methods <- c(
  "No Calibration",
  "platt", "isotonic", "beta",
  "locfit_0", "locfit_1", "locfit_2"
)
methods_labs <- c(
  "No Calibration", "Platt", "Isotonic", "Beta",
  "Locfit (deg=0)", "Locfit (deg=1)", "Locfit (deg=2)"
)

calib_curves <-
  calib_curves |>
  filter(sample %in% c("calibration", "test")) |>
  mutate(
    sample = factor(
      sample,
      levels = c("calibration", "test"),
      labels = c("Calibration", "Test")
    ),
    type = factor(
      type,
      levels = c("regression", "classification"),
      labels = c("Regression", "Classification")
    ),
    method = factor(
      method,
      levels = methods,
      labels = methods_labs
    )
  )

## Extract Uncalibrated Scores----

### Uncalibrated scores----

#### Train set----
scores_no_calib_reg_train <-
  map(metrics_rf_tuned_reg, "scores") |>
  map(~.x$"scores_train")
scores_no_calib_class_train <-
  map(metrics_rf_tuned_class, "scores") |>
  map(~.x$"scores_train")

#### Calibration set----
scores_no_calib_reg_calib <-
  map(metrics_rf_tuned_reg, "scores") |>
  map(~.x$"scores_calib")
scores_no_calib_class_calib <-
  map(metrics_rf_tuned_class, "scores") |>
  map(~.x$"scores_calib")

#### Test set----
scores_no_calib_reg_test <-
  map(metrics_rf_tuned_reg, "scores") |>
  map(~.x$"scores_test")
scores_no_calib_class_test <-
  map(metrics_rf_tuned_class, "scores") |>
  map(~.x$"scores_test")

### Recalibrated Scores----
scores_no_calib_reg_train <-
  map(metrics_rf_tuned_reg, "scores") |>
  map(~.x$"scores_train")


## Regression----
### No Calibration----
n_bins_no_calib_reg <-
  count_scores_simul(scores_no_calib_reg_train) |>
  mutate(sample = "train", method = "No Calibration", type = "Regression") |>
  bind_rows(
    count_scores_simul(scores_no_calib_reg_calib) |>
      mutate(sample = "calibration", method = "No Calibration", type = "Regression")
  ) |>
  bind_rows(
    count_scores_simul(scores_no_calib_reg_test) |>
      mutate(sample = "test", method = "No Calibration", type = "Regression")
  )
### With Recalibration methods----
n_bins_recalib_reg <- count_scores_simul_method(
  scores_simul_methods = metrics_rf_tuned_reg) |>
  mutate(type = "Regression")

## Extract Recalibrated Scores----

### Regression----
#### No Calibration----
n_bins_no_calib_class <-
  count_scores_simul(scores_no_calib_class_train) |>
  mutate(sample = "train", method = "No Calibration", type = "Classification") |>
  bind_rows(
    count_scores_simul(scores_no_calib_class_calib) |>
      mutate(sample = "calibration", method = "No Calibration", type = "Classification")
  ) |>
  bind_rows(
    count_scores_simul(scores_no_calib_class_test) |>
      mutate(sample = "test", method = "No Calibration", type = "Classification")
  )
#### With Recalibration methods----
n_bins_recalib_class <- count_scores_simul_method(
  scores_simul_methods = metrics_rf_tuned_class) |>
  mutate(type = "Classification")

### Merge all simulations----
n_bins <-
  n_bins_no_calib_reg |>
  bind_rows(n_bins_recalib_reg) |>
  bind_rows(n_bins_no_calib_class) |>
  bind_rows(n_bins_recalib_class)

n_bins <-
  n_bins |>
  filter(sample %in% c("calibration", "test")) |>
  mutate(
    sample = factor(
      sample,
      levels = c("calibration", "test"),
      labels = c("Calibration", "Test")
    ),
    type = factor(
      type,
      levels = c("Regression", "Classification"),
      labels = c("Regression", "Classification")
    ),
    method = factor(
      method,
      levels = !!methods,
      labels = !!methods_labs
    )
  )

## Visualizations----

colours_samples <- c(
  "Train" = "#0072B2", "Calibration" = "#D55E00", "Test" = "#009E73"
)
colours_calib <- c(
  # "#332288",
  "#117733",
  "#44AA99", "#88CCEE",
  "#DDCC77", "#CC6677", "#AA4499", "#882255")

plot_calib_locfit_simuls(
  calib_curve = calib_curves, type = "Regression"
)

plot_calib_locfit_simuls(
  calib_curve = calib_curves, type = "Classification"
)


# 5. Calibration vs. Performance----

# Load results (see `rf-simuls.R`)
load("output/compare_cal_gof_class.rda")
load("output/compare_cal_gof_reg.rda")

compare_class <- compare_cal_gof_class |> arrange(err_oob) |>
  mutate(
    AUC_train = AUC_train,
    AUC_rest = AUC_rest
  )
compare_reg <- compare_cal_gof_reg |> arrange(mse_oob)|>
  mutate(
    AUC_train = AUC_train,
    AUC_rest = AUC_rest
  )

samples <- c("train", "test")
sample_labels <- c("Train", "Test")
colours_samples <- c("Train" = "#0072B2", "Test" = "#009E73")

plot_compare(gof_metric = "AUC", calib_metric = "LCS")
plot_compare(gof_metric = "AUC", calib_metric = "brier")
plot_compare(gof_metric = "AUC", calib_metric = "ece")
