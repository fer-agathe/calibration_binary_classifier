# Measuring and visualizing recalibration on synthetic data

library(tidyverse)
library(locfit)
library(binom)
library(future)
library(betacal)

# Data grenerative process
source("functions/synthetic-data.R")

# Simulation + Calib + Recalib functions (specific to the simulations)
source("functions/simulations-synthetic.R")

# Standard performance metrics and calibration metrics
source("functions/metrics.R")

# 1. Simulations----

alphas <- gammas <- c(1/3, 1, 3)
n_repl <- 10 # number of replications (200 to get the results of the article)
n_obs <- 2000 # number of observations to draw
grid_alpha <- expand_grid(alpha = alphas, seed = 1:n_repl)
grid_gamma <- expand_grid(gamma = gammas, seed = 1:n_repl)

## Alpha----

nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_alpha))
  simul_recalib_alpha <- furrr::future_map(
    .x = 1:nrow(grid_alpha),
    .f = ~{
      p()
      f_simul_recalib(
        i = .x,
        grid = grid_alpha,
        n_obs = n_obs,
        type = "alpha",
        linspace = NULL)
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

## Gamma----
library(future)
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_gamma))
  simul_recalib_gamma <- furrr::future_map(
    .x = 1:nrow(grid_gamma),
    .f = ~{
      p()
      f_simul_recalib(
        i = .x,
        grid = grid_gamma,
        n_obs = n_obs,
        type = "gamma",
        linspace = NULL)
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

# 2. Goodness-of-fit----

## Alpha----
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_alpha))
  recalib_metrics_alpha <- furrr::future_map(
    .x = 1:nrow(grid_alpha),
    .f = ~{
      p()
      compute_gof_simul_recalib(
        i = .x,
        grid = grid_alpha,
        n_obs = n_obs,
        type = "alpha"
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

recalib_metrics_alpha <- list_rbind(recalib_metrics_alpha)

## Gamma----
library(future)
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_gamma))
  recalib_metrics_gamma <- furrr::future_map(
    .x = 1:nrow(grid_gamma),
    .f = ~{
      p()
      compute_gof_simul(
        i = .x,
        grid = grid_gamma,
        n_obs = n_obs,
        type = "gamma"
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

recalib_metrics_gamma <- list_rbind(recalib_metrics_gamma)

# c("mse", "accuracy", "sensitivity", "specificity", "roc", "auc")
metrics <- c("mse", "accuracy", "auc")
methods <- c("platt", "isotonic", "beta", "locfit_0", "locfit_1", "locfit_2")

# For alpha
for (method in methods) {
  current_recalib_metrics <- recalib_metrics_alpha |>
    filter(method == method)
  boxplot_simuls_metrics(
    tb_metrics = current_recalib_metrics,
    type = "alpha", metrics = metrics
  ) |> print()
}

# For gamma
for (method in methods) {
  current_recalib_metrics <- recalib_metrics_gamma |>
    filter(method == method)
  boxplot_simuls_metrics(
    tb_metrics = current_recalib_metrics,
    type = "gamma", metrics = metrics
  ) |> print()
}

# 3. Calibration Curves----

## 3.1 Quantile-based----

grid_count_alpha <-
  expand_grid(
    i = 1:nrow(grid_alpha),
    method = c(
      "True Prob.", "No Calibration",
      "platt", "isotonic", "beta", "locfit_0", "locfit_1", "locfit_2")
  )
grid_count_gamma <-
  expand_grid(
    i = 1:nrow(grid_gamma),
    method = c(
      "True Prob.", "No Calibration",
      "platt", "isotonic", "beta", "locfit_0", "locfit_1", "locfit_2")
  )

# Number of observation in each bin separating the [0,1] segment
# with uncalibrated and recalibrated scores
# (both on the calibration and the recalibration sets),
# for all the simulations.
count_scores_alpha <- map(
  .x = 1:nrow(grid_count_alpha),
  .f = ~get_count_simul_recalib(
    i = grid_count_alpha$i[.x],
    type = "alpha",
    method = grid_count_alpha$method[.x]
  ),
  .progress = TRUE
)
count_scores_gamma <- map(
  .x = 1:nrow(grid_count_gamma),
  .f = ~get_count_simul_recalib(
    i = grid_count_gamma$i[.x],
    type = "gamma",
    method = grid_count_gamma$method[.x]
  ),
  .progress = TRUE
)

# Average count per bin over the simulations.
count_scores_alpha <-
  count_scores_alpha |>
  list_rbind() |>
  group_by(method, type, bins, transform_scale) |>
  summarise(
    n_bins_calib = mean(n_bins_calib, na.rm = TRUE),
    n_bins_test = mean(n_bins_test, na.rm = TRUE),
    n_bins_c_calib = mean(n_bins_c_calib, na.rm = TRUE),
    n_bins_c_test = mean(n_bins_c_test, na.rm = TRUE),
    .groups = "drop"
  )

count_scores_gamma <-
  count_scores_gamma |>
  list_rbind() |>
  group_by(method,type, bins, transform_scale) |>
  summarise(
    n_bins_calib = mean(n_bins_calib, na.rm = TRUE),
    n_bins_test = mean(n_bins_test, na.rm = TRUE),
    n_bins_c_calib = mean(n_bins_c_calib, na.rm = TRUE),
    n_bins_c_test = mean(n_bins_c_test, na.rm = TRUE),
    .groups = "drop"
  )

### Alpha----
# For alpha
methods <- names(simul_recalib_alpha[[1]]$res_recalibration)
methods <- c("True Prob.", "No Calibration", methods)
summary_bins_simuls_alpha <- vector(mode = "list", length = length(methods))
names(summary_bins_simuls_alpha) <- methods
library(future)
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)

for (i_method in 1:length(methods)) {
  progressr::with_progress({
    p <- progressr::progressor(steps = nrow(grid_alpha))
    summary_bins_simuls_m <- furrr::future_map(
      .x = simul_recalib_alpha,
      .f = ~{
        p()
        get_summary_bins_simul(simul = .x, method = methods[i_method], k = 10)
      },
      .options = furrr::furrr_options(seed = NULL)
    )
  })
  # Add value for alpha
  for (j in 1:length(summary_bins_simuls_m)) {
    summary_bins_simuls_m[[j]]$scale_parameter <- grid_alpha$alpha[j]
  }
  summary_bins_simuls_alpha[[i_method]] <- summary_bins_simuls_m |>
    list_rbind(names_to = "i_row")
}
summary_bins_simuls_alpha <- list_rbind(
  summary_bins_simuls_alpha, names_to = "method"
) |>
  mutate(type = "alpha")

### Gamma----
# For gamma
methods <- names(simul_recalib_gamma[[1]]$res_recalibration)
methods <- c("True Prob.", "No Calibration", methods)
summary_bins_simuls_gamma <- vector(mode = "list", length = length(methods))
names(summary_bins_simuls_gamma) <- methods
library(future)
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)

for (i_method in 1:length(methods)) {
  progressr::with_progress({
    p <- progressr::progressor(steps = nrow(grid_gamma))
    summary_bins_simuls_m <- furrr::future_map(
      .x = simul_recalib_gamma,
      .f = ~{
        p()
        get_summary_bins_simul(simul = .x, method = methods[i_method], k = 10)
      },
      .options = furrr::furrr_options(seed = NULL)
    )
  })
  # Add value for alpha
  for (j in 1:length(summary_bins_simuls_m)) {
    summary_bins_simuls_m[[j]]$scale_parameter <- grid_gamma$gamma[j]
  }
  summary_bins_simuls_gamma[[i_method]] <- summary_bins_simuls_m |>
    list_rbind(names_to = "i_row")
}
summary_bins_simuls_gamma <- list_rbind(
  summary_bins_simuls_gamma, names_to = "method"
) |>
  mutate(type = "gamma")

### Merge----

summary_bins_simuls <- summary_bins_simuls_alpha |>
  bind_rows(summary_bins_simuls_gamma)

### Plots----

for (method in methods) {
  plot_calibration_quant_simul_recalib(method = method, type = "alpha") |>
    print()
}

for (method in methods) {
  plot_calibration_quant_simul_recalib(method = method, type = "gamma") |>
    print()
}


methods <- c(
  "True Prob.",
  "No Calibration",
  "platt", "isotonic", "beta",
  "locfit_0", "locfit_1", "locfit_2"
)
methods_labs <- c(
  "True Prob.",
  "No Calibration", "Platt", "Isotonic", "Beta",
  "Locfit (deg=0)", "Locfit (deg=1)", "Locfit (deg=2)"
)

for (method in methods) {
  plot_calibration_quant_simul_recalib_2(
    calib_curve = summary_bins_simuls,
    method = method
  ) |>
    print()
}

## 3.2 Local Regression----

### Alpha----
# For alpha
methods <- names(simul_recalib_alpha[[1]]$res_recalibration)
methods <- c("True Prob.", "No Calibration", methods)
calib_curve_locfit_simuls_alpha <- vector(mode = "list", length = length(methods))
names(calib_curve_locfit_simuls_alpha) <- methods

nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)

for (i_method in 1:length(methods)) {
  progressr::with_progress({
    p <- progressr::progressor(steps = nrow(grid_alpha))
    calib_curve_locfit_simuls_m <- furrr::future_map(
      .x = simul_recalib_alpha,
      .f = ~{
        p()
        calibration_curve_locfit_simul_recalib(
          simul = .x, method = methods[i_method], k = 10
        )
      },
      .options = furrr::furrr_options(seed = NULL)
    )
  })
  # Add value for alpha
  for (j in 1:length(calib_curve_locfit_simuls_m)) {
    calib_curve_locfit_simuls_m[[j]]$scale_parameter <- grid_alpha$alpha[j]
  }
  calib_curve_locfit_simuls_alpha[[i_method]] <- calib_curve_locfit_simuls_m |>
    list_rbind(names_to = "i_row")
}

calib_curve_locfit_simuls_alpha <- list_rbind(
  calib_curve_locfit_simuls_alpha, names_to = "method"
) |>
  mutate(type = "alpha")

### Gamma-----
# For gamma
methods <- names(simul_recalib_gamma[[1]]$res_recalibration)
# Remove isotonic which makes the R session crash...
methods <- c("True Prob.", "No Calibration", methods)
calib_curve_locfit_simuls_gamma <- vector(mode = "list", length = length(methods))
names(calib_curve_locfit_simuls_gamma) <- methods

nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)

for (i_method in 1:length(methods)) {
  progressr::with_progress({
    p <- progressr::progressor(steps = nrow(grid_gamma))
    calib_curve_locfit_simuls_m <- furrr::future_map(
      .x = simul_recalib_gamma,
      .f = ~{
        p()
        calibration_curve_locfit_simul_recalib(
          simul = .x, method = methods[i_method], k = 10
        )
      },
      .options = furrr::furrr_options(seed = NULL)
    )
  })
  # Add value for alpha
  for (j in 1:length(calib_curve_locfit_simuls_m)) {
    calib_curve_locfit_simuls_m[[j]]$scale_parameter <- grid_gamma$gamma[j]
  }
  calib_curve_locfit_simuls_gamma[[i_method]] <- calib_curve_locfit_simuls_m |>
    list_rbind(names_to = "i_row")
}
calib_curve_locfit_simuls_gamma <- list_rbind(
  calib_curve_locfit_simuls_gamma, names_to = "method"
) |>
  mutate(type = "gamma")

### Merge----

calib_curve_locfit_simuls <-
  calib_curve_locfit_simuls_alpha |>
  bind_rows(calib_curve_locfit_simuls_gamma) |>
  mutate(
    sample = factor(
      sample,
      levels = c("calibration", "test"), labels = c("Calibration", "Test")
    )
  )

### Plots----

for (method in methods) {
  plot_calibration_locfit_simuls_recalib(
    calib_curve = calib_curve_locfit_simuls,
    method = method, type = "alpha"
  ) |>
    print()
}

for (method in methods) {
  plot_calibration_locfit_simuls_recalib(
    calib_curve = calib_curve_locfit_simuls,
    method = method, type = "gamma"
  ) |>
    print()
}

# With both alpha / gamma on the same figure
for (method in methods) {
  plot_calib_locfit_simuls_recalib_2(
    calib_curve = calib_curve_locfit_simuls,
    method = method
  )
}


# 4. Boxplots of Metrics----

calib_metrics_simul_alpha <- map(simul_recalib_alpha, "calib_metrics") |>
  list_rbind()
calib_metrics_simul_gamma <- map(simul_recalib_gamma, "calib_metrics") |>
  list_rbind()

# In a single tibble
calib_metrics_simul <- calib_metrics_simul_alpha |>
  bind_rows(calib_metrics_simul_gamma) |>
  pivot_longer(
    cols = c(mse, brier, ece, qmse, wmse, lcs),
    names_to = "metric", values_to = "value"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("mse", "brier", "ece", "qmse", "wmse", "lcs"),
      labels = c("True MSE", "Brier Score", "ECE", "QMSE", "WMSE", "LCS")
    ),
    method = factor(
      method,
      levels = c(
        "True Prob.", "No Calibration",
        "platt", "isotonic", "beta", "locfit_0", "locfit_1", "locfit_2"),
      labels = c(
        "True Prob.", "No Calibration",
        "Platt Scaling", "Isotonic Reg.", "Beta Calib.",
        "Local Reg. (deg = 0)", "Local Reg. (deg = 1)", "Local Reg (deg = 2)"
      )
    ),
    sample = factor(
      sample,
      levels = c("Calibration", "Test")
    )
  )

metrics <- c("True MSE", "Brier Score", "ECE", "QMSE", "WMSE", "LCS")

for (metric in metrics) {
  plot_boxplot_metric_recalib(
    metric = metric,
    calib_metrics_simul = calib_metrics_simul,
    type = "alpha"
  ) |> print()
}

for (metric in metrics) {
  plot_boxplot_metric_recalib(
    metric = metric,
    calib_metrics_simul = calib_metrics_simul,
    type = "gamma"
  ) |> print()
}


# With another type of viz, as in the article

# Load metrics computed on uncalibrated scores (see `sim-calibration.R`)
load("output/standard_metrics_alpha.rda")
load("output/standard_metrics_gamma.rda")

# Standard Metrics
standard_metrics <- recalib_metrics_alpha |>
  bind_rows(recalib_metrics_gamma) |>
  # Without recalibration
  bind_rows(
    metrics_alpha |>
      mutate(method = "No Calibration")
  ) |>
  bind_rows(
    metrics_gamma |>
      mutate(method = "No Calibration")
  ) |>
  filter(threshold == .5) |>
  rename(transform_scale = scale_parameter) |>
  select(
    sample, seed, transform_scale, type, method,
    mse, accuracy, sensitivity, specificity, auc
  ) |>
  pivot_longer(
    cols = c(mse, accuracy, sensitivity, specificity, auc),
    names_to = "metric",
    values_to = "value"
  )

# Calibration metrics
calib_metrics_simul_alpha <- map(simul_recalib_alpha, "calib_metrics") |>
  list_rbind()
calib_metrics_simul_gamma <- map(simul_recalib_gamma, "calib_metrics") |>
  list_rbind()

calib_metrics_simul <- calib_metrics_simul_alpha |>
  bind_rows(calib_metrics_simul_gamma) |>
  pivot_longer(
    cols = c(mse, brier, ece, qmse, wmse, lcs),
    names_to = "metric", values_to = "value"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("mse", "brier", "ece", "qmse", "wmse", "lcs"),
      labels = c("True MSE", "Brier Score", "ECE", "QMSE", "WMSE", "LCS")
    ),
    sample = case_when(
      sample == "Calibration"~"calibration",
      sample == "Test"~"test",
      TRUE~NA_character_
    )
  )

metrics_all <-
  calib_metrics_simul |>
  bind_rows(standard_metrics) |>
  mutate(
    method = factor(
      method,
      levels = c("True Prob.", "No Calibration",
                 "platt", "isotonic", "beta", "locfit_0", "locfit_1", "locfit_2"),
      labels = c(
        "True Prob.",
        "No Calibration", "Platt", "Isotonic", "Beta",
        "Locfit (deg=0)", "Locfit (deg=1)", "Locfit (deg=2)")
    ),
    sample = factor(
      sample, levels = c("calibration", "test"), labels = c("Calibration", "Test")
    )
  )


metrics <- c("True MSE", "Brier Score", "ECE", "QMSE", "WMSE", "LCS")

for (metric in metrics) {
  boxplot_std_metrics_calib(
    tb_calib_metrics = metrics_all,
    metric = metric
  )
}



