# Measuring calibration, random forest
# Using real-world data
# Here: measuring calibration of both regression and classification RF
# Note: results presented in `rf-simuls-results.R`

library(tidyverse)
library(locfit)
library(randomForest)
library(binom)
library(future)
library(betacal)

# The functions used to compute calibration for given scores:
source("functions/recalibration-methods.R")

# Those to compute performance metrics and calibration metrics:
source("functions/metrics.R")

# Some additional functions, specific to the simulation exercise:
source("functions/rf-functions.R")


# 1. Load Data ----

# Data
# See `rf-grid-search.R`
load("output/data_credit_smote_train.rda")
load("output/data_credit_smote_rest.rda")

# Hyperparameters for each type of forest, for the grid
load("output/best_params_rf_reg.rda")
load("output/best_params_rf_classif.rda")

# 'Best' hyperparameters:

# For regression RF
best_params_rf_reg <-
  best_params_rf_reg |>
  arrange(mse_oob) |>
  slice(1)
best_params_rf_reg

# For classification RF
best_params_rf_classif <-
  best_params_rf_classif |>
  arrange(err_oob) |>
  slice(1)
best_params_rf_classif


# 2. Simulations----

# We consider 200 different splits of the data on which a recalibration
# technique will be applied.

# Number of replications
n_repl <- 5 # set to 200 to replicate the results in the paper
seed <- 1:n_repl

## Regression----

# Note: takes a few minutes on a standard computer.

nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = n_repl)
  metrics_rf_tuned_reg <- furrr::future_map(
    .x = 1:n_repl,
    .f = ~{
      resul <- simul_calib_rf(
        seed = .x,
        type = "regression",
        tuning = TRUE
      )
      p()
      resul
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

save(
  metrics_rf_tuned_reg,
  file = "output/metrics_rf_tuned_reg.rda"
)

## Classification----
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = n_repl)
  metrics_rf_tuned_class <- furrr::future_map(
    .x = 1:n_repl,
    .f = ~{
      resul <- simul_calib_rf(
        seed = .x,
        type = "classification",
        tuning = TRUE
      )
      p()
      resul
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

save(
  metrics_rf_tuned_class,
  file = "output/metrics_rf_tuned_class.rda"
)

# 3. Calibration vs Goodness-of-fit----

# Training forests for different hyperparameters
# For each set, computes calibration metrics and goodness-of-fit metrics
# The forest is trained on the train set
# The metrics are computed on both the train and the test set
# Note: no recalibration done here.


quick_example <- TRUE #set to FALSE to replicate results from the paper

if (quick_example) {
  grid_params <-
    expand_grid(
      num_trees = c(100),
      mtry = c(4, 8),
      nodesize = c(5, 10)
    )
} else {
  grid_params <-
    grid_params <-
    expand_grid(
      num_trees = c(100, 300, 500),
      mtry = seq(1, (ncol(tb_train) / 2)),
      nodesize = c(5, 10, 15, 20)
    )
}


## Regression----
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_params))
  compare_cal_gof_reg <- furrr::future_map(
    .x = 1:nrow(grid_params),
    .f = ~{
      # Estim random forest and get the evaluation metric
      rf <- randomForest(
        d ~ .,
        data = tb_train,
        mtry = grid_params$mtry[.x],
        nodesize = grid_params$nodesize[.x],
        ntree = grid_params$num_trees[.x],
        keep.inbag = TRUE
      )

      num_trees <- grid_params$num_trees[.x]

      # Identify out of bag observations in each tree
      out_of_bag <- map(
        .x = 1:nrow(tb_train),
        .f = ~which(rf[["inbag"]][.x,] == 0)
      )
      rf_pred_all <- predict(
        rf, tb_train,
        predict.all = TRUE,
        type = "response")$individual
      rf_pred <- unlist(
        map(
          .x = 1:nrow(tb_train),
          .f = ~mean(rf_pred_all[.x,out_of_bag[[.x]]])
        )
      )

      oob_err <- mse_function(pred = rf_pred, obs = tb_train |> pull(d))
      mse_oob <- oob_err

      # Predict RF on train/rest dataset
      scores_train <- predict(rf, newdata = tb_train, type = "response")
      scores_rest <- predict(rf, newdata = tb_rest, type = "response")

      # Calibration metrics (Brier Score and LCS) on train/rest dataset
      calib_metrics_train <- compute_metrics(obs = tb_train$d, scores = scores_train)
      calib_metrics_rest <- compute_metrics(obs = tb_rest$d, scores = scores_rest)

      # GOF metrics on train/rest dataset
      gof_metrics_train <- compute_gof(obs = tb_train$d, pred = scores_train)
      gof_metrics_rest <- compute_gof(obs = tb_rest$d, pred = scores_rest)
      # Update progressbar
      p()

      # Return object:
      tibble(
        mtry = grid_params$mtry[.x],
        nodesize = grid_params$nodesize[.x],
        num_trees = grid_params$num_trees[.x],
        mse_oob = mse_oob,
        brier_train = calib_metrics_train$brier,
        LCS_train = calib_metrics_train$lcs,
        ece_train = calib_metrics_train$ece,
        qmse_train = calib_metrics_train$qmse,
        wmse_train = calib_metrics_train$wmse,
        sensitivity_train = gof_metrics_train$sensitivity,
        specificity_train = gof_metrics_train$specificity,
        AUC_train = gof_metrics_train$AUC,
        accuracy_train = gof_metrics_train$accuracy,
        brier_rest = calib_metrics_rest$brier,
        LCS_rest = calib_metrics_rest$lcs,
        ece_rest = calib_metrics_rest$ece,
        qmse_rest = calib_metrics_rest$qmse,
        wmse_rest = calib_metrics_rest$wmse,
        sensitivity_rest = gof_metrics_rest$sensitivity,
        specificity_rest = gof_metrics_rest$specificity,
        AUC_rest = gof_metrics_rest$AUC,
        accuracy_rest = gof_metrics_rest$accuracy
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})
compare_cal_gof_reg <- list_rbind(compare_cal_gof_reg)

save(compare_cal_gof_reg,  file = "output/compare_cal_gof_reg.rda")

## Classification----

nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_params))
  compare_cal_gof_class <- furrr::future_map(
    .x = 1:nrow(grid_params),
    .f = ~{
      # Estim random forest and get the evaluation metric
      rf <- randomForest(
        as.factor(d) ~ .,
        data = tb_train,
        mtry = grid_params$mtry[.x],
        nodesize = grid_params$nodesize[.x],
        ntree = grid_params$num_trees[.x],
        keep.inbag = TRUE
      )

      num_trees <- grid_params$num_trees[.x]

      # Identify out of bag observations in each tree
      err_oob <- rf$err.rate[num_trees,1]

      # Predict RF on train/rest dataset
      scores_train <- predict(rf, newdata = tb_train, type = "vote")[, "1"]
      scores_rest <- predict(rf, newdata = tb_rest, type = "vote")[, "1"]

      # Calibration metrics (Brier Score and LCS) on train/rest dataset
      calib_metrics_train <- compute_metrics(obs = tb_train$d, scores = scores_train)
      calib_metrics_rest <- compute_metrics(obs = tb_rest$d, scores = scores_rest)

      # GOF metrics on train/rest dataset
      gof_metrics_train <- compute_gof(obs = tb_train$d, pred = scores_train)
      gof_metrics_rest <- compute_gof(obs = tb_rest$d, pred = scores_rest)

      # Update progressbar
      p()

      # Return object:
      tibble(
        mtry = grid_params$mtry[.x],
        nodesize = grid_params$nodesize[.x],
        num_trees = grid_params$num_trees[.x],
        err_oob = err_oob,
        brier_train = calib_metrics_train$brier,
        LCS_train = calib_metrics_train$lcs,
        ece_train = calib_metrics_train$ece,
        qmse_train = calib_metrics_train$qmse,
        wmse_train = calib_metrics_train$wmse,
        sensitivity_train = gof_metrics_train$sensitivity,
        specificity_train = gof_metrics_train$specificity,
        AUC_train = gof_metrics_train$AUC,
        accuracy_train = gof_metrics_train$accuracy,
        brier_rest = calib_metrics_rest$brier,
        LCS_rest = calib_metrics_rest$lcs,
        ece_rest = calib_metrics_rest$ece,
        qmse_rest = calib_metrics_rest$qmse,
        wmse_rest = calib_metrics_rest$wmse,
        sensitivity_rest = gof_metrics_rest$sensitivity,
        specificity_rest = gof_metrics_rest$specificity,
        AUC_rest = gof_metrics_rest$AUC,
        accuracy_rest = gof_metrics_rest$accuracy
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

compare_cal_gof_class <- list_rbind(compare_cal_gof_class)

save(compare_cal_gof_class,  file = "output/compare_cal_gof_class.rda")

# Results can be visualized in the next file: `rf-simuls-results.R`
