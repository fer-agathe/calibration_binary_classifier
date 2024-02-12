# Measuring calibration, random forest
# Using real-world data
# Here: grid search

library(tidyverse)
library(locfit)
library(randomForest)
library(future)

source("functions/rf-functions.R")


# 1. Data----
data_credit <- read.csv2(
  "data/default_credit.csv",
  header = TRUE, skip = 1
)
colnames(data_credit)[25] <- "d"
data_credit <- data_credit |> select(-ID)
data_credit <- data_credit |> as_tibble()
data_credit <- data_credit |>
  mutate(
    across(
      c("SEX", "EDUCATION", "MARRIAGE", "PAY_0", "PAY_2", "PAY_3",
        "PAY_4", "PAY_5", "PAY_6"),
      ~as.factor(.x)
    )
  )

## SMOTE----
library(performanceEstimation)

# need of a factor variable to apply SMOTE function
set.seed(123)
data_credit$d <- as.factor(data_credit$d)
new_df <- performanceEstimation::smote(d ~ ., data_credit, perc.over = 2)

data_credit <- new_df
data_credit$d <- as.numeric(data_credit$d) # change in the labels
data_credit <- data_credit |> mutate(d = ifelse(d == 1, 1, 2))
data_credit <- data_credit |> mutate(d = d-1)

# Save for later use
save(data_credit, file = "output/data_credit_smote.rda")

# Split dataset into 2 parts: train and rest
set.seed(1234)
ind_train <- sample(1:nrow(data_credit), size = .5 * nrow(data_credit))
tb_train <- data_credit |> slice(ind_train)
tb_rest <- data_credit |> slice(-ind_train)

# Save for later use
save(tb_train, file = "output/data_credit_smote_train.rda")
save(tb_rest, file = "output/data_credit_smote_rest.rda")


# 2. Grid Search----

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
    expand_grid(
      num_trees = c(100,300, 500),
      mtry = seq(1,(ncol(tb_train)/2)),
      nodesize = c(5, 10, 15, 20)
    )
}

## 2.1 Regression----

# Note: takes a few minutes
# To test the code, you can use a subset of data

nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_params))
  mse_oob_rf_reg <- furrr::future_map(
    .x = 1:nrow(grid_params),
    .f = ~{
      # Estim random forest and get the evaluation metric
      rf <- randomForest(
        d ~ .,
        data = tb_train,
        mtry = grid_params$mtry[.x],
        nodesize = grid_params$nodesize[.x],
        ntree = grid_params$num_trees[.x]
        keep.inbag = TRUE
      )

      num_trees <- grid_params$num_trees[.x]

      # Identify out of bag observations in each tree
      out_of_bag <- map(.x = 1:nrow(tb_train), .f = ~which(rf[["inbag"]][.x,] == 0))
      rf_pred_all <- predict(rf, tb_train,
                             predict.all = TRUE,
                             type = "response")$individual
      rf_pred <- unlist(map(.x = 1:nrow(tb_train), .f = ~mean(rf_pred_all[.x,out_of_bag[[.x]]])))

      oob_err <- mse_function(pred = rf_pred, obs = tb_train |> pull(d))
      mse_oob <- oob_err

      # Progress bar
      p()

      # Return object:
      tibble(
        mtry = grid_params$mtry[.x],
        nodesize = grid_params$nodesize[.x],
        num_trees = grid_params$num_trees[.x]
        mse_oob = mse_oob
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

# Ordering the sets of hyperparameters by increasing values of oob mse
best_params_rf_reg <-
  mse_oob_rf_reg |> list_rbind() |>
  arrange(mse_oob)

# Save for later use
save(best_params_rf_reg, file = "output/best_params_rf_reg.rda")

## 2.2 Classification----

nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_params))
  mse_oob_rf_classif <- furrr::future_map(
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

      p()

      # Return object:
      tibble(
        mtry = grid_params$mtry[.x],
        nodesize = grid_params$nodesize[.x],
        num_trees = grid_params$num_trees[.x],
        err_oob = err_oob
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

# Arrange by increasing values of oob error rate
best_params_rf_classif <-
  mse_oob_rf_classif |> list_rbind() |>
  arrange(err_oob)

save(best_params_rf_classif, file = "output/best_params_rf_classif.rda")
