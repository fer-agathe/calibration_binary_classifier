# Goodness-of-fit----

#' Mean Squared Error
mse_function <- function(pred, obs) mean((pred - obs)^2)

#' Accuracy with threshold 0.5
accuracy_function <- function(pred, obs, threshold = 0.5) {
  mean((as.numeric(as.character(pred)) > 0.5) == obs)
}

#' AUC (no threshold)
auc_function <- function(pred, obs){
  auc(obs, pred)
}

# Simulations----

#' Get calibration/test samples from a dataset
#'
#' @param seed seed to use to generate the data
#' @param n_obs number of desired observations
get_samples <- function(seed,
                        data) {
  set.seed(seed)
  n_obs <- nrow(data)

  # Calibration/test sets----
  calib_index <- sample(1:nrow(data), size = .5 * nrow(data), replace = FALSE)
  tb_calib <- data |> slice(calib_index)
  tb_test <- data |> slice(-calib_index)

  list(
    data = data,
    tb_calib = tb_calib,
    tb_test = tb_test,
    calib_index = calib_index,
    seed = seed,
    n_obs = n_obs
  )
}


#' Apply Random Forest algorithm (regression task). Predicted scores on each
#' sample are returned.
#'
#' @param train_data train dataset
#' @param calib_data calibration dataset
#' @param test_data test dataset
#'
#' @return predicted scores on train, calibration, and test sets
apply_rf <- function(train_data,
                     calib_data,
                     test_data,
                     mtry = max(floor(ncol(train_data)/3), 1),
                     nodesize = 1,
                     ntree = 500,
                     splitrule = "gini") {

  rf <- randomForest(
    d ~ .,
    data = train_data,
    mtry = mtry,
    nodesize = nodesize,
    ntree = ntree,
    splitrule =  splitrule
  )

  scores_train <- predict(rf, newdata = train_data, type = "response")
  scores_calib <- predict(rf, newdata = calib_data, type = "response")
  scores_test <- predict(rf, newdata = test_data, type = "response")

  list(
    scores_train = scores_train,
    scores_calib = scores_calib,
    scores_test = scores_test
  )
}

#' Apply Random Forest algorithm (classification task)
#'
#' @param train_data train dataset
#' @param calib_data calibration dataset
#' @param test_data test dataset
#'
#' @return predicted scores on train, calibration, and test sets
apply_rf_vote <- function(train_data,
                          calib_data,
                          test_data,
                          mtry = 1,
                          nodesize = 0.1 * nrow(train_data),
                          ntree = 500,
                          splitrule = "gini") {
  rf <- randomForest(
    d ~ .,
    data = train_data |> mutate(d = factor(d)),
    mtry = mtry,
    nodesize = nodesize,
    ntree = ntree,
    splitrule =  splitrule,
  )

  scores_train <- predict(rf, newdata = train_data, type = "vote")[, "1"]
  scores_calib <- predict(rf, newdata = calib_data, type = "vote")[, "1"]
  scores_test <- predict(rf, newdata = test_data, type = "vote")[, "1"]
  list(
    scores_train = scores_train,
    scores_calib = scores_calib,
    scores_test = scores_test
  )
}

## Wrapper function

simul_calib_rf <- function(seed,
                           type = c("regression", "classification"),
                           tuning = TRUE){

  set.seed(seed)
  linspace <- seq(0, 1, length.out = 101)

  # 1. Get and Split the dataset----

  data <- get_samples(seed = seed, data = tb_rest)
  tb_train <- tb_train
  tb_calib <- data$tb_calib
  tb_test <- data$tb_test

  # 2. Train RF----

  if (type == "regression") {
    # RF regressor parameters tuned
    if (tuning == TRUE) {
      # Use hyperparameters selected in `rf-grid-search.R`
      nodesize <- best_params_rf_reg$nodesize
      mtry <- best_params_rf_reg$mtry
      ntree <- best_params_rf_reg$num_trees
      splitrule <- best_params_rf_reg$splitrule
    } else {
      # Use default values from `randomForest()`
      nodesize <- 5
      # nodesize <- 0.1 * nrow(tb_train)
      mtry <- max(floor(ncol(tb_train)/3), 1)
      ntree <- 500
      splitrule <- "gini"
    }

    # Regression
    scores_reg <- apply_rf(
      train_data = tb_train,
      calib_data = tb_calib,
      test_data = tb_test,
      mtry = mtry,
      nodesize = nodesize,
      ntree = ntree,
      splitrule = splitrule
    )
    scores <- scores_reg
  } else {
    if (tuning == TRUE) {
      # Use hyperparameters selected in `rf-grid-search.R`
      nodesize <- best_params_rf_classif$nodesize
      mtry <- best_params_rf_classif$mtry
      ntree <- best_params_rf_classif$num_trees
      splitrule <- best_params_rf_classif$splitrule
    } else {
      # Use default values from `randomForest()`
      nodesize <- 0.1 * nrow(tb_train)
      mtry <- max(floor(ncol(tb_train)/3), 1)
      ntree <- 500
      splitrule <- "gini"
    }
    # Classification
    scores_classif <- apply_rf_vote(
      train_data = tb_train,
      calib_data = tb_calib,
      test_data = tb_test,
      mtry = mtry,
      nodesize = nodesize,
      ntree = ntree,
      splitrule = splitrule
    )
    scores <- scores_classif
  }

  # 3. Calibration----

  ## 3.1. Metrics----

  # For regression
  calibration_train <- compute_metrics(
    obs = tb_train$d,
    scores = scores$scores_train,
    k = 10
  ) |>
    mutate(sample = "train")

  calibration_calib <- compute_metrics(
    obs = tb_calib$d,
    scores = scores$scores_calib,
    k = 5
  ) |>
    mutate(sample = "calibration")

  calibration_test <- compute_metrics(
    obs = tb_test$d,
    scores = scores$scores_test,
    k = 5
  ) |>
    mutate(sample = "test")

  calibration <-
    calibration_train |>
    bind_rows(calibration_calib) |>
    bind_rows(calibration_test) |>
    mutate(model = type)

  # Comparison with standard metrics
  gof_train <- compute_gof(
    obs = tb_train$d,
    pred = scores$scores_train
  ) |> mutate(sample = "train")

  gof_calib <- compute_gof(
    obs = tb_calib$d,
    pred = scores$scores_calib
  ) |> mutate(sample = "calibration")

  gof_test <- compute_gof(
    obs = tb_test$d,
    pred = scores$scores_test
  ) |> mutate(sample = "test")

  gof <-
    gof_train |>
    bind_rows(gof_calib) |>
    bind_rows(gof_test) |>
    mutate(model = type)

  summary_metrics <-
    gof |>
    left_join(calibration, by = c("mse", "sample", "model")) |>
    relocate(sample, model, .before = "mse") |>
    mutate(seed = seed)

  # 5. Calibration curves with locfit ----

  scores_train <- scores$scores_train
  scores_calib <- scores$scores_calib
  scores_test <- scores$scores_test
  # Add a little noise to the scores, to avoir crashing R
  scores_train <- scores_train + rnorm(length(scores_train), 0, .001)
  scores_calib <- scores_calib + rnorm(length(scores_calib), 0, .001)
  scores_test <- scores_test + rnorm(length(scores_test), 0, .001)

  locfit_0_train <- locfit(
    formula = d ~ lp(score, nn = 0.15, deg = 0),
    kern = "rect", maxk = 200,
    data = tibble(d = tb_train$d, score = scores_train)
  )
  locfit_0_calib <- locfit(
    formula = d ~ lp(score, nn = 0.15, deg = 0),
    kern = "rect", maxk = 200,
    data = tibble(d = tb_calib$d, score = scores_calib)
  )
  locfit_0_test <- locfit(
    formula = d ~ lp(score, nn = 0.15, deg = 0),
    kern = "rect", maxk = 200,
    data = tibble(d = tb_test$d, score = scores_test)
  )

  score_locfit_0_train <- predict(locfit_0_train, newdata = linspace)
  score_locfit_0_calib <- predict(locfit_0_calib, newdata = linspace)
  score_locfit_0_test <- predict(locfit_0_test, newdata = linspace)
  # Make sure to have values in [0,1]
  score_locfit_0_train[score_locfit_0_train > 1] <- 1
  score_locfit_0_train[score_locfit_0_train < 0] <- 0

  score_locfit_0_calib[score_locfit_0_calib > 1] <- 1
  score_locfit_0_calib[score_locfit_0_calib < 0] <- 0

  score_locfit_0_test[score_locfit_0_test > 1] <- 1
  score_locfit_0_test[score_locfit_0_test < 0] <- 0

  res_train <- tibble(
    xlim = linspace,
    locfit_pred = score_locfit_0_train,
    sample = "train"
  )
  res_calib <- tibble(
    xlim = linspace,
    locfit_pred = score_locfit_0_calib,
    sample = "calibration"
  )
  res_test <- tibble(
    xlim = linspace,
    locfit_pred = score_locfit_0_test,
    sample = "test"
  )

  calibration_curves <-
    res_train |>
    bind_rows(
      res_calib
    ) |>
    bind_rows(
      res_test
    ) |>
    mutate(
      seed = seed,
      type = type,
      method = "No Calibration"
    )

  # 6. Recalibration----

  methods <- c("platt", "isotonic", "beta", "locfit", "locfit", "locfit")
  params <- list(
    NULL, NULL, NULL,
    list(nn = .15, deg = 0), list(nn = .15, deg = 1), list(nn = .15, deg = 2)
  )
  method_names <- c(
    "platt", "isotonic", "beta", "locfit_0", "locfit_1", "locfit_2"
  )

  scores_c <- vector(mode = "list", length = length(methods))
  names(scores_c) <- method_names
  res_c <- tibble()

  for (i_method in 1:length(methods)) {
    method <- methods[i_method]
    params_current <- params[[i_method]]

    ## 6.1 Recalibrated Scores----
    scores_c[[i_method]] <- recalibrate(
      obs_train = tb_train$d,
      pred_train = scores$scores_train,
      obs_calib = tb_calib$d,
      pred_calib = scores$scores_calib,
      obs_test = tb_test$d,
      pred_test = scores$scores_test,
      method = method,
      params = params_current
    )
    ## 6.2 Recalibration Curves----
    scores_c_train <- scores_c[[i_method]]$tb_score_c_train$p_c
    scores_c_calib <- scores_c[[i_method]]$tb_score_c_calib$p_c
    scores_c_test <- scores_c[[i_method]]$tb_score_c_test$p_c
    # Add a little noise to the scores, to avoir crashing R
    scores_c_train <- scores_c_train + rnorm(length(scores_c_train), 0, .001)
    scores_c_calib <- scores_c_calib + rnorm(length(scores_c_calib), 0, .001)
    scores_c_test <- scores_c_test + rnorm(length(scores_c_test), 0, .001)

    locfit_0_train <- locfit(
      formula = d ~ lp(score, nn = 0.15, deg = 0),
      kern = "rect", maxk = 200,
      data = tibble(d = tb_train$d, score = scores_c_train)
    )
    locfit_0_calib <- locfit(
      formula = d ~ lp(score, nn = 0.15, deg = 0),
      kern = "rect", maxk = 200,
      data = tibble(d = tb_calib$d, score = scores_c_calib)
    )
    locfit_0_test <- locfit(
      formula = d ~ lp(score, nn = 0.15, deg = 0),
      kern = "rect", maxk = 200,
      data = tibble(d = tb_test$d, score = scores_c_test)
    )

    score_c_locfit_0_train <- predict(locfit_0_train, newdata = linspace)
    score_c_locfit_0_calib <- predict(locfit_0_calib, newdata = linspace)
    score_c_locfit_0_test <- predict(locfit_0_test, newdata = linspace)

    # Make sure to have values in [0,1]
    score_c_locfit_0_train[score_c_locfit_0_train > 1] <- 1
    score_c_locfit_0_train[score_c_locfit_0_train < 0] <- 0

    score_c_locfit_0_calib[score_c_locfit_0_calib > 1] <- 1
    score_c_locfit_0_calib[score_c_locfit_0_calib < 0] <- 0

    score_c_locfit_0_test[score_c_locfit_0_test > 1] <- 1
    score_c_locfit_0_test[score_c_locfit_0_test < 0] <- 0

    res_c_train <- tibble(
      xlim = linspace,
      locfit_pred = score_c_locfit_0_train,
      sample = "train",
      method = method_names[i_method]
    )
    res_c_calib <- tibble(
      xlim = linspace,
      locfit_pred = score_c_locfit_0_calib,
      sample = "calibration",
      method = method_names[i_method]
    )
    res_c_test <- tibble(
      xlim = linspace,
      locfit_pred = score_c_locfit_0_test,
      sample = "test",
      method = method_names[i_method]
    )

    res_c <- res_c  |>
      bind_rows(
        res_c_train) |>
      bind_rows(
        res_c_calib
      ) |>
      bind_rows(
        res_c_test
      )
  }

  res_c <- res_c |>
    mutate(
      seed = seed,
      type = type
    )

  recalibration_curves <- calibration_curves |> bind_rows(res_c)

  # 7. Recalibration Metrics----

  ## 7.1. Standard Metrics----

  gof_c_train <- map(
    .x = scores_c,
    .f = ~compute_gof(
      obs = tb_train$d,
      pred = .x$tb_score_c_train$p_c
    )
  ) |>
    list_rbind(names_to = "method") |>
    mutate(sample = "train")

  gof_c_calib <- map(
    .x = scores_c,
    .f = ~compute_gof(
      obs = tb_calib$d,
      pred = .x$tb_score_c_calib$p_c
    )
  ) |>
    list_rbind(names_to = "method") |>
    mutate(sample = "calibration")

  gof_c_test <- map(
    .x = scores_c,
    .f = ~compute_gof(
      obs = tb_test$d,
      pred = .x$tb_score_c_test$p_c
    )
  ) |>
    list_rbind(names_to = "method") |>
    mutate(sample = "test")

  gof_c <-
    gof_c_train |>
    bind_rows(gof_c_calib) |>
    bind_rows(gof_c_test) |>
    mutate(model = type)

  ## 7.2. Calibration Metrics----

  calibration_c_train <- map(
    .x = scores_c,
    .f = ~compute_metrics(
      obs = tb_train$d,
      scores = .x$tb_score_c_train$p_c,
      k = 10
    )
  ) |>
    list_rbind(names_to = "method") |>
    mutate(sample = "train")

  calibration_c_calib <- map(
    .x = scores_c,
    .f = ~compute_metrics(
      obs = tb_calib$d,
      scores = .x$tb_score_c_calib$p_c,
      k = 5
    )
  ) |>
    list_rbind(names_to = "method") |>
    mutate(sample = "calibration")

  calibration_c_test <- map(
    .x = scores_c,
    .f = ~compute_metrics(
      obs = tb_test$d,
      scores = .x$tb_score_c_test$p_c,
      k = 5
    )
  ) |>
    list_rbind(names_to = "method") |>
    mutate(sample = "test")


  calibration_c <-
    calibration_c_train |>
    bind_rows(calibration_c_calib) |>
    bind_rows(calibration_c_test) |>
    mutate(model = type)

  summary_metrics_c <-
    gof_c |>
    left_join(calibration_c, by = c("mse", "sample", "model", "method")) |>
    relocate(sample, model, method, .before = "mse") |>
    mutate(seed = seed)

  summary_metrics <- summary_metrics |> mutate(method = "No Calibration")

  list(
    recalibration_curves = recalibration_curves,
    summary_metrics = list_rbind(list(summary_metrics, summary_metrics_c)),
    seed = seed,
    scores = scores,
    scores_c = scores_c,
    type = type,
    tuning = tuning
  )
}

# Visualizations----

## Boxplots Metrics----

#' Boxplot of performance metrics (Accuracy, Sensitivity, Specificity)
#'
boxplot_metrics_gof <- function(data_metrics){

  gof_metrics <- c("accuracy", "sensitivity", "specificity")
  gof_metrics_lab <- c("Accuracy", "Sensitivity", "Specificity")

  data_metrics_plot <-
    data_metrics |>
    select(-mse) |>
    select(sample, model, method, seed, !!!syms(gof_metrics)) |>
    pivot_longer(
      cols = !!gof_metrics,
      names_to = "metric", values_to = "value"
    )

  mat <- matrix(
    c(1:(6), rep(7,3)), ncol = 3, byrow = TRUE
  )
  layout(mat, heights = c(rep(3, 4),1))

  range_val <- data_metrics_plot |>
    group_by(metric) |>
    summarise(
      min_val = min(value),
      max_val = max(value)
    )

  # Goodness of Fit----
  for (i_model in 1:length(levels(data_metrics_plot$model))) {
    model <- levels(data_metrics_plot$model)[i_model]
    data_metrics_gof <- data_metrics_plot |>
      filter(model == !!model)

    for (i_metric in 1:length(gof_metrics)) {

      metric <- gof_metrics[i_metric]
      metric_lab <- gof_metrics_lab[i_metric]
      val_min <- range_val |> filter(metric == !!metric) |> pull(min_val)
      val_max <- range_val |> filter(metric == !!metric) |> pull(max_val)
      title <- ""
      if (i_metric == 2) {
        title <- str_c(model, "\n", metric_lab)
      } else {
        title <- str_c("\n", metric_lab)
      }

      nb_methods <- data_metrics_gof$method |> levels() |> length()

      data_metrics_gof_current <-
        data_metrics_gof |>
        filter(metric == !!metric)

      par(mar = c(2.5, 1, 4.1, 1))
      boxplot(
        value ~ sample + method,
        data = data_metrics_gof_current,
        col = colours,
        horizontal = TRUE,
        las = 1, xlab = "", ylab = "",
        main = title,
        border = c("black", adjustcolor("black", alpha.f = .5)),
        yaxt = "n",
        ylim = c(val_min, val_max)
      )
      # Horizontal lines
      for (i in seq(3, nb_methods * 2, by = 2) - .5) {
        abline(h = i, lty = 1, col = "gray")
      }
    }

  }

  # Legend----
  par(mar = c(0, 4.3, 0, 4.3))
  plot.new()
  legend(
    "center",
    legend = rev(levels(data_metrics_plot$method)),
    fill = rev(colours_calib),
    # lwd = 2,
    xpd = TRUE, ncol = 4
  )
}

#' Boxplots of calibration metrics (Brier Score, ECE, LSC)
#'
boxplot_metrics_calib <- function(data_metrics){

  calib_metrics <- c("brier", "ece", "lcs")
  calib_metrics_lab <- c("Brier Score", "ECE", "LCS")

  data_metrics_plot <-
    data_metrics |>
    select(-mse) |>
    select(sample, model, method, seed, !!!syms(calib_metrics)) |>
    pivot_longer(
      cols = !!calib_metrics,
      names_to = "metric", values_to = "value"
    )

  range_val <- data_metrics_plot |>
    group_by(metric) |>
    summarise(
      min_val = min(value),
      max_val = max(value)
    )

  mat <- matrix(
    c(1:(6), rep(7,3)), ncol = 3, byrow = TRUE
  )
  layout(mat, heights = c(rep(3, 4),1))

  # Calibration Metrics----
  for (i_model in 1:length(levels(data_metrics_plot$model))) {
    model <- levels(data_metrics_plot$model)[i_model]
    data_metrics_calib <- data_metrics_plot |>
      filter(model == !!model)

    for (i_metric in 1:length(calib_metrics)) {
      metric <- calib_metrics[i_metric]
      metric_lab <- calib_metrics_lab[i_metric]
      val_min <- range_val |> filter(metric == !!metric) |> pull(min_val)
      val_max <- range_val |> filter(metric == !!metric) |> pull(max_val)
      if (i_metric == 2) {
        title <- str_c(model, "\n", metric_lab)
      } else {
        title <- str_c("\n", metric_lab)
      }

      nb_methods <- data_metrics_calib$method |> levels() |> length()

      data_metrics_calib_current <-
        data_metrics_calib |>
        filter(metric == !!metric)

      par(mar = c(2.5, 1, 4.1, 1))
      boxplot(
        value ~ sample + method,
        data = data_metrics_calib_current,
        col = colours,
        horizontal = TRUE,
        las = 1, xlab = "", ylab = "",
        main = title,
        border = c("black", adjustcolor("black", alpha.f = .5)),
        yaxt = "n",
        ylim = c(val_min, val_max)
      )
      # Horizontal lines
      for (i in seq(3, nb_methods * 2, by = 2) - .5) {
        abline(h = i, lty = 1, col = "gray")
      }
    }

  }

  # Legend----
  par(mar = c(0, 4.3, 0, 4.3))
  plot.new()
  legend(
    "center",
    legend = rev(levels(data_metrics_plot$method)),
    fill = rev(colours_calib),
    # lwd = 2,
    xpd = TRUE, ncol = 4
  )
}


## Comparison Calibration / Performances----

plot_compare <- function(gof_metric, calib_metric) {

  gof_metrics <- c("sensitivity", "specificity", "AUC", "accuracy")
  calib_metrics <- c("brier", "LCS", "ece", "qmse", "wmse")
  types <- c("regression", "classification")
  types_labs <- c("Regression", "Classification")

  mat <- c(1,2,3,3) |> matrix(ncol=2, byrow = TRUE)
  layout(mat, heights = c(3,.5))

  for (i_type in 1:length(types)) {
    type <- types[i_type]
    type_lab <- types_labs[i_type]

    if (type == "regression"){
      data <- compare_reg
    } else if(type == "classification") {
      data <- compare_class
    } else {
      stop("wrong type")
    }

    train_gof <- paste(gof_metric, "train", sep = "_")
    rest_gof <- paste(gof_metric, "rest", sep = "_")
    train_calib <- paste(calib_metric, "train", sep = "_")
    rest_calib <- paste(calib_metric, "rest", sep = "_")

    # Plot boundaries
    xmin <- min(data[[train_gof]], data[[rest_gof]])
    xmax <- max(data[[train_gof]], data[[rest_gof]])
    ymin <- min(data[[train_calib]], data[[rest_calib]])
    ymax <- max(data[[train_calib]], data[[rest_calib]])


    par(mar = c(4.1, 4.1, 2.1, 1))

    # From train to test
    plot(
      xmin, xmax,
      xlim = c(xmin, xmax), ylim = c(ymin, ymax),
      xlab = gof_metric, ylab = calib_metric,
      main = type_lab,
    )
    segments(
      x0 = data[[train_gof]],
      y0 = data[[train_calib]],
      x1 = data[[rest_gof]],
      y1 = data[[rest_calib]],
      col = adjustcolor("gray", alpha.f = .3)
    )

    # Train set
    points(
      x = data[[train_gof]],
      y = data[[train_calib]],
      pch = 19, cex = .8,
      col = adjustcolor(colours_samples[["Train"]], alpha.f = .3)
    )
    # Loess
    lw_train <- loess(data[[train_calib]] ~ data[[train_gof]])
    ord <- order(data[[train_gof]])
    lines(data[[train_gof]][ord],lw_train$fitted[ord], col = colours_samples[["Train"]],lwd = 3, lty = 3)
    # Test set
    points(
      x = data[[rest_gof]],
      y = data[[rest_calib]],
      col = adjustcolor(colours_samples[["Test"]], alpha.f = .3),
      pch = 19, cex = .8
    )
    # Loess
    lw_rest <- loess(data[[rest_calib]] ~ data[[rest_gof]])
    ord <- order(data[[rest_gof]])
    lines(data[[rest_gof]][ord],lw_rest$fitted[ord], col = colours_samples[["Test"]],lwd = 3, lty = 3)

  }
  # Legend
  par(mar = c(0, 4.3, 0, 4.3))
  plot.new()
  legend(
    "center",
    legend = c("Train", "Test"),
    col = colours_samples,
    pch = 19,
    xpd = TRUE, ncol = 2
  )

}

# Calibration Curves----

#' For a vector of scores, compute the number of obs in each bin
#' The bins are defined on evenly separated values in [0,1]
#'
count_scores <- function(scores) {
  breaks <- seq(0, 1, by = .05)

  if (!is.null(scores)) {
    n_bins <- table(cut(scores, breaks = breaks))
  } else {
    n_bins <- NA_integer_
  }
  tibble(
    bins = names(table(cut(breaks, breaks = breaks))),
    n_bins = as.vector(n_bins)
  )
}

#' Applies the count_scores() function to a list of simulations
#'
count_scores_simul <- function(scores_simul) {
  map(scores_simul, count_scores) |>
    list_rbind() |>
    group_by(bins) |>
    summarise(
      n_bins = mean(n_bins)
    )
}

#' Extract the recalibrated scores from a list of simulations and computes
#' the number of scores in bins (see count_scores())
#' For each simulation, there are multiple correction techniques
#'
count_scores_simul_method <- function(scores_simul_methods) {
  map(scores_simul_methods, "scores_c") |>
    map(
      .f = function(simul){
        map(
          simul,
          .f = function(sim_method) {
            count_scores_simul(sim_method$tb_score_c_train) |>
              mutate(sample = "train") |>
              bind_rows(
                count_scores_simul(sim_method$tb_score_c_calib) |>
                  mutate(sample = "calibration")
              ) |>
              bind_rows(
                count_scores_simul(sim_method$tb_score_c_test) |>
                  mutate(sample = "test")
              )
          }
        ) |>
          list_rbind(names_to = "method")
      },
      .progress = TRUE
    ) |>
    list_rbind() |>
    group_by(method, sample, bins) |>
    summarise(
      n_bins = mean(n_bins),
      .groups = "drop"
    )
}

plot_calib_locfit_simuls <- function(calib_curve, type) {
  tb_calib_curve <- calib_curves |>filter(type == !!type)
  tb_n_bins <- n_bins |> filter(type == !!type)

  nb_methods <- tb_calib_curve$method |> unique() |> length()

  mat_init <- c(1, 1, 2, 4, 3, 5) |> matrix(ncol = 2, byrow = TRUE)
  mat <- mat_init
  for (j in 2:(ceiling(nb_methods/2))) {
    mat <- rbind(mat, mat_init + (5*(j-1)))
  }
  mat <- cbind(mat, mat + max(mat))

  layout(mat, height = rep(c(.35, 1, 3), nb_methods))

  y_lim <- c(0, 1)

  samples <- c("Calibration", "Test")

  for (i_method in 1:length(methods)) {
    method <- methods_labs[i_method]
    # Title for the models
    par(mar = c(0, 0.5, 0, 0.5))
    plot(
      0:1, 0:1,
      col = NULL,
      type="n",
      xlim = c(0,2), ylim = 0:1,
      xlab = "",
      ylab = "",
      main = "",
      xaxt = "n",
      yaxt = "n",
      new = TRUE,
      frame.plot = FALSE
    )
    # Get the center coordinates of the plot
    plot_center <- c(mean(par("usr")[1:2]), mean(par("usr")[3:4]))
    # Get the width and height of the text
    text_dim <- strwidth(method, cex = 1.2, font = 2)
    # Calculate the coordinates to center the text
    text_x <- plot_center[1] - text_dim / 2
    text_y <- plot_center[2]

    # Add text to the plot at the centered coordinates
    text(
      1, text_y,
      method, col = colours_calib[i_method],
      cex = 1, font = 2
    )
    for(i_sample in 1:length(samples)) {
      sample <- samples[i_sample]
      tb_plot <- tb_calib_curve |>
        filter(
          sample == !!sample,
          method == !!method,
          type == !!type
        )
      n_bins_current <- tb_n_bins |>
        filter(
          sample == !!sample,
          method == !!method,
          type == !!type
        )
      n_bins_no_calib_current <- tb_n_bins |>
        filter(
          sample == !!sample,
          method == "No Calibration",
          type == !!type
        )
      # Barplots on top
      par(mar = c(0.5, 4.3, 1.0, 0.5))
      y_lim_bp <- range(
        c(n_bins_no_calib_current$n_bins,
          n_bins_current$n_bins, na.rm = TRUE)
      )
      barplot(
        n_bins_no_calib_current$n_bins,
        col = adjustcolor("gray", alpha.f = .3),
        ylim = y_lim_bp,
        border = "white",
        axes = FALSE,
        xlab = "", ylab = "", main = ""
      )
      barplot(
        n_bins_current$n_bins,
        col = adjustcolor("#0072B2", alpha.f = .3),
        ylim = y_lim_bp,
        border = "white",
        axes = FALSE,
        xlab = "", ylab = "", main = "",
        add = TRUE
      )

      # Calib curve
      par(mar = c(4.1, 4.3, 0.5, 0.5), mgp = c(2, 1, 0))
      plot(
        tb_plot$xlim, tb_plot$mean,
        type = "l", col = colours_samples[sample],
        xlim = 0:1, ylim = 0:1,
        xlab = latex2exp::TeX("$p^u$"),
        ylab = latex2exp::TeX("$E(D | p^u = p^c)$"),
        main = ""
      )
      polygon(
        c(tb_plot$xlim, rev(tb_plot$xlim)),
        c(tb_plot$lower, rev(tb_plot$upper)),
        col = adjustcolor(col = colours_samples[sample], alpha.f = .4),
        border = NA
      )
      segments(0, 0, 1, 1, col = "black", lty = 2)
    }
  }

}
