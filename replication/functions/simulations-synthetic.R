# Calibration Metrics----
#' Performs one replication for a simulation
#'
#' @param i row number of the grid to use for the simulation
#' @param grid grid tibble with the seed number (column `seed`) and the deformations value (either `alpha` or `gamma`)
#' @param n_obs desired number of observation
#' @param type deformation probability type (either `alpha` or `gamma`); the
#' name should match with the `grid` tibble
#' @param linspace values at which to compute the mean observed event when computing the WMSE
f_simul <- function(i,
                    grid,
                    n_obs,
                    type = c("alpha", "gamma"),
                    linspace = NULL) {

  if (is.null(linspace)) linspace <- seq(0, 1, length.out = 100)

  current_seed <- grid$seed[i]
  if (type == "alpha") {
    transform_scale <- grid$alpha[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = transform_scale, gamma = 1
    )
  } else if (type == "gamma") {
    transform_scale <- grid$gamma[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = 1, gamma = transform_scale
    )
  } else {
    stop("Transform type should be either alpha or gamma.")
  }

  # Get the calib/test datasets with true probabilities
  data_all_calib <- current_data$data_all |>
    slice(current_data$calib_index)

  data_all_test <- current_data$data_all |>
    slice(-current_data$calib_index)

  # Transformed probabilities
  p_u_calib <- data_all_calib$p_u
  p_u_test <- data_all_test$p_u
  # Observed events
  d_calib <- data_all_calib$d
  d_test <- data_all_test$d

  # Mean observed events
  expected_events_calib <- map(
    .x = linspace,
    .f = ~local_ci_scores(
      obs = data_all_calib$d,
      scores = data_all_calib$p_u,
      tau = .x,
      nn = .15, prob = .5, method = "probit")
  ) |>
    bind_rows()

  expected_events_test <- map(
    .x = linspace,
    .f = ~local_ci_scores(
      obs = data_all_test$d,
      scores = data_all_test$p_u,
      tau = .x,
      nn = .15, prob = .5, method = "probit")
  ) |>
    bind_rows()

  # Compute Metrics
  ## Calibration set
  mse_calib <- mean((data_all_calib$p - data_all_calib$p_u)^2)
  brier_calib <- brier_score(obs = d_calib, score = p_u_calib)
  ece_calib <- e_calib_error(
    obs = d_calib, scores = p_u_calib, k = 10, threshold = .5
  )
  qmse_calib <- qmse_error(
    obs = d_calib, score = p_u_calib, k = 10, threshold = .5
  )
  wmse_calib <- weighted_mse(
    local_scores = expected_events_calib, scores = p_u_calib
  )
  lcs_calib <- local_calib_score(obs = d_calib, scores = p_u_calib)

  ## Test Set
  mse_test <- mean((data_all_test$p - data_all_test$p_u)^2)
  brier_test <- brier_score(obs = d_test, score = p_u_test)
  ece_test <- e_calib_error(
    obs = d_test, scores = p_u_test, k = 10, threshold = .5
  )
  qmse_test <- qmse_error(
    obs = d_test, score = p_u_test, k = 10, threshold = .5
  )
  wmse_test <- weighted_mse(
    local_scores = expected_events_test, scores = p_u_test
  )
  lcs_test <- local_calib_score(obs = d_test, scores = p_u_test)

  tibble(
    seed = grid$seed[i],
    scale_parameter = transform_scale,
    type = type,
    sample = "calibration",
    mse = mse_test,
    brier = brier_test,
    ece = ece_test,
    qmse = qmse_test,
    wmse = wmse_test,
    lcs = lcs_test
  ) |>
    bind_rows(
      tibble(
        seed = grid$seed[i],
        scale_parameter = transform_scale,
        type = type,
        sample = "test",
        mse = mse_test,
        brier = brier_test,
        ece = ece_test,
        qmse = qmse_test,
        wmse = wmse_test,
        lcs = lcs_test
      )
    )
}

#' Recalibrates scores using a calibration
#'
#' @param obs_calib vector of observed events in the calibration set
#' @param scores_calib vector of predicted probabilities in the calibration set
#' #' @param obs_test vector of observed events in the test set
#' @param scores_test vector of predicted probabilities in the test set
#' @param method recalibration method (`"platt"` for Platt-Scaling,
#'   `"isotonic"` for isotonic regression, `"beta"` for beta calibration,
#'   `"locfit"` for local regression)
#' @param iso_params list of named parameters to use in the local regression
#'   (`nn` for fraction of nearest neighbors to use, `deg` for degree)
#' @param linspace vector of alues at which to compute the recalibrated scores
#' @returns list of three elements: recalibrated scores on the calibration set,
#'   recalibrated scores on the test set, and recalibrated scores on a segment
#'   of values
recalibrate_simul <- function(obs_calib,
                              scores_calib,
                              obs_test,
                              scores_test,
                              method = c("platt", "isotonic", "beta", "locfit"),
                              iso_params = NULL,
                              linspace = NULL) {

  if (is.null(linspace)) linspace <- seq(0, 1, length.out = 100)

  data_calib <- tibble(d = obs_calib, p_u = scores_calib)
  data_test <- tibble(d = obs_test, p_u = scores_test)

  if (method == "platt") {
    lr <- glm(d ~ p_u, family = binomial(link = 'logit'), data = data_calib)
    # Recalibrated scores on calibration and test set
    score_c_calib <- predict(lr, newdata = data_calib, type = "response")
    score_c_test <- predict(lr, newdata = data_test, type = "response")
    # Recalibrated values along a segment
    score_c_linspace <- predict(
      lr,
      newdata = tibble(p_u = linspace),
      type = "response"
    )
  } else if (method == "isotonic") {
    iso <- isoreg(x = data_calib$p_u, y = data_calib$d)
    fit_iso <- as.stepfun(iso)
    # Recalibrated scores on calibration and test set
    score_c_calib <- fit_iso(data_calib$p_u)
    score_c_test <- fit_iso(data_test$p_u)
    # Recalibrated values along a segment
    score_c_linspace <- fit_iso(linspace)
  } else if (method == "beta") {
    capture.output({
      bc <- beta_calibration(
        p = data_calib$p_u,
        y = data_calib$d,
        parameters = "abm" # 3 parameters a, b & m
      )
    })
    # Recalibrated scores on calibration and test set
    score_c_calib <- beta_predict(p = data_calib$p_u, bc)
    score_c_test <- beta_predict(p = data_test$p_u, bc)
    # Recalibrated values along a segment
    score_c_linspace <- beta_predict(linspace, bc)
  } else if (method == "locfit") {
    # Deg 0
    locfit_reg <- locfit(
      formula = d ~ lp(p_u, nn = iso_params$nn, deg = iso_params$deg),
      kern = "rect", maxk = 200, data = data_calib
    )
    # Recalibrated scores on calibration and test set
    score_c_calib <- predict(locfit_reg, newdata = data_calib)
    score_c_calib[score_c_calib < 0] <- 0
    score_c_calib[score_c_calib > 1] <- 1

    score_c_test <- predict(locfit_reg, newdata = data_test)
    score_c_test[score_c_test < 0] <- 0
    score_c_test[score_c_test > 1] <- 1

    # Recalibrated values along a segment
    score_c_linspace <- predict(locfit_reg, newdata = linspace)
    score_c_linspace[score_c_linspace < 0] <- 0
    score_c_linspace[score_c_linspace > 1] <- 1
  } else {
    stop(str_c(
      'Wrong method. Use one of the following:',
      '"platt", "isotonic", "beta", "locfit"'
    ))
  }

  # Format results in tibbles:
  # For calibration set
  tb_score_c_calib <- tibble(
    d = obs_calib,
    p_u = scores_calib,
    p_c = score_c_calib
  )
  # For test set
  tb_score_c_test <- tibble(
    d = obs_test,
    p_u = scores_test,
    p_c = score_c_test
  )
  # For linear space
  tb_score_c_linspace <- tibble(
    linspace = linspace,
    p_c = score_c_linspace
  )

  list(
    tb_score_c_calib = tb_score_c_calib,
    tb_score_c_test = tb_score_c_test,
    tb_score_c_linspace = tb_score_c_linspace
  )
}

#' Performs one replication for a simulation (with recalibration)
#'
#' @param i row number of the grid to use for the simulation
#' @param grid grid tibble with the seed number (column `seed`) and the deformations value (either `alpha` or `gamma`)
#' @param n_obs desired number of observation
#' @param type deformation probability type (either `alpha` or `gamma`); the
#' name should match with the `grid` tibble
#' @param linspace values at which to compute the mean observed event when computing the WMSE
f_simul_recalib <- function(i,
                            grid,
                            n_obs,
                            type = c("alpha", "gamma"),
                            linspace = NULL) {

  if (is.null(linspace)) linspace <- seq(0, 1, length.out = 100)

  ## 1. Generate Data----
  current_seed <- grid$seed[i]
  if (type == "alpha") {
    transform_scale <- grid$alpha[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = transform_scale, gamma = 1
    )
  } else if (type == "gamma") {
    transform_scale <- grid$gamma[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = 1, gamma = transform_scale
    )
  } else {
    stop("Transform type should be either alpha or gamma.")
  }

  ## 2. Calibration/Test sets----
  # Datasets with true probabilities
  data_all_calib <- current_data$data_all |>
    slice(current_data$calib_index)

  data_all_test <- current_data$data_all |>
    slice(-current_data$calib_index)

  ## 3. Recalibration----
  methods <- c("platt", "isotonic", "beta", "locfit", "locfit", "locfit")
  params <- list(
    NULL, NULL, NULL,
    list(nn = .15, deg = 0), list(nn = .15, deg = 1), list(nn = .15, deg = 2)
  )
  method_names <- c(
    "platt", "isotonic", "beta", "locfit_0", "locfit_1", "locfit_2"
  )
  res_recalibration <- map2(
    .x = methods,
    .y = params,
    .f = ~recalibrate_simul(
      obs_calib = data_all_calib$d,
      scores_calib = data_all_calib$p_u,
      obs_test = data_all_test$d,
      scores_test = data_all_test$p_u,
      method = .x,
      iso_params = .y,
      linspace = linspace
    )
  )
  names(res_recalibration) <- method_names

  ## 4. Calibration metrics----

  ### Using True Probabilities
  #### Calibration Set
  calib_metrics_true_calib <- compute_metrics(
    obs = data_all_calib$d,
    scores = data_all_calib$p,
    true_probas = data_all_calib$p,
    linspace = linspace) |>
    mutate(method = "True Prob.", sample = "Calibration")
  #### Test Set
  calib_metrics_true_test <- compute_metrics(
    obs = data_all_test$d,
    scores = data_all_test$p,
    true_probas = data_all_test$p,
    linspace = linspace) |>
    mutate(method = "True Prob.", sample = "Test")

  ### Without Recalibration
  #### Calibration Set
  calib_metrics_without_calib <- compute_metrics(
    obs = data_all_calib$d,
    scores = data_all_calib$p_u,
    true_probas = data_all_calib$p,
    linspace = linspace) |>
    mutate(method = "No Calibration", sample = "Calibration")
  #### Test Set
  calib_metrics_without_test <- compute_metrics(
    obs = data_all_test$d,
    scores = data_all_test$p_u,
    true_probas = data_all_test$p,
    linspace = linspace) |>
    mutate(method = "No Calibration", sample = "Test")

  calib_metrics <-
    calib_metrics_true_calib |>
    bind_rows(calib_metrics_true_test) |>
    bind_rows(calib_metrics_without_calib) |>
    bind_rows(calib_metrics_without_test)

  ### With Recalibration: loop on methods
  for (method in method_names) {
    res_recalibration_current <- res_recalibration[[method]]
    #### Calibration Set
    calib_metrics_without_calib <- compute_metrics(
      obs = data_all_calib$d,
      scores = res_recalibration_current$tb_score_c_calib$p_c,
      true_probas = data_all_calib$p,
      linspace = linspace) |>
      mutate(method = method, sample = "Calibration")
    #### Test Set
    calib_metrics_without_test <- compute_metrics(
      obs = data_all_test$d,
      scores = res_recalibration_current$tb_score_c_test$p_c,
      true_probas = data_all_test$p,
      linspace = linspace) |>
      mutate(method = method, sample = "Test")

    calib_metrics <-
      calib_metrics |>
      bind_rows(calib_metrics_without_calib) |>
      bind_rows(calib_metrics_without_test)
  }

  calib_metrics <-
    calib_metrics |>
    mutate(
      seed = current_seed,
      transform_scale = transform_scale,
      type = type
    )

  list(
    res_recalibration = res_recalibration,
    linspace = linspace,
    calib_metrics = calib_metrics,
    data_all_calib = data_all_calib,
    data_all_test = data_all_test,
    seed = current_seed
  )
}

# Goodness-of-fit----

#' Computes goodness of fit metrics for a replication
#'
#' @param i row number of the grid to use for the simulation
#' @param grid grid tibble with the seed number (column `seed`) and the deformations value (either `alpha` or `gamma`)
#' @param n_obs desired number of observation
#' @param type deformation probability type (either `alpha` or `gamma`); the
#' name should match with the `grid` tibble
compute_gof_simul <- function(i,
                              grid,
                              n_obs,
                              type = c("alpha", "gamma")) {
  current_seed <- grid$seed[i]
  if (type == "alpha") {
    transform_scale <- grid$alpha[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = transform_scale, gamma = 1
    )
  } else if (type == "gamma") {
    transform_scale <- grid$gamma[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = 1, gamma = transform_scale
    )
  } else {
    stop("Transform type should be either alpha or gamma.")
  }


  # Get the calib/test datasets with true probabilities
  data_all_calib <- current_data$data_all |>
    slice(current_data$calib_index)

  data_all_test <- current_data$data_all |>
    slice(-current_data$calib_index)

  # Calibration set
  true_prob_calib <- data_all_calib$p_u
  obs_calib <- data_all_calib$d
  pred_calib <- data_all_calib$p
  # Test set
  true_prob_test <- data_all_test$p_u
  obs_test <- data_all_test$d
  pred_test <- data_all_test$p

  metrics_simul_calib <- map(
    .x = seq(0, 1, by = .01),
    .f = ~compute_gof(
      true_prob = true_prob_calib,
      obs = obs_calib,
      pred = pred_calib,
      threshold = .x
    )
  ) |>
    list_rbind()

  metrics_simul_test <- map(
    .x = seq(0, 1, by = .01),
    .f = ~compute_gof(
      true_prob = true_prob_test,
      obs = obs_test,
      pred = pred_test,
      threshold = .x
    )
  ) |>
    list_rbind()

  roc_calib <- pROC::roc(obs_calib, pred_calib)
  auc_calib <- as.numeric(pROC::auc(roc_calib))
  roc_test <- pROC::roc(obs_test, pred_test)
  auc_test <- as.numeric(pROC::auc(roc_test))

  metrics_simul_calib |>
    mutate(
      auc = auc_calib,
      seed = current_seed,
      scale_parameter = transform_scale,
      type = type,
      sample = "calibration"
    ) |>
    bind_rows(
      metrics_simul_test |>
        mutate(
          auc = auc_test,
          seed = current_seed,
          scale_parameter = transform_scale,
          type = type,
          sample = "test"
        )
    )
}

#' Computes goodness of fit metrics for a replication (recalibration)
#'
#' @param i row number of the grid to use for the simulation
#' @param grid grid tibble with the seed number (column `seed`) and the deformations value (either `alpha` or `gamma`)
#' @param n_obs desired number of observation
#' @param type deformation probability type (either `alpha` or `gamma`); the
#' name should match with the `grid` tibble
compute_gof_simul_recalib <- function(i,
                                      grid,
                                      n_obs,
                                      type = c("alpha", "gamma")) {
  current_seed <- grid$seed[i]
  if (type == "alpha") {
    transform_scale <- grid$alpha[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = transform_scale, gamma = 1
    )
  } else if (type == "gamma") {
    transform_scale <- grid$gamma[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = 1, gamma = transform_scale
    )
  } else {
    stop("Transform type should be either alpha or gamma.")
  }


  # Get the calib/test datasets with true probabilities
  data_all_calib <- current_data$data_all |>
    slice(current_data$calib_index)

  data_all_test <- current_data$data_all |>
    slice(-current_data$calib_index)

  # Calibration set
  true_prob_calib <- data_all_calib$p_u
  obs_calib <- data_all_calib$d
  pred_calib <- data_all_calib$p

  # Test set
  true_prob_test <- data_all_test$p_u
  obs_test <- data_all_test$d
  pred_test <- data_all_test$p

  # Recalibration
  methods <- c("platt", "isotonic", "beta", "locfit", "locfit", "locfit")
  params <- list(
    NULL, NULL, NULL,
    list(nn = .15, deg = 0), list(nn = .15, deg = 1), list(nn = .15, deg = 2)
  )
  method_names <- c(
    "platt", "isotonic", "beta", "locfit_0", "locfit_1", "locfit_2"
  )
  res_recalibration <- map2(
    .x = methods,
    .y = params,
    .f = ~recalibrate_simul(
      obs_calib = data_all_calib$d,
      scores_calib = data_all_calib$p_u,
      obs_test = data_all_test$d,
      scores_test = data_all_test$p_u,
      method = .x,
      iso_params = .y,
      linspace = NULL
    )
  )
  names(res_recalibration) <- method_names

  # Initialisation
  gof_metrics_simul_calib <- tibble()
  gof_metrics_simul_test <- tibble()

  # Calculate standard metrics
  ## With Recalibration: loop on methods
  for (method in method_names) {
    res_recalibration_current <- res_recalibration[[method]]
    ### Computation of metrics on the calibration set
    metrics_simul_calib <- map(
      .x = seq(0, 1, by = .01), # we vary the probability threshold
      .f = ~compute_gof(
        true_prob = true_prob_calib,
        obs = obs_calib,
        #### the predictions are now recalibrated:
        pred = res_recalibration_current$tb_score_c_calib$p_c,
        threshold = .x
      )
    ) |>
      list_rbind()

    ### Computation of metricson the test set
    metrics_simul_test <- map(
      .x = seq(0, 1, by = .01), # we vary the probability threshold
      .f = ~compute_gof(
        true_prob = true_prob_test,
        obs = obs_test,
        #### the predictions are now recalibrated:
        pred = res_recalibration_current$tb_score_c_test$p_c,
        threshold = .x
      )
    ) |>
      list_rbind()

    roc_calib <- pROC::roc(
      obs_calib,
      res_recalibration_current$tb_score_c_calib$p_c
    )
    auc_calib <- as.numeric(pROC::auc(roc_calib))

    roc_test <- pROC::roc(
      obs_test,
      res_recalibration_current$tb_score_c_test$p_c
    )
    auc_test <- as.numeric(pROC::auc(roc_test))

    metrics_simul_calib <- metrics_simul_calib |>
      mutate(
        auc = auc_calib,
        seed = current_seed,
        scale_parameter = transform_scale,
        type = type,
        method = method,
        sample = "calibration"
      )

    metrics_simul_test <- metrics_simul_test |>
      mutate(
        auc = auc_test,
        seed = current_seed,
        scale_parameter = transform_scale,
        type = type,
        method = method,
        sample = "test"
      )

    gof_metrics_simul_calib <- gof_metrics_simul_calib |>
      bind_rows(metrics_simul_calib)
    gof_metrics_simul_test <- gof_metrics_simul_test |>
      bind_rows(metrics_simul_test)
  }

  gof_metrics_simul_calib |>
    bind_rows(gof_metrics_simul_test)
}

# Visualization Tools----



#' Boxplots for the simulations to visualize the distribution of some
#' traditional metrics as a function of the probability threshold.
#' And, ROC curves
#' The resulting figure is a panel of graphs, with vayring values for the
#' transformation applied to the probabilities (in columns) and different
#' metrics (in rows).
#'
#' @param tb_metrics tibble with computed metrics for the simulations
#' @param type type of transformation: `"alpha"` or `"gamma"`
#' @param metrics names of the metrics computed
boxplot_simuls_metrics <- function(tb_metrics,
                                   type = c("alpha", "gamma"),
                                   metrics) {
  scale_parameters <- unique(tb_metrics$scale_parameter)

  par(mfrow = c(length(metrics), length(scale_parameters)))
  for (i_metric in 1:length(metrics)) {
    metric <- metrics[i_metric]
    for (i_scale_parameter in 1:length(scale_parameters)) {
      scale_parameter <- scale_parameters[i_scale_parameter]

      tb_metrics_current <- tb_metrics |>
        filter(scale_parameter == !!scale_parameter)

      if (metric == "roc") {
        seeds <- unique(tb_metrics_current$seed)
        if (i_metric == 1) {
          # first row
          title <- latex2exp::TeX(
            str_c("$\\", type, " = ", round(scale_parameter, 2), "$")
          )
          size_top <- 2.1
        } else if (i_metric == length(metrics)) {
          # Last row
          title <- ""
          size_top <- 1.1
        } else {
          title <- ""
          size_top <- 1.1
        }

        if (i_scale_parameter == 1) {
          # first column
          y_lab <- str_c(metric, "\n True Positive Rate")
          size_left <- 5.1
        } else {
          y_lab <- ""
          size_left <- 4.1
        }

        par(mar = c(4.5, size_left, size_top, 2.1))
        plot(
          0:1, 0:1,
          type = "l", col = NULL,
          xlim = 0:1, ylim = 0:1,
          xlab = "False Positive Rate",
          ylab = y_lab,
          main = ""
        )
        for (i_seed in 1:length(seeds)) {
          tb_metrics_current_seed <-
            tb_metrics_current |>
            filter(seed == seeds[i_seed])
          lines(
            x = tb_metrics_current_seed$FPR,
            y = tb_metrics_current_seed$sensitivity,
            lwd = 2, col = adjustcolor("black", alpha.f = .04)
          )
        }
        segments(0, 0, 1, 1, col = "black", lty = 2)

      } else {
        # not ROC
        tb_metrics_current <-
          tb_metrics_current |>
          filter(threshold %in% seq(0, 1, by = .1))
        form <- str_c(metric, "~threshold")
        if (i_metric == 1) {
          # first row
          title <- latex2exp::TeX(
            str_c("$\\", type, " = ", round(scale_parameter, 2), "$")
          )
          size_top <- 2.1
        } else if (i_metric == length(metrics)) {
          # Last row
          title <- ""
          size_top <- 1.1
        } else {
          title <- ""
          size_top <- 1.1
        }

        if (i_scale_parameter == 1) {
          # first column
          y_lab <- metric
        } else {
          y_lab <- ""
        }

        par(mar = c(4.5, 4.1, size_top, 2.1))
        boxplot(
          formula(form), data = tb_metrics_current,
          xlab = "Threshold", ylab = y_lab,
          main = title
        )
      }
    }
  }
}

#' @param current_metric name of the metric to display (MSE, Brier Score,
#'   ECE, QMSE, or WMSE)
#' @param calib_metrics tibble with the metrics computed on the test set
#' @param type deformation probability type (either `alpha` or `gamma`); the
#' name should match with the `grid` tibble
plot_boxplot_metric <- function(current_metric,
                                calib_metrics,
                                type = c("alpha", "gamma"),
                                sample = "test") {

  data_plot <- calib_metrics |>
    filter(sample == !!sample) |>
    filter(metric == current_metric, type == !!type) |>
    mutate(
      label = str_c("$\\", type, "$=", scale_parameter, "$")
    )

  labels_y <- rep(unique(data_plot$label))
  par(mar = c(4.1, 4.1, 2.1, 2.1))
  boxplot(
    value ~ scale_parameter,
    data = data_plot,
    xlab = latex2exp::TeX(str_c("$\\", type, "$")),
    ylab = current_metric
  )
}

plot_boxplot_metric_2 <- function(tb_metrics, type, metrics, metrics_labs) {

  df_plot <- tb_metrics |> filter(threshold == 0.5)
  colours <- RColorBrewer::brewer.pal(
    5+1, name = "Blues"
  )
  colours <- colours[-1]
  colours[3] <- "orange"
  if (length(unique(df_plot$scale_parameter)) == 3) {
    colours <- colours[c(1, 3, 5)]
  }



  for (i_metric in 1:length(metrics)) {
    y_lab <- str_c("$", type, "=", sort(round(unique(df_plot$scale_parameter), 2)), "$")
    metric <- metrics[i_metric]
    metric_lab <- metrics_labs[i_metric]
    boxplot(
      formula(str_c(metric, "~ scale_parameter")),
      data = df_plot,
      las = 1,
      xlab = "", ylab = "",
      main = metric_lab,
      col = rev(colours),
      horizontal = TRUE,
      yaxt = "n"
    )
    axis(
      side = 2, at = 1:length(y_lab),
      labels = latex2exp::TeX(y_lab), las = 2
    )
  }
}


#' @param current_metric name of the metric to display (MSE, Brier Score,
#'   ECE, QMSE, or WMSE)
#' @param calib_metrics tibble with the metrics computed on the test set
#' @param type deformation probability type (either `alpha` or `gamma`); the
#' name should match with the `grid` tibble
plot_boxplot_metric_3 <- function(current_metric,
                                  calib_metrics,
                                  type = c("alpha", "gamma"),
                                  sample = "test") {

  data_plot <- calib_metrics |>
    filter(metric == current_metric, type == !!type) |>
    filter(sample == !!sample) |>
    mutate(
      label = str_c("$\\", type, "$=", round(scale_parameter, 2), "$")
    ) |>
    mutate(scale_parameter = round(scale_parameter, 2))


  if (as.character(current_metric) == "QMSE") {
    title <- "QBS"
  } else if (as.character(current_metric) == "MSE") {
    title <- "True MSE"
  } else {
    title <- current_metric
  }

  labels_y <- rep(unique(data_plot$label))
  colours <- RColorBrewer::brewer.pal(
    5+1, name = "Blues"
  )
  colours <- colours[-1]
  colours[3] <- "orange"
  if (length(labels_y) == 3) {
    colours <- colours[c(1, 3, 5)]
  }
  boxplot(
    value ~ scale_parameter,
    data = data_plot,
    ylab = "",
    xlab = "",
    las = 1,
    main = title,
    col = colours,
    horizontal = TRUE,
    yaxt = "n"
  )
  axis(
    side = 2, at = 1:length(labels_y),
    labels = latex2exp::TeX(labels_y), las = 2
  )
}

plot_boxplot_metric_recalib <- function(metric,
                                        calib_metrics_simul,
                                        type) {
  data_plot <- calib_metrics_simul |>
    filter(metric == !!metric, type == !!type) |>
    arrange(transform_scale)

  methods <- levels(data_plot$method)
  labels_y <- unique(data_plot$transform_scale) |> round(2)


  par(mfrow = c(4,2))
  for (method in methods) {
    data_plot_current <- data_plot |> filter(method == !!method)
    # par(mar = c(2.1, 12.1, 2.1, 2.1))
    par(mar = c(3.1, 4.1, 2.1, 2.1))
    boxplot(
      value ~ sample + transform_scale,
      data = data_plot_current,
      col = c("#D55E00", "#009E73"),
      horizontal = FALSE,
      main = method,
      las = 1, xlab = "", ylab = "",
      xaxt = "n"
    )
    # ind_benchmark <- which(labels_y == 1)
    labs_y <- str_c("$\\", type, "=", labels_y, "$")
    # labs_y[ind_benchmark] <- str_c(labs_y[ind_benchmark], " (benchmark)")
    axis(
      side = 1, at = seq(1, 2*length(labels_y), by = 2) + .5,
      labels = latex2exp::TeX(labs_y),
      las = 1,
      # col.axis = "black"
    )
    # # Horizontal lines
    # for (i in seq(1, 2*(length(labels_y)-1), by = 2) + 1.5) {
    #   abline(h = i, lty = 1, col = "gray")
    # }
    # Vertical lines
    for (i in seq(1, 2*(length(labels_y)-1), by = 2) + 1.5) {
      abline(v = i, lty = 1, col = "gray")
    }
  }
}

boxplot_std_metrics_calib <- function(tb_calib_metrics, metric, x_lim = NULL) {
  scale_parameters <- unique(tb_calib_metrics$transform_scale)
  nb <- length(scale_parameters)
  mat <- mat_init <- matrix(c(1:(2*nb), rep(2*nb+1, nb)), ncol = nb, byrow = TRUE)

  # colours_calib <- c(
  #   "#332288",
  #   # "#117733",
  #   "#44AA99", "#88CCEE",
  #   "#DDCC77", "#CC6677", "#AA4499", "#882255") |> rev()

  colours_names <- c(
    "True Prob.",
    "No Calibration", "Platt", "Isotonic", "Beta",
    "Locfit (deg=0)", "Locfit (deg=1)", "Locfit (deg=2)"
  )
  colours_calib <- colours_legend <- c(
    "#332288", "#117733", "#44AA99", "#88CCEE",
    "#DDCC77", "#CC6677", "#AA4499", "#882255"
  )
  colours_test <- adjustcolor(colours_calib, alpha.f = .5)

  colours <- NULL
  for (k in 1:length(colours_calib))
    colours <- c(colours, colours_calib[k], colours_test[k])


  layout(mat, heights = c(.45, .45, .15))

  for (type in unique(tb_calib_metrics$type)) {
    for (i_scale in 1:length(scale_parameters)) {
      scale_parameter <- scale_parameters[i_scale]
      tb_metrics_current <- tb_calib_metrics |>
        filter(
          transform_scale == !!scale_parameter,
          type == !!type,
          metric == !!metric
        ) |>
        mutate(
          method = fct_rev(fct_drop(method))
        )
      title <- latex2exp::TeX(
        str_c("$\\", type, " = ", round(scale_parameter, 2), "$")
      )
      methods_bp <- tb_metrics_current$method |> levels()

      ind_colours <- match(methods_bp, colours_names)

      colours <- NULL
      for (k in 1:length(colours_calib[ind_colours]))
        colours <- c(colours, colours_calib[ind_colours][k],
                     colours_test[ind_colours][k])


      form <- str_c("value~sample + method")
      par(mar = c(1.5, 1.5, 3.1, 1))
      boxplot(
        formula(form),
        data = tb_metrics_current |>
          mutate(
            sample = fct_rev(sample)
            # method = fct_rev(method)
          ),
        xlab = "",
        ylab = "",
        main = title,
        horizontal = TRUE,
        las = 1,
        col = colours,
        ylim = x_lim,
        border = c("black", adjustcolor("black", alpha.f = .5)),
        # sep = ", "
        yaxt = "n"
      )
      # Horizontal lines
      for (i in seq(3, length(methods_bp) * 2, by = 2) - .5) {
        abline(h = i, lty = 1, col = "gray")
      }
    }
  }
  par(mar = c(0, 4.3, 0, 4.3))
  plot.new()
  legend(
    "center",
    legend = colours_names,
    fill = colours_legend,
    # lwd = 2,
    xpd = TRUE, ncol = 4
  )
}

# Calibration Curves----

## Quantile-based----

#' Get the calibration curve for one simulation
#' using the quantile-based approach
#'
#' @param i row number of the grid to use for the simulation
#' @param grid grid tibble with the seed number (column `seed`) and the deformations value (either `alpha` or `gamma`)
#' @param n_obs desired number of observation
#' @param type deformation probability type (either `alpha` or `gamma`); the
#' name should match with the `grid` tibble
#' @param linspace values at which to compute the mean observed event when computing the WMSE
calibration_curve_quant_simul <- function(i,
                                          grid,
                                          n_obs,
                                          type = c("alpha", "gamma"),
                                          linspace = NULL) {

  if (is.null(linspace)) linspace <- seq(0, 1, length.out = 100)

  current_seed <- grid$seed[i]
  if (type == "alpha") {
    transform_scale <- grid$alpha[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = transform_scale, gamma = 1
    )
  } else if (type == "gamma") {
    transform_scale <- grid$gamma[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = 1, gamma = transform_scale
    )
  } else {
    stop("Transform type should be either alpha or gamma.")
  }

  # Get the test dataset with true probabilities
  data_all_test <- current_data$data_all |>
    slice(-current_data$calib_index)

  summary_bins <- get_summary_bins(
    obs = data_all_test$d,
    scores = data_all_test$p_u,
    k = 10, threshold = .5)
  summary_bins |>
    select(score_class, mean_score, mean_obs) |>
    mutate(
      seed = grid$seed[i],
      scale_parameter = transform_scale,
      type = type
    )
}

#' Counts the number of observation in evenly-spaced bind in [0,1]
#'
get_count <- function(seed, alpha = 1, gamma = 1) {
  breaks <- seq(0, 1, by = .05)
  current_data <- get_samples(seed = seed, n_obs = 2000, alpha = alpha, gamma = gamma)
  data_all_test <- current_data$data_all |>
    slice(-current_data$calib_index)

  nb_0_bins <- table(cut(data_all_test |> filter(d == 0) |> pull(p_u), breaks = breaks))
  nb_1_bins <- table(cut(data_all_test |> filter(d == 1) |> pull(p_u), breaks = breaks))
  tibble(
    bins = names(nb_0_bins),
    nb_0_bins = as.vector(nb_0_bins),
    nb_1_bins = as.vector(nb_1_bins),
    seed = seed,
    alpha = alpha,
    gamma = gamma
  )
}

#' @param i index of the simulation to use (in `simul_recalib_alpha` or
#'   `simul_recalib_gamma`)
#' @param type type of transformed probabilities (made on `alpha` or `gamma`)
#' @param method name of the recalibration method to focus on
get_count_simul_recalib <- function(i, type, method) {
  if (type == "alpha") {
    simul <- simul_recalib_alpha[[i]]
    transform_scale <- grid_alpha$alpha[i]
  } else if (type == "gamma") {
    simul <- simul_recalib_gamma[[i]]
    transform_scale <- grid_gamma$gamma[i]
  } else {
    stop("Wrong value for argument `type`.")
  }

  # Counting number of obs in bins defined over [0,1]
  breaks <- seq(0, 1, by = .05)
  if (method == "True Prob.") {
    scores_calib <- simul$data_all_calib$p
    scores_test <- simul$data_all_test$p
    scores_c_calib <- scores_c_test <- NULL
  } else if (method == "No Calibration") {
    scores_calib <- simul$data_all_calib$p_u
    scores_test <- simul$data_all_test$p_u
    scores_c_calib <- scores_c_test <- NULL
  } else {
    tb_score_c_calib <- simul$res_recalibration[[method]]$tb_score_c_calib
    tb_score_c_test <- simul$res_recalibration[[method]]$tb_score_c_test
    scores_calib <- tb_score_c_calib$p_u
    scores_test <- tb_score_c_test$p_u
    scores_c_calib <- tb_score_c_calib$p_c
    scores_c_test <- tb_score_c_test$p_c
  }

  n_bins_calib <- table(cut(scores_calib, breaks = breaks))
  n_bins_test <- table(cut(scores_test, breaks = breaks))
  if (!is.null(scores_c_calib)) {
    n_bins_c_calib <- table(cut(scores_c_calib, breaks = breaks))
  } else {
    n_bins_c_calib <- NA_integer_
  }
  if (!is.null(scores_c_test)) {
    n_bins_c_test <- table(cut(scores_c_test, breaks = breaks))
  } else {
    n_bins_c_test <- NA_integer_
  }

  n_bins_scores <- tibble(
    bins = names(table(cut(breaks, breaks = breaks))),
    n_bins_calib = as.vector(n_bins_calib),
    n_bins_test = as.vector(n_bins_test),
    n_bins_c_calib = as.vector(n_bins_c_calib),
    n_bins_c_test = as.vector(n_bins_c_test),
    method = method,
    seed = simul$seed,
    type = type,
    transform_scale = transform_scale
  )
  n_bins_scores
}

#' @param simul a single replication result
#' @param method name of the method used to recalibrate for which to compute the calibration curve
#' @param k number of bins to create (quantiles, default to `10`)
get_summary_bins_simul <- function(simul, method, k = 10) {
  obs_calib <- simul$data_all_calib$d
  obs_test <- simul$data_all_test$d

  if (method == "True Prob.") {
    scores_calib <- simul$data_all_calib$p
    scores_test <- simul$data_all_test$p
  } else if (method == "No Calibration") {
    scores_calib <- simul$data_all_calib$p_u
    scores_test <- simul$data_all_test$p_u
  } else {
    tb_score_c_calib <- simul$res_recalibration[[method]]$tb_score_c_calib
    tb_score_c_test <- simul$res_recalibration[[method]]$tb_score_c_test
    scores_calib <- tb_score_c_calib$p_c
    scores_test <- tb_score_c_test$p_c
  }

  summary_bins_calib <- get_summary_bins(
    obs = obs_calib, scores = scores_calib, k = k)
  summary_bins_test <- get_summary_bins(
    obs = obs_test, scores = scores_test, k = k)

  summary_bins_calib |> mutate(sample = "Calibration") |>
    bind_rows(summary_bins_test |> mutate(sample = "Test")) |>
    mutate(method = method, seed = simul$seed)
}

#' Plot the calibration curve for a set of simulations, using quantile-based
#' bins to obtain the estimations
#'
#' @param tb_calibration_curve
#' @param counts_samples tibble with the average number of observation with
#'   true probability in each bin defined along the [0,1] segment.
#' @param type deformation probability type (either `alpha` or `gamma`)
plot_calibration_quant_simul <- function(tb_calibration_curve,
                                         counts_samples,
                                         type = c("alpha", "gamma")) {
  tb_calibration_curve_curr <-  tb_calibration_curve |>
    filter(type == !!type)

  scale_params <- unique(tb_calibration_curve_curr$scale_parameter)
  seeds <- unique(tb_calibration_curve_curr$seed)

  for (scale_param in scale_params) {
    title <- str_c("$\\", type, " = $", round(scale_param, 2))

    # Histogram
    ## Calculate the heights for stacking
    if (type == "alpha") {
      counts_samples_current <-
        counts_samples |>
        filter(gamma == 1, alpha == !!scale_param)
    } else {
      counts_samples_current <-
        counts_samples |>
        filter(alpha == 1, gamma == !!scale_param)
    }
    heights <- rbind(
      counts_samples_current$nb_0_bins,
      counts_samples_current$nb_1_bins
    )
    col_bars <- c("#CC79A7", "#E69F00")
    par(mar = c(0.5, 4.5, 3.0, 0.5))
    barplot(
      heights,
      col = col_bars,
      border = "white",
      space = 0,
      xlab = "", ylab = "", main = latex2exp::TeX(title),
      axes = FALSE,
    )
    par(mar = c(4.1, 4.5, 0.5, 0.5))
    colour <- ifelse(scale_param == 1, "#D55E00", "#0072B2")
    plot(
      0:1, 0:1,
      type = "l", col = NULL,
      xlim = 0:1, ylim = 0:1,
      xlab = latex2exp::TeX("$p^u$"),
      ylab = latex2exp::TeX("$\\hat{E}(D | p^u = p^c)$"),
      main = ""
    )
    for (i_simul in seeds) {
      tb_current <- tb_calibration_curve_curr |>
        filter(
          scale_parameter == scale_param,
          seed == i_simul
        )
      lines(
        tb_current$mean_score, tb_current$mean_obs,
        lwd = 1, cex = .1,
        col = adjustcolor(colour, alpha.f = 0.1), t = "b",
      )
      segments(0, 0, 1, 1, col = "black", lty = 2)
    }
  }
}

plot_calibration_quant_simul_recalib <- function(method, type) {
  tb_calibration_curve <- summary_bins_simuls |>
    filter(
      method == !!method,
      type == !!type
    )

  if (type == "alpha") {
    count_scores <- count_scores_alpha |> filter(method == !!method)
  } else if (type == "gamma") {
    count_scores <- count_scores_gamma |> filter(method == !!method)
  } else {
    stop("Type should be either \"alpha\" or \"gamma\"")
  }

  scale_params <- unique(tb_calibration_curve$scale_parameter)
  seeds <- unique(tb_calibration_curve$seed)
  colours <- c("Calibration" = "#D55E00", "Test" = "#009E73")

  nb <- length(scale_params)
  mat <- mat_init <- matrix(c(1:4), ncol = 2)
  for (j in 1:(nb-1)) {
    mat <- rbind(mat, mat_init + j * 4)
  }
  layout(mat, heights = rep(c(1, 3), nb))

  y_lim <- c(0, 1)

  for (scale_param in scale_params) {
    title <- str_c("$\\", type, " = $", round(scale_param, 2))

    for (sample in c("Calibration", "Test")) {
      if (sample == "Calibration"){
        n_bins_current <- count_scores |>
          filter(transform_scale == !!scale_param) |> pull("n_bins_calib")
        n_bins_c_current <- count_scores |>
          filter(transform_scale == !!scale_param) |> pull("n_bins_c_calib")
      } else {
        n_bins_current <- count_scores |>
          filter(transform_scale == !!scale_param) |> pull("n_bins_test")
        n_bins_c_current <- count_scores |>
          filter(transform_scale == !!scale_param) |> pull("n_bins_c_test")
      }
      par(mar = c(0.5, 4.3, 1.0, 0.5))
      y_lim_bp <- range(c(n_bins_current, n_bins_c_current), na.rm = TRUE)
      barplot(
        n_bins_current,
        col = adjustcolor("#000000", alpha.f = .3),
        ylim = y_lim_bp,
        border = "white",
        axes = FALSE,
        xlab = "", ylab = "",
        main = latex2exp::TeX(title)
      )
      barplot(
        n_bins_c_current,
        col = adjustcolor("#0072B2", alpha.f = .3),
        ylim = y_lim_bp,
        border = "white",
        axes = FALSE,
        xlab = "", ylab = "", main = "",
        add = TRUE
      )
      par(mar = c(4.1, 4.3, 0.5, 0.5))
      plot(
        0:1, 0:1,
        type = "l", col = NULL,
        xlim = 0:1, ylim = 0:1,
        xlab = latex2exp::TeX("$p^u$"),
        ylab = latex2exp::TeX("$\\hat{E}(D | p^u = p^c)$"),
        main = ""
      )
      for (i_simul in seeds) {
        tb_current <- tb_calibration_curve |>
          filter(
            scale_parameter == scale_param,
            seed == i_simul,
            sample == !!sample
          )
        lines(
          tb_current$mean_score, tb_current$mean_obs,
          lwd = 2, col = adjustcolor(colours[sample], alpha.f = 0.1), t = "b",
          cex = .1, pch = 19
        )
      }
      segments(0, 0, 1, 1, col = "black", lty = 2)
    }
  }
}

plot_calibration_quant_simul_recalib_2 <- function(calib_curve, method, type) {
  tb_calibration_curve_both <- calib_curve |>
    filter(
      method == !!method
    )

  mat <- c(
    1, 3, 13, 15,
    2, 4, 14, 16,
    5, 7, 17, 19,
    6, 8, 18, 20,
    9, 11, 21, 23,
    10, 12, 22, 24
  ) |>
    matrix(ncol = 4, byrow = TRUE)

  layout(mat, height = rep(c(1, 3), 3))

  y_lim <- c(0, 1)

  for (type in c("alpha", "gamma")) {

    if (type == "alpha") {
      count_scores <- count_scores_alpha |> filter(method == !!method)
    } else if (type == "gamma") {
      count_scores <- count_scores_gamma |> filter(method == !!method)
    } else {
      stop("Type should be either \"alpha\" or \"gamma\"")
    }

    tb_calibration_curve <-
      tb_calibration_curve_both |>
      filter(type == !!type)

    scale_params <- unique(tb_calibration_curve$scale_parameter)
    seeds <- unique(tb_calibration_curve$seed)
    colours <- c("Calibration" = "#D55E00", "Test" = "#009E73")

    for (scale_param in scale_params) {
      title <- str_c("$\\", type, " = $", round(scale_param, 2))

      for (sample in c("Calibration", "Test")) {

        if (sample == "Calibration"){
          n_bins_current <- count_scores |>
            filter(transform_scale == !!scale_param) |> pull("n_bins_calib")
          n_bins_c_current <- count_scores |>
            filter(transform_scale == !!scale_param) |> pull("n_bins_c_calib")
        } else {
          n_bins_current <- count_scores |>
            filter(transform_scale == !!scale_param) |> pull("n_bins_test")
          n_bins_c_current <- count_scores |>
            filter(transform_scale == !!scale_param) |> pull("n_bins_c_test")
        }

        par(mar = c(0.5, 4.3, 1.0, 0.5))
        y_lim_bp <- range(c(n_bins_current, n_bins_c_current), na.rm = TRUE)
        barplot(
          n_bins_current,
          col = adjustcolor("#000000", alpha.f = .3),
          ylim = y_lim_bp,
          border = "white",
          axes = FALSE,
          xlab = "", ylab = "",
          main = latex2exp::TeX(title)
        )
        barplot(
          n_bins_c_current,
          col = adjustcolor("#0072B2", alpha.f = .3),
          ylim = y_lim_bp,
          border = "white",
          axes = FALSE,
          xlab = "", ylab = "", main = "",
          add = TRUE
        )

        tb_current <-
          tb_calibration_curve |>
          filter(
            scale_parameter == !!scale_param,
            sample == !!sample
          )
        par(mar = c(4.1, 4.3, 0.5, 0.5), mgp = c(2, 1, 0))
        plot(
          0:1, 0:1,
          type = "l", col = NULL,
          xlim = 0:1, ylim = 0:1,
          xlab = latex2exp::TeX("$p^u$"),
          ylab = latex2exp::TeX("$E(D | p^u = p^c)$"),
          main = ""
        )
        for (i_simul in seeds) {
          tb_current <- tb_calibration_curve |>
            filter(
              scale_parameter == scale_param,
              seed == i_simul,
              sample == !!sample
            )
          lines(
            tb_current$mean_score, tb_current$mean_obs,
            lwd = 2, col = adjustcolor(colours[sample], alpha.f = 0.1), t = "b",
            cex = .1, pch = 19
          )
          segments(0, 0, 1, 1, col = "black", lty = 2)
        }
      }
    }
  }
}

## Local Regression----

#' Get the calibration curve for one simulation
#' using the local regression-based approach
#'
#' @param i row number of the grid to use for the simulation
#' @param grid grid tibble with the seed number (column `seed`) and the
#'   deformations value (either `alpha` or `gamma`)
#' @param n_obs desired number of observation
#' @param type deformation probability type (either `alpha` or `gamma`); the
#' name should match with the `grid` tibble
calibration_curve_locfit_simul <- function(i,
                                           grid,
                                           n_obs,
                                           type = c("alpha", "gamma")) {


  current_seed <- grid$seed[i]

  if (type == "alpha") {
    transform_scale <- grid$alpha[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = transform_scale, gamma = 1
    )
  } else if (type == "gamma") {
    transform_scale <- grid$gamma[i]
    current_data <- get_samples(
      seed = current_seed, n_obs = n_obs, alpha = 1, gamma = transform_scale
    )
  } else {
    stop("Transform type should be either alpha or gamma.")
  }

  # Get the test dataset with true probabilities
  data_all_test <- current_data$data_all |>
    slice(-current_data$calib_index)

  locfit_0 <- locfit(
    formula = d ~ lp(p_u, nn = 0.15, deg = 0),
    kern = "rect", maxk = 200, data = data_all_test
  )

  scores <- data_all_test$p_u
  # Predictions on [0,1]
  linspace_raw <- seq(0, 1, length.out = 100)
  # Restricting this space to the range of observed scores
  keep_linspace <- which(linspace_raw >= min(scores) & linspace_raw <= max(scores))
  linspace <- linspace_raw[keep_linspace]

  score_c_locfit_0 <- predict(locfit_0, newdata = linspace)
  score_c_locfit_0[score_c_locfit_0 > 1] <- 1
  score_c_locfit_0[score_c_locfit_0 < 0] <- 0

  tibble(
    seed = grid$seed[i],
    scale_parameter = transform_scale,
    type = type,
    xlim = linspace,
    locfit_pred = score_c_locfit_0
  )
}

calibration_curve_locfit_simul_recalib <- function(simul,
                                                   method,
                                                   k = 10) {


  obs_calib <- simul$data_all_calib$d
  obs_test <- simul$data_all_test$d

  if (method == "True Prob.") {
    scores_calib <- simul$data_all_calib$p
    scores_test <- simul$data_all_test$p
  } else if (method == "No Calibration") {
    scores_calib <- simul$data_all_calib$p_u
    scores_test <- simul$data_all_test$p_u
  } else {
    tb_score_c_calib <- simul$res_recalibration[[method]]$tb_score_c_calib
    tb_score_c_test <- simul$res_recalibration[[method]]$tb_score_c_test
    scores_calib <- tb_score_c_calib$p_c
    scores_test <- tb_score_c_test$p_c
  }

  # Add a little noise (otherwise, R may crash...)
  scores_calib <- scores_calib + rnorm(length(scores_calib), 0, .001)
  scores_test <- scores_test + rnorm(length(scores_test), 0, .001)

  tb_calib <- tibble(
    obs = obs_calib,
    scores = scores_calib
  )

  tb_test <- tibble(
    obs = obs_test,
    scores = scores_test
  )

  locfit_0_calib <- locfit(
    formula = obs ~ lp(scores, nn = 0.15, deg = 0),
    kern = "rect", maxk = 200, data = tb_calib
  )

  # Predictions on [0,1]
  linspace_raw <- seq(0, 1, length.out = 100)

  # Restricting this space to the range of observed scores
  keep_linspace_calib <- which(
    linspace_raw >= min(scores_calib) & linspace_raw <= max(scores_calib)
  )
  linspace_calib <- linspace_raw[keep_linspace_calib]

  score_c_locfit_0_calib <- predict(locfit_0_calib, newdata = linspace_calib)
  score_c_locfit_0_calib[score_c_locfit_0_calib < 0] <- 0
  score_c_locfit_0_calib[score_c_locfit_0_calib > 1] <- 1

  locfit_0_test <- locfit(
    formula = obs ~ lp(scores, nn = 0.15, deg = 0),
    kern = "rect", maxk = 200, data = tb_test
  )

  keep_linspace_test <- which(
    linspace_raw >= min(scores_test) & linspace_raw <= max(scores_test)
  )
  linspace_test <- linspace_raw[keep_linspace_test]

  score_c_locfit_0_test <- predict(locfit_0_test, newdata = linspace_test)
  score_c_locfit_0_test[score_c_locfit_0_test < 0] <- 0
  score_c_locfit_0_test[score_c_locfit_0_test > 1] <- 1

  tb_calibration_curve_locfit <- tibble(
    xlim = linspace_calib,
    locfit_pred = score_c_locfit_0_calib,
    method = method,
    seed = simul$seed,
    sample = "calibration"
  ) |>
    bind_rows(
      tibble(
        xlim = linspace_test,
        locfit_pred = score_c_locfit_0_test,
        method = method,
        seed = simul$seed,
        sample = "test"
      )
    )

  tb_calibration_curve_locfit
}

#' Plot the calibration curve for a set of simulations, using local regression
#' to obtain the estimations
#'
#' @param tb_calibration_curve
#' @param type deformation probability type (either `alpha` or `gamma`)
#' @param counts_samples tibble with the average number of observation with
#'   true probability in each bin defined along the [0,1] segment.
plot_calibration_locfit_simuls <- function(tb_calibration_curve,
                                           type = c("alpha", "gamma"),
                                           counts_samples = counts_samples) {

  tb_calibration_curve_curr <-  tb_calibration_curve |>
    filter(type == !!type)

  scale_params <- unique(tb_calibration_curve_curr$scale_parameter)


  df_plot <- tb_calibration_curve_curr |>
    group_by(type, scale_parameter, xlim) |>
    summarise(
      mean = mean(locfit_pred),
      lower = quantile(locfit_pred, probs = .025),
      upper = quantile(locfit_pred, probs = .975),
      .groups = "drop"
    )

  for (scale_param in scale_params) {
    title <- str_c("$\\", type, " = $", round(scale_param, 2))
    # Histogram
    ## Calculate the heights for stacking
    if (type == "alpha") {
      counts_samples_current <-
        counts_samples |>
        filter(gamma == 1, alpha == !!scale_param)
    } else {
      counts_samples_current <-
        counts_samples |>
        filter(alpha == 1, gamma == !!scale_param)
    }
    heights <- rbind(
      counts_samples_current$nb_0_bins,
      counts_samples_current$nb_1_bins
    )
    col_bars <- c("#CC79A7", "#E69F00")
    par(mar = c(0.5, 4.5, 3.0, 0.5))
    barplot(
      heights,
      col = col_bars,
      border = "white",
      space = 0,
      xlab = "", ylab = "", main = latex2exp::TeX(title),
      axes = FALSE,
    )

    par(mar = c(4.1, 4.5, 0.5, 0.5))
    df_plot_current <- df_plot |> filter(scale_parameter == scale_param)
    colour <- ifelse(scale_param == 1, "#D55E00", "#0072B2")
    plot(
      df_plot_current$xlim, df_plot_current$mean,
      type = "l", col = colour,
      xlim = 0:1, ylim = 0:1,
      xlab = latex2exp::TeX("$p^u$"),
      ylab = latex2exp::TeX("$\\hat{E}(D | p^u = p^c)$"),
      main = ""
    )
    polygon(
      c(df_plot_current$xlim, rev(df_plot_current$xlim)),
      c(df_plot_current$lower, rev(df_plot_current$upper)),
      col = adjustcolor(col = colour, alpha.f = .4),
      border = NA
    )
    segments(0, 0, 1, 1, col = "black", lty = 2)
  }
}

plot_calibration_locfit_simuls_recalib <- function(calib_curve, method, type) {
  tb_calibration_curve <- calib_curve |>
    filter(
      method == !!method,
      type == !!type
    )
  if (type == "alpha") {
    count_scores <- count_scores_alpha |> filter(method == !!method)
  } else if (type == "gamma") {
    count_scores <- count_scores_gamma |> filter(method == !!method)
  } else {
    stop("Type should be either alpha or gamma")
  }

  scale_params <- unique(tb_calibration_curve$scale_parameter)
  seeds <- unique(tb_calibration_curve$seed)
  colours <- c("Calibration" = "#D55E00", "Test" = "#009E73")

  nb <- length(scale_params)
  mat <- mat_init <- matrix(c(1:4), ncol = 2)
  for (j in 1:(nb-1)) {
    mat <- rbind(mat, mat_init + j * 4)
  }
  layout(mat, heights = rep(c(1, 3), nb))

  y_lim <- c(0, 1)

  for (scale_param in scale_params) {
    title <- str_c("$\\", type, " = $", round(scale_param, 2))

    for (sample in c("Calibration", "Test")) {

      if (sample == "Calibration"){
        n_bins_current <- count_scores |>
          filter(transform_scale == !!scale_param) |> pull("n_bins_calib")
        n_bins_c_current <- count_scores |>
          filter(transform_scale == !!scale_param) |> pull("n_bins_c_calib")
      } else {
        n_bins_current <- count_scores |>
          filter(transform_scale == !!scale_param) |> pull("n_bins_test")
        n_bins_c_current <- count_scores |>
          filter(transform_scale == !!scale_param) |> pull("n_bins_c_test")
      }

      par(mar = c(0.5, 4.3, 1.0, 0.5))
      y_lim_bp <- range(c(n_bins_current, n_bins_c_current), na.rm = TRUE)
      barplot(
        n_bins_current,
        col = adjustcolor("#000000", alpha.f = .3),
        ylim = y_lim_bp,
        border = "white",
        axes = FALSE,
        xlab = "", ylab = "",
        main = latex2exp::TeX(title)
      )
      barplot(
        n_bins_c_current,
        col = adjustcolor("#0072B2", alpha.f = .3),
        ylim = y_lim_bp,
        border = "white",
        axes = FALSE,
        xlab = "", ylab = "", main = "",
        add = TRUE
      )

      tb_current <-
        tb_calibration_curve |>
        filter(
          scale_parameter == !!scale_param,
          sample == !!sample
        ) |>
        group_by(type, scale_parameter, xlim) |>
        summarise(
          mean = mean(locfit_pred),
          lower = quantile(locfit_pred, probs = .025),
          upper = quantile(locfit_pred, probs = .975),
          .groups = "drop"
        )
      par(mar = c(4.1, 4.3, 0.5, 0.5))
      plot(
        tb_current$xlim, tb_current$mean,
        type = "l", col = colours[sample],
        xlim = 0:1, ylim = 0:1,
        xlab = "Predicted probability", ylab = "Mean predicted probability",
        main = ""
      )
      polygon(
        c(tb_current$xlim, rev(tb_current$xlim)),
        c(tb_current$lower, rev(tb_current$upper)),
        col = adjustcolor(col = colours[sample], alpha.f = .4),
        border = NA
      )
      segments(0, 0, 1, 1, col = "black", lty = 2)
    }
  }
}

plot_calib_locfit_simuls_recalib_2 <- function(calib_curve, method) {
  tb_calibration_curve_both <- calib_curve |>
    filter(
      method == !!method
    )

  mat <- c(
    1, 3, 13, 15,
    2, 4, 14, 16,
    5, 7, 17, 19,
    6, 8, 18, 20,
    9, 11, 21, 23,
    10, 12, 22, 24
  ) |>
    matrix(ncol = 4, byrow = TRUE)

  layout(mat, height = rep(c(1, 3), 3))

  y_lim <- c(0, 1)

  for (type in c("alpha", "gamma")) {

    tb_calibration_curve <-
      tb_calibration_curve_both |>
      filter(type == !!type)

    if (type == "alpha") {
      count_scores <- count_scores_alpha |> filter(method == !!method)
    } else if (type == "gamma") {
      count_scores <- count_scores_gamma |> filter(method == !!method)
    }

    scale_params <- unique(tb_calibration_curve$scale_parameter)
    seeds <- unique(tb_calibration_curve$seed)
    colours <- c("Calibration" = "#D55E00", "Test" = "#009E73")

    for (scale_param in scale_params) {
      title <- str_c("$\\", type, " = $", round(scale_param, 2))

      for (sample in c("Calibration", "Test")) {

        if (sample == "Calibration"){
          n_bins_current <- count_scores |>
            filter(transform_scale == !!scale_param) |> pull("n_bins_calib")
          n_bins_c_current <- count_scores |>
            filter(transform_scale == !!scale_param) |> pull("n_bins_c_calib")
        } else {
          n_bins_current <- count_scores |>
            filter(transform_scale == !!scale_param) |> pull("n_bins_test")
          n_bins_c_current <- count_scores |>
            filter(transform_scale == !!scale_param) |> pull("n_bins_c_test")
        }

        par(mar = c(0.5, 4.3, 1.0, 0.5))
        y_lim_bp <- range(c(n_bins_current, n_bins_c_current), na.rm = TRUE)
        barplot(
          n_bins_current,
          col = adjustcolor("#000000", alpha.f = .3),
          ylim = y_lim_bp,
          border = "white",
          axes = FALSE,
          xlab = "", ylab = "",
          main = latex2exp::TeX(title)
        )
        barplot(
          n_bins_c_current,
          col = adjustcolor("#0072B2", alpha.f = .3),
          ylim = y_lim_bp,
          border = "white",
          axes = FALSE,
          xlab = "", ylab = "", main = "",
          add = TRUE
        )

        tb_current <-
          tb_calibration_curve |>
          filter(
            scale_parameter == !!scale_param,
            sample == !!sample
          ) |>
          group_by(type, scale_parameter, xlim, sample) |>
          summarise(
            mean = mean(locfit_pred),
            lower = quantile(locfit_pred, probs = .025),
            upper = quantile(locfit_pred, probs = .975),
            .groups = "drop"
          )
        par(mar = c(4.1, 4.3, 0.5, 0.5), mgp = c(2, 1, 0))
        plot(
          tb_current$xlim, tb_current$mean,
          type = "l", col = colours[sample],
          xlim = 0:1, ylim = 0:1,
          xlab = latex2exp::TeX("$p^u$"),
          ylab = latex2exp::TeX("$E(D | p^u = p^c)$"),
          main = ""
        )
        polygon(
          c(tb_current$xlim, rev(tb_current$xlim)),
          c(tb_current$lower, rev(tb_current$upper)),
          col = adjustcolor(col = colours[sample], alpha.f = .4),
          border = NA
        )
        segments(0, 0, 1, 1, col = "black", lty = 2)
      }
    }
  }
}
