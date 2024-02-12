# Measuring and visualizing calibration on synthetic data

library(tidyverse)
library(locfit)
library(binom)

# Functions to simulate data
source("functions/synthetic-data.R")

op <- par() # Initial graphical parameters

# 1. Data Generation----

alphas <- c(1/3, 1, 3)
data_alphas <- map(
  .x = alphas,
  .f = ~sim_data(
    n_obs = 2000, seed = 1, alpha = .x, gamma = 1,
    a = c(-0.1, 0.05, 0.2, -0.05)
  )
)

gammas <- c(c(1/3, 1, 3))
data_gammas <-map(
  .x = gammas,
  .f = ~sim_data(
    n_obs = 2000, seed = 1, alpha = 1, gamma = .x,
    a = c(-0.1, 0.05, 0.2, -0.05)
  )
)


## Viz of the transformed probabilities----

### Alpha----
true_prob_alphas <- map(data_alphas, "p")
transformed_prob_alphas <- map(data_alphas, "p_u")
# Ordering by increasing values of the true proba
order_true_alphas <- map(true_prob_alphas, order)
true_prob_alphas <- map2(
  .x = true_prob_alphas, .y = order_true_alphas, .f = ~.x[.y]
)
transformed_prob_alphas <- map2(
  .x = transformed_prob_alphas, .y = order_true_alphas, .f = ~.x[.y]
)

colours <- RColorBrewer::brewer.pal(
  length(true_prob_alphas)+1, name = "Blues"
)
colours <- colours[-1]
colours[alphas == 1] <- "orange"

par(mar = c(4.1, 4.3, 2.1, 0.5))
plot(
  x = true_prob_alphas[[1]],
  y = transformed_prob_alphas[[1]], type = "l",
  xlab = latex2exp::TeX("$p$"),
  ylab = latex2exp::TeX("$p^u$"),
  col = colours[1],
  ylim = c(0, 1)
)
for (i in 2:length(true_prob_alphas)) {
  lines(
    x = true_prob_alphas[[i]],
    y = transformed_prob_alphas[[i]],
    col = colours[i]
  )
}
legend(
  "bottomright", col = colours, lty = 1,
  legend = latex2exp::TeX(str_c("$\\alpha = ",round(alphas, 2), "$"))
)

### Gamma----
true_prob_gammas <- map(data_gammas, "p")
transformed_prob_gammas <- map(data_gammas, "p_u")
# Ordering by increasing values of the true proba
order_true_gammas <- map(true_prob_gammas, order)
true_prob_gammas <- map2(
  .x = true_prob_gammas, .y = order_true_gammas, .f = ~.x[.y])
transformed_prob_gammas <- map2(
  .x = transformed_prob_gammas, .y = order_true_gammas, .f = ~.x[.y])

colours <- RColorBrewer::brewer.pal(
  length(true_prob_gammas)+1, name = "Blues"
)
colours <- colours[-1]
colours[gammas == 1] <- "orange"

par(mar = c(4.1, 4.3, 2.1, 0.5))
plot(
  x = true_prob_gammas[[1]],
  y = transformed_prob_gammas[[1]], type = "l",
  xlab = latex2exp::TeX("$p$"),
  ylab = latex2exp::TeX("$p^u$"),
  col = colours[1],
  ylim = c(0, 1)
)
for (i in 2:length(true_prob_gammas)) {
  lines(
    x = true_prob_gammas[[i]],
    y = transformed_prob_gammas[[i]],
    col = colours[i]
  )
}
legend(
  "bottomright", col = colours, lty = 1,
  legend = latex2exp::TeX(str_c("$\\gamma = ",round(gammas, 2), "$"))
)



# Restore initial graphical parameters
par(op)

## With a histogram----

### Alpha----
p_alphas <- map(data_alphas, ~sort(.x$p_u))
colours <- RColorBrewer::brewer.pal(
  length(p_alphas)+1, name = "Blues"
)
colours <- colours[-1]
colours[alphas == 1] <- "orange"

par(mar = c(4.1, 4.3, 2.1, 0.5), mfrow = c(3,1))
hist(
  p_alphas[[1]], breaks = seq(0, 1, by = .05),
  col = adjustcolor(colours[1], alpha.f = .4),
  xlab = "p", ylab = "Freq.",
  xlim = c(0, 1),
  main = latex2exp::TeX(str_c("$\\alpha = ", round(alphas[1],2), "$"))
)
for (i in 2:length(p_alphas)) {
  hist(
    p_alphas[[i]], breaks = seq(0, 1, by = .05),
    col = adjustcolor(colours[i], alpha.f = .4),
    xlab = "p", ylab = "Freq.",
    main = latex2exp::TeX(str_c("$\\alpha = ", round(alphas[i], 2), "$"))
  )
}

### Gamma----
p_gammas <- map(data_gammas, ~sort(.x$p_u))
colours <- RColorBrewer::brewer.pal(
  length(p_gammas)+1, name = "Blues"
)
colours <- colours[-1]
colours[gammas == 1] <- "orange"

par(mar = c(4.1, 4.3, 2.1, 0.5), mfrow = c(3,1))
hist(
  p_gammas[[1]], breaks = seq(0, 1, by = .05),
  col = adjustcolor(colours[1], alpha.f = .4),
  xlab = "p", ylab = "Freq.",
  xlim = c(0, 1),
  main = latex2exp::TeX(str_c("$\\gamma = ", round(gammas[1],2), "$"))
)
for (i in 2:length(p_gammas)) {
  hist(
    p_gammas[[i]], breaks = seq(0, 1, by = .05),
    col = adjustcolor(colours[i], alpha.f = .4),
    xlab = "p", ylab = "Freq.",
    main = latex2exp::TeX(str_c("$\\gamma = ", round(gammas[i], 2), "$"))
  )
}

# Restore initial graphical parameters
par(op)

# 2. Measuring Standard Metrics----

# Functions to compute Goodness of fit for the simulated data
source("functions/metrics.R")
source("functions/simulations-synthetic.R")

# CHANGE TO 200 to get the results of the paper
n_repl <- 10 # number of replications

n_obs <- 2000 # number of observations to draw
grid_alpha <- expand_grid(alpha = c(1/3, 1, 3), seed = 1:n_repl)
grid_gamma <- expand_grid(gamma = c(1/3, 1, 3), seed = 1:n_repl)

## Varying Alpha----
library(future)
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_alpha))
  metrics_alpha <- furrr::future_map(
    .x = 1:nrow(grid_alpha),
    .f = ~{
      p()
      compute_gof_simul(
        i = .x,
        grid = grid_alpha,
        n_obs = n_obs,
        type = "alpha"
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

metrics_alpha <- list_rbind(metrics_alpha)
save(
  metrics_alpha,
  file = "output/standard_metrics_alpha.rda"
)


## Varying gamma----
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_gamma))
  metrics_gamma <- furrr::future_map(
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

metrics_gamma <- list_rbind(metrics_gamma)
save(
  metrics_gamma,
  file = "output/standard_metrics_gamma.rda"
)

# "accuracy", "sensitivity", "specificity", "roc", "auc"
metrics <- c("mse", "auc", "roc")
boxplot_simuls_metrics(
  tb_metrics = metrics_alpha |> filter(sample == "test"),
  type = "alpha", metrics = metrics
)

metrics <- c("accuracy", "sensitivity", "specificity", "auc")
metrics_labs <- c("Accuracy", "Sensitivity", "Specificity", "AUC")

par(mar = c(4.1, 4.1, 2.1, 2.1), mfrow = c(2,4))
plot_boxplot_metric_2(
  tb_metrics = metrics_alpha |> filter(scale_parameter %in% c(1/3, 1, 3)) |>
    filter(sample == "test"),
  type = "alpha",
  metrics = metrics, metrics_labs = metrics_labs
)
plot_boxplot_metric_2(
  tb_metrics = metrics_gamma |> filter(scale_parameter %in% c(1/3, 1, 3)) |>
    filter(sample == "test"),
  type = "gamma",
  metrics = metrics, metrics_labs = metrics_labs
)

# 3. Measuring Calibration----

## Alpha----
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_alpha))
  simul_alpha <- furrr::future_map(
    .x = 1:nrow(grid_alpha),
    .f = ~{
      p()
      f_simul(
        i = .x,
        grid = grid_alpha,
        n_obs = n_obs,
        type = "alpha",
        linspace = NULL
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

simul_alpha <- list_rbind(simul_alpha)

## Gamma----
nb_cores <- future::availableCores()-1
plan(multisession, workers = nb_cores)
progressr::with_progress({
  p <- progressr::progressor(steps = nrow(grid_gamma))
  simul_gamma <- furrr::future_map(
    .x = 1:nrow(grid_gamma),
    .f = ~{
      p()
      f_simul(
        i = .x,
        grid = grid_gamma,
        n_obs = n_obs,
        type = "gamma",
        linspace = NULL
      )
    },
    .options = furrr::furrr_options(seed = NULL)
  )
})

simul_gamma <- list_rbind(simul_gamma)

## Binding Results----
calib_metrics <-
  simul_alpha |>
  bind_rows(simul_gamma) |>
  pivot_longer(
    cols = c("mse", "brier", "ece", "qmse", "wmse", "lcs"),
    names_to = "metric", values_to = "value"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("mse", "brier", "ece", "qmse", "wmse", "lcs"),
      labels = c("MSE", "Brier Score", "ECE", "QMSE", "WMSE", "LCS")
    )
  )

par(op)
plot_boxplot_metric(
  current_metric = "MSE",
  calib_metrics = calib_metrics,
  type = "alpha"
)

plot_boxplot_metric(
  current_metric = "LCS",
  calib_metrics = calib_metrics,
  type = "alpha"
)

tb_lab <- expand_grid(
  metric = c("MSE", "Brier Score", "ECE", "LCS"),
  type = c("alpha", "gamma")
) |>
  mutate(
    metric = factor(metric, levels = c("MSE", "Brier Score", "ECE", "LCS")),
    type = factor(type, levels = c("alpha", "gamma"))
  ) |>
  arrange(type, metric) |>
  mutate(
    metric = as.character(metric)
  )
par(mfrow = c(2,4), mar = c(4.1, 2.5, 2.1, 2.1))
for (i in 1:nrow(tb_lab)) {
  plot_boxplot_metric_3(
    current_metric = tb_lab$metric[i],
    calib_metrics = calib_metrics,
    type = tb_lab$type[i]
  )
}

# 4. Calibration Curves----

## 4.1 Quantile-based----

### Alpha----
tb_calibration_curve_quant_alphas <- map(
  .x = 1:nrow(grid_alpha),
  .f = ~calibration_curve_quant_simul(
    i = .x,
    grid = grid_alpha,
    n_obs = n_obs,
    type = "alpha"
  )
) |>
  list_rbind()

### Gamma----
tb_calibration_curve_quant_gammas <- map(
  .x = 1:nrow(grid_gamma),
  .f = ~calibration_curve_quant_simul(
    i = .x,
    grid = grid_gamma,
    n_obs = n_obs,
    type = "gamma"
  )
) |>
  list_rbind()

### Merge----
tb_calibration_curve_quant <-
  tb_calibration_curve_quant_alphas |>
  bind_rows(tb_calibration_curve_quant_gammas)

### Counting scores in bins----
grid <- expand_grid(seed = 1:n_repl, alpha = alphas) |>
  mutate(gamma = 1) |>
  bind_rows(
    expand_grid(seed = 1:n_repl, gamma = gammas) |>
      mutate(alpha = 1)
  )

counts_samples <- map(
  .x = 1:nrow(grid),
  ~get_count(
    seed = grid$seed[.x],
    alpha = grid$alpha[.x],
    gamma = grid$gamma[.x]
  ),
  .progress = TRUE
) |>
  list_rbind() |>
  group_by(bins, alpha, gamma) |>
  summarise(
    nb_0_bins = mean(nb_0_bins),
    nb_1_bins = mean(nb_1_bins),
    .groups = "drop"
  )

# Save for later
save(counts_samples, file = "output/simul-calib-counts_samples.rda")



### Plots----
par(op)
mat <- matrix(1:6, nrow = 2)
layout(mat, heights = c(1,3))
plot_calibration_quant_simul(
  tb_calibration_curve = tb_calibration_curve_quant,
  type = "alpha",
  counts_samples = counts_samples
)

mat <- matrix(1:6, nrow = 2)
layout(mat, heights = c(1,3))
plot_calibration_quant_simul(
  tb_calibration_curve = tb_calibration_curve_quant,
  type = "gamma", counts_samples = counts_samples
)


## 4.2 Local Regression----

### Alpha----
tb_calibration_curve_locfit_alphas <- map(
  .x = 1:nrow(grid_alpha),
  .f = ~calibration_curve_locfit_simul(
    i = .x,
    grid = grid_alpha,
    n_obs = n_obs,
    type = "alpha"
  )
) |>
  list_rbind()

### Gamma----
tb_calibration_curve_locfit_gammas <- map(
  .x = 1:nrow(grid_gamma),
  .f = ~calibration_curve_locfit_simul(
    i = .x,
    grid = grid_gamma,
    n_obs = n_obs,
    type = "gamma"
  )
) |>
  list_rbind()

### Merge----
tb_calibration_curve_locfit <-
  tb_calibration_curve_locfit_alphas |>
  bind_rows(tb_calibration_curve_locfit_gammas)

### Plots----
par(op)

mat <- matrix(1:6, nrow = 2)
layout(mat, heights = c(1,3, 1,3))
plot_calibration_locfit_simuls(
  tb_calibration_curve = tb_calibration_curve_locfit,
  type = "alpha", counts_samples = counts_samples
)

mat <- matrix(1:6, nrow = 2)
layout(mat, heights = c(1,3, 1,3))
plot_calibration_locfit_simuls(
  tb_calibration_curve = tb_calibration_curve_locfit,
  type = "gamma", counts_samples = counts_samples
)

## 4.3 Moving Average----

### Alpha----
calib_curve_alpha_ci <- map(
  .x = alphas,
  .f = function(alpha) {
    linspace_raw <- seq(0, 1, length.out = 100)
    scores <- data_alphas[[which(alphas == alpha)]]$p_u
    keep_linspace <- which(
      linspace_raw >= min(scores) & linspace_raw <= max(scores)
    )
    linspace <- linspace_raw[keep_linspace]
    map(
      .x = linspace,
      .f = ~local_ci_scores(
        obs = data_alphas[[which(alphas == alpha)]]$d,
        scores = data_alphas[[which(alphas == alpha)]]$p_u,
        tau = .x,
        nn = .15, prob = .5, method = "probit")
    ) |>
      bind_rows() |>
      mutate(alpha = alpha)
  }
)

### Gamma----

calib_curve_gamma_ci <- map(
  .x = gammas,
  .f = function(gamma) {
    map(
      .x = seq(0, 1, length.out = 100),
      .f = ~local_ci_scores(
        obs = data_gammas[[which(gammas == gamma)]]$d,
        scores = data_gammas[[which(gammas == gamma)]]$p_u,
        tau = .x,
        nn = .15, prob = .5, method = "probit")
    ) |>
      bind_rows() |>
      mutate(gamma = gamma)
  }
)

### Plot----

par(op)

mat <- matrix(1:6, nrow = 2)
layout(mat, heights = c(1,3, 1,3))

for (i in 1:length(calib_curve_alpha_ci)) {
  calib_curve_alpha_ci_curr <- calib_curve_alpha_ci[[i]]
  alpha <- unique(calib_curve_alpha_ci_curr$alpha)
  title <- str_c("$\\alpha = $", round(alpha, 2))

  # Histogram
  ## Calculate the heights for stacking
  counts_samples_current <-
    counts_samples |>
    filter(gamma == 1, alpha == !!alpha)
  heights <- rbind(
    counts_samples_current$nb_0_bins,
    counts_samples_current$nb_1_bins
  )
  col_bars <- c("#CC79A7", "#E69F00")
  par(mar = c(0.5, 4.3, 3.0, 0.5))
  barplot(
    heights,
    col = col_bars,
    border = "white",
    space = 0,
    xlab = "", ylab = "", main = latex2exp::TeX(title),
    axes = FALSE,
  )

  par(mar = c(4.1, 4.3, 0.5, 0.5))

  plot(
    calib_curve_alpha_ci_curr$xlim, calib_curve_alpha_ci_curr$mean,
    type = "l", main = "",
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = latex2exp::TeX("$p^u$"),
    ylab = latex2exp::TeX("$\\hat{E}(D | p^u = p^c)$"),
  )
  col_ic <- ifelse(alpha == 1, "#D55E00", "#0072B2")
  polygon(
    c(calib_curve_alpha_ci_curr$xlim, rev(calib_curve_alpha_ci_curr$xlim)),
    c(calib_curve_alpha_ci_curr$lower, rev(calib_curve_alpha_ci_curr$upper)),
    col = adjustcolor(col = col_ic, alpha.f = .4),
    border = NA
  )
  segments(0, 0, 1, 1, col = "black", lty = 2)
}

par(op)
calib_curve_gamma_ci <- map(
  .x = gammas,
  .f = function(gamma) {
    map(
      .x = seq(0, 1, length.out = 100),
      .f = ~local_ci_scores(
        obs = data_gammas[[which(gammas == gamma)]]$d,
        scores = data_gammas[[which(gammas == gamma)]]$p_u,
        tau = .x,
        nn = .15, prob = .5, method = "probit")
    ) |>
      bind_rows() |>
      mutate(gamma = gamma)
  }
)

par(op)
mat <- matrix(1:6, nrow = 2)
layout(mat, heights = c(1,3, 1,3))

for (i in 1:length(calib_curve_gamma_ci)) {
  calib_curve_gamma_ci_curr <- calib_curve_gamma_ci[[i]]
  gamma <- unique(calib_curve_gamma_ci_curr$gamma)
  title <- str_c("$\\gamma = $", round(gamma, 2))

  # Histogram
  ## Calculate the heights for stacking
  counts_samples_current <-
    counts_samples |>
    filter(alpha == 1, gamma == !!gamma)
  heights <- rbind(
    counts_samples_current$nb_0_bins,
    counts_samples_current$nb_1_bins
  )
  col_bars <- c("#CC79A7", "#E69F00")
  par(mar = c(0.5, 4.3, 3.0, 0.5))
  barplot(
    heights,
    col = col_bars,
    border = "white",
    space = 0,
    xlab = "", ylab = "", main = latex2exp::TeX(title),
    axes = FALSE,
  )

  par(mar = c(4.1, 4.3, 0.5, 0.5))

  plot(
    calib_curve_gamma_ci_curr$xlim, calib_curve_gamma_ci_curr$mean,
    type = "l", main = "",
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = latex2exp::TeX("$p^u$"),
    ylab = latex2exp::TeX("$\\hat{E}(D | p^u = p^c)$"),
  )
  col_ic <- ifelse(gamma == 1, "#D55E00", "#0072B2")
  polygon(
    c(calib_curve_gamma_ci_curr$xlim, rev(calib_curve_gamma_ci_curr$xlim)),
    c(calib_curve_gamma_ci_curr$lower, rev(calib_curve_gamma_ci_curr$upper)),
    col = adjustcolor(col = col_ic, alpha.f = .4),
    border = NA
  )
  segments(0, 0, 1, 1, col = "black", lty = 2)
}

