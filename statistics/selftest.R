#!/usr/bin/env Rscript
# 功能：统计工程自检；数据全部为合成，不可作为临床结果。
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_directory <- dirname(normalizePath(script_path))
for (file in c("common.R", "imputation.R", "inference.R", "patterns.R"))
  source(file.path(script_directory, file), encoding = "UTF-8")
arguments <- commandArgs(TRUE)
if (length(arguments) != 1L) stop("请指定 archive 内的新测试目录")
directory <- normalizePath(arguments[1], mustWork = FALSE)
if (!grepl("/archive/", directory, fixed = TRUE)) stop("测试只允许写入 archive")
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
set.seed(6127)
sample_size <- 180L
data <- data.frame(subject_id = sprintf("test_%04d", seq_len(sample_size)),
  age = rnorm(sample_size, 63, 9), sex = factor(sample(c("1", "2"), sample_size, TRUE)),
  log_volume = rnorm(sample_size, 7, 1))
features <- c(sprintf("g_%03d", 1:3), sprintf("w_%03d", 1:3))
for (feature in features) data[[feature]] <- runif(sample_size)
latent <- 2.5 * data$g_001 + 1.8 * data$w_001 + .03 * data$age + rlogis(sample_size)
data$mrs_3m <- as.integer(cut(latent, breaks = c(-Inf, 2, 3, 4, 5, 6, 7, Inf))) - 1L
data$moca_3m <- 28 - 3 * data$g_001 - 2 * data$w_001 + rnorm(sample_size)
data$age[sample(seq_len(sample_size), 18)] <- NA
data$sex[sample(seq_len(sample_size), 12)] <- NA
variables <- list(age = list(type = "continuous"), sex = list(type = "binary", levels = c("1", "2")),
  log_volume = list(type = "continuous"))
outcomes <- list(list(name = "mrs_3m", type = "ordinal", family = "mrs_single", levels = 0:6,
  covariates = c("age", "sex", "log_volume")), list(name = "moca_3m", type = "continuous",
    family = "moca", maximum = 30, covariates = c("age", "sex", "log_volume")))
job <- list(output_dir = file.path(directory, "association"), feature_columns = features,
  mode_pairs = lapply(1:3, function(mode) list(mode = mode, g = features[mode], w = features[mode + 3])),
  variables = variables, outcomes = outcomes, reference_id = "synthetic",
  config = list(statistics = list(mi = list(m = 3L, maxit = 3L, increment = 2L,
    max_iterations = 5L, seed = 6127L, psrf_threshold = 1.1), threshold_sensitivity = TRUE),
    patterns = list(targets = "mrs_3m", folds = 2L, n_lambda = 3L, lambda_min = .02,
      lambda_max = .2, seed = 6127L, maxiter = 300L, tolerance = 1e-5)))

# ---------- Rubin 与 D1 独立核对 ----------
estimates <- rbind(g = c(.4, .5, .45), w = c(.2, .25, .3))
covariances <- array(rep(c(.03, .005, .005, .04), 3), c(2, 2, 3),
  dimnames = list(c("g", "w"), c("g", "w"), NULL))
pooled <- pool_parameters(estimates, covariances)
# 独立计算 D1 的 Li-Raghunathan-Rubin F 与自由度，而不是再次调用被测库。
within <- apply(covariances, c(1, 2), mean)
between <- cov(t(estimates))
relative <- (1 + 1 / 3) * sum(diag(between %*% solve(within))) / 2
f_value <- drop(t(rowMeans(estimates)) %*% solve(within, rowMeans(estimates))) / (2 * (1 + relative))
t_value <- 2 * (3 - 1)
degrees <- t_value * (1 + 1 / 2) * (1 + 1 / relative)^2 / 2
expected_p <- pf(f_value, 2, degrees, lower.tail = FALSE)
stopifnot(abs(pooled$p_joint - expected_p) < 1e-8) # mitml 使用数值 Delta 导数。
expected <- apply(covariances, c(1, 2), mean) + (1 + 1 / 3) * cov(t(estimates))
stopifnot(max(abs(pooled$total_covariance - expected)) < 1e-12)
adjusted <- family_adjust(data.frame(family = c("nihss", "nihss", "moca"), p_joint = c(.01, NA, .02)))
stopifnot(abs(adjusted$q[1] - .02) < 1e-12, is.na(adjusted$q[2]), adjusted$q[3] == .02)

# ---------- MI 不修改影像、已观察结局及已观察协变量 ----------
selected <- endpoint_data(data, outcomes[[1]], job)
imputed <- build_imputations(selected, outcomes[[1]], job, file.path(directory, "mi_contract"))
for (index in seq_len(imputed$m)) {
  completed <- complete_covariates(selected, imputed, index)
  stopifnot(identical(completed[, features], selected[, features]),
    identical(completed$mrs_3m, selected$mrs_3m), !anyNA(completed[, outcomes[[1]]$covariates]))
  for (variable in outcomes[[1]]$covariates) {
    observed <- !is.na(selected[[variable]])
    stopifnot(identical(completed[[variable]][observed], selected[[variable]][observed]))
  }
}
association <- association_stage(data, job)
stopifnot(nrow(association) == 6L, all(is.finite(association$p_joint)),
  all(association$planned_family_tests == 3L))
diagnostics <- utils::read.csv(file.path(job$output_dir, "model_diagnostics.csv"))
stopifnot(any(is.finite(diagnostics$po_min_p[diagnostics$outcome == "mrs_3m"])))

# ---------- 训练折预处理与验证标签隔离 ----------
job$output_dir <- file.path(directory, "patterns")
job$outcomes <- outcomes[1]
patterns <- pattern_stage(data, job)
stopifnot(patterns$status[1] == "frozen_pending_validation")
model <- jsonlite::read_json(file.path(job$output_dir, "mrs_3m", "model.json"))
scores <- frozen_score(data, model, "synthetic")
altered <- data
altered$mrs_3m <- rev(altered$mrs_3m)
stopifnot(identical(scores, frozen_score(altered, model, "synthetic")))
train <- seq_len(nrow(data)) <= 120L
scaler <- image_scaler(data, features, train)
altered[!train, features] <- altered[!train, features] + 1000
stopifnot(identical(scaler, image_scaler(altered, features, train)))
job$frozen_models_dir <- job$output_dir
job$output_dir <- file.path(directory, "validation")
job$outcomes <- outcomes
external <- data
external$subject_id <- paste0("external_", external$subject_id)
validation <- validation_stage(external, job)
stopifnot(nrow(validation) == 2L, all(is.finite(validation$p_joint)))
write_json(list(status = "passed", synthetic_only = TRUE,
  checks = c("rubin", "d1", "planned_family_size", "mi_only_covariates",
    "ordinal_and_hc3", "pattern_learning", "frozen_score_label_invariance",
    "training_scaler_isolation", "external_association")), file.path(directory, "selftest.json"))
cat("统计工程自检通过\n")
