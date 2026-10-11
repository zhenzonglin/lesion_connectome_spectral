#!/usr/bin/env Rscript
# 功能：258 特征、129 模态、20 组合的合成工程测试；不代表临床发现。
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_directory <- dirname(normalizePath(script_path))
for (file in c("common.R", "imputation.R", "inference.R", "patterns.R"))
  source(file.path(script_directory, file), encoding = "UTF-8")
arguments <- commandArgs(TRUE)
if (length(arguments) != 1L || !grepl("/archive/", arguments[1], fixed = TRUE)) stop("请指定 archive 内的新测试目录")
directory <- arguments[1]
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
set.seed(51623)
sample_size <- 320L
features <- c(sprintf("g_%03d", 1:129), sprintf("w_%03d", 1:129))
data <- data.frame(subject_id = sprintf("synthetic_%04d", seq_len(sample_size)),
  age = rnorm(sample_size, 63, 9), sex = factor(sample(c("1", "2"), sample_size, TRUE)),
  log_lesion_volume = rnorm(sample_size, 7, 1),
  prestroke_mrs = factor(sample(0:2, sample_size, TRUE)), acute_nihss = rpois(sample_size, 8),
  education = factor(sample(1:5, sample_size, TRUE)))
for (feature in features) data[[feature]] <- runif(sample_size)
targets <- c("mrs_3m", paste0("nihss_", c("1a", "1b", "1c", "2", "3", "4", "5a", "5b", "6a", "6b", "7", "8", "9", "10", "11", "consciousness", "motor_left", "motor_right", "motor_total")))
variables <- list(age = list(type = "continuous"), sex = list(type = "binary", levels = c("1", "2")),
  log_lesion_volume = list(type = "continuous"), prestroke_mrs = list(type = "ordinal", levels = 0:2),
  acute_nihss = list(type = "continuous"), education = list(type = "nominal", levels = 1:5))
base_covariates <- c("age", "sex", "log_lesion_volume")
outcomes <- list()
for (index in seq_along(targets)) {
  target <- targets[index]
  latent <- 6 * data$g_001 + 4 * data$w_001 + .02 * data$age + rlogis(sample_size)
  data[[target]] <- as.integer(cut(latent, breaks = c(-Inf, 5, 7, Inf))) - 1L
  outcomes[[index]] <- list(name = target, type = "ordinal", levels = 0:2,
    family = if (index == 1L) "mrs_single" else "nihss",
    covariates = c(base_covariates, if (index == 1L) c("prestroke_mrs", "acute_nihss")))
}
data$moca_3m <- 29 - 4 * data$g_001 - 3 * data$w_001 + rnorm(sample_size)
moca <- list(name = "moca_3m", type = "continuous", family = "moca", maximum = 30,
  covariates = c(base_covariates, "education"))
job <- list(output_dir = file.path(directory, "association"), feature_columns = features,
  variables = variables, outcomes = list(outcomes[[1]], outcomes[[2]], moca), reference_id = "synthetic_258",
  config = list(statistics = list(mi = list(m = 2L, maxit = 2L, increment = 2L,
    max_iterations = 2L, seed = 51623L), threshold_sensitivity = TRUE),
    patterns = list(targets = targets, folds = 2L, alpha = .5, n_lambda = 3L,
      lambda_min = .02, lambda_max = .2, seed = 51623L, maxiter = 300L,
      tolerance = 1e-5, selection = "one_standard_error")))
write_table(data, file.path(directory, "synthetic_input.csv"))
association <- association_stage(data, job)
stopifnot(nrow(association) == 387L, all(is.finite(association$p_joint)))
job$output_dir <- file.path(directory, "patterns")
job$outcomes <- outcomes
patterns <- pattern_stage(data, job)
stopifnot(nrow(patterns) == 20L, all(patterns$status == "frozen_pending_validation"))
for (target in targets) {
  model <- jsonlite::read_json(file.path(job$output_dir, target, "model.json"))
  stopifnot(length(model$feature_columns) == 258L, length(model$coefficients) == 258L)
}
# 全独立合成外部行；这里只验证入口及冻结模型，绝不称作真实外部复现。
validation <- data
validation$subject_id <- sub("synthetic", "external", validation$subject_id)
set.seed(926)
for (feature in features) validation[[feature]] <- runif(sample_size)
for (target in targets) {
  latent <- 6 * validation$g_001 + 4 * validation$w_001 + .02 * validation$age + rlogis(sample_size)
  validation[[target]] <- as.integer(cut(latent, breaks = c(-Inf, 5, 7, Inf))) - 1L
}
job$frozen_models_dir <- job$output_dir
job$output_dir <- file.path(directory, "validation")
result <- validation_stage(validation, job)
stopifnot(nrow(result) == 39L, all(is.finite(result$p_joint)))
write_json(list(status = "passed", synthetic_only = TRUE, patients = sample_size,
  features = 258L, association_rows = nrow(association), combinations = 20L,
  validation_tests = nrow(result), engineering_overrides = job$config), file.path(directory, "fulltest.json"))
cat("258 特征、129 模态、20 组合测试通过\n")
