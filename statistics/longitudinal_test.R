#!/usr/bin/env Rscript
# 功能：合成纵向 mRS、有缺失访视、随机截距和次级时间交互的工程验证。
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_directory <- dirname(normalizePath(script_path))
for (file in c("common.R", "imputation.R", "inference.R", "patterns.R"))
  source(file.path(script_directory, file), encoding = "UTF-8")
directory <- commandArgs(TRUE)[1]
if (!grepl("/archive/", directory, fixed = TRUE)) stop("必须使用 archive 测试目录")
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
set.seed(581)
size <- 160L
data <- data.frame(subject_id = sprintf("long_%04d", 1:size), age = rnorm(size, 60, 8),
  sex = factor(sample(1:2, size, TRUE)), log_lesion_volume = rnorm(size, 8),
  g_001 = runif(size), w_001 = runif(size))
intercept <- rnorm(size)
for (index in seq_along(c(3, 6, 12))) {
  latent <- 2 * data$g_001 + 2 * data$w_001 + intercept - index * .3 + rlogis(size)
  data[[paste0("mrs_", c(3, 6, 12)[index], "m")]] <- as.integer(cut(latent, c(-Inf, 1, 3, Inf))) - 1L
}
data$mrs_6m[sample(1:size, 20)] <- NA
data$mrs_12m[sample(1:size, 15)] <- NA
data$age[sample(1:size, 8)] <- NA
outcome <- list(name = "mrs_longitudinal", type = "longitudinal", family = "mrs_longitudinal",
  levels = 0:2, visit_columns = c("mrs_3m", "mrs_6m", "mrs_12m"), visit_times = c(3, 6, 12),
  covariates = c("age", "sex", "log_lesion_volume"))
job <- list(output_dir = directory, feature_columns = c("g_001", "w_001"),
  mode_pairs = list(list(mode = 1, g = "g_001", w = "w_001")), outcomes = list(outcome),
  variables = list(age = list(type = "continuous"), sex = list(type = "binary", levels = 1:2),
    log_lesion_volume = list(type = "continuous")),
  config = list(statistics = list(jmax_sensitivity = TRUE,
    mi = list(m = 3L, maxit = 4L, increment = 2L, max_iterations = 6L, seed = 581L),
    models = list(longitudinal = list(n_agq = 1L, time_interaction = TRUE, threshold_sensitivity = TRUE)))))
result <- association_stage(data, job)
stopifnot(nrow(result) == 1L, is.finite(result$p_joint), result$n == size,
  result$observations == size * 3 - 35L)
jmax <- utils::read.csv(file.path(directory, "jmax_sensitivity.csv"))
interaction <- utils::read.csv(file.path(directory, "time_interaction.csv"))
sensitivity <- utils::read.csv(file.path(directory, "threshold_sensitivity.csv"))
stopifnot(is.finite(jmax$p_joint), is.finite(interaction$p_joint), all(is.finite(sensitivity$p_joint)))
write_json(list(status = "passed", synthetic_only = TRUE, patients = size,
  observations = result$observations, secondary_time_interaction = interaction$p_joint,
  threshold_models_keep_random_intercept = TRUE, jmax_separate_sensitivity = TRUE),
  file.path(directory, "longitudinal_test.json"))
cat("纵向随机截距、时间交互和 Jmax 敏感性测试通过\n")
