#!/usr/bin/env Rscript
# 功能：固定训练样本后，改变验证侧的标签和数值不能改变训练侧填补结果。
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_directory <- dirname(normalizePath(script_path))
for (file in c("common.R", "imputation.R", "inference.R", "patterns.R"))
  source(file.path(script_directory, file), encoding = "UTF-8")
directory <- commandArgs(TRUE)[1]
if (!grepl("/archive/", directory, fixed = TRUE)) stop("请使用 archive 测试目录")
set.seed(319)
size <- 150L
data <- data.frame(subject_id = sprintf("mi_%04d", 1:size), age = rnorm(size, 60, 10),
  sex = factor(sample(1:2, size, TRUE)), log_lesion_volume = rnorm(size, 8),
  g_001 = runif(size), w_001 = runif(size), mrs_3m = sample(0:3, size, TRUE))
data$age[seq(1, size, by = 7)] <- NA
data$sex[seq(2, size, by = 11)] <- NA
outcome <- list(name = "mrs_3m", type = "ordinal", levels = 0:3,
  covariates = c("age", "sex", "log_lesion_volume"))
job <- list(feature_columns = c("g_001", "w_001"),
  variables = list(age = list(type = "continuous"), sex = list(type = "binary", levels = 1:2)),
  config = list(statistics = list(mi = list(m = 3L, maxit = 4L, increment = 2L,
    max_iterations = 6L, seed = 319L))))
train <- seq_len(size) <= 100L
changed <- data
changed$g_001[!train] <- changed$g_001[!train] + 8
changed$age[!train] <- changed$age[!train] + 90
changed$age[c(120, 125, 130)] <- NA
changed$sex[c(123, 127)] <- NA
changed$mrs_3m[!train] <- sample(0:3, sum(!train), TRUE)
original <- fold_imputations(data, outcome, job, file.path(directory, "original"), train)
altered <- fold_imputations(changed, outcome, job, file.path(directory, "altered"), train)
for (index in seq_len(original$m)) {
  one <- complete_covariates(data, original, index)
  two <- complete_covariates(changed, altered, index)
  stopifnot(identical(one[train, outcome$covariates], two[train, outcome$covariates]))
}
write_json(list(status = "passed", changed_validation_values_and_labels = TRUE,
  changed_validation_missingness = TRUE,
  training_imputations_identical = TRUE), file.path(directory, "isolation.json"))
cat("训练侧插补独立性测试通过\n")
