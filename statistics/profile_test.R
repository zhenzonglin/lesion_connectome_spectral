#!/usr/bin/env Rscript
# 功能：七时点 mRS、19 项 NIHSS 与三个 MoCA 的冻结组合外部关联接口测试。
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_directory <- dirname(normalizePath(script_path))
for (file in c("common.R", "imputation.R", "inference.R", "patterns.R"))
  source(file.path(script_directory, file), encoding = "UTF-8")
arguments <- commandArgs(TRUE)
if (length(arguments) != 2L || !grepl("/archive/", arguments[1], fixed = TRUE))
  stop("使用方式：profile_test.R 新archive目录 已通过fulltest的目录")
directory <- arguments[1]
prior <- arguments[2]
data <- utils::read.csv(file.path(prior, "synthetic_input.csv"))
data$subject_id <- paste0("profile_external_", data$subject_id)
features <- c(sprintf("g_%03d", 1:129), sprintf("w_%03d", 1:129))
set.seed(976)
for (feature in features) data[[feature]] <- runif(nrow(data))
base_covariates <- c("age", "sex", "log_lesion_volume")
mrs_covariates <- c(base_covariates, "prestroke_mrs", "acute_nihss")
targets <- c("mrs_3m", paste0("nihss_", c("1a", "1b", "1c", "2", "3", "4", "5a", "5b", "6a", "6b", "7", "8", "9", "10", "11", "consciousness", "motor_left", "motor_right", "motor_total")))
outcomes <- list()
for (target in c(paste0("mrs_", c(3, 6, 12, 24, 36, 48, 60), "m"), targets[-1])) {
  latent <- 6 * data$g_001 + 4 * data$w_001 + .02 * data$age + rlogis(nrow(data))
  data[[target]] <- as.integer(cut(latent, c(-Inf, 5, 7, Inf))) - 1L
  is_mrs <- startsWith(target, "mrs_")
  outcomes[[length(outcomes) + 1L]] <- list(name = target, type = "ordinal", levels = 0:2,
    family = if (is_mrs) "mrs_single" else "nihss",
    covariates = if (is_mrs) mrs_covariates else base_covariates)
}
for (month in c(3, 6, 12)) {
  name <- paste0("moca_", month, "m")
  data[[name]] <- 29 - 4 * data$g_001 - 3 * data$w_001 + rnorm(nrow(data))
  outcomes[[length(outcomes) + 1L]] <- list(name = name, type = "continuous", family = "moca",
    maximum = 30, covariates = c(base_covariates, "education"))
}
job <- list(output_dir = directory, feature_columns = features, outcomes = outcomes,
  frozen_models_dir = file.path(prior, "patterns"), reference_id = "synthetic_258",
  variables = list(age = list(type = "continuous"), sex = list(type = "binary", levels = 1:2),
    log_lesion_volume = list(type = "continuous"), prestroke_mrs = list(type = "ordinal", levels = 0:2),
    acute_nihss = list(type = "continuous"), education = list(type = "nominal", levels = 1:5)),
  config = list(statistics = list(mi = list(m = 2, maxit = 2, max_iterations = 2)),
    patterns = list(targets = targets), validation = list(mrs_pattern_profile = TRUE)))
result <- validation_stage(data, job)
stopifnot(nrow(result) == 48L, all(is.finite(result$p_joint)),
  all(result$planned_family_tests[result$family == "nihss"] == 38L),
  all(result$planned_family_tests[result$family == "mrs_single"] == 7L),
  all(result$planned_family_tests[result$family == "moca"] == 3L))
# 四个家族分别校正；三个 MoCA 时点共享 moca 家族，失败项继续占据分母。
counts <- c(nihss = 19 * 129, mrs_single = 7 * 129, mrs_longitudinal = 129,
            moca = 3 * 129)
fixture <- do.call(rbind, lapply(names(counts), function(family) data.frame(family = family,
  p_joint = c(NA, seq(.001, .9, length.out = counts[[family]] - 1)))))
adjusted <- family_adjust(fixture)
for (family in names(counts)) {
  selected <- adjusted$family == family
  stopifnot(all(adjusted$planned_family_tests[selected] == counts[[family]]))
  finite <- selected & is.finite(adjusted$p_joint)
  stopifnot(max(abs(adjusted$q[finite] - p.adjust(adjusted$p_joint[finite], "BH", n = counts[[family]]))) < 1e-12)
}
write_json(list(status = "passed", synthetic_only = TRUE, validation_tests = 48,
  mrs_timepoints = 7, nihss_items = 19, moca_timepoints = 3,
  four_family_planned_tests = as.list(counts)), file.path(directory, "profile_test.json"))
cat("48 个冻结组合临床关联与四家族边界测试通过\n")
