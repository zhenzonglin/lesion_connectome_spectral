#!/usr/bin/env Rscript
# 功能：统一统计入口。使用方式：env/bin/Rscript statistics/run.R job.json
# 输入 JSON 只声明字段、类型和参数，不执行任意 R 表达式。

# ---------- 模块加载 ----------
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_directory <- dirname(normalizePath(script_path))
for (file in c("common.R", "imputation.R", "inference.R", "patterns.R"))
  source(file.path(script_directory, file), encoding = "UTF-8")

# ---------- 参数与来源冻结 ----------
arguments <- commandArgs(TRUE)
if (length(arguments) != 1L) stop("使用方式: Rscript statistics/run.R job.json")
job <- jsonlite::read_json(arguments[1], simplifyVector = FALSE)
required <- c("stage", "input_csv", "output_dir", "outcomes", "config")
if (!all(required %in% names(job))) stop("统计作业缺少字段: ", paste(setdiff(required, names(job)), collapse = ","))
if (!job$stage %in% c("imputation", "association", "patterns", "validation")) stop("未知统计阶段")
packages <- c("jsonlite", "mice", "mitml", "ordinal", "sandwich", "ordinalNet", "lme4")
for (package in packages) if (!requireNamespace(package, quietly = TRUE)) stop("缺少 R 包: ", package)
dir.create(job$output_dir, recursive = TRUE, showWarnings = FALSE)
signature <- object_hash(list(job = job[setdiff(names(job), c("stop_file"))],
  input = unname(tools::md5sum(job$input_csv)),
  code = unname(tools::md5sum(list.files(script_directory, pattern = "\\.R$", full.names = TRUE)))))
manifest_path <- file.path(job$output_dir, "job_manifest.json")
if (file.exists(manifest_path)) {
  previous <- jsonlite::read_json(manifest_path)
  if (!identical(previous$signature, signature)) stop("统计阶段输入、配置或代码已改变，请使用新的阶段目录")
} else write_json(list(signature = signature, stage = job$stage,
  input_checksum = unname(tools::md5sum(job$input_csv)),
  package_versions = setNames(lapply(packages, function(package) as.character(utils::packageVersion(package))), packages)), manifest_path)

# ---------- 可恢复阶段执行 ----------
status <- tryCatch({
  data <- read_input(job)
  check_stop(job)
  result <- switch(job$stage,
    imputation = {
      records <- lapply(job$outcomes, function(outcome) {
        check_stop(job)
        tryCatch({
          selected <- endpoint_data(data, outcome, job)
          imputed <- build_imputations(selected, outcome, job, file.path(job$output_dir, outcome$name))
          data.frame(outcome = outcome$name, status = imputed$status, n = nrow(selected), reason = "")
        }, safe_stop = function(error) stop(error), error = function(error)
          data.frame(outcome = outcome$name, status = "not_estimable", n = NA_integer_, reason = conditionMessage(error)))
      })
      result <- do.call(rbind, records)
      write_table(result, file.path(job$output_dir, "imputation_summary.csv"))
      result
    },
    association = association_stage(data, job),
    patterns = pattern_stage(data, job),
    validation = validation_stage(data, job))
  failed <- if ("status" %in% names(result)) sum(result$status == "not_estimable") else 0L
  write_json(list(status = if (failed > 0) "completed_with_unestimable_tests" else "completed",
                   stage = job$stage, rows = nrow(result), unestimable = failed,
                   signature = signature), file.path(job$output_dir, "status.json"))
  0L
}, safe_stop = function(error) {
  write_json(list(status = "stopped", stage = job$stage, reason = conditionMessage(error), signature = signature),
              file.path(job$output_dir, "status.json"))
  75L
}, error = function(error) {
  write_json(list(status = "failed", stage = job$stage, reason = conditionMessage(error), signature = signature),
              file.path(job$output_dir, "status.json"))
  message(conditionMessage(error))
  2L
})
quit(status = status)
