# 功能：统计后端共享的输入合同、原子写入和安全检查。

# ---------- 简单数据与路径工具 ----------
`%or%` <- function(value, fallback) if (is.null(value)) fallback else value
as_chars <- function(value) as.character(unlist(value, use.names = FALSE))

write_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".part")
  jsonlite::write_json(value, temporary, auto_unbox = TRUE, pretty = TRUE,
                       null = "null", na = "null", digits = 16)
  if (!file.rename(temporary, path)) stop("无法提交文件: ", path)
}

write_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".part")
  saveRDS(value, temporary, compress = "xz")
  if (!file.rename(temporary, path)) stop("无法提交检查点: ", path)
}

write_table <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".part")
  utils::write.csv(value, temporary, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  if (!file.rename(temporary, path)) stop("无法提交表格: ", path)
}

check_stop <- function(job) {
  if (!is.null(job$stop_file) && file.exists(job$stop_file)) {
    condition <- structure(list(message = "已在统计检查点安全停止", call = NULL),
                           class = c("safe_stop", "error", "condition"))
    stop(condition)
  }
}

object_hash <- function(value) {
  temporary <- tempfile()
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary, version = 3, compress = FALSE)
  unname(tools::md5sum(temporary))
}

safe_name <- function(value) {
  if (length(value) != 1L || !grepl("^[a-zA-Z][a-zA-Z0-9_]*$", value))
    stop("统计变量名必须是字母起始的安全名称: ", paste(value, collapse = ","))
  value
}

# ---------- 只读标准临床宽表；不在 R 中修改原始 SAS ----------
read_input <- function(job) {
  separator <- job$sep %or% if (grepl("\\.tsv(\\.gz)?$", job$input_csv)) "\t" else ","
  result <- utils::read.table(job$input_csv, header = TRUE, sep = separator,
                             check.names = FALSE, stringsAsFactors = FALSE,
                             na.strings = c("", "NA"), quote = '"', comment.char = "",
                             colClasses = c(subject_id = "character"))
  # 临床 NA 编码不能吞掉字面被试 ID "NA"；单独按字符串恢复 ID。
  if ("subject_id" %in% names(result)) {
    classes <- rep("NULL", ncol(result))
    classes[match("subject_id", names(result))] <- "character"
    identifiers <- utils::read.table(job$input_csv, header = TRUE, sep = separator,
      check.names = FALSE, stringsAsFactors = FALSE, na.strings = "", quote = '"',
      comment.char = "", colClasses = classes)
    result$subject_id <- identifiers$subject_id
  }
  if (!"subject_id" %in% names(result) || anyNA(result$subject_id) ||
      anyDuplicated(result$subject_id)) stop("临床宽表 subject_id 缺失或重复")
  invisible(lapply(setdiff(names(result), "subject_id"), safe_name))
  result
}

feature_names <- function(job) {
  as_chars(job$feature_columns %or% c(sprintf("g_%03d", 1:129), sprintf("w_%03d", 1:129)))
}

mode_pairs <- function(job) {
  job$mode_pairs %or% lapply(seq_len(129), function(mode) {
    list(mode = mode, g = sprintf("g_%03d", mode), w = sprintf("w_%03d", mode))
  })
}

mi_settings <- function(job) {
  defaults <- list(m = 30L, initial_iterations = 20L, iteration_increment = 20L,
                   maximum_iterations = 100L, seed = 72631L, psrf_threshold = 1.1)
  utils::modifyList(defaults, job$config$statistics$mi %or% job$config$mi %or% list())
}

statistics_settings <- function(job) {
  settings <- utils::modifyList(list(q = 0.05, po_alpha = 0.05, minimum_category_count = 5L,
                         information_condition_limit = 1e12, n_agq = 5L,
                         threshold_sensitivity = TRUE, gradient_tolerance = 1e-4,
                         time_interaction = FALSE, ceiling_fraction_warning = .15),
                         job$config$statistics %or% list())
  models <- job$config$statistics$models
  settings$q <- job$config$statistics$fdr$q %or% settings$q
  settings$po_alpha <- models$ordinal$po_alpha %or% settings$po_alpha
  settings$minimum_category_count <- models$ordinal$sparse_level_minimum %or% settings$minimum_category_count
  settings$information_condition_limit <- models$ordinal$information_condition_limit %or% settings$information_condition_limit
  settings$gradient_tolerance <- models$ordinal$gradient_tolerance %or% settings$gradient_tolerance
  settings$threshold_sensitivity <- models$ordinal$threshold_sensitivity %or% settings$threshold_sensitivity
  settings$longitudinal_threshold_sensitivity <- models$longitudinal$threshold_sensitivity %or% settings$threshold_sensitivity
  settings$n_agq <- models$longitudinal$n_agq %or% settings$n_agq
  settings$longitudinal_optimizer_iterations <- models$longitudinal$optimizer_iterations %or% 300L
  settings$longitudinal_relative_tolerance <- models$longitudinal$optimizer_relative_tolerance %or% 1e-10
  settings$longitudinal_optimizer <- models$longitudinal$optimizer %or% "ucminf"
  settings$time_interaction <- models$longitudinal$time_interaction %or% settings$time_interaction
  settings$ceiling_fraction_warning <- models$linear$ceiling_fraction_warning %or% settings$ceiling_fraction_warning
  settings
}

prepare_types <- function(data, variables) {
  for (name in intersect(names(variables), names(data))) {
    definition <- variables[[name]]
    type <- definition$type %or% "continuous"
    if (type %in% c("binary", "nominal", "ordinal", "categorical")) {
      levels <- as_chars(definition$levels)
      if (!length(levels)) stop("分类变量缺少明确等级: ", name)
      invalid <- !is.na(data[[name]]) & !as.character(data[[name]]) %in% levels
      if (any(invalid)) stop("分类变量出现未声明等级: ", name)
      # 协变量一律采用分类虚拟编码；有序标记仅供选择插补方法。
      data[[name]] <- factor(as.character(data[[name]]), levels = levels, ordered = FALSE)
    } else {
      old <- data[[name]]
      data[[name]] <- suppressWarnings(as.numeric(old))
      if (any(is.na(data[[name]]) & !is.na(old))) stop("连续变量包含非数值: ", name)
      if (any(!is.finite(data[[name]]) & !is.na(data[[name]]))) stop("连续变量包含无限值: ", name)
    }
  }
  data
}

# 固定家族总检验数包含失败模型；失败 p 保留缺失，不被记为无效应。
family_adjust <- function(table, p_column = "p_joint") {
  table$q <- NA_real_
  table$planned_family_tests <- NA_integer_
  for (family in unique(table$family)) {
    selected <- which(table$family == family)
    observed <- selected[is.finite(table[[p_column]][selected])]
    table$planned_family_tests[selected] <- length(selected)
    if (length(observed)) table$q[observed] <- p.adjust(table[[p_column]][observed],
                                                       method = "BH", n = length(selected))
  }
  table
}
