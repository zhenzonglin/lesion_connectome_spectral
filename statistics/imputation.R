# 功能：只填补缺失协变量；保存紧凑填充值与可恢复的链检查点。

# ---------- 端点纳入名单：结局已观察且影像完整 ----------
endpoint_data <- function(data, outcome, job) {
  features <- feature_names(job)
  covariates <- as_chars(outcome$covariates)
  columns <- unique(c("subject_id", features, covariates,
                      if (!identical(outcome$type, "longitudinal")) outcome$name,
                      as_chars(outcome$visit_columns)))
  missing <- setdiff(columns, names(data))
  if (length(missing)) stop("临床宽表缺少列: ", paste(missing, collapse = ","))
  complete_images <- stats::complete.cases(data[, features, drop = FALSE]) &
    apply(data[, features, drop = FALSE], 1, function(values) all(is.finite(values)))
  if (identical(outcome$type, "longitudinal")) {
    observed <- rowSums(!is.na(data[, as_chars(outcome$visit_columns), drop = FALSE])) > 0
  } else observed <- !is.na(data[[outcome$name]])
  result <- data[complete_images & observed, columns, drop = FALSE]
  if (nrow(result) < 3L) stop("有效结局和影像样本不足: ", outcome$name)
  result <- prepare_types(result, job$variables %or% list())
  immutable <- unique(c("lesion_volume_mm3", "log_lesion_volume", "log_volume", as_chars(job$immutable_columns)))
  if (anyNA(result[, intersect(immutable, names(result)), drop = FALSE]))
    stop("病灶体积及其变换来自影像，不能插补")
  for (name in covariates) {
    if (all(is.na(result[[name]]))) stop("必要协变量整列缺失: ", name)
    if (length(unique(stats::na.omit(result[[name]]))) < 2L)
      stop("必要协变量在该结局样本中为常量: ", name)
  }
  if (!identical(outcome$type, "longitudinal") && identical(outcome$type, "ordinal")) {
    levels <- as_chars(outcome$levels)
    if (any(!as.character(result[[outcome$name]]) %in% levels)) stop("结局含未定义等级")
    result[[outcome$name]] <- ordered(result[[outcome$name]], levels = levels)
  }
  rownames(result) <- NULL
  result
}

# ---------- 构造 MI 输入：预测端和拟合端清楚分离 ----------
imputation_frame <- function(data, covariates, features, outcome, inference) {
  result <- data[, unique(c(covariates, features)), drop = FALSE]
  if (inference && identical(outcome$type, "longitudinal")) {
    for (name in as_chars(outcome$visit_columns)) {
      observed <- !is.na(data[[name]])
      # 已观察等级和未观察状态作为完整辅助变量，不填补任何访视结局。
      value <- ifelse(observed, as.character(data[[name]]), "not_observed")
      result[[paste0("aux_", name)]] <- factor(value,
        levels = c(as_chars(outcome$levels), "not_observed"))
      result[[paste0("observed_", name)]] <- as.integer(observed)
    }
  } else if (inference) result$observed_outcome <- data[[outcome$name]]
  result
}

mi_method <- function(definition, value, methods = list()) {
  type <- definition$type %or% if (is.factor(value)) "nominal" else "continuous"
  method <- methods[[type]] %or% switch(type, continuous = "pmm", binary = "logreg", nominal = "polyreg",
         categorical = "polyreg", ordinal = "polr", stop("不支持的插补类型: ", type))
  if (!method %in% c("pmm", "logreg", "polyreg", "polr")) stop("不允许的插补方法: ", method)
  method
}

# Vehtari 等（2021）的秩归一化、折叠及分链 PSRF；不依赖额外的 Stan 运行库。
rank_split_psrf <- function(chains) {
  if (any(!is.finite(chains)) || nrow(chains) < 4L || ncol(chains) < 2L) return(NA_real_)
  half <- floor(nrow(chains) / 2L)
  split <- cbind(chains[seq_len(half), , drop = FALSE],
                  chains[nrow(chains) - rev(seq_len(half)) + 1L, , drop = FALSE])
  calculate <- function(values) {
    normalized <- stats::qnorm((rank(as.vector(values), ties.method = "average") - .375) /
                                 (length(values) + .25))
    normalized <- matrix(normalized, nrow(values), ncol(values))
    within <- mean(apply(normalized, 2, stats::var))
    between <- nrow(values) * stats::var(colMeans(normalized))
    if (within <= 0) return(NA_real_)
    sqrt(((nrow(values) - 1) / nrow(values) * within + between / nrow(values)) / within)
  }
  max(calculate(split), calculate(abs(split - stats::median(split))))
}

mi_diagnostics <- function(imputed, covariates, threshold) {
  records <- list()
  for (name in covariates) {
    for (parameter in c("mean", "variance")) {
      array <- if (parameter == "mean") imputed$chainMean else imputed$chainVar
      for (iteration in seq_len(imputed$iteration)) {
        chains <- matrix(array[name, seq_len(iteration), ], nrow = iteration, ncol = imputed$m)
        records[[length(records) + 1L]] <- data.frame(.it = iteration, vrb = name,
          parameter = parameter, psrf = rank_split_psrf(chains))
      }
    }
  }
  diagnostics <- do.call(rbind, records)
  current <- diagnostics[diagnostics$.it == max(diagnostics$.it), , drop = FALSE]
  finite <- current$psrf[is.finite(current$psrf)]
  maximum <- if (length(finite)) max(finite) else NA_real_
  # 未定义诊断不视为通过；保留可核查的状态与链轨迹。
  status <- if (!nrow(current) || any(!is.finite(current$psrf))) "diagnostic_incomplete" else
    if (maximum > threshold) "psrf_warning" else "diagnostic_passed"
  list(status = status, maximum_psrf = maximum, table = diagnostics)
}

build_imputations <- function(data, outcome, job, directory, inference = TRUE,
                              ignored = rep(FALSE, nrow(data)), validation_only = FALSE) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  settings <- mi_settings(job)
  count <- as.integer(settings$m)
  if (count < 2L && !validation_only) stop("正式多重插补要求 m >= 2")
  initial <- as.integer(settings$maxit %or% settings$initial_iterations)
  increment <- as.integer(settings$increment %or% settings$iteration_increment)
  maximum <- as.integer(settings$max_iterations %or% settings$maximum_iterations)
  if (initial < 2L || increment < 1L || maximum < initial) stop("插补迭代参数无效")
  covariates <- as_chars(outcome$covariates)
  features <- feature_names(job)
  frame <- imputation_frame(data, covariates, features, outcome, inference)
  signature <- object_hash(list(frame = frame, ids = data$subject_id, ignored = ignored,
                                inference = inference, validation_only = validation_only, settings = settings,
                                variables = job$variables))
  final_path <- file.path(directory, "imputations.rds")
  if (file.exists(final_path)) {
    cached <- readRDS(final_path)
    if (!identical(cached$signature, signature)) stop("插补检查点指纹不一致，必须新建运行")
    return(cached)
  }
  missing <- is.na(frame)
  where <- matrix(FALSE, nrow(frame), ncol(frame), dimnames = dimnames(missing))
  where[, covariates] <- missing[, covariates, drop = FALSE]
  needs <- covariates[colSums(where[, covariates, drop = FALSE]) > 0]
  methods <- setNames(rep("", ncol(frame)), names(frame))
  for (name in needs) methods[name] <- mi_method(job$variables[[name]] %or% list(), frame[[name]], settings$methods %or% list())
  predictor <- matrix(0L, ncol(frame), ncol(frame), dimnames = list(names(frame), names(frame)))
  predictor[needs, ] <- 1L
  diag(predictor) <- 0L
  # 训练折中的常量不能成为预测变量；不查看验证折分布来决定列选择。
  constants <- vapply(frame, function(column) length(unique(stats::na.omit(column[!ignored]))) < 2L, logical(1))
  predictor[, constants] <- 0L
  if (any(constants[needs])) stop("训练行中待插补协变量无有效变异: ", paste(names(constants)[constants & names(constants) %in% needs], collapse = ","))
  write_table(data.frame(variable = names(frame), method = methods, constant_in_training = constants),
              file.path(directory, "methods.csv"))
  write_table(data.frame(variable = rownames(predictor), predictor, check.names = FALSE),
              file.path(directory, "predictor_matrix_requested.csv"))
  write_table(data.frame(subject_id = data$subject_id, ignored_for_estimation = ignored),
              file.path(directory, "participants.csv"))
  if (!length(needs)) {
    result <- list(signature = signature, m = count, fills = list(),
                   ids = data$subject_id, status = "no_missing_covariates", iterations = 0L)
    write_rds(result, final_path)
    write_json(result[c("signature", "m", "status", "iterations")], file.path(directory, "diagnostics.json"))
    return(result)
  }
  state_path <- file.path(directory, "mi_state.rds")
  if (file.exists(state_path)) {
    saved <- readRDS(state_path)
    if (!identical(saved$signature, signature)) stop("插补链指纹不一致")
    imputed <- saved$imputed
  } else {
    # 初始化同样只从训练行取值；防止默认随机初值借用验证行分布。
    set.seed(as.integer(settings$seed))
    initialized <- frame
    for (name in needs) {
      donors <- frame[[name]][!ignored & !is.na(frame[[name]])]
      if (!length(donors)) stop("训练行没有可用插补供者: ", name)
      initialized[[name]][is.na(initialized[[name]])] <- sample(donors, sum(is.na(initialized[[name]])), replace = TRUE)
    }
    imputed <- mice::mice(frame, m = count, maxit = initial, method = methods,
                          predictorMatrix = predictor, where = where, ignore = ignored,
                          data.init = initialized, seed = as.integer(settings$seed),
                          remove.constant = FALSE, remove.collinear = FALSE,
                          printFlag = FALSE)
    write_rds(list(signature = signature, imputed = imputed), state_path)
  }
  repeat {
    diagnostic <- mi_diagnostics(imputed, needs, settings$psrf_threshold)
    # 只有训练行的链可使用诊断提前停止；验证条件抽样固定迭代上限。
    if ((!validation_only && diagnostic$status == "diagnostic_passed") || imputed$iteration >= maximum) break
    check_stop(job)
    iterations <- min(increment, maximum - imputed$iteration)
    imputed <- mice::mice.mids(imputed, maxit = iterations, printFlag = FALSE)
    write_rds(list(signature = signature, imputed = imputed), state_path)
  }
  # 仅储存缺失协变量的 m 份填充值，不复制 G/W 或已观察结局。
  fills <- lapply(needs, function(name) {
    positions <- which(is.na(frame[[name]]))
    values <- lapply(seq_len(count), function(index) {
      as.character(mice::complete(imputed, index)[[name]][positions])
    })
    list(variable = name, rows = positions, values = values)
  })
  result <- list(signature = signature, m = count, fills = fills,
                 ids = data$subject_id, status = diagnostic$status,
                 iterations = imputed$iteration)
  write_rds(result, final_path)
  write_rds(list(mean = imputed$chainMean[needs, , , drop = FALSE],
                 variance = imputed$chainVar[needs, , , drop = FALSE]),
            file.path(directory, "chain_traces.rds"))
  if (!is.null(diagnostic$table)) write_table(diagnostic$table, file.path(directory, "chain_diagnostics.csv"))
  if (!is.null(imputed$loggedEvents)) write_table(imputed$loggedEvents, file.path(directory, "logged_events.csv"))
  write_table(data.frame(variable = rownames(imputed$predictorMatrix), imputed$predictorMatrix,
                         check.names = FALSE), file.path(directory, "predictor_matrix_actual.csv"))
  write_json(list(status = diagnostic$status, maximum_psrf = diagnostic$maximum_psrf,
                   iterations = imputed$iteration, m = count, signature = signature,
                   inference = inference, only_covariates_imputed = TRUE,
                   diagnostic_method = "rank_normalized_split_folded_psrf_chain_mean_and_variance",
                   events = if (is.null(imputed$loggedEvents)) 0 else nrow(imputed$loggedEvents)),
              file.path(directory, "diagnostics.json"))
  unlink(state_path) # 最终填充值已原子提交，删除含原始宽表的临时链对象。
  check_stop(job)
  result
}

# ---------- 严格训练/验证隔离：先训练 MI，再对每套固定训练供者做条件填补 ----------
fold_imputations <- function(data, outcome, job, directory, train) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  signature <- object_hash(list(data = data, outcome = outcome, settings = mi_settings(job),
                                variables = job$variables, train = train))
  path <- file.path(directory, "fold_imputations.rds")
  if (file.exists(path)) {
    cached <- readRDS(path)
    if (!identical(cached$signature, signature)) stop("训练/验证插补检查点指纹变化")
    return(cached)
  }
  # 此调用完全没有验证侧的行、标签或数值。
  training_data <- data[train, , drop = FALSE]
  rownames(training_data) <- NULL
  training <- build_imputations(training_data, outcome, job, file.path(directory, "training"), FALSE)
  covariates <- as_chars(outcome$covariates)
  needs <- covariates[colSums(is.na(data[, covariates, drop = FALSE])) > 0]
  missing_rows <- setNames(lapply(needs, function(name) which(is.na(data[[name]]))), needs)
  all_values <- setNames(lapply(needs, function(name) vector("list", training$m)), needs)
  validation_traces <- list()
  for (index in seq_len(training$m)) {
    check_stop(job)
    # 每套已完成的训练协变量成为固定的拟合行；where 永远不能改动训练行。
    combined <- data
    filled_training <- complete_covariates(training_data, training, index)
    for (name in covariates) combined[[name]][train] <- filled_training[[name]]
    validation_job <- job
    validation_job$config$statistics$mi <- mi_settings(job)
    validation_job$config$statistics$mi$m <- 1L
    validation_job$config$statistics$mi$seed <- mi_settings(job)$seed + index
    draw_directory <- file.path(directory, sprintf("validation_draw_%03d", index))
    conditional <- build_imputations(combined, outcome, validation_job, draw_directory,
                                     inference = FALSE, ignored = !train, validation_only = TRUE)
    completed <- complete_covariates(combined, conditional, 1L)
    for (name in needs) all_values[[name]][[index]] <- as.character(completed[[name]][missing_rows[[name]]])
    trace_path <- file.path(draw_directory, "chain_traces.rds")
    if (file.exists(trace_path)) validation_traces[[index]] <- readRDS(trace_path)
  }
  validation_status <- "no_missing_covariates"
  if (length(validation_traces)) {
    first <- validation_traces[[1]]$mean
    dimensions <- c(dim(first)[1:2], training$m)
    combined_chains <- list(m = training$m, iteration = dim(first)[2],
      chainMean = array(unlist(lapply(validation_traces, `[[`, "mean")), dimensions,
                        dimnames = list(dimnames(first)[[1]], NULL, NULL)),
      chainVar = array(unlist(lapply(validation_traces, `[[`, "variance")), dimensions,
                       dimnames = list(dimnames(first)[[1]], NULL, NULL)))
    diagnostic <- mi_diagnostics(combined_chains, dimnames(first)[[1]], mi_settings(job)$psrf_threshold)
    validation_status <- diagnostic$status
    write_table(diagnostic$table, file.path(directory, "validation_chain_diagnostics.csv"))
  }
  clean <- c("diagnostic_passed", "no_missing_covariates")
  status <- if (!training$status %in% clean) "training_mi_warning" else
    if (!validation_status %in% clean) "validation_conditional_mi_warning" else "diagnostic_passed"
  if (!length(needs)) status <- "no_missing_covariates"
  result <- list(signature = signature, m = training$m, ids = data$subject_id,
    status = status, training_status = training$status, validation_status = validation_status,
    fills = lapply(needs, function(name) list(variable = name, rows = missing_rows[[name]], values = all_values[[name]])))
  write_rds(result, path)
  write_json(list(status = status, training_status = training$status, validation_status = validation_status,
    m = training$m, method = "training_only_mi_then_fixed_training_conditional_validation",
    training_rows = sum(train), validation_rows = sum(!train), signature = signature),
    file.path(directory, "diagnostics.json"))
  result
}

complete_covariates <- function(data, imputations, index) {
  if (!identical(data$subject_id, imputations$ids)) stop("插补病例顺序不一致")
  result <- data
  for (fill in imputations$fills) {
    value <- fill$values[[index]]
    if (is.factor(result[[fill$variable]])) {
      result[[fill$variable]][fill$rows] <- value
    } else result[[fill$variable]][fill$rows] <- as.numeric(value)
  }
  result
}
