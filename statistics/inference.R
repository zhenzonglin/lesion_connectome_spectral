# 功能：多重插补后的联合 D1 检验、Rubin 区间与模型诊断。

# ---------- 合并统计量；绝不平均 p 值、标准误或似然比 ----------
pool_parameters <- function(estimates, covariances, complete_df = Inf) {
  count <- ncol(estimates)
  parameters <- nrow(estimates)
  if (count < 2L || any(!is.finite(estimates)) || any(!is.finite(covariances)))
    stop("合并输入无效或有效插补少于两套")
  average <- rowMeans(estimates)
  within <- apply(covariances, c(1, 2), mean)
  between <- stats::cov(t(estimates))
  if (parameters == 1L) {
    within <- matrix(within, 1, 1)
    between <- matrix(between, 1, 1)
  }
  total <- within + (1 + 1 / count) * between
  if (any(eigen(within, symmetric = TRUE, only.values = TRUE)$values <= 0))
    stop("联合参数协方差非正定，不能合并")
  if (parameters == 1L) {
    statistic <- degrees <- p_joint <- NA_real_ # 后面采用精确的一维 Rubin t。
  } else if (max(abs(between)) < 1e-20) {
    statistic <- drop(t(average) %*% solve(within, average)) / parameters
    degrees <- complete_df
    p_joint <- stats::pf(statistic, parameters, degrees, lower.tail = FALSE)
  } else {
    pooled <- mitml::testConstraints(qhat = estimates, uhat = covariances,
      constraints = paste0("`", rownames(estimates), "`"), method = "D1",
      df.com = if (is.finite(complete_df) && parameters * (count - 1) > 4) complete_df else NULL)
    values <- unname(pooled$test)
    statistic <- values[1L]
    degrees <- values[3L]
    p_joint <- values[4L]
  }
  coefficients <- lapply(seq_len(parameters), function(index) {
    variance <- total[index, index]
    relative <- (1 + 1 / count) * between[index, index] / within[index, index]
    old_df <- if (relative < 1e-14) Inf else (count - 1) * (1 + 1 / relative)^2
    # Barnard-Rubin 小样本自由度；无穷完整数据自由度时退化为普通 Rubin。
    fraction <- (1 + 1 / count) * between[index, index] / variance
    observed_df <- if (is.finite(complete_df))
      (complete_df + 1) / (complete_df + 3) * complete_df * (1 - fraction) else Inf
    degrees_one <- 1 / (1 / old_df + 1 / observed_df)
    critical <- stats::qt(0.975, degrees_one)
    standard_error <- sqrt(variance)
    data.frame(term = rownames(estimates)[index], estimate = average[index],
      standard_error = standard_error, lower = average[index] - critical * standard_error,
      upper = average[index] + critical * standard_error, df = degrees_one,
      p = 2 * stats::pt(-abs(average[index] / standard_error), degrees_one))
  })
  coefficient_table <- do.call(rbind, coefficients)
  if (parameters == 1L) {
    # 单参数检验直接采用 Rubin t；平方后是同一 Wald F，避免小 m 下
    # mitml 的多参数小样本自由度公式在一维约束时返回 NaN。
    statistic <- (coefficient_table$estimate / coefficient_table$standard_error)^2
    degrees <- coefficient_table$df
    p_joint <- coefficient_table$p
  }
  if (!is.finite(p_joint)) stop("合并检验自由度或 p 值不可计算")
  list(p_joint = p_joint, statistic = statistic, df1 = parameters, df2 = degrees,
       coefficients = coefficient_table, total_covariance = total)
}

formula_for <- function(terms, longitudinal = FALSE, interactions = character()) {
  invisible(lapply(terms, safe_name))
  right <- paste(terms, collapse = " + ")
  if (!nzchar(right)) right <- "1"
  if (longitudinal) right <- paste(right, "+ time + (1 | subject_id)")
  if (length(interactions)) right <- paste(right, "+", paste(paste0(interactions, ":time"), collapse = " + "))
  stats::as.formula(paste("response ~", right), env = parent.frame())
}

make_response <- function(data, outcome) {
  if (identical(outcome$type, "longitudinal")) {
    columns <- as_chars(outcome$visit_columns)
    visits <- as_chars(outcome$visit_times %or% seq_along(columns))
    records <- lapply(seq_along(columns), function(index) {
      selected <- !is.na(data[[columns[index]]])
      visit <- data[selected, , drop = FALSE]
      visit$response <- ordered(as.character(visit[[columns[index]]]), levels = as_chars(outcome$levels))
      visit$time <- factor(visits[index], levels = visits)
      visit
    })
    result <- do.call(rbind, records)
    result$subject_id <- factor(result$subject_id)
    result$response <- droplevels(result$response)
  } else {
    result <- data
    result$response <- if (identical(outcome$type, "ordinal"))
      droplevels(ordered(as.character(data[[outcome$name]]), levels = as_chars(outcome$levels))) else
        as.numeric(data[[outcome$name]])
  }
  result
}

# ---------- 每套插补模型：数值失败和假设警报必须单列 ----------
fit_model <- function(data, outcome, tested, job, threshold = NULL, interaction = FALSE) {
  settings <- statistics_settings(job)
  terms <- c(tested, as_chars(outcome$covariates))
  longitudinal <- identical(outcome$type, "longitudinal")
  frame <- make_response(data, outcome)
  # 连续协变量只做等价的数值重参数化，改善年龄与微小损伤量纲并存时的条件数。
  # G/W 不缩放，因此输出仍是原始损伤单位；同一设计矩阵列空间保持不变。
  for (name in as_chars(outcome$covariates)) {
    if (is.numeric(frame[[name]])) {
      scale <- stats::sd(frame[[name]])
      if (is.finite(scale) && scale > 0) frame[[name]] <- (frame[[name]] - mean(frame[[name]])) / scale
    }
  }
  ordinal <- outcome$type %in% c("ordinal", "longitudinal")
  warnings <- character()
  if (!is.null(threshold)) {
    frame$response <- as.integer(as.integer(frame$response) > threshold)
    if (length(unique(frame$response)) < 2L) stop("阈值敏感性中只有一个结局等级")
  } else if (ordinal && nlevels(frame$response) < 2L) stop("有序结局只有一个等级")
  fitted <- withCallingHandlers({
    if (!is.null(threshold) && longitudinal) {
      lme4::glmer(formula_for(terms, TRUE), data = frame, family = stats::binomial(),
        nAGQ = as.integer(settings$n_agq), control = lme4::glmerControl(optimizer = "bobyqa"))
    } else if (!is.null(threshold)) {
      stats::glm(formula_for(terms), data = frame, family = stats::binomial())
    } else if (longitudinal) {
      optimizer_control <- if (settings$longitudinal_optimizer == "ucminf")
        ordinal::clmm.control(method = "ucminf", maxIter = settings$longitudinal_optimizer_iterations,
          gradTol = min(1e-6, settings$gradient_tolerance / 100), innerCtrl = "warnOnly",
          maxeval = settings$longitudinal_optimizer_iterations * 4L, grtol = 1e-7) else
        ordinal::clmm.control(method = "nlminb",
          maxIter = settings$longitudinal_optimizer_iterations,
          gradTol = min(1e-6, settings$gradient_tolerance / 100), innerCtrl = "warnOnly",
          iter.max = settings$longitudinal_optimizer_iterations,
          eval.max = settings$longitudinal_optimizer_iterations * 4L,
          rel.tol = settings$longitudinal_relative_tolerance)
      ordinal::clmm(formula_for(terms, TRUE, if (interaction) tested else character()), data = frame, link = "logit", Hess = TRUE,
        nAGQ = as.integer(settings$n_agq), control = optimizer_control)
    } else if (ordinal) {
      ordinal::clm(formula_for(terms), data = frame, link = "logit", Hess = TRUE)
    } else stats::lm(formula_for(terms), data = frame)
  }, warning = function(warning) {
    warnings <<- c(warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  })
  coefficients <- if (inherits(fitted, "merMod")) lme4::fixef(fitted) else stats::coef(fitted)
  if (interaction) {
    selected <- unlist(lapply(tested, function(term) grep(paste0("^", term, ":time"), names(coefficients), value = TRUE)))
    if (!length(selected)) stop("时间交互参数不可估计")
    tested <- selected
  }
  covariance <- if (!ordinal && is.null(threshold))
    sandwich::vcovHC(fitted, type = "HC3") else as.matrix(stats::vcov(fitted))
  if (!all(tested %in% names(coefficients)) || !all(tested %in% rownames(covariance)))
    stop("目标损伤参数不可估计（常量或完全共线）")
  if (any(!is.finite(coefficients)) || any(!is.finite(covariance))) stop("模型估计或协方差非有限值")
  eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  condition <- max(eigenvalues) / min(eigenvalues)
  if (min(eigenvalues) <= 0 || condition > settings$information_condition_limit)
    stop("信息矩阵非正定或严重病态")
  gradient_vector <- if (inherits(fitted, "merMod")) NULL else fitted$gradient
  raw_gradient <- if (is.null(gradient_vector)) NA_real_ else max(abs(gradient_vector))
  # 使用信息矩阵尺度上的梯度，避免绝对梯度随样本量、变量单位而改变。
  gradient <- if (length(gradient_vector) == nrow(covariance))
    max(abs(forwardsolve(t(chol(solve(covariance))), as.numeric(gradient_vector)))) else raw_gradient
  if (is.finite(gradient) && gradient > settings$gradient_tolerance)
    stop("有序模型未满足标准化梯度收敛标准: ", format(gradient, digits = 4))
  convergence <- if (inherits(fitted, "clm")) fitted$convergence$code else
    if (inherits(fitted, "clmm") && settings$longitudinal_optimizer == "nlminb") fitted$optRes$convergence else NULL
  if (length(convergence) && any(convergence != 0)) stop("优化器未收敛")
  if (inherits(fitted, "glm") && !isTRUE(fitted$converged)) stop("二分类敏感性模型未收敛")
  if (inherits(fitted, "merMod")) {
    messages <- fitted@optinfo$conv$lme4$messages
    if (length(messages)) stop("随机效应敏感性模型未收敛: ", paste(messages, collapse = ";"))
  }
  diagnostic <- list(condition = condition, gradient = gradient, raw_gradient = raw_gradient, warnings = unique(warnings),
                     po_status = "not_applicable", po_min_p = NA_real_)
  if (ordinal && is.null(threshold)) {
    diagnostic$category_counts <- as.list(table(frame$response))
    diagnostic$sparse_category <- any(table(frame$response) < settings$minimum_category_count)
    if (longitudinal) diagnostic$po_status <- "requires_random_effect_threshold_sensitivity" else {
      # nominal_test 会 update 原调用；保存已展开公式，不能让其重新执行 formula_for(terms)。
      fitted$call$formula <- stats::formula(fitted)
      # update 在 ordinal 的命名空间求值，不能在那里重新寻找本函数局部的 frame。
      fitted$call$data <- frame[, c("response", terms), drop = FALSE]
      nominal <- tryCatch(suppressWarnings(ordinal::nominal_test(fitted)), error = function(error) NULL)
      probabilities <- if (is.null(nominal)) numeric() else nominal[["Pr(>Chi)"]]
      probabilities <- probabilities[is.finite(probabilities)]
      diagnostic$po_min_p <- if (length(probabilities)) min(probabilities) else NA_real_
      diagnostic$po_status <- if (!length(probabilities)) {
        if (nlevels(frame$response) == 2L) "binary_no_po_constraint" else "test_unavailable"
      } else if (min(probabilities) < settings$po_alpha) "assumption_warning" else "no_warning_detected"
      # 每套插补的诊断 p 单独记录，不把它们平均成正式合并检验。
      diagnostic$po_tests <- if (is.null(nominal)) NULL else as.data.frame(nominal)
    }
  }
  if (!ordinal && is.null(threshold)) {
    residuals <- stats::residuals(fitted)
    fitted_values <- stats::fitted(fitted)
    variance_test <- stats::lm(I(residuals^2) ~ fitted_values)
    diagnostic$heteroscedasticity_p <- stats::pf(summary(variance_test)$fstatistic[1],
      summary(variance_test)$fstatistic[2], summary(variance_test)$fstatistic[3], lower.tail = FALSE)
    diagnostic$residual_quantiles <- as.list(stats::quantile(residuals, c(0, .05, .25, .5, .75, .95, 1)))
    diagnostic$ceiling_fraction <- if (is.null(outcome$maximum)) NA_real_ else mean(frame$response == outcome$maximum)
    diagnostic$ceiling_warning <- is.finite(diagnostic$ceiling_fraction) &&
      diagnostic$ceiling_fraction > settings$ceiling_fraction_warning
    diagnostic$covariance <- "HC3"
  }
  list(estimate = coefficients[tested], covariance = covariance[tested, tested, drop = FALSE],
       diagnostic = diagnostic, n = length(unique(as.character(frame$subject_id))),
       observations = nrow(frame), complete_df = if (!ordinal) stats::df.residual(fitted) else Inf)
}

pool_fits <- function(fits) {
  estimates <- do.call(cbind, lapply(fits, `[[`, "estimate"))
  covariances <- array(unlist(lapply(fits, `[[`, "covariance")),
    dim = c(nrow(estimates), nrow(estimates), length(fits)),
    dimnames = list(rownames(estimates), rownames(estimates), NULL))
  pool_parameters(estimates, covariances, min(vapply(fits, `[[`, numeric(1), "complete_df")))
}

analyze_terms <- function(data, outcome, tested, imputations, job) {
  fits <- lapply(seq_len(imputations$m), function(index) {
    fit_model(complete_covariates(data, imputations, index), outcome, tested, job)
  })
  pooled <- pool_fits(fits)
  diagnostics <- lapply(fits, `[[`, "diagnostic")
  po_warning <- any(vapply(diagnostics, function(value)
    value$po_status %in% c("assumption_warning", "test_unavailable"), logical(1)))
  sparse <- any(vapply(diagnostics, function(value) isTRUE(value$sparse_category), logical(1)))
  sensitivities <- list()
  run_sensitivity <- outcome$type %in% c("ordinal", "longitudinal") &&
    isTRUE(if (outcome$type == "longitudinal") statistics_settings(job)$longitudinal_threshold_sensitivity else
             statistics_settings(job)$threshold_sensitivity) &&
    (po_warning || outcome$type == "longitudinal")
  if (run_sensitivity) {
    levels <- levels(make_response(data, outcome)$response)
    for (threshold in seq_len(length(levels) - 1L)) {
      sensitivity <- tryCatch({
        threshold_fits <- lapply(seq_len(imputations$m), function(index)
          fit_model(complete_covariates(data, imputations, index), outcome, tested, job, threshold))
        result <- pool_fits(threshold_fits)
        list(threshold = levels[threshold], status = "estimated", pooled = result)
      }, error = function(error) list(threshold = levels[threshold], status = "not_estimable", reason = conditionMessage(error)))
      sensitivities[[length(sensitivities) + 1L]] <- sensitivity
    }
  }
  interaction_result <- NULL
  if (outcome$type == "longitudinal" && isTRUE(statistics_settings(job)$time_interaction)) {
    interaction_result <- tryCatch({
      interaction_fits <- lapply(seq_len(imputations$m), function(index)
        fit_model(complete_covariates(data, imputations, index), outcome, tested, job, interaction = TRUE))
      list(status = "estimated", pooled = pool_fits(interaction_fits))
    }, error = function(error) list(status = "not_estimable", reason = conditionMessage(error)))
  }
  list(status = if (po_warning || sparse) "assumption_limited" else
    if (!imputations$status %in% c("diagnostic_passed", "no_missing_covariates")) "mi_diagnostic_warning" else "estimated",
    pooled = pooled, diagnostics = diagnostics, sensitivities = sensitivities,
    interaction = interaction_result,
    n = fits[[1]]$n, observations = fits[[1]]$observations, mi_status = imputations$status)
}

# ---------- 按结局、按模态保存检查点 ----------
association_stage <- function(data, job) {
  rows <- list()
  sensitivity_rows <- list()
  interaction_rows <- list()
  diagnostic_rows <- list()
  pairs <- mode_pairs(job)
  analyses <- lapply(pairs, function(pair) list(pair = pair, analysis = "gw_primary"))
  if (isTRUE(job$config$statistics$jmax_sensitivity)) analyses <- c(analyses,
    lapply(pairs, function(pair) list(pair = pair, analysis = "jmax_sensitivity")))
  for (outcome in job$outcomes) {
    safe_name(outcome$name)
    directory <- file.path(job$output_dir, outcome$name)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    preparation <- tryCatch({
      selected <- endpoint_data(data, outcome, job)
      imputation_directory <- file.path(job$upstream_imputation_dir %or% file.path(job$output_dir, "imputation"), outcome$name)
      imputed <- build_imputations(selected, outcome, job, imputation_directory)
      list(data = selected, imputed = imputed)
    }, safe_stop = function(error) stop(error), error = function(error) list(error = conditionMessage(error)))
    for (specification in analyses) {
      check_stop(job)
      pair <- specification$pair
      analysis <- specification$analysis
      tested <- if (analysis == "gw_primary") c(pair$g, pair$w) else "jmax"
      analysis_data <- preparation$data
      if (!is.null(analysis_data) && analysis == "jmax_sensitivity")
        analysis_data$jmax <- pmax(analysis_data[[pair$g]], analysis_data[[pair$w]])
      checkpoint <- file.path(directory, sprintf(if (analysis == "gw_primary") "mode_%03d.rds" else
                                                   "jmax_%03d.rds", as.integer(pair$mode)))
      signature <- object_hash(list(outcome = outcome, pair = pair, config = job$config,
                                    input = preparation$data, imputation = preparation$imputed$signature))
      if (file.exists(checkpoint)) {
        result <- readRDS(checkpoint)
        if (!identical(signature, result$signature)) stop("模型检查点指纹变化: ", checkpoint)
      } else {
        result <- if (!is.null(preparation$error)) list(status = "not_estimable", reason = preparation$error) else
          tryCatch(analyze_terms(analysis_data, outcome, tested, preparation$imputed, job),
                   error = function(error) list(status = "not_estimable", reason = conditionMessage(error)))
        result$signature <- signature
        write_rds(result, checkpoint)
      }
      row <- data.frame(outcome = outcome$name, family = outcome$family, mode = pair$mode, analysis = analysis,
        status = result$status, reason = result$reason %or% "", n = result$n %or% NA_integer_,
        observations = result$observations %or% NA_integer_, p_joint = result$pooled$p_joint %or% NA_real_,
        mi_status = result$mi_status %or% "not_run")
      for (prefix in c("g", "w", "jmax")) {
        term <- if (prefix == "g") pair$g else if (prefix == "w") pair$w else "jmax"
        coefficients <- result$pooled$coefficients
        value <- if (is.null(coefficients)) NULL else coefficients[coefficients$term == term, , drop = FALSE]
        for (field in c("estimate", "standard_error", "lower", "upper", "p"))
          row[[paste0(prefix, "_", field)]] <- if (is.null(value) || !nrow(value)) NA_real_ else value[[field]]
      }
      rows[[length(rows) + 1L]] <- row
      for (index in seq_along(result$diagnostics)) {
        diagnostic <- result$diagnostics[[index]]
        diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(outcome = outcome$name,
          family = outcome$family, mode = pair$mode, analysis = analysis, imputation = index,
          po_status = diagnostic$po_status, po_min_p = diagnostic$po_min_p,
          sparse_category = diagnostic$sparse_category %or% FALSE,
          information_condition = diagnostic$condition, gradient = diagnostic$gradient,
          heteroscedasticity_p = diagnostic$heteroscedasticity_p %or% NA_real_,
          ceiling_fraction = diagnostic$ceiling_fraction %or% NA_real_,
          warnings = paste(diagnostic$warnings, collapse = ";"))
      }
      for (sensitivity in result$sensitivities) {
        terms <- sensitivity$pooled$coefficients
        if (is.null(terms)) terms <- data.frame(term = NA_character_, estimate = NA_real_,
          standard_error = NA_real_, lower = NA_real_, upper = NA_real_, df = NA_real_, p = NA_real_)
        sensitivity_rows[[length(sensitivity_rows) + 1L]] <- cbind(data.frame(outcome = outcome$name,
          family = outcome$family, mode = pair$mode, analysis = analysis, threshold = sensitivity$threshold,
          status = sensitivity$status, reason = sensitivity$reason %or% "",
          p_joint = sensitivity$pooled$p_joint %or% NA_real_), terms)
      }
      if (!is.null(result$interaction)) interaction_rows[[length(interaction_rows) + 1L]] <-
        data.frame(outcome = outcome$name, family = outcome$family, mode = pair$mode,
          analysis = analysis, status = result$interaction$status,
          reason = result$interaction$reason %or% "", p_joint = result$interaction$pooled$p_joint %or% NA_real_)
      write_json(list(outcome = outcome$name, mode = pair$mode, status = result$status,
                       completed = length(rows), planned = length(job$outcomes) * length(analyses)),
                 file.path(job$output_dir, "progress.json"))
    }
  }
  all_rows <- do.call(rbind, rows)
  result <- family_adjust(all_rows[all_rows$analysis == "gw_primary", , drop = FALSE])
  write_table(result, file.path(job$output_dir, "associations.csv"))
  if (isTRUE(job$config$statistics$jmax_sensitivity)) write_table(
    family_adjust(all_rows[all_rows$analysis == "jmax_sensitivity", , drop = FALSE]),
    file.path(job$output_dir, "jmax_sensitivity.csv"))
  if (length(diagnostic_rows)) write_table(do.call(rbind, diagnostic_rows), file.path(job$output_dir, "model_diagnostics.csv"))
  if (length(sensitivity_rows)) write_table(do.call(rbind, sensitivity_rows), file.path(job$output_dir, "threshold_sensitivity.csv"))
  if (length(interaction_rows)) {
    interactions <- do.call(rbind, interaction_rows)
    primary <- family_adjust(interactions[interactions$analysis == "gw_primary", , drop = FALSE])
    write_table(primary, file.path(job$output_dir, "time_interaction.csv"))
  }
  result
}
