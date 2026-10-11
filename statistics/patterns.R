# 功能：协变量调整后的有序弹性网组合，训练折隔离预处理和插补。

pattern_settings <- function(job) {
  utils::modifyList(list(folds = 5L, alpha = 0.5, lambda_min = 1e-4,
    lambda_max = 1, n_lambda = 60L, seed = 72631L, selection = "one_se",
    maxiter = 1000L, tolerance = 1e-6), job$config$patterns %or% list())
}

# ---------- 训练折内的标准化；验证折不参与均值、尺度或常量判断 ----------
image_scaler <- function(data, features, train) {
  values <- as.matrix(data[train, features, drop = FALSE])
  means <- colMeans(values)
  scales <- apply(values, 2, stats::sd)
  retained <- is.finite(scales) & scales > 1e-12
  scales[!retained] <- 1
  list(mean = means, scale = scales, retained = retained)
}

pattern_matrix <- function(data, features, covariates, scaler, train) {
  image <- sweep(sweep(as.matrix(data[, features, drop = FALSE]), 2, scaler$mean, "-"),
                  2, scaler$scale, "/")[, scaler$retained, drop = FALSE]
  formula <- stats::reformulate(covariates)
  covariate <- stats::model.matrix(formula, data = data)
  covariate <- covariate[, colnames(covariate) != "(Intercept)", drop = FALSE]
  # 未惩罚协变量使用训练行均值尺度，改善年龄/体积与虚拟变量的数值稳定性。
  if (ncol(covariate)) {
    scales <- apply(covariate[train, , drop = FALSE], 2, stats::sd)
    keep <- is.finite(scales) & scales > 1e-12
    covariate <- covariate[, keep, drop = FALSE]
    scales <- scales[keep]
    means <- colMeans(covariate[train, , drop = FALSE])
    covariate <- sweep(sweep(covariate, 2, means, "-"), 2, scales, "/")
  }
  # 影像在前，协变量在后；临床项 penaltyFactors=0。
  list(matrix = cbind(image, covariate), penalty = c(rep(1, ncol(image)), rep(0, ncol(covariate))))
}

stratified_folds <- function(outcome, folds, seed) {
  if (folds < 2L) stop("交叉验证至少两折")
  if (any(table(outcome) < 2L)) stop("组合学习某一结局等级不足两例，不能稳定交叉验证")
  set.seed(seed)
  result <- integer(length(outcome))
  for (level in levels(outcome)) {
    indices <- which(outcome == level)
    result[indices] <- sample(rep(seq_len(folds), length.out = length(indices)))
  }
  if (any(tabulate(result, folds) == 0L)) stop("交叉验证出现空折")
  result
}

fit_net <- function(matrix, outcome, penalty, lambdas, settings) {
  if (length(unique(outcome)) < 2L || any(table(outcome) == 0)) stop("训练折缺少结局等级")
  ordinalNet::ordinalNet(matrix, outcome, family = "cumulative", link = "logit",
    reverse = TRUE, parallelTerms = TRUE, nonparallelTerms = FALSE,
    alpha = settings$alpha, standardize = FALSE, penaltyFactors = penalty,
    lambdaVals = lambdas, maxiterOut = as.integer(settings$maxiter),
    maxiterIn = as.integer(settings$maxiter), threshOut = settings$tolerance,
    threshIn = settings$tolerance, keepTrainingData = FALSE, printIter = FALSE)
}

net_converged <- function(fitted, settings) {
  is.finite(fitted$loglik) & is.finite(fitted$dif) &
    !(fitted$iterOut >= settings$maxiter & abs(fitted$dif) > settings$tolerance * 10)
}

net_image_coefficients <- function(fitted, features, scaler, index = 1L) {
  coefficient <- stats::coef(fitted, whichLambda = index, matrix = TRUE)
  retained <- features[scaler$retained]
  if (!all(retained %in% rownames(coefficient))) stop("ordinalNet 系数名称未与输入影像列一致")
  values <- coefficient[retained, , drop = FALSE]
  if (ncol(values) > 1L && any(abs(values - values[, 1]) > 1e-7))
    stop("组合不是比例优势的平行系数模型")
  result <- setNames(rep(0, length(features)), features)
  result[retained] <- values[, 1]
  result
}

# ---------- 每折检查点保存平均验证概率，MI 不作为独立样本 ----------
fit_pattern <- function(data, outcome, job, directory) {
  settings <- pattern_settings(job)
  if (!settings$selection %in% c("one_se", "one_standard_error", "minimum")) stop("未知的惩罚参数选择规则")
  features <- feature_names(job)
  covariates <- as_chars(outcome$covariates)
  selected <- endpoint_data(data, outcome, job)
  if (!identical(outcome$type, "ordinal")) stop("当前正式组合学习仅支持有序结局")
  selected[[outcome$name]] <- droplevels(selected[[outcome$name]])
  levels <- levels(selected[[outcome$name]])
  folds <- stratified_folds(selected[[outcome$name]], as.integer(settings$folds), as.integer(settings$seed))
  lambdas <- exp(seq(log(settings$lambda_max), log(settings$lambda_min), length.out = settings$n_lambda))
  signature <- object_hash(list(data = selected, outcome = outcome, settings = settings,
                                mi = mi_settings(job), variables = job$variables,
                                reference_id = job$reference_id))
  write_table(data.frame(subject_id = selected$subject_id, fold = folds), file.path(directory, "folds.csv"))
  losses <- matrix(NA_real_, settings$folds, length(lambdas))
  fold_mi_status <- rep(NA_character_, settings$folds)
  for (fold in seq_len(settings$folds)) {
    check_stop(job)
    fold_directory <- file.path(directory, paste0("fold_", fold))
    dir.create(fold_directory, recursive = TRUE, showWarnings = FALSE)
    checkpoint <- file.path(fold_directory, "result.rds")
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      if (!identical(saved$signature, signature)) stop("交叉验证检查点指纹变化")
      losses[fold, ] <- saved$loss
      fold_mi_status[fold] <- saved$mi_status
      next
    }
    train <- folds != fold
    validation <- !train
    scaler <- image_scaler(selected, features, train)
    imputed <- fold_imputations(selected, outcome, job, file.path(fold_directory, "imputation"), train)
    fold_mi_status[fold] <- imputed$status
    average <- array(0, dim = c(sum(validation), length(levels), length(lambdas)))
    converged <- rep(TRUE, length(lambdas))
    for (index in seq_len(imputed$m)) {
      check_stop(job)
      # 每次模型拟合是一个可恢复边界，避免整折重算已完成的插补模型。
      fit_checkpoint <- file.path(fold_directory, sprintf("imputation_%03d.rds", index))
      if (file.exists(fit_checkpoint)) {
        partial <- readRDS(fit_checkpoint)
        if (!identical(partial$signature, signature)) stop("弹性网拟合检查点指纹变化")
      } else {
        completed <- complete_covariates(selected, imputed, index)
        design <- pattern_matrix(completed, features, covariates, scaler, train)
        fitted <- fit_net(design$matrix[train, , drop = FALSE], selected[[outcome$name]][train],
                           design$penalty, lambdas, settings)
        probabilities <- array(NA_real_, dim = dim(average))
        valid <- net_converged(fitted, settings)
        for (lambda in seq_along(lambdas)) {
          probabilities[, , lambda] <- stats::predict(fitted,
            newx = design$matrix[validation, , drop = FALSE], whichLambda = lambda, type = "response")
        }
        partial <- list(signature = signature, probabilities = probabilities, converged = valid)
        write_rds(partial, fit_checkpoint)
      }
      average <- average + partial$probabilities / imputed$m
      converged <- converged & partial$converged
    }
    observed <- as.integer(selected[[outcome$name]][validation])
    for (lambda in seq_along(lambdas)) {
      probability <- average[, , lambda]
      if (converged[lambda] && all(is.finite(probability)) && all(probability >= 0) &&
          max(abs(rowSums(probability) - 1)) < 1e-6)
        losses[fold, lambda] <- -mean(log(pmax(probability[cbind(seq_along(observed), observed)], 1e-15)))
    }
    write_rds(list(signature = signature, loss = losses[fold, ], scaler = scaler,
                    converged = converged, mi_status = imputed$status), checkpoint)
    # 验证概率仅供拟合恢复使用；折结束后删除逐插补的重复概率矩阵。
    unlink(list.files(fold_directory, pattern = "^imputation_[0-9]+\\.rds$", full.names = TRUE))
    write_json(list(target = outcome$name, fold = fold, folds = settings$folds,
                     stage = "cross_validation"), file.path(job$output_dir, "progress.json"))
  }
  means <- colMeans(losses)
  standard_errors <- apply(losses, 2, stats::sd) / sqrt(settings$folds)
  if (!any(is.finite(means))) stop("全部惩罚参数存在未收敛折，无法选择组合")
  best <- which.min(replace(means, !is.finite(means), Inf))
  eligible <- which(is.finite(means) & means <= means[best] + standard_errors[best])
  chosen <- if (settings$selection %in% c("one_se", "one_standard_error")) eligible[which.max(lambdas[eligible])] else best
  write_table(data.frame(lambda = lambdas, mean_log_loss = means, standard_error = standard_errors,
                         selected = seq_along(lambdas) == chosen, converged_all_folds = is.finite(means)),
              file.path(directory, "cross_validation.csv"))
  train <- rep(TRUE, nrow(selected))
  scaler <- image_scaler(selected, features, train)
  imputed <- build_imputations(selected, outcome, job, file.path(directory, "full_imputation"), inference = FALSE)
  coefficients <- matrix(0, length(features), imputed$m, dimnames = list(features, NULL))
  for (index in seq_len(imputed$m)) {
    check_stop(job)
    checkpoint <- file.path(directory, sprintf("refit_%03d.rds", index))
    if (file.exists(checkpoint)) {
      saved <- readRDS(checkpoint)
      if (!identical(saved$signature, signature)) stop("全发现集重拟合检查点指纹变化")
      coefficients[, index] <- saved$coefficients
      next
    }
    completed <- complete_covariates(selected, imputed, index)
    design <- pattern_matrix(completed, features, covariates, scaler, train)
    fitted <- fit_net(design$matrix, selected[[outcome$name]], design$penalty, lambdas[chosen], settings)
    if (!all(net_converged(fitted, settings))) stop("所选惩罚参数在全发现集重拟合未收敛")
    coefficients[, index] <- net_image_coefficients(fitted, features, scaler)
    write_rds(list(signature = signature, coefficients = coefficients[, index]), checkpoint)
  }
  average <- rowMeans(coefficients)
  raw <- average / scaler$scale
  similarities <- vapply(seq_len(imputed$m), function(index) {
    denominator <- sqrt(sum(coefficients[, index]^2) * sum(average^2))
    if (denominator > 0) sum(coefficients[, index] * average) / denominator else NA_real_
  }, numeric(1))
  training_mi_passed <- all(c(fold_mi_status, imputed$status) %in%
                              c("diagnostic_passed", "no_missing_covariates"))
  direction_stable <- all(is.finite(similarities) & similarities > 0)
  pattern <- list(target = outcome$name, family = outcome$family, feature_columns = as.list(features),
    feature_mean = as.list(unname(scaler$mean)), feature_scale = as.list(unname(scaler$scale)),
    coefficients = as.list(unname(average)), raw_coefficients = as.list(unname(raw)),
    centering_offset = sum(raw * scaler$mean), lambda = lambdas[chosen], alpha = settings$alpha,
    lambda_boundary = chosen %in% c(1L, length(lambdas)), seed = settings$seed,
    covariates = as.list(covariates), levels = as.list(levels), mi_m = imputed$m,
    training_ids_sha = object_hash(selected$subject_id), reference_id = job$reference_id,
    training_signature = signature, outcome_direction = "higher_score_higher_severity",
    status = if (all(average == 0)) "no_image_pattern" else "frozen_pending_validation",
    mi_status = imputed$status, fold_mi_status = as.list(fold_mi_status),
    training_diagnostics_passed = training_mi_passed,
    coefficient_direction_stable = direction_stable,
    coefficient_cosine_to_ensemble = as.list(similarities),
    coefficients_are_mi_ensemble = TRUE)
  write_json(pattern, file.path(directory, "model.json"))
  write_rds(coefficients, file.path(directory, "imputation_coefficients.rds"))
  write_table(data.frame(feature = features, mean = scaler$mean, scale = scaler$scale,
    coefficient = average, raw_coefficient = raw,
    selection_frequency = rowMeans(abs(coefficients) > 1e-10),
    sign_agreement_with_ensemble = rowMeans(sign(coefficients) == sign(average)),
    coefficient_sd = apply(coefficients, 1, stats::sd)), file.path(directory, "coefficients.csv"))
  write_table(data.frame(subject_id = selected$subject_id,
    score = drop(as.matrix(selected[, features, drop = FALSE]) %*% raw) - pattern$centering_offset),
    file.path(directory, "discovery_scores.csv"))
  pattern
}

pattern_stage <- function(data, job) {
  targets <- as_chars(job$config$patterns$targets %or% vapply(job$outcomes, `[[`, character(1), "name"))
  summary <- list()
  for (target in targets) {
    check_stop(job)
    outcome <- Filter(function(value) identical(value$name, target), job$outcomes)
    directory <- file.path(job$output_dir, target)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    result <- tryCatch({
      if (length(outcome) != 1L) stop("组合目标缺少唯一结局定义: ", target)
      fit_pattern(data, outcome[[1]], job, directory)
    }, safe_stop = function(error) stop(error), error = function(error)
      list(target = target, status = "not_estimable", reason = conditionMessage(error)))
    write_json(result, file.path(directory, "status.json"))
    summary[[length(summary) + 1L]] <- data.frame(target = target, status = result$status,
      reason = result$reason %or% "", lambda = result$lambda %or% NA_real_,
      lambda_boundary = result$lambda_boundary %or% NA)
  }
  table <- do.call(rbind, summary)
  write_table(table, file.path(job$output_dir, "patterns.csv"))
  table
}

# ---------- 外部队列只计算已冻结的影像得分 ----------
frozen_score <- function(data, model, reference_id = NULL) {
  if (!is.null(reference_id) && !identical(reference_id, model$reference_id)) stop("外部队列健康参考版本不一致")
  features <- as_chars(model$feature_columns)
  if (!all(features %in% names(data))) stop("外部验证缺少已冻结的影像特征")
  coefficients <- as.numeric(unlist(model$raw_coefficients))
  if (length(coefficients) != length(features)) stop("冻结权重长度错误")
  drop(as.matrix(data[, features, drop = FALSE]) %*% coefficients) - model$centering_offset
}

validation_stage <- function(data, job) {
  targets <- as_chars(job$config$patterns$targets)
  primary <- job$config$patterns$primary_target %or% "mrs_3m"
  plans <- list()
  for (target in targets) {
    outcomes <- if (identical(target, primary) && !identical(job$config$validation$mrs_pattern_profile, FALSE)) job$outcomes else
      Filter(function(value) identical(value$name, target), job$outcomes)
    for (outcome in outcomes) plans[[length(plans) + 1L]] <- list(target = target, outcome = outcome)
  }
  # 预先保存全部检验，包括不存在/未复现的模式；校正分母不可事后缩水。
  write_table(do.call(rbind, lapply(plans, function(plan) data.frame(pattern = plan$target,
    outcome = plan$outcome$name, family = plan$outcome$family))), file.path(job$output_dir, "planned_tests.csv"))
  rows <- list()
  for (plan in plans) {
    check_stop(job)
    outcome <- plan$outcome
    target <- plan$target
    directory <- file.path(job$output_dir, target, outcome$name)
    model <- NULL
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    result <- tryCatch({
      model_path <- file.path(job$frozen_models_dir, target, "model.json")
      if (!file.exists(model_path)) stop("没有可用的冻结组合模型")
      model <- jsonlite::read_json(model_path, simplifyVector = FALSE)
      membership <- file.path(job$frozen_models_dir, target, "folds.csv")
      if (!file.exists(membership)) stop("冻结模型缺少发现集成员清单，无法核查队列独立性")
      training_ids <- utils::read.csv(membership, colClasses = c(subject_id = "character"))$subject_id
      if (length(intersect(training_ids, data$subject_id))) stop("发现集与外部验证集存在重叠 subject_id")
      if (model$status == "no_image_pattern") stop("发现阶段没有选中影像损伤组合")
      selected <- endpoint_data(data, outcome, job)
      selected$fixed_score <- frozen_score(selected, model, job$reference_id)
      if (stats::sd(selected$fixed_score) < 1e-12) stop("外部队列固定得分没有变异")
      checkpoint <- file.path(directory, "result.rds")
      signature <- object_hash(list(model = model, data = selected, outcome = outcome, config = job$config))
      if (file.exists(checkpoint)) {
        value <- readRDS(checkpoint)
        if (!identical(value$signature, signature)) stop("外部验证检查点指纹变化")
      } else {
        imputed <- build_imputations(selected, outcome, job, file.path(directory, "imputation"), inference = TRUE)
        value <- analyze_terms(selected, outcome, "fixed_score", imputed, job)
        value$signature <- signature
        write_rds(value, checkpoint)
      }
      value
    }, safe_stop = function(error) stop(error), error = function(error)
      list(status = "not_estimable", reason = conditionMessage(error)))
    coefficients <- result$pooled$coefficients
    rows[[length(rows) + 1L]] <- data.frame(pattern = target, outcome = outcome$name,
      family = outcome$family, status = result$status, reason = result$reason %or% "",
      n = result$n %or% NA_integer_, p_joint = result$pooled$p_joint %or% NA_real_,
      estimate = if (is.null(coefficients)) NA_real_ else coefficients$estimate,
      lower = if (is.null(coefficients)) NA_real_ else coefficients$lower,
      upper = if (is.null(coefficients)) NA_real_ else coefficients$upper)
    rows[[length(rows)]]$training_diagnostics_passed <- isTRUE(model$training_diagnostics_passed)
    rows[[length(rows)]]$coefficient_direction_stable <- isTRUE(model$coefficient_direction_stable)
  }
  table <- family_adjust(do.call(rbind, rows))
  write_table(table, file.path(job$output_dir, "validation.csv"))
  own <- table[table$pattern == table$outcome, , drop = FALSE]
  own$annotation_eligible <- own$status == "estimated" & is.finite(own$q) &
    own$q <= statistics_settings(job)$q & is.finite(own$estimate) & own$estimate > 0 &
    own$training_diagnostics_passed & own$coefficient_direction_stable
  own$annotation_eligible[is.na(own$annotation_eligible)] <- FALSE
  own$annotation_reason <- vapply(seq_len(nrow(own)), function(index) {
    if (own$annotation_eligible[index]) return("replicated_with_passed_training_and_validation_diagnostics")
    reasons <- character()
    if (!own$training_diagnostics_passed[index]) reasons <- c(reasons, "training_mi_warning_or_unverified")
    if (!own$coefficient_direction_stable[index]) reasons <- c(reasons, "mi_coefficient_direction_unstable_or_unverified")
    if (own$status[index] != "estimated") reasons <- c(reasons, "validation_model_diagnostic_warning")
    if (!is.finite(own$q[index]) || own$q[index] > statistics_settings(job)$q) reasons <- c(reasons, "not_replicated")
    if (!is.finite(own$estimate[index]) || own$estimate[index] <= 0) reasons <- c(reasons, "direction_not_replicated")
    paste(reasons, collapse = ";")
  }, character(1))
  write_table(own, file.path(job$output_dir, "pattern_eligibility.csv"))
  table
}
