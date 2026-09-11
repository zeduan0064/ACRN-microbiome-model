###############################################################################
# Six machine-learning algorithms comparison
#
# 1. Logistic regression
# 2. Neural network
# 3. XGBoost
# 4. Gradient boosting machine
# 5. Support vector machine
# 6. AdaBoost
#
# 包含：
# - 1–47菌属逐步累加 AUC 比较
# - repeated 10-fold CV × 5
# - 最终47菌六模型 ROC/AUC 叠加图
# - 六模型混淆矩阵
# - Logistic vs 其他模型 DeLong test
###############################################################################

library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(patchwork)
library(gbm)
library(kernlab)
library(xgboost)
library(adabag)
library(nnet)

set.seed(13)

# =============================================================================
# 1. 文件
# =============================================================================

input_file <- file.path("data", "model_input_47genera.csv")
output_dir <- file.path("results", "machine_learning")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. 数据
# =============================================================================

dat <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"Result" %in% colnames(dat)) {
  stop("输入文件必须包含 Result 列。")
}

if (all(as.character(dat$Result) %in% c("0", "1"))) {
  dat$Result <- factor(
    dat$Result,
    levels = c(0, 1),
    labels = c("No", "Yes")
  )
} else {
  dat$Result <- factor(
    dat$Result,
    levels = c("No", "Yes")
  )
}

# 47个菌属按输入文件中的列顺序进行逐步累加
ordered_features <- setdiff(
  colnames(dat),
  c("SampleID", "sample-id", "Result")
)

if (length(ordered_features) != 47) {
  stop("模型输入文件必须包含47个菌属特征。")
}

analysis_data <- dat[
  complete.cases(dat[, c("Result", ordered_features)]),
  c("Result", ordered_features),
  drop = FALSE
]

for (x in ordered_features) {
  analysis_data[[x]] <- as.numeric(analysis_data[[x]])
}

if (any(!is.finite(as.matrix(analysis_data[, ordered_features])))) {
  stop("菌属数据中存在 NA、Inf 或 -Inf。")
}

Y <- analysis_data$Result

# =============================================================================
# 3. 六种模型
# =============================================================================

model_names <- c(
  "Logistic",
  "NeuralNetwork",
  "Xgboost",
  "GBM",
  "SVM",
  "Adaboost"
)

model_full_names <- c(
  Logistic = "Logistic regression",
  NeuralNetwork = "Neural network",
  Xgboost = "Extreme gradient boosting",
  GBM = "Gradient boosting machine",
  SVM = "Support vector machine",
  Adaboost = "Adaptive boosting"
)

# 相同CV折叠
set.seed(13)

cv_index <- createMultiFolds(
  Y,
  k = 10,
  times = 5
)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  index = cv_index,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  allowParallel = FALSE
)

# =============================================================================
# 4. 训练函数
# =============================================================================

train_one_model <- function(model_name, X) {

  if (model_name == "Logistic") {
    return(
      train(
        x = X, y = Y,
        method = "glm",
        family = binomial(),
        metric = "ROC",
        trControl = ctrl
      )
    )
  }

  if (model_name == "NeuralNetwork") {
    return(
      train(
        x = X, y = Y,
        method = "nnet",
        tuneGrid = expand.grid(
          size = 1,
          decay = 1.0
        ),
        trace = FALSE,
        MaxNWts = 10000,
        metric = "ROC",
        trControl = ctrl
      )
    )
  }

  if (model_name == "Xgboost") {
    return(
      train(
        x = X, y = Y,
        method = "xgbTree",
        tuneGrid = expand.grid(
          nrounds = 300,
          max_depth = 2,
          eta = 0.005,
          gamma = 2,
          colsample_bytree = 0.7,
          min_child_weight = 20,
          subsample = 0.7
        ),
        verbose = FALSE,
        metric = "ROC",
        trControl = ctrl
      )
    )
  }

  if (model_name == "GBM") {
    return(
      train(
        x = X, y = Y,
        method = "gbm",
        tuneGrid = expand.grid(
          n.trees = 90,
          interaction.depth = 2,
          shrinkage = 0.005,
          n.minobsinnode = 50
        ),
        verbose = FALSE,
        metric = "ROC",
        trControl = ctrl
      )
    )
  }

  if (model_name == "SVM") {
    return(
      train(
        x = X, y = Y,
        method = "svmRadial",
        tuneGrid = expand.grid(
          sigma = 0.001,
          C = 0.09
        ),
        metric = "ROC",
        trControl = ctrl
      )
    )
  }

  if (model_name == "Adaboost") {
    return(
      train(
        x = X, y = Y,
        method = "AdaBoost.M1",
        tuneGrid = expand.grid(
          mfinal = 2,
          maxdepth = 2,
          coeflearn = "Zhu"
        ),
        metric = "ROC",
        trControl = ctrl
      )
    )
  }
}

# =============================================================================
# 5. 1–47菌属逐步累加：六模型AUC
# =============================================================================

incremental_results <- data.frame()

for (step in seq_along(ordered_features)) {

  X_step <- analysis_data[
    ,
    ordered_features[1:step],
    drop = FALSE
  ]

  for (i in seq_along(model_names)) {

    model_name <- model_names[i]

    set.seed(13 + step * 100 + i)

    fit <- train_one_model(
      model_name,
      X_step
    )

    prob <- predict(
      fit,
      newdata = X_step,
      type = "prob"
    )[, "Yes"]

    roc_obj <- roc(
      Y,
      prob,
      levels = c("No", "Yes"),
      direction = "<",
      quiet = TRUE
    )

    auc_ci <- ci.auc(
      roc_obj,
      method = "delong"
    )

    cv_auc <- fit$resample$ROC
    cv_auc <- cv_auc[is.finite(cv_auc)]

    incremental_results <- rbind(
      incremental_results,
      data.frame(
        Step = step,
        Added_Genus = ordered_features[step],
        Model = model_name,
        Apparent_AUC = as.numeric(auc(roc_obj)),
        Lower_95CI = as.numeric(auc_ci[1]),
        Upper_95CI = as.numeric(auc_ci[3]),
        CV_Mean_AUC = mean(cv_auc),
        CV_SD_AUC = sd(cv_auc)
      )
    )
  }
}

write.csv(
  incremental_results,
  file.path(
    output_dir,
    "Incremental_AUC_6Models_AllResults.csv"
  ),
  row.names = FALSE
)

# AUC累加叠加图
p_incremental <- ggplot(
  incremental_results,
  aes(
    x = Step,
    y = Apparent_AUC,
    color = Model,
    group = Model
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 47)
  ) +
  scale_y_continuous(
    limits = c(0.50, 1.00),
    breaks = seq(0.50, 1.00, 0.05)
  ) +
  labs(
    x = "Number of genera",
    y = "Area under the ROC curve (AUC)",
    color = NULL
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "Incremental_AUC_6Models.pdf"
  ),
  p_incremental,
  width = 9,
  height = 6.5
)

# CV AUC累加叠加图
p_cv <- ggplot(
  incremental_results,
  aes(
    x = Step,
    y = CV_Mean_AUC,
    color = Model,
    group = Model
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5) +
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 47)
  ) +
  scale_y_continuous(
    limits = c(0.50, 1.00),
    breaks = seq(0.50, 1.00, 0.05)
  ) +
  labs(
    x = "Number of genera",
    y = "Cross-validated AUC",
    color = NULL
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "Incremental_AUC_6Models_CV.pdf"
  ),
  p_cv,
  width = 9,
  height = 6.5
)

# =============================================================================
# 6. 最终47菌：重新训练六模型
# =============================================================================

X_final <- analysis_data[
  ,
  ordered_features,
  drop = FALSE
]

final_fits <- list()
final_predictions <- data.frame(Result = Y)
final_performance <- data.frame()

for (i in seq_along(model_names)) {

  model_name <- model_names[i]

  set.seed(1000 + i)

  fit <- train_one_model(
    model_name,
    X_final
  )

  final_fits[[model_name]] <- fit

  prob <- predict(
    fit,
    newdata = X_final,
    type = "prob"
  )[, "Yes"]

  final_predictions[[model_name]] <- prob

  roc_obj <- roc(
    Y,
    prob,
    levels = c("No", "Yes"),
    direction = "<",
    quiet = TRUE
  )

  auc_ci <- ci.auc(
    roc_obj,
    method = "delong"
  )

  cv_auc <- fit$resample$ROC

  final_performance <- rbind(
    final_performance,
    data.frame(
      Model = model_name,
      AUC = as.numeric(auc(roc_obj)),
      AUC_CI_Low = as.numeric(auc_ci[1]),
      AUC_CI_High = as.numeric(auc_ci[3]),
      CV_Mean_AUC = mean(cv_auc, na.rm = TRUE),
      CV_SD_AUC = sd(cv_auc, na.rm = TRUE)
    )
  )
}

write.csv(
  final_predictions,
  file.path(
    output_dir,
    "Final_47genera_predictions.csv"
  ),
  row.names = FALSE
)

write.csv(
  final_performance,
  file.path(
    output_dir,
    "Final_47genera_6Models_performance.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 7. 最终47菌：六模型 ROC/AUC 叠加图
# =============================================================================

roc_list <- lapply(
  model_names,
  function(m) {
    roc(
      Y,
      final_predictions[[m]],
      levels = c("No", "Yes"),
      direction = "<",
      quiet = TRUE
    )
  }
)

names(roc_list) <- model_names

roc_labels <- sapply(
  model_names,
  function(m) {
    paste0(
      model_full_names[m],
      " (AUC=",
      sprintf("%.3f", auc(roc_list[[m]])),
      ")"
    )
  }
)

p_roc_overlay <- ggroc(
  roc_list,
  legacy.axes = TRUE,
  linewidth = 1.2
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  scale_color_discrete(
    breaks = model_names,
    labels = roc_labels
  ) +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    color = NULL
  ) +
  theme_bw()

ggsave(
  file.path(
    output_dir,
    "Final_47genera_ROC_AUC_overlay.pdf"
  ),
  p_roc_overlay,
  width = 8,
  height = 7
)

# =============================================================================
# 8. 混淆矩阵 + Sensitivity/Specificity/Precision/F1
# =============================================================================

metric_table <- data.frame()
cm_plots <- list()

for (i in seq_along(model_names)) {

  m <- model_names[i]

  roc_obj <- roc_list[[m]]

  best_threshold <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = "threshold",
    transpose = FALSE
  )$threshold[1]

  pred <- factor(
    ifelse(
      final_predictions[[m]] > best_threshold,
      "Yes",
      "No"
    ),
    levels = c("No", "Yes")
  )

  ref <- factor(
    Y,
    levels = c("No", "Yes")
  )

  cm <- confusionMatrix(
    pred,
    ref,
    positive = "Yes"
  )

  sensitivity <- unname(cm$byClass["Sensitivity"])
  specificity <- unname(cm$byClass["Specificity"])
  precision <- unname(cm$byClass["Pos Pred Value"])
  f1 <- unname(cm$byClass["F1"])

  metric_table <- rbind(
    metric_table,
    data.frame(
      Model = m,
      Threshold = best_threshold,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Precision = precision,
      F1 = f1
    )
  )

  cm_df <- as.data.frame(
    table(
      Predicted = pred,
      Actual = ref
    )
  )

  p_cm <- ggplot(
    cm_df,
    aes(
      x = Predicted,
      y = Actual,
      fill = Freq
    )
  ) +
    geom_tile(color = "white") +
    geom_text(
      aes(label = Freq),
      size = 6
    ) +
    labs(
      title = model_full_names[m],
      subtitle = paste0(
        "Sensitivity = ", sprintf("%.3f", sensitivity),
        "\nSpecificity = ", sprintf("%.3f", specificity),
        "\nPrecision = ", sprintf("%.3f", precision),
        "\nF1 = ", sprintf("%.3f", f1)
      ),
      x = "Predicted",
      y = "Actual"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )

  cm_plots[[m]] <- p_cm
}

write.csv(
  metric_table,
  file.path(
    output_dir,
    "Final_47genera_confusion_matrix_metrics.csv"
  ),
  row.names = FALSE
)

p_cm_all <- wrap_plots(
  cm_plots,
  ncol = 3
)

ggsave(
  file.path(
    output_dir,
    "Final_47genera_6Models_ConfusionMatrices.pdf"
  ),
  p_cm_all,
  width = 15,
  height = 10
)

# =============================================================================
# 9. Logistic vs 其他5模型 DeLong test
# =============================================================================

delong_results <- lapply(
  setdiff(model_names, "Logistic"),
  function(m) {

    test <- roc.test(
      roc_list[["Logistic"]],
      roc_list[[m]],
      method = "delong"
    )

    data.frame(
      Baseline = "Logistic",
      Comparison = m,
      AUC_Logistic =
        as.numeric(auc(roc_list[["Logistic"]])),
      AUC_Model =
        as.numeric(auc(roc_list[[m]])),
      P_value =
        as.numeric(test$p.value)
    )
  }
)

delong_results <- bind_rows(
  delong_results
)

write.csv(
  delong_results,
  file.path(
    output_dir,
    "Logistic_vs_other_models_DeLong.csv"
  ),
  row.names = FALSE
)

cat("\n完成：六模型比较、1-47菌累加AUC、ROC叠加及混淆矩阵。\n")