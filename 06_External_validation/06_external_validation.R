###############################################################################
# External validation in two independent cohorts
#
# Core analyses:
# 1. Data preparation
# 2. Z-score standardization
# 3. Multivariable logistic regression
# 4. ROC / AUC
# 5. Sensitivity / Specificity with bootstrap 95% CI
# 6. Confusion matrix
# 7. Bootstrap logistic feature importance
#
# Public version:
# - no real local paths
# - no real sample identifiers
# - no patient-level information
###############################################################################

library(readxl)
library(dplyr)
library(ggplot2)
library(pROC)
library(patchwork)

set.seed(123)

# =============================================================================
# 1. Input / output
# =============================================================================

cohort_files <- c(
  Cohort1 = file.path("data", "external_cohort_1.xlsx"),
  Cohort2 = file.path("data", "external_cohort_2.xlsx")
)

output_dir <- file.path("results", "external_validation")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

importance_boot <- 1000
metric_boot <- 1000

# =============================================================================
# 2. Data preparation
# =============================================================================

prepare_data <- function(file_path) {

  df <- read_excel(file_path)

  if (!"Result" %in% colnames(df)) {
    stop("External validation file must contain a Result column.")
  }

  # 0 = HC
  # 1 = AA
  # 2 = CRC
  # AA + CRC = ACRN
  result_raw <- trimws(
    as.character(df$Result)
  )

  df$Group <- case_when(
    result_raw %in%
      c("0", "HC", "No", "Control", "Healthy") ~ "HC",

    result_raw %in%
      c("1", "2", "AA", "CRC", "ACRN", "Yes") ~ "ACRN",

    TRUE ~ NA_character_
  )

  if (any(is.na(df$Group))) {
    stop("Unrecognized values were found in Result.")
  }

  df$Group <- factor(
    df$Group,
    levels = c("HC", "ACRN")
  )

  # External files contain only the genera available
  # from the predefined microbiome panel.
  features <- grep(
    "^g__",
    colnames(df),
    value = TRUE
  )

  if (length(features) < 2) {
    stop("Fewer than two genus features were detected.")
  }

  for (f in features) {

    df[[f]] <- suppressWarnings(
      as.numeric(df[[f]])
    )

    # Median imputation
    if (anyNA(df[[f]])) {

      med <- median(
        df[[f]],
        na.rm = TRUE
      )

      if (!is.finite(med)) {
        med <- 0
      }

      df[[f]][is.na(df[[f]])] <- med
    }
  }

  # Remove zero-variance features
  feature_sd <- sapply(
    df[, features, drop = FALSE],
    sd,
    na.rm = TRUE
  )

  features <- features[
    is.finite(feature_sd) &
      feature_sd > 0
  ]

  list(
    data = df,
    features = features
  )
}

# =============================================================================
# 3. Z-score standardization
# =============================================================================

zscore_matrix <- function(df, features) {

  X <- as.matrix(
    df[, features, drop = FALSE]
  )

  storage.mode(X) <- "double"

  mu <- colMeans(
    X,
    na.rm = TRUE
  )

  sigma <- apply(
    X,
    2,
    sd,
    na.rm = TRUE
  )

  sigma[
    !is.finite(sigma) |
      sigma == 0
  ] <- 1

  X <- sweep(
    X,
    2,
    mu,
    "-"
  )

  X <- sweep(
    X,
    2,
    sigma,
    "/"
  )

  colnames(X) <- features

  X
}

# =============================================================================
# 4. Logistic regression
# =============================================================================

fit_logistic <- function(X, y) {

  X_design <- cbind(
    "(Intercept)" = 1,
    X
  )

  suppressWarnings(
    glm.fit(
      x = X_design,
      y = y,
      family = binomial()
    )
  )
}

# =============================================================================
# 5. Bootstrap feature importance
#
# Importance = median(|standardized beta|)
# Each bootstrap resample is standardized again before model fitting.
# =============================================================================

bootstrap_importance <- function(
    df,
    features,
    n_boot = 1000,
    seed = 123
) {

  set.seed(seed)

  n <- nrow(df)

  beta_matrix <- matrix(
    NA_real_,
    nrow = n_boot,
    ncol = length(features),
    dimnames = list(NULL, features)
  )

  y <- ifelse(
    df$Group == "ACRN",
    1,
    0
  )

  for (b in seq_len(n_boot)) {

    idx <- sample(
      seq_len(n),
      size = n,
      replace = TRUE
    )

    y_b <- y[idx]

    if (length(unique(y_b)) < 2) {
      next
    }

    df_b <- df[
      idx,
      ,
      drop = FALSE
    ]

    X_b <- zscore_matrix(
      df_b,
      features
    )

    fit_b <- tryCatch(
      fit_logistic(
        X_b,
        y_b
      ),
      error = function(e) NULL
    )

    if (is.null(fit_b)) {
      next
    }

    beta_b <- fit_b$coefficients[-1]

    good <- is.finite(beta_b)

    beta_matrix[
      b,
      good
    ] <- beta_b[good]
  }

  importance <- data.frame(
    Genus = features,

    Median_Importance =
      apply(
        abs(beta_matrix),
        2,
        median,
        na.rm = TRUE
      ),

    Lower_95CI =
      apply(
        abs(beta_matrix),
        2,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),

    Upper_95CI =
      apply(
        abs(beta_matrix),
        2,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      ),

    Median_Signed_Beta =
      apply(
        beta_matrix,
        2,
        median,
        na.rm = TRUE
      ),

    stringsAsFactors = FALSE
  )

  importance %>%
    filter(
      is.finite(Median_Importance)
    ) %>%
    arrange(
      desc(Median_Importance)
    ) %>%
    mutate(
      Rank = row_number(),
      Genus_Label = gsub(
        "^g__",
        "",
        Genus
      )
    )
}

# =============================================================================
# 6. ROC / Sensitivity / Specificity
# =============================================================================

calculate_metrics <- function(
    actual,
    probability,
    n_boot = 1000,
    seed = 123
) {

  roc_obj <- roc(
    actual,
    probability,
    levels = c("No", "Yes"),
    direction = "<",
    quiet = TRUE
  )

  best_threshold <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = "threshold",
    transpose = FALSE
  )$threshold[1]

  calc_stats <- function(act, prob) {

    pred <- factor(
      ifelse(
        prob > best_threshold,
        "Yes",
        "No"
      ),
      levels = c("No", "Yes")
    )

    ref <- factor(
      act,
      levels = c("No", "Yes")
    )

    cm <- table(
      Prediction = pred,
      Reference = ref
    )

    tn <- cm["No", "No"]
    fp <- cm["Yes", "No"]
    fn <- cm["No", "Yes"]
    tp <- cm["Yes", "Yes"]

    c(
      Sensitivity =
        tp / (tp + fn),

      Specificity =
        tn / (tn + fp)
    )
  }

  point <- calc_stats(
    actual,
    probability
  )

  set.seed(seed)

  boot_stats <- replicate(
    n_boot,
    {

      idx <- sample(
        seq_along(actual),
        replace = TRUE
      )

      if (
        length(
          unique(
            actual[idx]
          )
        ) < 2
      ) {
        return(
          c(
            Sensitivity = NA_real_,
            Specificity = NA_real_
          )
        )
      }

      calc_stats(
        actual[idx],
        probability[idx]
      )
    }
  )

  ci <- apply(
    boot_stats,
    1,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )

  auc_ci <- ci.auc(
    roc_obj
  )

  list(
    ROC = roc_obj,
    Threshold = best_threshold,

    AUC =
      as.numeric(
        auc(roc_obj)
      ),

    AUC_Lower =
      as.numeric(
        auc_ci[1]
      ),

    AUC_Upper =
      as.numeric(
        auc_ci[3]
      ),

    Sensitivity =
      unname(
        point["Sensitivity"]
      ),

    Sensitivity_Lower =
      ci[1, "Sensitivity"],

    Sensitivity_Upper =
      ci[2, "Sensitivity"],

    Specificity =
      unname(
        point["Specificity"]
      ),

    Specificity_Lower =
      ci[1, "Specificity"],

    Specificity_Upper =
      ci[2, "Specificity"]
  )
}

# =============================================================================
# 7. Run one external cohort
# =============================================================================

run_external_cohort <- function(
    file_path,
    cohort_name,
    seed
) {

  prepared <- prepare_data(
    file_path
  )

  df <- prepared$data
  features <- prepared$features

  X <- zscore_matrix(
    df,
    features
  )

  y <- ifelse(
    df$Group == "ACRN",
    1,
    0
  )

  fit <- fit_logistic(
    X,
    y
  )

  probability <- fit$fitted.values

  actual <- factor(
    ifelse(
      y == 1,
      "Yes",
      "No"
    ),
    levels = c("No", "Yes")
  )

  metrics <- calculate_metrics(
    actual,
    probability,
    n_boot = metric_boot,
    seed = seed
  )

  importance <- bootstrap_importance(
    df,
    features,
    n_boot = importance_boot,
    seed = seed
  )

  prediction <- factor(
    ifelse(
      probability >
        metrics$Threshold,
      "Yes",
      "No"
    ),
    levels = c("No", "Yes")
  )

  cm <- as.data.frame(
    table(
      Predicted = prediction,
      Actual = actual
    )
  )

  cm$Cohort <- cohort_name

  list(
    cohort = cohort_name,
    n = nrow(df),
    features = features,
    metrics = metrics,
    importance = importance,
    confusion = cm
  )
}

# =============================================================================
# 8. Run two cohorts
# =============================================================================

results <- list()

results[["Cohort1"]] <- run_external_cohort(
  cohort_files["Cohort1"],
  "Cohort1",
  123
)

results[["Cohort2"]] <- run_external_cohort(
  cohort_files["Cohort2"],
  "Cohort2",
  456
)

# =============================================================================
# 9. Export results
# =============================================================================

performance <- bind_rows(
  lapply(
    results,
    function(x) {

      data.frame(
        Cohort = x$cohort,
        N = x$n,
        Number_of_Genera =
          length(
            x$features
          ),

        AUC =
          x$metrics$AUC,

        AUC_Lower_95CI =
          x$metrics$AUC_Lower,

        AUC_Upper_95CI =
          x$metrics$AUC_Upper,

        Sensitivity =
          x$metrics$Sensitivity,

        Sensitivity_Lower_95CI =
          x$metrics$Sensitivity_Lower,

        Sensitivity_Upper_95CI =
          x$metrics$Sensitivity_Upper,

        Specificity =
          x$metrics$Specificity,

        Specificity_Lower_95CI =
          x$metrics$Specificity_Lower,

        Specificity_Upper_95CI =
          x$metrics$Specificity_Upper,

        Threshold =
          x$metrics$Threshold
      )
    }
  )
)

write.csv(
  performance,
  file.path(
    output_dir,
    "External_validation_performance.csv"
  ),
  row.names = FALSE
)

for (cohort_name in names(results)) {

  write.csv(
    results[[cohort_name]]$importance,
    file.path(
      output_dir,
      paste0(
        cohort_name,
        "_bootstrap_importance.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    results[[cohort_name]]$confusion,
    file.path(
      output_dir,
      paste0(
        cohort_name,
        "_confusion_matrix.csv"
      )
    ),
    row.names = FALSE
  )
}

# =============================================================================
# 10. Simple six-panel figure
# =============================================================================

make_importance_plot <- function(x) {

  d <- x$importance %>%
    arrange(
      Median_Importance
    )

  d$Genus_Label <- factor(
    d$Genus_Label,
    levels = d$Genus_Label
  )

  ggplot(
    d,
    aes(
      x = Median_Importance,
      y = Genus_Label
    )
  ) +
    geom_col() +
    geom_errorbarh(
      aes(
        xmin = Lower_95CI,
        xmax = Upper_95CI
      ),
      height = 0.2
    ) +
    labs(
      x = "Median |standardized beta|",
      y = NULL
    ) +
    theme_classic()
}

make_roc_plot <- function(x) {

  r <- x$metrics$ROC

  roc_df <- data.frame(
    FPR =
      1 - r$specificities,
    Sensitivity =
      r$sensitivities
  )

  ggplot(
    roc_df,
    aes(
      x = FPR,
      y = Sensitivity
    )
  ) +
    geom_line(
      linewidth = 1.2
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +
    labs(
      subtitle = paste0(
        "AUC = ",
        sprintf(
          "%.2f",
          x$metrics$AUC
        ),
        " (95% CI ",
        sprintf(
          "%.2f",
          x$metrics$AUC_Lower
        ),
        "-",
        sprintf(
          "%.2f",
          x$metrics$AUC_Upper
        ),
        ")"
      ),
      x = "1 - Specificity",
      y = "Sensitivity"
    ) +
    theme_bw() +
    coord_equal()
}

make_cm_plot <- function(x) {

  d <- x$confusion

  ggplot(
    d,
    aes(
      x = Predicted,
      y = Actual,
      fill = Freq
    )
  ) +
    geom_tile() +
    geom_text(
      aes(label = Freq),
      size = 6
    ) +
    labs(
      subtitle = paste0(
        "Sensitivity = ",
        sprintf(
          "%.2f",
          x$metrics$Sensitivity
        ),
        "\nSpecificity = ",
        sprintf(
          "%.2f",
          x$metrics$Specificity
        )
      ),
      x = "Predicted",
      y = "Actual"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none"
    )
}

plot1 <- make_importance_plot(
  results[["Cohort1"]]
)

plot2 <- make_roc_plot(
  results[["Cohort1"]]
)

plot3 <- make_cm_plot(
  results[["Cohort1"]]
)

plot4 <- make_importance_plot(
  results[["Cohort2"]]
)

plot5 <- make_roc_plot(
  results[["Cohort2"]]
)

plot6 <- make_cm_plot(
  results[["Cohort2"]]
)

final_plot <-
  (
    plot1 |
      plot2 |
      plot3
  ) /
  (
    plot4 |
      plot5 |
      plot6
  )

ggsave(
  file.path(
    output_dir,
    "External_validation_summary.pdf"
  ),
  final_plot,
  width = 18,
  height = 10
)

cat(
  "\nExternal validation completed.\n"
)