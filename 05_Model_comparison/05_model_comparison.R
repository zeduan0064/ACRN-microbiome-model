###############################################################################
# Model comparison
#
# Includes:
# 1. PRS / ERS distributions
# 2. ROC / AUC comparison
# 3. Decision curve analysis (DCA)
# 4. Continuous NRI and IDI
# 5. PPV / NPV across assumed ACRN prevalence
#
# Models:
# - Microbiome
# - PRS-ERS
# - PRS-ERS & Microbiome
#
# Public version: paths and sample identifiers are anonymized.
###############################################################################

library(tidyverse)
library(pROC)
library(rmda)
library(patchwork)

set.seed(123)

# =============================================================================
# 1. Input / output
# =============================================================================

input_file <- file.path("data", "model_comparison_input.csv")
output_dir <- file.path("results", "model_comparison")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters retained from the original analysis
dca_bootstraps <- 50
nri_idi_bootstraps <- 100

# PPV / NPV parameters
sensitivity <- 0.783
specificity <- 0.847
prevalence_range <- seq(3, 15, by = 0.1)
prevalence_points <- c(3, 6, 9, 12, 15)

# =============================================================================
# 2. Read data
# =============================================================================

dat <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required <- c("Result", "PRS", "ERS")

if (!all(required %in% colnames(dat))) {
  stop("输入文件必须包含 Result、PRS 和 ERS。")
}

# Microbiome features = all columns except identifiers / outcome / PRS / ERS
micro_features <- setdiff(
  colnames(dat),
  c("SampleID", "sample-id", "Result", "PRS", "ERS")
)

if (length(micro_features) != 47) {
  stop(
    paste0(
      "当前识别到 ",
      length(micro_features),
      " 个菌属特征；本分析要求47个。"
    )
  )
}

analysis_vars <- c(
  "Result",
  "PRS",
  "ERS",
  micro_features
)

dat <- dat[
  complete.cases(dat[, analysis_vars]),
  ,
  drop = FALSE
]

# Outcome: 0 = HC, 1 = ACRN
if (all(as.character(dat$Result) %in% c("0", "1"))) {

  dat$Outcome <- as.numeric(as.character(dat$Result))

} else if (all(as.character(dat$Result) %in% c("HC", "ACRN"))) {

  dat$Outcome <- ifelse(
    dat$Result == "ACRN",
    1,
    0
  )

} else {

  stop("Result 必须为 0/1 或 HC/ACRN。")
}

dat$Group <- factor(
  ifelse(dat$Outcome == 1, "ACRN", "HC"),
  levels = c("ACRN", "HC")
)

dat$PRS <- as.numeric(dat$PRS)
dat$ERS <- as.numeric(dat$ERS)

# =============================================================================
# 3. PRS / ERS distributions
# =============================================================================

prs_test <- t.test(
  PRS ~ Group,
  data = dat
)

ers_test <- t.test(
  ERS ~ Group,
  data = dat
)

distribution_stats <- bind_rows(
  dat %>%
    group_by(Group) %>%
    summarise(
      Variable = "PRS",
      Mean = mean(PRS),
      SD = sd(PRS),
      .groups = "drop"
    ),
  dat %>%
    group_by(Group) %>%
    summarise(
      Variable = "ERS",
      Mean = mean(ERS),
      SD = sd(ERS),
      .groups = "drop"
    )
)

write.csv(
  distribution_stats,
  file.path(output_dir, "PRS_ERS_summary.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    Variable = c("PRS", "ERS"),
    P_value = c(
      prs_test$p.value,
      ers_test$p.value
    )
  ),
  file.path(output_dir, "PRS_ERS_t_test.csv"),
  row.names = FALSE
)

plot_score <- function(variable, p_value) {

  ggplot(
    dat,
    aes(
      x = Group,
      y = .data[[variable]]
    )
  ) +
    geom_jitter(
      aes(color = Group),
      width = 0.15,
      size = 1.8,
      alpha = 0.7
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      size = 3
    ) +
    stat_summary(
      fun.data = mean_sdl,
      fun.args = list(mult = 1),
      geom = "errorbar",
      width = 0.18
    ) +
    labs(
      title = variable,
      subtitle = paste0(
        "P = ",
        format.pval(p_value, digits = 3)
      ),
      x = NULL,
      y = "Score"
    ) +
    theme_classic() +
    theme(
      legend.position = "none"
    )
}

p_prs <- plot_score(
  "PRS",
  prs_test$p.value
)

p_ers <- plot_score(
  "ERS",
  ers_test$p.value
)

ggsave(
  file.path(
    output_dir,
    "PRS_ERS_distribution.pdf"
  ),
  p_prs + p_ers,
  width = 9,
  height = 4.5
)

# =============================================================================
# 4. Three logistic models
# =============================================================================

model_features <- list(
  Microbiome = micro_features,
  `PRS-ERS` = c("PRS", "ERS"),
  `PRS-ERS & Microbiome` = c(
    "PRS",
    "ERS",
    micro_features
  )
)

model_names <- names(model_features)

fit_list <- list()
prob_list <- list()
roc_list <- list()

for (model_name in model_names) {

  xvars <- model_features[[model_name]]

  formula_obj <- as.formula(
    paste(
      "Outcome ~",
      paste(
        paste0("`", xvars, "`"),
        collapse = " + "
      )
    )
  )

  fit <- glm(
    formula_obj,
    data = dat,
    family = binomial()
  )

  prob <- predict(
    fit,
    type = "response"
  )

  fit_list[[model_name]] <- fit
  prob_list[[model_name]] <- prob

  roc_list[[model_name]] <- roc(
    dat$Outcome,
    prob,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
}

# =============================================================================
# 5. AUC + 95% CI
# =============================================================================

auc_summary <- bind_rows(
  lapply(
    model_names,
    function(model_name) {

      r <- roc_list[[model_name]]
      ci <- ci.auc(r)

      data.frame(
        Model = model_name,
        AUC = as.numeric(auc(r)),
        Lower_95CI = as.numeric(ci[1]),
        Upper_95CI = as.numeric(ci[3])
      )
    }
  )
)

write.csv(
  auc_summary,
  file.path(
    output_dir,
    "AUC_3Models.csv"
  ),
  row.names = FALSE
)

auc_labels <- setNames(
  paste0(
    auc_summary$Model,
    " (AUC ",
    sprintf("%.2f", auc_summary$AUC),
    ", 95% CI ",
    sprintf("%.2f", auc_summary$Lower_95CI),
    "-",
    sprintf("%.2f", auc_summary$Upper_95CI),
    ")"
  ),
  auc_summary$Model
)

p_roc <- ggroc(
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
    labels = auc_labels
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
    "ROC_AUC_3Models.pdf"
  ),
  p_roc,
  width = 7,
  height = 7
)

# =============================================================================
# 6. Decision curve analysis
# =============================================================================

dca_data <- data.frame(
  Outcome = dat$Outcome
)

for (model_name in model_names) {
  dca_data[[model_name]] <- prob_list[[model_name]]
}

dca_list <- list()

for (model_name in model_names) {

  dca_formula <- as.formula(
    paste0(
      "Outcome ~ `",
      model_name,
      "`"
    )
  )

  dca_list[[model_name]] <- decision_curve(
    dca_formula,
    data = dca_data,
    study.design = "cohort",
    bootstraps = dca_bootstraps
  )
}

pdf(
  file.path(
    output_dir,
    "DCA_3Models.pdf"
  ),
  width = 7,
  height = 7
)

plot_decision_curve(
  dca_list,
  curve.names = model_names,
  confidence.intervals = "none",
  lwd = 2,
  legend.position = "topright"
)

dev.off()

# =============================================================================
# 7. Continuous NRI / IDI
# Reference model = Microbiome
# =============================================================================

calc_nri_idi <- function(
    p_ref,
    p_new,
    outcome,
    comparison,
    n_boot = 100
) {

  calc_once <- function(ref, new, y) {

    nri_event <-
      mean(new[y == 1] > ref[y == 1]) -
      mean(new[y == 1] < ref[y == 1])

    nri_nonevent <-
      mean(new[y == 0] < ref[y == 0]) -
      mean(new[y == 0] > ref[y == 0])

    nri <- nri_event + nri_nonevent

    idi <-
      (
        mean(new[y == 1]) -
        mean(ref[y == 1])
      ) -
      (
        mean(new[y == 0]) -
        mean(ref[y == 0])
      )

    c(
      NRI = nri,
      IDI = idi
    )
  }

  estimate <- calc_once(
    p_ref,
    p_new,
    outcome
  )

  set.seed(123)

  boot_result <- replicate(
    n_boot,
    {
      idx <- sample(
        seq_along(outcome),
        replace = TRUE
      )

      calc_once(
        p_ref[idx],
        p_new[idx],
        outcome[idx]
      )
    }
  )

  summarize_boot <- function(
      values,
      estimate_value
  ) {

    ci <- quantile(
      values,
      c(0.025, 0.975),
      na.rm = TRUE
    )

    se <- sd(
      values,
      na.rm = TRUE
    )

    p_value <- ifelse(
      se > 0,
      2 * (
        1 -
        pnorm(
          abs(
            estimate_value / se
          )
        )
      ),
      NA_real_
    )

    c(
      Estimate = estimate_value,
      Lower_95CI = ci[1],
      Upper_95CI = ci[2],
      P_value = p_value
    )
  }

  nri_result <- summarize_boot(
    boot_result["NRI", ],
    estimate["NRI"]
  )

  idi_result <- summarize_boot(
    boot_result["IDI", ],
    estimate["IDI"]
  )

  data.frame(
    Comparison = comparison,

    NRI = nri_result["Estimate"],
    NRI_Lower_95CI = nri_result["Lower_95CI"],
    NRI_Upper_95CI = nri_result["Upper_95CI"],
    NRI_P = nri_result["P_value"],

    IDI = idi_result["Estimate"],
    IDI_Lower_95CI = idi_result["Lower_95CI"],
    IDI_Upper_95CI = idi_result["Upper_95CI"],
    IDI_P = idi_result["P_value"]
  )
}

nri_idi_results <- bind_rows(
  calc_nri_idi(
    prob_list[["Microbiome"]],
    prob_list[["PRS-ERS"]],
    dat$Outcome,
    "PRS-ERS vs. Microbiome",
    nri_idi_bootstraps
  ),

  calc_nri_idi(
    prob_list[["Microbiome"]],
    prob_list[["PRS-ERS & Microbiome"]],
    dat$Outcome,
    "PRS-ERS & Microbiome vs. Microbiome",
    nri_idi_bootstraps
  )
)

write.csv(
  nri_idi_results,
  file.path(
    output_dir,
    "NRI_IDI_results.csv"
  ),
  row.names = FALSE
)

# =============================================================================
# 8. PPV / NPV across assumed ACRN prevalence
# =============================================================================

calc_ppv <- function(
    prevalence,
    sensitivity,
    specificity
) {

  (
    prevalence * sensitivity
  ) /
    (
      prevalence * sensitivity +
      (1 - prevalence) *
        (1 - specificity)
    )
}

calc_npv <- function(
    prevalence,
    sensitivity,
    specificity
) {

  (
    (1 - prevalence) *
      specificity
  ) /
    (
      (1 - prevalence) *
        specificity +
      prevalence *
        (1 - sensitivity)
    )
}

prev_decimal <- prevalence_range / 100

predictive_value_curve <- data.frame(
  Prevalence = prevalence_range,
  PPV = calc_ppv(
    prev_decimal,
    sensitivity,
    specificity
  ) * 100,
  NPV = calc_npv(
    prev_decimal,
    sensitivity,
    specificity
  ) * 100
)

predictive_value_points <- data.frame(
  Prevalence = prevalence_points
) %>%
  mutate(
    PPV = calc_ppv(
      Prevalence / 100,
      sensitivity,
      specificity
    ) * 100,
    NPV = calc_npv(
      Prevalence / 100,
      sensitivity,
      specificity
    ) * 100
  )

write.csv(
  predictive_value_points,
  file.path(
    output_dir,
    "PPV_NPV_selected_prevalence.csv"
  ),
  row.names = FALSE
)

pv_long <- predictive_value_curve %>%
  pivot_longer(
    cols = c("PPV", "NPV"),
    names_to = "Metric",
    values_to = "Value"
  )

p_pv <- ggplot(
  pv_long,
  aes(
    x = Prevalence,
    y = Value,
    color = Metric
  )
) +
  geom_line(
    linewidth = 1.2
  ) +
  geom_point(
    data = predictive_value_points %>%
      pivot_longer(
        cols = c("PPV", "NPV"),
        names_to = "Metric",
        values_to = "Value"
      ),
    size = 2
  ) +
  scale_x_continuous(
    breaks = prevalence_points
  ) +
  scale_y_continuous(
    limits = c(0, 100)
  ) +
  labs(
    x = "ACRN prevalence (%)",
    y = "Predictive value (%)",
    color = NULL
  ) +
  theme_classic()

ggsave(
  file.path(
    output_dir,
    "PPV_NPV_curve.pdf"
  ),
  p_pv,
  width = 7,
  height = 5
)

cat(
  "\n完成：PRS/ERS分布、三模型AUC、DCA、NRI/IDI和PPV/NPV分析。\n"
)