###############################################################################
# Iterative LASSO feature selection
#
# 输入文件：
#   data/genus_abundance.tsv
#
# 数据格式：
#   第1列：匿名样本编号
#   中间列：菌属丰度
#   最后1列：二分类标签（0/1）
#
# 公开版本已删除真实路径、真实样本编号和其他敏感信息。
###############################################################################

library(glmnet)
library(compositions)

set.seed(26)

# 参数
input_file <- file.path("data", "genus_abundance.tsv")
output_dir <- file.path("results", "feature_selection")

n_iterations <- 50
train_fraction <- 0.70
n_folds <- 10
top_fraction <- 0.20

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 读取数据
dat <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = TRUE
)

X <- as.matrix(dat[, 2:(ncol(dat) - 1)])
y <- as.numeric(dat[, ncol(dat)])

# CLR 转换
X_clr <- clr(X + 1e-6)
feature_names <- colnames(X)

# 保存每个特征的入选次数和系数
selection_count <- setNames(integer(length(feature_names)), feature_names)

coef_matrix <- matrix(
  0,
  nrow = n_iterations,
  ncol = length(feature_names),
  dimnames = list(NULL, feature_names)
)

# 50次迭代 LASSO
for (i in seq_len(n_iterations)) {

  train_id <- sample(
    seq_len(nrow(X_clr)),
    size = floor(train_fraction * nrow(X_clr))
  )

  X_train <- scale(X_clr[train_id, , drop = FALSE])

  cv_fit <- cv.glmnet(
    x = X_train,
    y = y[train_id],
    alpha = 1,
    family = "binomial",
    nfolds = n_folds
  )

  coef_i <- as.matrix(
    coef(cv_fit, s = "lambda.min")
  )[-1, 1]

  coef_matrix[i, names(coef_i)] <- coef_i

  selected <- names(coef_i)[coef_i != 0]

  if (length(selected) > 0) {
    selection_count[selected] <- selection_count[selected] + 1L
  }
}

# 汇总 LASSO 结果
summary_df <- data.frame(
  Feature = feature_names,
  LASSO_Count = as.integer(selection_count[feature_names]),
  LASSO_Frequency = as.integer(selection_count[feature_names]) / n_iterations,
  Mean_Coefficient = colMeans(coef_matrix),
  stringsAsFactors = FALSE
)

summary_df <- summary_df[
  order(-summary_df$LASSO_Count),
  ,
  drop = FALSE
]

# 仅保留至少被选中过一次的特征
selected_df <- summary_df[
  summary_df$LASSO_Count > 0,
  ,
  drop = FALSE
]

# 按 LASSO 入选频率排序后取前20%
n_top <- ceiling(nrow(selected_df) * top_fraction)

top20_df <- selected_df[
  seq_len(n_top),
  ,
  drop = FALSE
]

# 输出
write.csv(
  summary_df,
  file.path(output_dir, "LASSO_50x_summary.csv"),
  row.names = FALSE
)

write.csv(
  top20_df,
  file.path(output_dir, "LASSO_top20_percent.csv"),
  row.names = FALSE
)

cat(
  "\n完成：50次迭代 LASSO。\n",
  "候选特征按入选频率排序后保留前20%。\n",
  sep = ""
)