###############################################################################
# Genus-pathway correlation and microbial co-occurrence network
#
# 公开版说明：
# - 删除真实样本编号、真实本地/服务器路径及其他敏感信息。
# - 分组信息统一从匿名化 metadata.tsv 读取。
# - 仅保留论文分析与可视化所需的核心步骤。
###############################################################################

library(tidyverse)
library(psych)
library(compositions)
library(igraph)
library(ggraph)

set.seed(123)

# =============================================================================
# 1. 文件与参数
# =============================================================================

genus_file   <- file.path("data", "genus_47_relative_abundance.csv")
pathway_file <- file.path("data", "pathway_abundance.csv")
metadata_file <- file.path("data", "metadata.tsv")

output_dir <- file.path("results", "correlation_network")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05
rho_cutoff <- 0.20
prevalence_cutoff <- 0.05
n_boot <- 1000

# =============================================================================
# 2. 读取并匹配数据
# =============================================================================

metadata <- read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

id_col <- intersect(c("SampleID", "sample-id"), colnames(metadata))[1]
if (is.na(id_col) || !"Group" %in% colnames(metadata)) {
  stop("metadata.tsv 需包含样本ID列（SampleID或sample-id）和 Group 列。")
}
metadata$SampleID <- metadata[[id_col]]

genus <- read.csv(
  genus_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

genus_id <- intersect(c("SampleID", "sample-id"), colnames(genus))[1]
if (is.na(genus_id)) {
  stop("菌属丰度文件需包含样本ID列。")
}
genus$SampleID <- genus[[genus_id]]

# 文件中除样本ID外的列均作为47个目标菌属
genera <- setdiff(colnames(genus), c("SampleID", "sample-id"))

genus_matrix <- as.matrix(
  genus[, genera, drop = FALSE]
)
storage.mode(genus_matrix) <- "numeric"
rownames(genus_matrix) <- genus$SampleID

# Pathway文件：行 = pathway，列 = sample
pathway_raw <- read.csv(
  pathway_file,
  row.names = 1,
  check.names = FALSE
)

pathway_matrix <- t(as.matrix(pathway_raw))
storage.mode(pathway_matrix) <- "numeric"

# Pathway CLR transformation
pathway_clr <- clr(pathway_matrix + 1e-6)

common_samples <- Reduce(
  intersect,
  list(
    rownames(genus_matrix),
    rownames(pathway_clr),
    metadata$SampleID
  )
)

if (length(common_samples) == 0) {
  stop("菌属、通路和metadata之间没有匹配到共同样本。")
}

genus_matrix <- genus_matrix[common_samples, , drop = FALSE]
pathway_clr <- pathway_clr[common_samples, , drop = FALSE]

group <- metadata$Group[
  match(common_samples, metadata$SampleID)
]

if (!all(c("HC", "ACRN") %in% group)) {
  stop("metadata中的Group需包含HC和ACRN。")
}

# =============================================================================
# 3. HC / ACRN 菌属 × Pathway Spearman相关性
# =============================================================================

run_group_correlation <- function(group_name) {

  idx <- group == group_name

  result <- psych::corr.test(
    pathway_clr[idx, , drop = FALSE],
    genus_matrix[idx, , drop = FALSE],
    method = "spearman",
    adjust = "fdr",
    ci = FALSE
  )

  write.csv(
    result$r,
    file.path(output_dir, paste0(group_name, "_genus_pathway_rho.csv"))
  )

  write.csv(
    result$p,
    file.path(output_dir, paste0(group_name, "_genus_pathway_FDR.csv"))
  )

  list(r = result$r, fdr = result$p)
}

cor_ACRN <- run_group_correlation("ACRN")
cor_HC   <- run_group_correlation("HC")

# 统一Pathway顺序
mean_cor <- (cor_ACRN$r + cor_HC$r) / 2

if (nrow(mean_cor) > 1) {
  pathway_order <- rownames(mean_cor)[
    hclust(dist(mean_cor))$order
  ]
} else {
  pathway_order <- rownames(mean_cor)
}

to_long <- function(res, group_name) {

  rho_long <- as.data.frame(res$r) %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(
      -Pathway,
      names_to = "Genus",
      values_to = "Rho"
    )

  fdr_long <- as.data.frame(res$fdr) %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(
      -Pathway,
      names_to = "Genus",
      values_to = "FDR"
    )

  left_join(
    rho_long,
    fdr_long,
    by = c("Pathway", "Genus")
  ) %>%
    mutate(
      Group = group_name,
      Significance = case_when(
        FDR < 0.001 ~ "***",
        FDR < 0.01  ~ "**",
        FDR < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
}

heatmap_df <- bind_rows(
  to_long(cor_ACRN, "ACRN"),
  to_long(cor_HC, "HC")
) %>%
  mutate(
    Group = factor(Group, levels = c("ACRN", "HC")),
    Pathway = factor(Pathway, levels = rev(pathway_order)),
    Genus = factor(
      gsub("^g__", "", Genus),
      levels = gsub("^g__", "", genera)
    )
  )

max_r <- max(abs(heatmap_df$Rho), na.rm = TRUE)

heatmap_plot <- ggplot(
  heatmap_df,
  aes(x = Genus, y = Pathway, fill = Rho)
) +
  geom_tile() +
  geom_text(aes(label = Significance), size = 2.5) +
  facet_grid(Group ~ .) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#2166AC",
    midpoint = 0,
    limits = c(-max_r, max_r),
    name = "Spearman r"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "italic"
    )
  )

ggsave(
  file.path(output_dir, "Genus_Pathway_Spearman_heatmap.pdf"),
  heatmap_plot,
  width = 13,
  height = 10
)

# =============================================================================
# 4. 47菌属共现网络
# Relative abundance -> Spearman -> BH-FDR -> 1000 bootstrap
# =============================================================================

mean_abundance <- colMeans(genus_matrix, na.rm = TRUE)
prevalence <- colMeans(genus_matrix > 0, na.rm = TRUE)

eligible_genera <- names(
  prevalence[prevalence >= prevalence_cutoff]
)

if (length(eligible_genera) < 2) {
  stop("满足prevalence阈值的菌属不足2个。")
}

pair_index <- combn(
  seq_along(eligible_genera),
  2
)

# 原始Spearman与P值
edge_list <- lapply(
  seq_len(ncol(pair_index)),
  function(i) {

    g1 <- eligible_genera[pair_index[1, i]]
    g2 <- eligible_genera[pair_index[2, i]]

    test <- suppressWarnings(
      cor.test(
        genus_matrix[, g1],
        genus_matrix[, g2],
        method = "spearman",
        exact = FALSE
      )
    )

    data.frame(
      From = g1,
      To = g2,
      Original_Rho = unname(test$estimate),
      P_value = test$p.value
    )
  }
)

edges <- bind_rows(edge_list)
edges$FDR <- p.adjust(edges$P_value, method = "BH")

# 1000次Bootstrap
eligible_matrix <- genus_matrix[
  ,
  eligible_genera,
  drop = FALSE
]

boot_rho <- matrix(
  NA_real_,
  nrow = n_boot,
  ncol = ncol(pair_index)
)

for (b in seq_len(n_boot)) {

  boot_id <- sample.int(
    nrow(eligible_matrix),
    replace = TRUE
  )

  boot_cor <- suppressWarnings(
    cor(
      eligible_matrix[boot_id, , drop = FALSE],
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )

  boot_rho[b, ] <- boot_cor[
    cbind(
      pair_index[1, ],
      pair_index[2, ]
    )
  ]
}

edges$Bootstrap_Median_Rho <- apply(
  boot_rho,
  2,
  median,
  na.rm = TRUE
)

edges$Bootstrap_CI_Low <- apply(
  boot_rho,
  2,
  quantile,
  probs = 0.025,
  na.rm = TRUE
)

edges$Bootstrap_CI_High <- apply(
  boot_rho,
  2,
  quantile,
  probs = 0.975,
  na.rm = TRUE
)

edges$Abs_Rho <- abs(edges$Bootstrap_Median_Rho)

write.csv(
  edges,
  file.path(output_dir, "Network_all_pairs.csv"),
  row.names = FALSE
)

# 最终网络边
final_edges <- edges %>%
  filter(
    FDR < fdr_cutoff,
    Abs_Rho >= rho_cutoff,
    Bootstrap_CI_Low > 0 |
      Bootstrap_CI_High < 0
  ) %>%
  mutate(
    Correlation = ifelse(
      Bootstrap_Median_Rho > 0,
      "Positive",
      "Negative"
    )
  )

write.csv(
  final_edges,
  file.path(output_dir, "Network_final_edges.csv"),
  row.names = FALSE
)

# =============================================================================
# 5. 网络属性与简单网络图
# =============================================================================

if (nrow(final_edges) > 0) {

  connected <- unique(
    c(final_edges$From, final_edges$To)
  )

  nodes <- data.frame(
    name = connected,
    Label = gsub("^g__", "", connected),
    Mean_Relative_Abundance = mean_abundance[connected],
    Prevalence = prevalence[connected],
    stringsAsFactors = FALSE
  )

  graph <- graph_from_data_frame(
    final_edges %>%
      transmute(
        from = From,
        to = To,
        Rho = Bootstrap_Median_Rho,
        Abs_Rho = Abs_Rho,
        Correlation = Correlation
      ),
    directed = FALSE,
    vertices = nodes
  )

  V(graph)$Degree <- degree(graph)
  E(graph)$weight <- E(graph)$Abs_Rho
  V(graph)$Strength <- strength(
    graph,
    weights = E(graph)$weight
  )

  node_stats <- data.frame(
    Genus = V(graph)$name,
    Mean_Relative_Abundance =
      V(graph)$Mean_Relative_Abundance,
    Prevalence = V(graph)$Prevalence,
    Degree = V(graph)$Degree,
    Strength = V(graph)$Strength
  ) %>%
    arrange(desc(Degree), desc(Strength))

  write.csv(
    node_stats,
    file.path(output_dir, "Network_node_statistics.csv"),
    row.names = FALSE
  )

  network_summary <- data.frame(
    Nodes = vcount(graph),
    Edges = ecount(graph),
    Density = edge_density(graph),
    Mean_Degree = mean(degree(graph)),
    Mean_Strength = mean(
      strength(
        graph,
        weights = E(graph)$weight
      )
    )
  )

  write.csv(
    network_summary,
    file.path(output_dir, "Network_summary.csv"),
    row.names = FALSE
  )

  network_plot <- ggraph(
    graph,
    layout = "fr"
  ) +
    geom_edge_link(
      aes(
        edge_colour = Correlation,
        edge_width = Abs_Rho
      ),
      alpha = 0.7
    ) +
    scale_edge_colour_manual(
      values = c(
        Positive = "#B2182B",
        Negative = "#2166AC"
      )
    ) +
    geom_node_point(
      aes(size = Mean_Relative_Abundance),
      shape = 21,
      fill = "white"
    ) +
    geom_node_text(
      aes(label = Label),
      repel = TRUE,
      size = 3
    ) +
    theme_void()

  ggsave(
    file.path(output_dir, "Microbial_cooccurrence_network.pdf"),
    network_plot,
    width = 10,
    height = 8
  )
}

cat(
  "\n完成：菌属-通路相关性分析与菌属共现网络分析。\n"
)