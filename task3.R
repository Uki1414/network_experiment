#!/usr/bin/env Rscript
# 課題3: 保存済みグラフの中心性・ランキング・ノードサイズによる可視化。
library(igraph)
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else "output/graphs.rds"
out <- if (length(args) >= 2) args[2] else "output/task3"
if (!file.exists(input)) stop("先に Rscript task1.R を実行してください。")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
all_graphs <- readRDS(input)
ids <- c("er_p005", "ws_p005_k5", "ba_m2", "sbm_separated")
stopifnot(all(ids %in% names(all_graphs)))
graphs <- all_graphs[ids]
model_names <- c("ER", "WS", "BA", "SBM")
fields <- c("degree_normalized", "closeness_normalized", "betweenness_normalized", "eigenvector")
metric_names <- c("Degree", "Closeness", "Betweenness", "Eigenvector")
# 連続値の計算誤差を区別しないため、有効数字12桁で同順位を判定する。
ranking <- function(x) rank(-signif(x, 12), ties.method = "min")
tables <- list()
layouts <- list()
for (i in seq_along(graphs)) {
  g <- graphs[[i]]
  id <- ids[i]
  n <- vcount(g)
  stopifnot(!is_directed(g), is_simple(g), is_connected(g), n > 2)
  d <- degree(g)
  close_raw <- closeness(g, mode = "all", weights = NA, normalized = FALSE)
  close_norm <- closeness(g, mode = "all", weights = NA, normalized = TRUE)
  between <- betweenness(g, directed = FALSE, weights = NA, normalized = FALSE)
  between_norm <- betweenness(g, directed = FALSE, weights = NA, normalized = TRUE)
  set.seed(3100+i)
  eig <- eigen_centrality(g, directed = FALSE, weights = NA)
  # 元の位置は課題1で保存していないので、課題3用の配置seedを別途固定。
  set.seed(3200+i)
  layouts[[id]] <- layout_with_fr(g, niter = 1500)

  # 定義・正規化の検証。媒介中心性の総和は全ノード対の(d-1)の総和。
  D <- distances(g, mode = "all", weights = NA)
  A <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
  stopifnot(isTRUE(all.equal(close_raw, as.numeric(1/rowSums(D)), check.attributes = FALSE)),
    isTRUE(all.equal(close_norm, as.numeric((n-1)/rowSums(D)), check.attributes = FALSE)),
    isTRUE(all.equal(between_norm, between/choose(n-1, 2))),
    isTRUE(all.equal(sum(between), sum(D[upper.tri(D)]-1))),
    max(abs(A %*% eig$vector - eig$value*eig$vector)) < 1e-7)
  block <- if ("block" %in% vertex_attr_names(g)) V(g)$block else rep(NA_integer_, n)
  external <- if (all(is.na(block))) rep(NA_integer_, n) else
    vapply(seq_len(n), function(v) sum(block[as.integer(neighbors(g, v))] != block[v]), 1L)
  t <- data.frame(id, model = model_names[i], node = seq_len(n), block,
    external_degree = external, degree_raw = d, degree_normalized = d/(n-1),
    closeness_raw = close_raw, closeness_normalized = close_norm,
    betweenness_raw = between, betweenness_normalized = between_norm,
    eigenvector = eig$vector)
  for (f in fields) t[[paste0("rank_", f)]] <- ranking(t[[f]])
  tables[[id]] <- t
}
all_scores <- do.call(rbind, tables)
rownames(all_scores) <- NULL
write.csv(all_scores, file.path(out, "node_centralities.csv"), row.names = FALSE)
saveRDS(layouts, file.path(out, "layouts.rds"))

# 上位5位以内を全て出力。同率5位なら6ノード以上になることもある。
tops <- list()
pair_tables <- list()
for (id in ids) {
  t <- tables[[id]]
  for (j in seq_along(fields)) {
    f <- fields[j]
    ranks <- t[[paste0("rank_", f)]]
    pos <- order(ranks, t$node)
    pos <- pos[ranks[pos] <= 5]
    tops[[length(tops)+1]] <- data.frame(id, model = t$model[1], metric = metric_names[j],
      rank = ranks[pos], node = t$node[pos], score = t[[f]][pos])
  }
  for (pair in combn(seq_along(fields), 2, simplify = FALSE)) {
    a <- t[[paste0("rank_", fields[pair[1]])]]
    b <- t[[paste0("rank_", fields[pair[2]])]]
    sa <- t$node[a <= 5]; sb <- t$node[b <= 5]
    pair_tables[[length(pair_tables)+1]] <- data.frame(id,
      metric_a = metric_names[pair[1]], metric_b = metric_names[pair[2]],
      spearman = cor(a, b, method = "spearman"),
      top5_a_count = length(sa), top5_b_count = length(sb),
      overlap_count = length(intersect(sa, sb)),
      jaccard = length(intersect(sa, sb))/length(union(sa, sb)))
  }
}
top <- do.call(rbind, tops)
pairs <- do.call(rbind, pair_tables)
write.csv(top, file.path(out, "top_rankings.csv"), row.names = FALSE)
write.csv(pairs, file.path(out, "ranking_comparison.csv"), row.names = FALSE)

palette <- c("#2878B5", "#E49B29", "#499D77", "#B66BA3")
plot_centrality <- function(id, f, title, label_rank = 1, label_cex = 0.85, monochrome = FALSE) {
  g <- graphs[[id]]; t <- tables[[id]]; scores <- t[[f]]
  top5 <- t[[paste0("rank_", f)]] <= 5
  # 全点を見せるため最小サイズを確保。面積の差を抑えるため平方根を使う。
  sizes <- 2.5 + 8*sqrt(scores/max(scores))
  colors <- if (id == "sbm_separated") palette[t$block] else palette[match(id,ids)]
  if (monochrome) colors <- "white"
  labels <- ifelse(t[[paste0("rank_", f)]] <= label_rank, t$node, NA)
  plot(g, layout = layouts[[id]], vertex.size = sizes,
    vertex.color = colors, vertex.frame.color = if (monochrome) "black" else ifelse(top5,"#263238",NA),
    vertex.label = labels, vertex.label.cex = label_cex, vertex.label.dist = 0.7,
    vertex.label.color = if (monochrome) "black" else "#18252B",
    edge.color = if (monochrome) "#909090" else "#8A98A64D", edge.width = 0.6,
    main = title, margin = 0.12)
  mtext(sprintf("Max=%.4f | labels: rank <= %d", max(scores), label_rank),
        side = 1, line = -0.7, cex = 0.68)
}
save_figure <- function(name, width, height, draw, res = 120) {
  png(file.path(out, "figures", paste0(name,".png")), width = width,
      height = height, units = "in", res = res, type = "cairo")
  draw(); dev.off()
  pdf(file.path(out, "figures", paste0(name,".pdf")), width = width, height = height)
  draw(); dev.off()
}
save_figure("task3_2_models", 12, 12, function() {
  par(mfrow = c(3,3), mar = c(1.7,0.3,2.5,0.3))
  for (i in 1:3) for (j in 1:3)
    plot_centrality(ids[i], fields[j], paste(model_names[i], "-", metric_names[j]))
})
save_figure("task3_2_sbm", 16, 5, function() {
  par(mfrow = c(1,4), mar = c(1.7,0.3,2.5,0.3))
  for (j in 1:4) plot_centrality("sbm_separated", fields[j], paste("SBM -",metric_names[j]))
})
save_figure("task3_3_eigenvector", 12, 5, function() {
  par(mfrow = c(1,3), mar = c(1.7,0.3,2.5,0.3))
  for (i in 1:3) plot_centrality(ids[i], "eigenvector", paste(model_names[i],"- Eigenvector"))
})
save_figure("task3_submission_sbm_betweenness", 7.5, 7.5, function() {
  par(mfrow = c(1,1), mar = c(2,0.4,3,0.4))
  plot_centrality("sbm_separated", "betweenness_normalized", "SBM - Betweenness centrality", label_rank = 5, label_cex = 1)
  legend("topleft", legend = paste("Given block", 1:4), pch = 19, col = palette, bty = "n", cex = 0.85)
}, res = 150)
save_figure("task3_submission_ba_comparison", 15, 5.8, function() {
  par(mfrow = c(1,3), mar = c(1.7,0.3,2.5,0.3), oma = c(1.2,0,2.2,0))
  for (j in 1:3) plot_centrality("ba_m2", fields[j], metric_names[j])
  mtext("BA model: N=200, m=2 (undirected)", side = 3, outer = TRUE, line = 0.3, font = 2, cex = 1.15)
  mtext("Same node positions in all panels. Larger nodes = higher centrality. Outlined nodes: top 5 ranks.",
        side = 1, outer = TRUE, line = 0, cex = 0.85)
}, res = 150)
save_figure("task3_submission_ba_betweenness", 7.5, 7.5, function() {
  par(mfrow = c(1,1), mar = c(2.5,0.4,3.5,0.4))
  plot_centrality("ba_m2", "betweenness_normalized", "BA - Betweenness centrality", label_rank = 5, label_cex = 1)
  mtext("N=200, m=2, undirected | Larger nodes = higher centrality", side = 3, line = 0.2, cex = 0.85)
}, res = 150)
save_figure("task3_submission_ba_betweenness_bw", 7.5, 7.5, function() {
  par(mfrow = c(1,1), mar = c(2.5,0.4,3.5,0.4))
  plot_centrality("ba_m2", "betweenness_normalized", "BA - Betweenness centrality",
                  label_rank = 5, label_cex = 1, monochrome = TRUE)
  mtext("N=200, m=2, undirected | Larger nodes = higher centrality", side = 3, line = 0.2, cex = 0.85)
}, res = 200)

# ランキングの概要も実測値から生成する。
lines <- c("# 中心性ランキング（実測値）", "", "ノードIDはRの1始まり。課題1の同一グラフ内で比較する。", "",
           "|モデル|中心性|上位5位以内のノードID（順位順）|", "|---|---|---|")
for (id in ids) for (m in metric_names) {
  s <- top[top$id == id & top$metric == m, ]
  lines <- c(lines, sprintf("|%s|%s|%s|", s$model[1], m,
              paste(sprintf("%d（%d位）",s$node,s$rank), collapse = ", ")))
}
writeLines(lines, file.path(out,"rankings.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out,"sessionInfo.txt"))
write.csv(data.frame(input = input, md5 = unname(tools::md5sum(input))),
          file.path(out,"input_fingerprint.csv"), row.names = FALSE)
print(top[top$rank == 1, ], row.names = FALSE)
print(pairs[pairs$metric_a == "Degree", ], row.names = FALSE)
cat("\nSaved: ", normalizePath(out), "\n", sep = "")
