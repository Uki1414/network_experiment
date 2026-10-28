#!/usr/bin/env Rscript
# 課題2。課題1の保存済みグラフを読み込み、同一のグラフを分析する。
# リポジトリのルートで Rscript task2.R を実行する。
library(igraph)
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else "output/graphs.rds"
out <- if (length(args) >= 2) args[2] else "output/task2"
if (!file.exists(input)) stop("入力がありません。先に Rscript task1.R を実行してください。")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
graphs <- readRDS(input)
meta <- read.csv(file.path(dirname(input), "metrics.csv"))
baseline <- c("er_p005", "ws_p005_k5", "ba_m2")
stopifnot(all(c(baseline, "ba_m5", "ws_p0_k5", "ws_p1_k5") %in% names(graphs)))

# 課題2-1: k=0を含む確率分布。度数と割合を区別する。
degree_tables <- lapply(names(graphs), function(id) {
  g <- graphs[[id]]
  k <- 0:max(degree(g))
  count <- tabulate(degree(g)+1L, nbins = length(k))
  stopifnot(sum(count) == vcount(g), sum(k*count) == 2*ecount(g))
  data.frame(id, degree = k, count, probability = count/vcount(g))
})
dd <- do.call(rbind, degree_tables)
write.csv(dd, file.path(out, "degree_distribution.csv"), row.names = FALSE)

# 課題2-2, 2-3, 2-5: 定義を明記し、非連結グラフの平均距離を区別する。
metrics <- do.call(rbind, lapply(names(graphs), function(id) {
  g <- graphs[[id]]
  stopifnot(!is_directed(g), is_simple(g), vcount(g) > 1)
  d <- degree(g)
  cc <- components(g)
  lcc <- induced_subgraph(g, which(cc$membership == which.max(cc$csize)))
  C <- transitivity(g, type = "global")
  L <- mean_distance(g, directed = FALSE, weights = NA, unconnected = FALSE)
  # igraphとは別の集計経路でも値を検証。自己間距離は除外。
  C_check <- sum(count_triangles(g))/sum(choose(d, 2))
  D <- distances(g, mode = "all", weights = NA)
  pairs <- D[upper.tri(D)]
  stopifnot(isTRUE(all.equal(C, C_check)), isTRUE(all.equal(L, mean(pairs))))
  data.frame(id, model = meta$model[match(id, meta$id)],
    parameters = meta$parameters[match(id, meta$id)],
    seed = meta$seed[match(id, meta$id)],
    nodes = vcount(g), edges = ecount(g), mean_degree = mean(d), max_degree = max(d),
    components = cc$no, largest_component = vcount(lcc),
    global_clustering = C,
    average_local_clustering = transitivity(g, type = "average", isolates = "NaN"),
    mean_path_all_pairs = L,
    mean_path_reachable_pairs = mean(pairs[is.finite(pairs)]),
    mean_path_lcc = mean_distance(lcc, directed = FALSE, weights = NA, unconnected = FALSE),
    unreachable_pair_fraction = mean(!is.finite(pairs)),
    degree_assortativity = assortativity_degree(g, directed = FALSE))
}))
rownames(metrics) <- NULL
write.csv(metrics, file.path(out, "metrics.csv"), row.names = FALSE)
nodes <- do.call(rbind, lapply(names(graphs), function(id) {
  g <- graphs[[id]]
  data.frame(id, node = seq_len(vcount(g)), degree = degree(g),
             local_clustering = transitivity(g, type = "local", isolates = "NaN"))
}))
write.csv(nodes, file.path(out, "node_metrics.csv"), row.names = FALSE)

cols <- c("#2878B5", "#E49B29", "#499D77", "#B66BA3")
save_figure <- function(name, width, height, draw) {
  png(file.path(out, "figures", paste0(name, ".png")), width = width,
      height = height, units = "in", res = 180, type = "cairo")
  draw(); dev.off()
  pdf(file.path(out, "figures", paste0(name, ".pdf")), width = width, height = height)
  draw(); dev.off()
}
save_figure("task2_1_degree_distribution", 13, 8, function() {
  par(mfrow = c(2,3), mar = c(4,4,2.7,1), oma = c(0,0,0,0))
  titles <- c("ER: p=0.05", "WS: nei=5, p=0.05", "BA: m=2")
  max_k <- max(dd$degree[dd$id %in% baseline])
  max_p <- max(dd$probability[dd$id %in% baseline])
  for (logarithmic in c(FALSE, TRUE)) for (i in seq_along(baseline)) {
    t <- dd[dd$id == baseline[i], ]
    if (logarithmic) {
      # 0の対数は定義できない。小さい正数への置換や再正規化はしない。
      t <- t[t$degree > 0 & t$probability > 0, ]
      plot(t$degree, t$probability, log = "xy", type = "p", pch = 19,
           col = cols[i], xlim = c(1,max_k), ylim = c(1/200, max_p*1.1),
           xlab = "Degree k (log scale)", ylab = "P(k) (log scale)",
           main = paste(titles[i], "- log-log"))
    } else {
      plot(t$degree, t$probability, type = "h", lwd = 2.5,
           col = cols[i], xlim = c(0,max_k), ylim = c(0,max_p*1.1),
           xlab = "Degree k", ylab = "P(k)", main = paste(titles[i], "- linear"))
    }
  }
})
bar_metric <- function(ids, labels, field, title, ylab, colors) {
  vals <- metrics[[field]][match(ids, metrics$id)]
  pos <- barplot(vals, names.arg = labels, col = colors, border = NA,
                 main = title, ylab = ylab, ylim = c(0, max(vals)*1.25), las = 1)
  text(pos, vals, labels = sprintf("%.3f", vals), pos = 3, cex = 0.9)
}
save_figure("task2_2_3_comparison", 11, 5, function() {
  par(mfrow = c(1,2), mar = c(4,4,3,1))
  ids <- c(baseline, "ba_m5")
  labels <- c("ER", "WS", "BA m=2", "BA m=5")
  bar_metric(ids, labels, "global_clustering", "Global clustering", "C", cols)
  bar_metric(ids, labels, "mean_path_all_pairs", "Mean shortest path length", "L", cols)
})
save_figure("task2_ws_rewiring", 10, 4.5, function() {
  par(mfrow = c(1,2), mar = c(4,4,3,1))
  ids <- c("ws_p0_k5", "ws_p005_k5", "ws_p1_k5")
  labels <- c("p=0", "p=0.05", "p=1")
  bar_metric(ids, labels, "global_clustering", "WS: global clustering", "C", cols[1:3])
  bar_metric(ids, labels, "mean_path_all_pairs", "WS: shortcuts reduce distance", "L", cols[1:3])
})

# 課題2-4: 実測値からメール本文の草稿を生成する。送信はしない。
r <- metrics[match(baseline, metrics$id), ]
mail <- c(
  "各200ノードのER（p=0.05）、WS（nei=5, p=0.05）、BA（m=2）について、次数分布・クラスタリング係数・平均経路長を比較した。",
  sprintf("ERは次数が平均%.2f付近に集中し、全体クラスタリング係数C=%.4f、平均経路長L=%.4fだった。", r$mean_degree[1], r$global_clustering[1], r$mean_path_all_pairs[1]),
  sprintf("WSは次数10付近への集中が強く、C=%.4fとERより高い一方、L=%.4fで、局所的なまとまりと短い経路が共存していた。", r$global_clustering[2], r$mean_path_all_pairs[2]),
  sprintf("BAは低次数のノードが多い一方、最大次数%dのハブが存在し、C=%.4f、L=%.4fだった。BAは平均次数%.2fとER・WSより疎なため、指標の差には辺数の違いも含まれる。", r$max_degree[3], r$global_clustering[3], r$mean_path_all_pairs[3], r$mean_degree[3]),
  "以上は各条件1試行の結果である。3グラフはすべて連結しており、Lは自己間を除く全ノード対の最短距離の平均、Cは全体クラスタリング係数を用いた。両対数図だけからべき乗則の成立は断定しない。")
writeLines(mail, file.path(out, "email_draft.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
write.csv(data.frame(input = input, md5 = unname(tools::md5sum(input))),
          file.path(out, "input_fingerprint.csv"), row.names = FALSE)
print(metrics[match(c(baseline, "ba_m5"), metrics$id),
              c("id", "mean_degree", "max_degree", "global_clustering", "mean_path_all_pairs", "degree_assortativity")], row.names = FALSE)
cat("\nSaved: ", normalizePath(out), "\n", sep = "")
