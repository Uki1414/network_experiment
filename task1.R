#!/usr/bin/env Rscript
# 情報メディア実験・課題1。Rscript task1.R [出力先] で実行する。
# 対話実行する場合: source("task1.R", encoding = "UTF-8")
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop('igraph が必要です。先に install.packages("igraph") を実行してください。')
}
library(igraph)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else "output"
for (d in c(out_dir, file.path(out_dir, "figures"), file.path(out_dir, "graphs"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 図の配置に使う乱数と、グラフ生成に使う乱数を分ける。
# 比較条件ごとに seed を固定。同じグラフは複数の図でも再利用する。
graphs <- list()
metadata <- list()
add_graph <- function(id, model, parameters, seed, generator) {
  set.seed(seed)
  g <- generator()
  stopifnot(vcount(g) == 200L, !is_directed(g), is_simple(g))
  graphs[[id]] <<- g
  metadata[[id]] <<- data.frame(id, model, parameters, seed)
}

# 課題1-1: 補足資料と同じ n, p, nei, m。
# 旧関数名 erdos.renyi.game / watts.strogatz.game / ba.game は
# 現在の sample_gnp / sample_smallworld / sample_pa に置き換える。
add_graph("er_p005", "ER", "n=200, p=0.05", 101,
          function() sample_gnp(200, 0.05, directed = FALSE, loops = FALSE))
add_graph("ws_p005_k5", "WS", "dim=1, size=200, nei=5, p=0.05", 201,
          function() sample_smallworld(1, 200, 5, 0.05, loops = FALSE, multiple = FALSE))
add_graph("ba_m2", "BA", "n=200, m=2, power=1, zero.appeal=1", 301,
          function() sample_pa(200, m = 2, directed = FALSE))

# 課題1-2: 一度に変更するパラメータは一つにする。
for (i in seq_along(c(0.01, 0.10))) {
  p <- c(0.01, 0.10)[i]
  add_graph(c("er_p001", "er_p010")[i], "ER", sprintf("n=200, p=%g", p), 110+i,
            function() sample_gnp(200, p, directed = FALSE, loops = FALSE))
}
for (i in seq_along(c(0, 1))) {
  p <- c(0, 1)[i]
  add_graph(c("ws_p0_k5", "ws_p1_k5")[i], "WS",
            sprintf("dim=1, size=200, nei=5, p=%g", p), 210+i,
            function() sample_smallworld(1, 200, 5, p, loops = FALSE, multiple = FALSE))
}
for (i in seq_along(c(2, 10))) {
  k <- c(2, 10)[i]
  add_graph(paste0("ws_p005_k", k), "WS",
            sprintf("dim=1, size=200, nei=%d, p=0.05", k), 220+i,
            function() sample_smallworld(1, 200, k, 0.05, loops = FALSE, multiple = FALSE))
}
for (i in seq_along(c(1, 5))) {
  m <- c(1, 5)[i]
  add_graph(paste0("ba_m", m), "BA", sprintf("n=200, m=%d, power=1, zero.appeal=1", m),
            310+i, function() sample_pa(200, m = m, directed = FALSE))
}

# 課題1-3: 確率的ブロックモデル (SBM)。4集団 x 各50ノード。
# 同集団の確率を0.15に固定し、異なる集団の確率のみ変更。
for (i in seq_along(c(0.002, 0.05))) {
  p_out <- c(0.002, 0.05)[i]
  pm <- matrix(p_out, 4, 4)
  diag(pm) <- 0.15
  add_graph(c("sbm_separated", "sbm_mixed")[i], "SBM",
            sprintf("blocks=4x50, p_in=0.15, p_out=%g", p_out), 400+i,
            function() {
              g <- sample_sbm(200, pm, rep(50, 4), directed = FALSE, loops = FALSE)
              V(g)$block <- rep(1:4, each = 50)
              g
            })
}

# 観察を補助する指標。平均距離は非連結時の誤解を避けて最大連結成分のみ。
metrics <- do.call(rbind, lapply(names(graphs), function(id) {
  g <- graphs[[id]]
  cc <- components(g)
  largest <- induced_subgraph(g, which(cc$membership == which.max(cc$csize)))
  data.frame(metadata[[id]], nodes = vcount(g), edges = ecount(g),
    mean_degree = mean(degree(g)), max_degree = max(degree(g)),
    density = edge_density(g), components = cc$no, largest_component = vcount(largest),
    global_clustering = transitivity(g, type = "global"),
    mean_distance_lcc = mean_distance(largest, directed = FALSE),
    stringsAsFactors = FALSE)
}))
rownames(metrics) <- NULL
write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)
saveRDS(graphs, file.path(out_dir, "graphs.rds"))
for (id in names(graphs)) {
  # DOTは論理型のグラフ属性を持てないので、書き出すコピーだけ0/1にする。
  dot_graph <- graphs[[id]]
  for (a in graph_attr_names(dot_graph)) {
    value <- graph_attr(dot_graph, a)
    if (is.logical(value)) dot_graph <- set_graph_attr(dot_graph, a, as.integer(value))
  }
  write_graph(dot_graph, file.path(out_dir, "graphs", paste0(id, ".dot")), "dot")
  write_graph(graphs[[id]], file.path(out_dir, "graphs", paste0(id, ".graphml")), "graphml")
}

# FR: ばねモデルによる配置。距離や位置自体に物理的な意味はない。
# WSのパラメータ比較には固定された円周配置も使い、長距離辺を見やすくする。
layouts <- lapply(seq_along(graphs), function(i) {
  set.seed(1000+i)
  layout_with_fr(graphs[[i]], niter = 1000)
})
names(layouts) <- names(graphs)
palette <- c("#2878B5", "#E49B29", "#499D77", "#B66BA3")
plot_graph <- function(id, title, circle = FALSE) {
  g <- graphs[[id]]
  cols <- if ("block" %in% vertex_attr_names(g)) palette[V(g)$block] else "#2878B5"
  xy <- if (circle) layout_in_circle(g) else layouts[[id]]
  plot(g, layout = xy, vertex.size = 3.5, vertex.label = NA,
       vertex.color = cols, vertex.frame.color = NA,
       edge.color = "#8A98A655", edge.width = 0.65, main = title,
       margin = 0.05)
  r <- metrics[metrics$id == id, ]
  mtext(sprintf("N=%d  E=%d  mean k=%.2f  max k=%d", r$nodes, r$edges,
                r$mean_degree, r$max_degree), side = 1, line = -0.7, cex = 0.76)
}
save_figure <- function(name, width, height, draw) {
  # PNGとPDFで全く同じ座標を使う。英語ラベルでフォント依存を避ける。
  png(file.path(out_dir, "figures", paste0(name, ".png")),
      width = width, height = height, units = "in", res = 180, type = "cairo")
  draw()
  dev.off()
  pdf(file.path(out_dir, "figures", paste0(name, ".pdf")), width = width, height = height)
  draw()
  dev.off()
}
three <- function(ids, titles, circle = FALSE) {
  par(mfrow = c(1,3), mar = c(1.8,0.4,2.7,0.4), oma = c(0,0,0,0))
  for (i in seq_along(ids)) plot_graph(ids[i], titles[i], circle)
}
save_figure("task1_1_models", 15, 5.4, function() three(
  c("er_p005", "ws_p005_k5", "ba_m2"),
  c("ER: p=0.05", "WS: nei=5, p=0.05", "BA: m=2")))
save_figure("task1_1_degree", 12, 4.5, function() {
  par(mfrow = c(1,3), mar = c(4,4,3,1))
  ids <- c("er_p005", "ws_p005_k5", "ba_m2")
  max_k <- max(metrics$max_degree[match(ids, metrics$id)])
  for (i in seq_along(ids)) {
    counts <- tabulate(degree(graphs[[ids[i]]])+1, nbins = max_k+1)
    plot(0:max_k, counts, type = "h", lwd = 2, col = palette[i],
         xlab = "Degree k", ylab = "Number of nodes", main = c("ER", "WS", "BA")[i],
         xlim = c(0, max_k),
         ylim = c(0, max(vapply(ids, function(id) max(table(degree(graphs[[id]]))), 1))))
  }
})
save_figure("task1_2_er", 15, 5.4, function() three(
  c("er_p001", "er_p005", "er_p010"), c("ER: p=0.01", "ER: p=0.05", "ER: p=0.10")))
save_figure("task1_2_ws_p", 15, 5.4, function() three(
  c("ws_p0_k5", "ws_p005_k5", "ws_p1_k5"),
  c("WS: p=0, nei=5", "WS: p=0.05, nei=5", "WS: p=1, nei=5"), circle = TRUE))
save_figure("task1_2_ws_nei", 15, 5.4, function() three(
  c("ws_p005_k2", "ws_p005_k5", "ws_p005_k10"),
  c("WS: nei=2, p=0.05", "WS: nei=5, p=0.05", "WS: nei=10, p=0.05"), circle = TRUE))
save_figure("task1_2_ba", 15, 5.4, function() three(
  c("ba_m1", "ba_m2", "ba_m5"), c("BA: m=1", "BA: m=2", "BA: m=5")))
save_figure("task1_3_sbm", 12, 5.8, function() {
  par(mfrow = c(1,2), mar = c(3.2,0.4,2.7,0.4))
  plot_graph("sbm_separated", "SBM: p_in=0.15, p_out=0.002")
  legend("bottomleft", legend = paste("Block", 1:4), col = palette, pch = 19,
         bty = "n", cex = 0.8, ncol = 2)
  plot_graph("sbm_mixed", "SBM: p_in=0.15, p_out=0.05")
})
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
print(metrics[, c("id", "edges", "mean_degree", "max_degree", "components", "global_clustering", "mean_distance_lcc")], row.names = FALSE)
cat("\nSaved figures, graphs and metrics to: ", normalizePath(out_dir), "\n", sep = "")
