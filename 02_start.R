# リポジトリのルートを作業ディレクトリにし、上から順に実行する。
library(igraph)
graphs <- readRDS("output/graphs.rds")
g <- graphs$er_p005 # ws_p005_k5 / ba_m2 に替えて比較する。

# 課題2-1: P(k) = 次数kのノード数 / 全ノード数。
dd <- degree_distribution(g)
k <- seq_along(dd)-1
plot(k, dd, type = "h", xlab = "Degree k", ylab = "P(k)")
keep <- k > 0 & dd > 0
plot(k[keep], dd[keep], log = "xy", pch = 19, xlab = "Degree k", ylab = "P(k)")

# 課題2-2: 三角形数をT、中心ノードを持つ三つ組の数をWとして C=3T/W。
transitivity(g, type = "global")

# 課題2-3: 自分自身を除くノード対の最短距離の平均。
# 非連結ならInfを返す指定。今回の基準3グラフはいずれも連結。
mean_distance(g, directed = FALSE, weights = NA, unconnected = FALSE)

# 課題2-5（発展）: 隣接ノードの次数の相関。
assortativity_degree(g, directed = FALSE)
