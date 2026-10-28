# 課題3の計算を数行ずつ確認する。作業ディレクトリはリポジトリのルート。
library(igraph)
graphs <- readRDS("output/graphs.rds")
g <- graphs$ws_p005_k5 # er_p005、ba_m2、sbm_separatedにも変更してみる。
stopifnot(is_connected(g), !is_directed(g))

# 次数中心性: 直接つながる相手の数。
d <- degree(g)
# 近接中心性: 全体への平均最短距離の逆数（正規化あり）。
c <- closeness(g, mode="all", weights=NA, normalized=TRUE)
# 媒介中心性: 他ノード対の最短経路を仲介する割合の合計（正規化あり）。
b <- betweenness(g, directed=FALSE, weights=NA, normalized=TRUE)
# 固有ベクトル中心性（発展）: 重要な相手とのつながりを評価。
e <- eigen_centrality(g, directed=FALSE, weights=NA)$vector

head(data.frame(node=seq_len(vcount(g)), degree=d, closeness=c, betweenness=b, eigenvector=e))
head(order(d, decreasing=TRUE), 5) # 表示上は同順位でも5件で切る簡易版。
head(order(c, decreasing=TRUE), 5)
head(order(b, decreasing=TRUE), 5)

# 配置は同じまま、中心性に応じてサイズを変える。
set.seed(3202)
xy <- layout_with_fr(g, niter=1500)
par(mfrow=c(1,3))
for (name in c("Degree", "Closeness", "Betweenness")) {
  s <- switch(name, Degree=d, Closeness=c, Betweenness=b)
  plot(g, layout=xy, vertex.size=2.5+8*sqrt(s/max(s)), vertex.label=NA, main=name)
}
par(mfrow=c(1,1))
