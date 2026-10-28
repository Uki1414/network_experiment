# 課題1-1を一つずつ確認するための短いスクリプト。
# Rで上から数行ずつ実行する。保存を含む全課題は task1.R を使う。
library(igraph)

# ER: すべてのノード対に、独立に確率pで辺を張る。
set.seed(101)
er <- sample_gnp(n = 200, p = 0.05, directed = FALSE, loops = FALSE)
plot(er, vertex.label = NA, vertex.size = 4, main = "ER: n=200, p=0.05")

# WS: 円環状に近傍を結んだグラフの辺を、確率pでつなぎ替える。
# dim=1では各ノードの初期次数は2*nei=10。
set.seed(201)
ws <- sample_smallworld(dim = 1, size = 200, nei = 5, p = 0.05)
plot(ws, layout = layout_in_circle(ws), vertex.label = NA,
     vertex.size = 4, main = "WS: nei=5, p=0.05")

# BA: ノードを追加しながら、既存ノードに優先的選択でつなぐ。
# mは新しいノードが張る辺数。初期段階では相手が足りない場合がある。
set.seed(301)
ba <- sample_pa(n = 200, m = 2, directed = FALSE)
plot(ba, vertex.label = NA, vertex.size = 4, main = "BA: n=200, m=2")

# 並べて見る。
par(mfrow = c(1, 3))
plot(er, vertex.label = NA, vertex.size = 3, main = "ER")
plot(ws, vertex.label = NA, vertex.size = 3, main = "WS")
plot(ba, vertex.label = NA, vertex.size = 3, main = "BA")
par(mfrow = c(1, 1))
