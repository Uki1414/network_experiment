# 情報メディア実験：ネットワーク分析

## 課題2

課題1と同じグラフを読み込み、次数分布（線形・両対数）、全体クラスタリング係数、平均経路長を求める。

```sh
Rscript task2.R
```

- `02_start.R`：数行ずつ試すための短いコード。
- `task2.R`：課題2-1〜2-5の集計、図、報告本文の草稿を生成。
- `task2_notes.md`：定義、結果、考察、実行方法。
- `output/task2/figures/`：次数分布と指標の比較図（PNG / PDF）。
- `output/task2/metrics.csv`：課題1の全13条件の指標。
- `output/task2/email_draft.txt`：課題2-4の報告本文案。送信は別途行う。

入力の `output/graphs.rds` がない場合は、先に `Rscript task1.R` を実行する。

## 課題1

課題1-1〜1-3の実行用コードと観察用の図。全グラフは200ノード、無向、自己ループ・多重辺なし。補足資料に合わせてRとigraphを使う。

## 最初に見るファイル

- `01_start.R`：課題1-1を上から少しずつ実行する短いコード。
- `task1.R`：課題1の全条件の生成・図とデータの保存。
- `observations.md`：今回の実行結果と、図を読むときの着眼点。
- `output/figures/`：比較図（PNG / PDF）。
- `output/metrics.csv`：条件・乱数seed・観察を補助する指標。
- `output/graphs/`：生成したグラフ（DOT / GraphML）。
- `output/graphs.rds`：Rで読み直せる全グラフ。
- `output/sessionInfo.txt`：実行したRとパッケージのバージョン。

## 実行

ターミナルから全課題を実行する。

```sh
cd network_experiment
Rscript task1.R
```

Rの対話環境を使う場合は、次のようにする。

```r
# 作業ディレクトリをこのリポジトリのルートにして実行する。
source("01_start.R", encoding = "UTF-8")
# 全条件を保存する場合
source("task1.R", encoding = "UTF-8")
```

再実行すると `output/` の同名ファイルを更新する。seedとR・igraphのバージョンが同じなら同じグラフを再現できる。配置用乱数はグラフ生成用と分けてあり、PNGとPDFは同じ座標を使用する。

別環境でigraphがない場合は、Rで `install.packages("igraph", repos="https://cloud.r-project.org")` を一度実行する。

## 課題1-1：3モデルを生成し、特徴を観察する

基準条件は補足資料と同じ。古いAPI名は現在の関数名に置き換えた。

|モデル|補足資料の関数|今回の関数|基準条件|
|---|---|---|---|
|ER|`erdos.renyi.game()`|`sample_gnp()`|`n=200, p=0.05`|
|WS|`watts.strogatz.game()`|`sample_smallworld()`|`dim=1, size=200, nei=5, p=0.05`|
|BA|`ba.game()`|`sample_pa()`|`n=200, m=2, directed=FALSE`|

1. `output/figures/task1_1_models.png` で全体構造を見る。
2. `task1_1_degree.png` で次数のばらつきを確認する。次数とは1ノードにつながる辺の数。
3. ERのランダムなつながり、WSの局所的なつながり、BAの高次数ノード（ハブ）に着目する。

ERの `p` は各ノード対に辺を張る確率で、実現した密度そのものではない。期待辺数は `n*(n-1)/2*p`、期待平均次数は `(n-1)*p`。WSの `nei=5` は1次元円環上で左右5個ずつにつなぐ設定なので初期次数は10。BAの `m=2` は追加ノードが既存ノードに張る辺数で、平均次数はおおむね4になる。

基準のER・WSは平均次数がおおむね10、BAはおおむね4となり、見た目の違いには密度の違いも含まれる。BAの `m=5` の図も参照すると、平均次数を近づけて比較できる。ただし厳密に辺数を一致させた比較ではない。

## 課題1-2：パラメータ変更

|図|固定する条件|変える条件|確認すること|
|---|---|---|---|
|`task1_2_er`|`n=200`|`p=0.01, 0.05, 0.10`|辺数、孤立・分断、密度|
|`task1_2_ws_p`|`size=200, nei=5`|`p=0, 0.05, 1`|規則的な近傍辺から長距離辺への変化|
|`task1_2_ws_nei`|`size=200, p=0.05`|`nei=2, 5, 10`|近傍の広がりと辺数|
|`task1_2_ba`|`n=200`|`m=1, 2, 5`|辺数、ハブ、木構造からの変化|

各比較では一度に一つのパラメータを変える。乱数seedは条件ごとに異なるため、個別の辺の差をパラメータ変更だけの効果と断定しない。この実行は各条件1回の観察用であり、モデルの一般的傾向を統計的に示すものではない。

WSの比較図はノードの順番を固定した円周配置。規則的な近傍辺とつなぎ替えられた辺の区別がつきやすい。実装はigraph版で、原論文のWSと異なり辺の両端をつなぎ替える。したがって `p=1` では同じノード数・辺数のランダムグラフになる。

## 課題1-3：確率的ブロックモデル（SBM）

所属する集団（ブロック）に応じて、ノード対を結ぶ確率を設定するモデル。例えばクラス内・クラス間の交流を単純化した構造を作れる。集団は生成時に与えており、データから推定した集団ではない。

今回は4集団各50ノードとし、同集団内の確率 `p_in=0.15` を固定する。集団間確率 `p_out=0.002` と `0.05` を比較する。`pref.matrix` の対角要素が集団内、それ以外が集団間の確率に対応する。

```r
pm <- matrix(0.002, nrow = 4, ncol = 4)
diag(pm) <- 0.15
g <- sample_sbm(n = 200, pref.matrix = pm,
                block.sizes = rep(50, 4), directed = FALSE, loops = FALSE)
```

`task1_3_sbm.png` で同色の集団がまとまって見えるか観察する。色は既知の所属を表し、ノード座標は全体に対するばねモデルで決めている。位置を集団ごとに手動で分離してはいない。`p_out` を大きくすると集団間の辺数と全体密度も増えるため、密度一定での比較ではない。

## 図・指標の読み方

- 円はノード、線は辺。ノードサイズは全て一定。
- ばねモデル（Fruchterman–Reingold法）で近くに配置されたことは、直接の辺があることを必ずしも意味しない。
- `components`：連結成分数。1なら全ノードがひとつながり。
- `global_clustering`：全体クラスタ係数。三つ組に対して三角形がどれくらいあるかを示す。
- `mean_distance_lcc`：**最大連結成分内**のノード対の平均最短距離。非連結グラフを含むため、全ノード対の平均として扱わない。
- 200ノードの図・次数分布だけから「べき乗則に従う」とは断定できない。

## 参照

- [授業ページ（課題1）](https://www.snlab.cs.tsukuba.ac.jp/s-tugawa/jikken/)
- [生成と可視化の補足資料](https://www.snlab.cs.tsukuba.ac.jp/s-tugawa/jikken/suppl1.html)
- [igraph：ER](https://r.igraph.org/reference/sample_gnp.html)
- [igraph：WSと原モデルとの差](https://r.igraph.org/reference/sample_smallworld.html)
- [igraph：優先的選択（BA系）と既定値](https://r.igraph.org/reference/sample_pa.html)
- [igraph：SBM](https://r.igraph.org/reference/sample_sbm.html)

BAは補足資料の `ba.game()` に合わせた `sample_pa()` の既定設定（`power=1, zero.appeal=1`）を使う。優先的選択モデルには設定の異なる変種があり、原論文のBAと全く同じ条件だとは扱わない。
