# FMO境界処理と矩形変換による混成軌道の除去

FMO法（Fragment Molecular Orbital法）でフラグメント間の共有結合を切断したときに生じる「切断方向の混成軌道」を、SCF計算からどのように取り除くかを解説した文書です。

- 解説文書（PDF）: [`fmo_matrix_projection_kaisetsu.pdf`](fmo_matrix_projection_kaisetsu.pdf)
- LaTeXソース: [`fmo_matrix_projection_kaisetsu.tex`](fmo_matrix_projection_kaisetsu.tex)

## 概要

FMO法では、大きな分子をフラグメントに分割する際に共有結合を切断することがあります。切断点の炭素（α炭素）には、切った結合の方向を向いた混成軌道（例: sp³混成軌道の1本）が残り、これをそのままにするとフラグメントの電子が本来相手側に属する方向へ入り込んでしまいます。

この文書では、その混成軌道を変分空間から除去する2つの方法を比較し、特に後者の「矩形変換によるdeflation」を詳しく解説します。

1. **従来法（レベルシフト射影 / hybrid orbital projection）**
   Fock行列に射影項 $`\mathbf{S A B A}^\dagger \mathbf{S}`$ を加え、大きなレベルシフト $`B`$ で不要軌道をエネルギー的に押し上げる。有限の $`B`$ では成分が $`O(B^{-1})`$ 程度残り、 $`B`$ を大きくしすぎると条件数が悪化する。

2. **矩形変換によるdeflation**
   条件 $`\mathbf{A}^\dagger \mathbf{S T} = \mathbf{0}`$ を満たす $`N \times (N-m)`$ の矩形変換行列 $`\mathbf{T}`$ を作り、 $`\mathbf{C} = \mathbf{T C}_\mathrm{r}`$ と置いて縮約空間でSCFを解く。不要な方向は最初から変分空間に存在しないため、残留成分は機械精度でゼロになり、レベルシフト $`B`$ の調整問題そのものが消える。矩形変換はレベルシフト法の $`B \to \infty`$ 極限に対応する。

## 主な内容

- FMO境界処理で何が問題になるか（切断結合方向の混成軌道）
- 除去条件 $`\mathbf{A}^\dagger \mathbf{S C} = \mathbf{0}`$ の導出（非直交AO基底で重なり行列 $`\mathbf{S}`$ が現れる理由）
- レベルシフト法の仕組みとその弱点（残留成分と条件数のトレードオフ）
- 矩形変換行列 $`\mathbf{T}`$ の構成（ $`\ker(\mathbf{A}^\dagger \mathbf{S})`$ の基底と $`S`$-正規直交化）
- 縮約空間でのSCF方程式 $`\mathbf{F}_\mathrm{r} \mathbf{C}_\mathrm{r} = \mathbf{S}_\mathrm{r} \mathbf{C}_\mathrm{r} \boldsymbol{\varepsilon}`$
- sp³炭素での具体例（4→3の次元削減）
- 結合切断時の原子核・電子数の扱い（核とAOは分割されず、減るのは変分自由度のみ）
- カーテシアン d 軌道の 6→5 変換との類似点と相違点
- 両手法の比較表と実装の流れ

## ビルド方法

upLaTeX でコンパイルします（UTF-8）。

```sh
uplatex fmo_matrix_projection_kaisetsu.tex
dvipdfmx fmo_matrix_projection_kaisetsu.dvi
```

`latexmk` を使う場合:

```sh
latexmk -pdfdvi -latex=uplatex fmo_matrix_projection_kaisetsu.tex
```
