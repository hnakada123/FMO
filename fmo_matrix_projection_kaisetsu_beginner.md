# FMO境界処理と矩形変換による混成軌道の除去

## この文書の目的

この文書は、`fmo_matrix_projection_kaisetsu.md` の内容を、量子化学を学び始めた人にも追えるように書き直したものです。

扱うテーマは、FMO 法で分子をフラグメントに分けたとき、切断した結合の方向を向いた「使ってはいけない混成軌道」をどうやって計算から取り除くか、です。

結論を先に言うと、従来のように大きなレベルシフト $B$ でその軌道を無理に押し上げる必要はありません。代わりに、最初からその方向を含まない小さい空間で SCF を解けばよい、というのが核心です。この操作をここでは「矩形変換による deflation」と呼びます。

元文書では、フラグメント $x$ に属する量として $\mathbf H^x$ や $\mathbf S^x$ のように上付き添字 $x$ を付けています。この初学者向け版では、式を見やすくするため、基本的にこの $x$ を省いて書きます。

---

## 1. まず何が問題なのか

FMO は Fragment Molecular Orbital の略で、大きな分子をいくつかの小さなフラグメントに分けて量子化学計算を行う方法です。タンパク質や高分子のように全体を一度に扱うのが重い系では、この考え方が重要になります。

ただし、分子をフラグメントに分けるとき、共有結合を切ることがあります。たとえば、炭素と炭素の単結合を途中で切って、左側をフラグメント $x$、右側を別のフラグメントにする、という状況です。

共有結合を切ると、切断点には「もともと結合を作っていた方向」が残ります。炭素が sp$^3$ 混成している場合なら、4本の sp$^3$ 混成軌道のうち、切った結合の方向を向いた1本が問題になります。

模式的には次のような状況です。

```text
フラグメント x 側                  相手フラグメント側

   残す方向 ── Cα ────── [ 切断方向の sp3 混成軌道 ]
                                  ↑
                         x 側の計算では使わせたくない
```

この「切断方向の混成軌道」を、この文書では

$$
\varphi_i^h
$$
と書きます。複数本ある場合は $i=1,2,\ldots$ と番号を付けます。

フラグメント $x$ の計算でこの軌道をそのまま使える状態にしておくと、フラグメント $x$ の電子が、本来は相手フラグメント側に属する方向へ入り込んでしまいます。これは物理的に望ましくありません。

したがって、FMO の境界処理では次のことを行う必要があります。

> 切断した結合の方向を向いた混成軌道を、フラグメントの SCF 変分空間から取り除く。

「変分空間」とは、SCF が分子軌道を探すときに動かしてよい自由度の集合です。そこに不要な方向が残っていると、SCF はその方向も使ってしまう可能性があります。

---

## 2. 最低限の量子化学の言葉

ここで使う記号を、先に直感的に整理します。

### AO と MO

量子化学計算では、分子軌道 MO を原子軌道 AO の線形結合として表します。

$$
\psi_p = \sum_\mu C_{\mu p}\chi_\mu
$$
ここで、

| 記号 | 意味 |
|---|---|
| $\chi_\mu$ | AO。原子中心の基底関数 |
| $\psi_p$ | MO。分子全体に広がる軌道 |
| $C_{\mu p}$ | AO をどれだけ混ぜて MO を作るかを表す係数 |

すべての MO 係数をまとめた行列を

$$
\mathbf C
$$
と書きます。SCF は、ざっくり言えば、この $\mathbf C$ を調整して電子のエネルギーを下げる手続きです。

### AO は普通、直交していない

AO 同士の重なりは

$$
S_{\mu\nu} = \langle \chi_\mu | \chi_\nu \rangle
$$
で表します。これを集めたものが重なり行列

$$
\mathbf S
$$
です。

もし AO が完全に直交していれば $\mathbf S=\mathbf I$ です。しかし実際の原子中心基底では、異なる原子上の AO も空間的に重なります。そのため、一般には

$$
\mathbf S \neq \mathbf I
$$
です。

この「AO が非直交である」という事実が、後で式に $\mathbf S$ が何度も出てくる理由です。

### 除去したい混成軌道

切断方向の混成軌道も AO の線形結合で表せます。

$$
\varphi_i^h = \sum_\mu A_{\mu i}\chi_\mu
$$
この係数を列に並べた行列を

$$
\mathbf A
$$
と書きます。AO の数を $N$、除去したい混成軌道の本数を $m$ とすると、$\mathbf A$ は

$$
N \times m
$$
の行列です。

たとえば、切断方向の sp$^3$ 混成軌道が1本だけなら $m=1$ で、$\mathbf A$ は1列のベクトルです。

---

## 3. 「使ってはいけない方向を使っていない」とはどういう意味か

ある MO $\psi_p$ が、除去したい混成軌道 $\varphi_i^h$ を含んでいるかどうかは、両者の内積で見ます。

$$
\langle \varphi_i^h | \psi_p \rangle
$$
これがゼロなら、$\psi_p$ は $\varphi_i^h$ 方向の成分を持っていません。

すべての除去対象 $i$ とすべての MO $p$ についてゼロにしたいので、行列で書くと

$$
\mathbf A^\dagger \mathbf S \mathbf C = \mathbf 0
$$
です。

ここに $\mathbf S$ が入る理由は、AO が非直交だからです。実際、

$$
\begin{aligned}
\langle \varphi_i^h | \psi_p \rangle
&=
\left\langle
\sum_\mu A_{\mu i}\chi_\mu
\middle|
\sum_\nu C_{\nu p}\chi_\nu
\right\rangle \\
&=
\sum_{\mu\nu} A_{\mu i}^*
\langle \chi_\mu | \chi_\nu \rangle
C_{\nu p} \\
&=
\sum_{\mu\nu} A_{\mu i}^* S_{\mu\nu} C_{\nu p}.
\end{aligned}
$$
これを行列でまとめると $\mathbf A^\dagger\mathbf S\mathbf C$ になります。

つまり、境界処理の数学的な目的は次の一行にまとめられます。

$$
\boxed{\mathbf A^\dagger \mathbf S \mathbf C = \mathbf 0}
$$
この式は、

> どの MO も、切断方向の混成軌道を含まない

という意味です。

---

## 4. 従来法: レベルシフトで押し上げる

従来の方法では、使ってはいけない軌道に非常に大きなエネルギーを与えます。

考え方は単純です。

> その方向を使うとエネルギー的に大損するようにして、SCF がその方向を選ばないようにする。

この方法をレベルシフト射影、または hybrid orbital projection と呼びます。

コアハミルトニアンを $\mathbf H$、環境静電ポテンシャルを $\mathbf V$ とすると、レベルシフト込みの行列は次の形になります。

$$
\begin{aligned}
\widetilde{\mathbf H}
&= \mathbf H + \mathbf V \\
&\quad + \mathbf S\mathbf A\mathbf B\mathbf A^\dagger\mathbf S
\end{aligned}
$$
ここで $\mathbf B$ は大きな正の値を並べた対角行列です。

$$
\mathbf B = \mathrm{diag}(B_1, B_2, \ldots, B_m)
$$
典型的には、$B$ として非常に大きな値を入れます。これにより、$\varphi_i^h$ 方向の軌道エネルギーが大きく上がり、被占有軌道から追い出されることを期待します。

### なぜ $\mathbf S\mathbf A\mathbf B\mathbf A^\dagger\mathbf S$ なのか

射影項の成分表示は

$$
\begin{aligned}
P_{\mu\nu}
&= \sum_i B_i
\langle \chi_\mu | \varphi_i^h \rangle
\langle \varphi_i^h | \chi_\nu \rangle
\end{aligned}
$$
です。

混成軌道 $\varphi_i^h$ は AO 係数 $\mathbf A$ で表されますが、AO は非直交です。そのため、

$$
\langle \chi_\mu | \varphi_i^h \rangle
= (\mathbf S\mathbf A)_{\mu i}
$$
となります。同様に右側の内積からも $\mathbf S$ が出ます。したがって、

$$
\mathbf P = \mathbf S\mathbf A\mathbf B\mathbf A^\dagger\mathbf S
$$
となります。

もし AO が最初から直交していて $\mathbf S=\mathbf I$ なら、射影項は単に

$$
\mathbf A\mathbf B\mathbf A^\dagger
$$
になります。つまり、左右に出てくる2つの $\mathbf S$ は、非直交 AO を使っていることの結果です。

---

## 5. レベルシフト法の弱点

レベルシフト法には、わかりやすい弱点があります。

$B$ を大きくすれば、不要な混成軌道の成分は小さくなります。しかし、有限の $B$ では完全にはゼロになりません。残る成分は大まかには

$$
O(B^{-1})
$$
の大きさです。

一方で、$B$ を大きくしすぎると、行列の中に極端に大きな値が入ります。すると条件数が悪化し、数値的に不安定になります。SCF の収束にも悪影響が出る可能性があります。

つまり、レベルシフト法には次のトレードオフがあります。

| $B$ の選び方 | 良い点 | 悪い点 |
|---|---|---|
| 小さめ | 数値的には比較的安定 | 不要成分が残りやすい |
| 大きめ | 不要成分は小さくなる | 条件数が悪化しやすい |
| 無限大 | 理想的には完全除去 | 実際の計算では使えない |

この問題を避けるには、そもそも $B$ を使わない方法が必要です。

---

## 6. 新しい考え方: 罰するのではなく、最初から選べなくする

レベルシフト法は「使ってはいけない方向を選んだら罰する」方法です。

一方、矩形変換による方法では、発想を変えます。

> 使ってはいけない方向を、変数の中から最初から消しておく。

SCF が動かす変数は MO 係数 $\mathbf C$ です。そこで、$\mathbf C$ を自由に直接動かすのではなく、より小さい係数行列 $\mathbf C_r$ を使って

$$
\mathbf C = \mathbf T\mathbf C_r
$$
と表します。

ここで $\mathbf T$ は

$$
N \times (N-m)
$$
の縦長の行列です。列数が $N$ ではなく $N-m$ なので、自由度が $m$ 個減っています。このため「矩形変換」と呼びます。

$\mathbf T$ は、次の条件を満たすように作ります。

$$
\boxed{\mathbf A^\dagger \mathbf S \mathbf T = \mathbf 0}
$$
すると、

$$
\begin{aligned}
\mathbf A^\dagger \mathbf S \mathbf C
&= \mathbf A^\dagger \mathbf S \mathbf T\mathbf C_r \\
&= \mathbf 0\mathbf C_r \\
&= \mathbf 0
\end{aligned}
$$
となります。

これは非常に重要です。$\mathbf C_r$ がどんな値になっても、できあがる $\mathbf C=\mathbf T\mathbf C_r$ は必ず

$$
\mathbf A^\dagger \mathbf S \mathbf C = \mathbf 0
$$
を満たします。

つまり、SCF がどのように進んでも、不要な混成軌道は最初から使えません。

---

## 7. T は何を表しているのか

$\mathbf T$ は、許された方向だけを並べた行列です。

除去したい方向を表す条件は

$$
\mathbf A^\dagger \mathbf S \mathbf v = \mathbf 0
$$
です。これを満たすベクトル $\mathbf v$ は、除去したい混成軌道と $\mathbf S$ の意味で直交しています。

このような $\mathbf v$ 全体の空間を

$$
\ker(\mathbf A^\dagger \mathbf S)
$$
と書きます。$\ker$ は kernel、つまり核または零空間です。

$\mathbf T$ の列は、この核を張る基底です。

言い換えると、$\mathbf T$ は「不要な混成軌道と直交する方向だけを並べた行列」です。

AO 空間全体の次元が $N$、取り除く混成軌道が $m$ 本なら、残る許容空間の次元は

$$
N - m
$$
です。そのため、$\mathbf T$ は $N\times(N-m)$ の行列になります。

---

## 8. 直感的なたとえ

2つの方法の違いは、次のように考えると分かりやすいです。

| 方法 | たとえ | 何が起きるか |
|---|---|---|
| レベルシフト | 入ってはいけない道に高い坂を作る | 坂が有限なら、少しは入り込む可能性がある |
| 矩形変換 | 地図からその道を消す | そもそもその方向へ進めない |

レベルシフトは「罰する」方法です。矩形変換は「選択肢から消す」方法です。

この違いにより、矩形変換では不要成分が原理的にゼロになります。

---

## 9. 矩形変換を使うとレベルシフト項は消える

レベルシフト込みの行列は

$$
\begin{aligned}
\widetilde{\mathbf H}
&= \mathbf H + \mathbf V \\
&\quad + \mathbf S\mathbf A\mathbf B\mathbf A^\dagger\mathbf S
\end{aligned}
$$
でした。

これを、許された空間へ変換します。つまり、左から $\mathbf T^\dagger$、右から $\mathbf T$ を掛けます。

$$
\begin{aligned}
\mathbf T^\dagger\widetilde{\mathbf H}\mathbf T
&= \mathbf T^\dagger(\mathbf H+\mathbf V)\mathbf T \\
&\quad +
\mathbf T^\dagger
\mathbf S\mathbf A\mathbf B\mathbf A^\dagger\mathbf S
\mathbf T
\end{aligned}
$$
第2項は次のように書き換えられます。

$$
\begin{aligned}
&\mathbf T^\dagger
\mathbf S\mathbf A\mathbf B\mathbf A^\dagger\mathbf S
\mathbf T \\
&\quad =
(\mathbf A^\dagger\mathbf S\mathbf T)^\dagger
\mathbf B
(\mathbf A^\dagger\mathbf S\mathbf T)
\end{aligned}
$$
しかし、$\mathbf T$ は

$$
\mathbf A^\dagger\mathbf S\mathbf T = \mathbf 0
$$
を満たすように作っています。したがって、

$$
\begin{aligned}
(\mathbf A^\dagger\mathbf S\mathbf T)^\dagger
\mathbf B
(\mathbf A^\dagger\mathbf S\mathbf T)
&= \mathbf 0
\end{aligned}
$$
です。

よって、縮約空間で使う行列は

$$
\boxed{\mathbf H_r
= \mathbf T^\dagger(\mathbf H+\mathbf V)\mathbf T}
$$
になります。

ここには、もう $B$ が出てきません。

これは偶然ではありません。$B$ が担当していた「不要な方向を排除する」という役割を、$\mathbf T$ の条件

$$
\mathbf A^\dagger\mathbf S\mathbf T = \mathbf 0
$$
が完全に引き受けたからです。

---

## 10. 縮約空間での SCF 方程式

通常の Roothaan-Hall 方程式は、AO が非直交なので一般化固有値問題になります。

矩形変換後は、縮約空間で次の式を解きます。

$$
\begin{aligned}
\mathbf F_r\mathbf C_r
&= \mathbf S_r\mathbf C_r\boldsymbol\varepsilon
\end{aligned}
$$
ここで、

$$
\begin{aligned}
\mathbf S_r
&= \mathbf T^\dagger\mathbf S\mathbf T
\end{aligned}
$$
$$
\begin{aligned}
\mathbf F_r
&= \mathbf T^\dagger
\left[
\mathbf H+\mathbf V+\mathbf G[\mathbf D]
\right]
\mathbf T
\end{aligned}
$$
です。

$\mathbf G[\mathbf D]$ は二電子項、つまり Coulomb 項と exchange 項を表します。

密度行列も縮約空間から元の AO 空間へ戻して使います。

$$
\begin{aligned}
\mathbf D
&= \mathbf T\mathbf D_r\mathbf T^\dagger
\end{aligned}
$$
このようにしておけば、二電子積分や Fock 行列の構築は元の AO 空間の量として扱いつつ、SCF の自由度だけを $N-m$ 次元に減らせます。

---

## 11. S-正規直交化

$\mathbf T$ の列は、除去したい軌道と $\mathbf S$ の意味で直交していれば十分です。ただし、そのままだと縮約空間の重なり行列

$$
\mathbf S_r = \mathbf T^\dagger\mathbf S\mathbf T
$$
が単位行列にならないことがあります。

そこで、$\mathbf T$ の列を $\mathbf S$ の計量で正規直交化して、

$$
\mathbf T^\dagger\mathbf S\mathbf T = \mathbf I
$$
となるようにできます。

この場合、縮約空間の SCF 方程式は

$$
\begin{aligned}
\mathbf F_r\mathbf C_r
&= \mathbf C_r\boldsymbol\varepsilon
\end{aligned}
$$
という普通の固有値問題になります。

実装では、まず

$$
\ker(\mathbf A^\dagger\mathbf S)
$$
の基底を求め、それに対して Löwdin 直交化や正準直交化を行って $\mathbf T$ を作ります。

---

## 12. sp3 炭素での具体例

炭素の価電子 AO を単純化して

```text
s, px, py, pz
```

の4つだけとします。この4つから、4本の sp$^3$ 混成軌道を作れます。

理想的な場合、たとえば次のように書けます。

$$
\begin{aligned}
h_1 &= \frac12(s+p_x+p_y+p_z), \\
h_2 &= \frac12(s+p_x-p_y-p_z), \\
h_3 &= \frac12(s-p_x+p_y-p_z), \\
h_4 &= \frac12(s-p_x-p_y+p_z).
\end{aligned}
$$
もし切断した結合の方向が $h_1$ なら、フラグメント $x$ の計算では $h_1$ を使わせたくありません。

このとき、許される方向は $h_2,h_3,h_4$ です。局所的には、

$$
4 \to 3
$$
の次元削減が起きます。

この単純な直交基底の例では、炭素上の変換行列は

$$
\mathbf T_C =
\begin{array}{c|ccc}
& h_2 & h_3 & h_4\\
\hline
s   & \frac12 & \frac12  & \frac12 \\
p_x & \frac12 & -\frac12 & -\frac12 \\
p_y & -\frac12 & \frac12 & -\frac12 \\
p_z & -\frac12 & -\frac12 & \frac12
\end{array}
$$
です。

この $\mathbf T_C$ は、$h_1$ の列を落として、残りの3本だけを並べたものだと考えられます。

実際の計算では、基底はもっと大きく、AO も非直交です。そのため、係数はこのようなきれいな $1/2$ にはなりません。現実には $\mathbf S$ を使って

$$
\mathbf A^\dagger\mathbf S\mathbf T = \mathbf 0
$$
を満たすように $\mathbf T$ を作ります。

それでも、考え方は同じです。

> 切断方向の1本を除き、残りの許された方向だけで MO を作る。

---

## 13. d 軌道の 6→5 変換との似ている点と違う点

量子化学では、カーテシアン d 関数6個から、球面調和 d 関数5個へ変換することがあります。

これも

$$
6 \to 5
$$
の次元削減なので、形だけ見ると矩形変換に似ています。

ただし、意味は違います。

| 比較項目 | カーテシアン d から球面 d | FMO 境界の矩形変換 |
|---|---|---|
| 次元 | 6 から 5 | $N$ から $N-m$ |
| 落とすもの | d 関数に混じる冗長な s 的成分 | 切断結合方向の混成軌道 |
| 行列の決まり方 | 角運動量で決まるほぼ普遍的な変換 | 結合方向、$\mathbf A$、$\mathbf S$ に依存 |
| 目的 | 基底関数の表現を整理する | フラグメント境界の不要成分を除去する |

どちらも「長方形の行列で次元を減らす」という点では似ています。しかし、FMO 境界処理の $\mathbf T$ は、系ごとに決まる物理的な拘束を入れるための行列です。

---

## 14. レベルシフト法との関係

矩形変換による方法は、レベルシフト法と無関係な別手法ではありません。

むしろ、次のように理解できます。

$$
\boxed{\text{deflation} \Longleftrightarrow B\to\infty}
$$

つまり、矩形変換はレベルシフト法の $B\to\infty$ 極限に対応します。
レベルシフト法では、$B$ を大きくするほど不要な軌道の成分が小さくなります。理想的には $B$ を無限大にすれば完全に除去できます。

しかし、実際の数値計算で無限大の $B$ は使えません。非常に大きな有限値を使うと、条件数が悪くなります。

矩形変換は、この「$B$ を無限大にした理想状態」を、実際には $B$ を使わずに実現します。

そのため、

> レベルシフト $B$ をどの値にするか

という調整問題そのものが消えます。

---

## 15. 両者の比較

| 項目 | レベルシフト射影 | 矩形変換による deflation |
|---|---|---|
| 基本発想 | 不要方向をエネルギー的に罰する | 不要方向を変数から消す |
| 使う条件 | Fock 行列に $\mathbf S\mathbf A\mathbf B\mathbf A^\dagger\mathbf S$ を足す | $\mathbf A^\dagger\mathbf S\mathbf T=\mathbf 0$ を満たす $\mathbf T$ を使う |
| $B$ | 必要 | 不要 |
| 除去の正確さ | 有限 $B$ では近似的 | 原理的に厳密 |
| 残留成分 | およそ $O(B^{-1})$ | 機械精度でゼロ |
| 次元 | $N$ のまま | $N-m$ |
| 数値安定性 | 大きな $B$ で条件数が悪化しうる | $B$ による条件数悪化がない |
| 見方 | $B$ を大きくして追い出す | $B\to\infty$ 極限を有限の変換で実現する |

---

## 16. 実装の流れ

実装するなら、流れは次のようになります。

1. 切断した結合方向の混成軌道を AO 係数として作る: $\mathbf A$。
2. AO 重なり行列を用意する: $\mathbf S$。
3. 次の核を張る基底を求める: $\ker(\mathbf A^\dagger\mathbf S)$。
4. その基底を必要に応じて $\mathbf S$-正規直交化し、$\mathbf T^\dagger\mathbf S\mathbf T=\mathbf I$ を満たす $\mathbf T$ を作る。
5. Fock 行列や重なり行列を縮約空間へ変換する: $\mathbf F_r=\mathbf T^\dagger\mathbf F\mathbf T$、$\mathbf S_r=\mathbf T^\dagger\mathbf S\mathbf T$。
6. 縮約空間で SCF 方程式 $\mathbf F_r\mathbf C_r=\mathbf S_r\mathbf C_r\boldsymbol\varepsilon$ を解く。
7. 必要に応じて元の AO 空間の MO 係数や密度に戻す: $\mathbf C=\mathbf T\mathbf C_r$、$\mathbf D=\mathbf T\mathbf D_r\mathbf T^\dagger$。

この流れでは、レベルシフト $B$ はどこにも登場しません。

---

## 17. まとめ

FMO では、共有結合を切ってフラグメントを作るため、切断方向の混成軌道をそのフラグメントの SCF から取り除く必要があります。

従来のレベルシフト法は、不要な混成軌道に大きなエネルギー $B$ を与えて、被占有軌道から追い出します。しかし、有限の $B$ では不要成分が完全には消えず、$B$ を大きくしすぎると数値的に不安定になります。

矩形変換では、不要な混成軌道と $\mathbf S$-直交する空間だけを張る行列 $\mathbf T$ を作り、

$$
\mathbf C = \mathbf T\mathbf C_r
$$
と置きます。このとき

$$
\mathbf A^\dagger\mathbf S\mathbf T = \mathbf 0
$$
なので、どんな $\mathbf C_r$ を使っても

$$
\mathbf A^\dagger\mathbf S\mathbf C = \mathbf 0
$$
が成り立ちます。したがって、不要な混成軌道は最初から MO に混ざりません。

この方法は、レベルシフト法の $B\to\infty$ 極限を、$B$ を使わずに有限の行列操作で実現するものです。結果として、レベルシフトの値を調整する必要がなくなり、残留成分と条件数悪化のトレードオフも避けられます。

---

## 記号一覧

| 記号 | 意味 |
|---|---|
| $\mathbf H$ | コアハミルトニアン |
| $\mathbf V$ | 他フラグメントなどから来る環境静電ポテンシャル |
| $\mathbf S$ | AO 重なり行列 |
| $\chi_\mu$ | AO 基底関数 |
| $\psi_p$ | MO |
| $\mathbf C$ | 元の AO 空間での MO 係数 |
| $\mathbf C_r$ | 縮約空間での MO 係数 |
| $\varphi_i^h$ | 除去したい混成軌道 |
| $\mathbf A$ | 除去したい混成軌道の AO 係数 |
| $m$ | 除去したい混成軌道の本数 |
| $N$ | AO 基底関数の数 |
| $\mathbf B$ | レベルシフト行列 |
| $\mathbf T$ | 許容空間だけを張る矩形変換行列 |
| $\mathbf F$ | Fock 行列 |
| $\mathbf G[\mathbf D]$ | 二電子項 |
| $\mathbf D$ | 密度行列 |
| $\boldsymbol\varepsilon$ | 軌道エネルギー |
