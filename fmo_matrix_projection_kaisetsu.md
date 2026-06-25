# FMO境界処理における不要混成軌道の除去：レベルシフト射影と矩形変換

> FMO のフラグメント境界処理において、切断結合方向の混成軌道を変分空間から除去する2つの方法
> ——(i) 巨大レベルシフトによる射影、(ii) 矩形変換による厳密 deflation——を対比し、両者が
> $B\to\infty$ 極限で一致することを示す。式番号は本ノート内で完結する通し番号 (1)–(8)。

---

## 1. 問題設定

FMO では共有結合を跨いでフラグメント分割するため、結合切断点（BDA: bond-detached atom）に、
切断結合方向を向いた混成軌道 $\{\varphi_i^{h}\}$ が現れる。これらをフラグメント $x$ のモノマー SCF の
変分空間に残すと、被占有 MO に非物理的な成分が混入する。したがって $\{\varphi_i^{h}\}$ を
**変分空間から排除**する必要がある。標準的手法は、これらに巨大なレベルシフトを与える
射影（hybrid orbital projection, HOP）である。

記号： $\mathbf H^x$ コアハミルトニアン、 $\mathbf V^x$ 環境静電ポテンシャル（他フラグメントの ESP）、
$\mathbf S^x$ AO 重なり行列、 $\mathbf A^x$ は $\varphi_i^{h}=\sum_\mu A^x_{\mu i}\chi_\mu$ の展開係数を列に持つ
$N\times m$ 行列（ $m$ は除去本数）、 $\mathbf B=\mathrm{diag}(B_i)$ レベルシフト。

---

## 2. レベルシフト射影（式1）

成分表示（ $\mu,\nu$ は AO 添字）：

```math
\widetilde{H}^{x}_{\mu\nu} = H^{x}_{\mu\nu} + V^{x}_{\mu\nu} + \sum_i B_i \langle \mu | \varphi_i^{h} \rangle \langle \varphi_i^{h} | \nu \rangle \qquad (1)
```

ここで $\langle\mu|\varphi_i^{h}\rangle = (\mathbf S^x\mathbf A^x)_{\mu i}$。第3項は射影演算子
$\hat P=\sum_i B_i\,|\varphi_i^{h}\rangle\langle\varphi_i^{h}|$ の AO 表現であり、行列でまとめると

```math
\widetilde{\mathbf H}^{x} = \mathbf H^{x} + \mathbf V^{x} + \mathbf S^{x}\mathbf A^{x}\mathbf B\,(\mathbf A^{x})^{\dagger}\mathbf S^{x} \qquad (1')
```

### 成分表示 (1) から行列形 (1′) へ — なぜ重なり行列 $\mathbf S^x$ が2個現れるか

行列形 (1′) の射影項には $\mathbf S^x$ が左右に2個現れる。これは **AO 基底が非直交であること**の帰結であり、射影項の左右の内積からそれぞれ1個ずつ湧く。順に追う。

射影項の成分は

```math
P^{x}_{\mu\nu} = \sum_i B_i\,\langle \mu | \varphi_i^{h}\rangle\,\langle \varphi_i^{h} | \nu\rangle
```

で、 $|\mu\rangle=|\chi_\mu\rangle$ は AO 基底関数、 $|\varphi_i^{h}\rangle$ は除去対象の混成軌道である。混成軌道は AO の線形結合（展開係数が $\mathbf A^x$）で表される：

```math
|\varphi_i^{h}\rangle = \sum_\lambda A^{x}_{\lambda i}\,|\chi_\lambda\rangle
```

これをブラ側の内積に代入すると、 AO 同士の内積がそのまま重なり行列 $\langle\chi_\mu|\chi_\lambda\rangle = S^{x}_{\mu\lambda}$ になり、 **1個目の $\mathbf S^x$** が現れる：

```math
\langle \mu | \varphi_i^{h}\rangle = \sum_\lambda \langle\chi_\mu|\chi_\lambda\rangle\,A^{x}_{\lambda i} = \sum_\lambda S^{x}_{\mu\lambda} A^{x}_{\lambda i} = (\mathbf S^{x}\mathbf A^{x})_{\mu i}
```

ケット側の内積も同様に、 AO 同士の重なりから **2個目の $\mathbf S^x$** が現れる：

```math
\langle \varphi_i^{h} | \nu\rangle = \sum_\sigma A^{x*}_{\sigma i}\,\langle\chi_\sigma|\chi_\nu\rangle = \big((\mathbf A^{x})^{\dagger}\mathbf S^{x}\big)_{i\nu}
```

両者を合成し、 $B_i$ を対角に持つ $\mathbf B$ を挟むと、

```math
P^{x}_{\mu\nu} = \sum_i (\mathbf S^{x}\mathbf A^{x})_{\mu i}\,B_i\,\big((\mathbf A^{x})^{\dagger}\mathbf S^{x}\big)_{i\nu} = \big(\mathbf S^{x}\mathbf A^{x}\,\mathbf B\,(\mathbf A^{x})^{\dagger}\mathbf S^{x}\big)_{\mu\nu}
```

となり、行列形 $\mathbf P^x = \mathbf S^{x}\mathbf A^{x}\mathbf B(\mathbf A^{x})^{\dagger}\mathbf S^{x}$ を得る。これが (1) → (1′) の対応で、 $\mathbf S^x$ が2個現れる理由である。

**本質**：演算子 $\hat P=\sum_i B_i|\varphi_i^{h}\rangle\langle\varphi_i^{h}|$ の AO 行列要素 $\langle\chi_\mu|\hat P|\chi_\nu\rangle$ を作るとき、左右の内積 $\langle\chi|\varphi^{h}\rangle$ が各1個の $\mathbf S^x$ を生む。もし AO 基底が正規直交（ $\mathbf S^x=\mathbf I$）であれば $\langle\mu|\varphi_i^{h}\rangle = A^{x}_{\mu i}$ となり、射影項は $\mathbf A^x\mathbf B(\mathbf A^x)^{\dagger}$ で $\mathbf S^x$ は出てこない。すなわち **2個の $\mathbf S^x$ は、 $\mathbf A^x$ が「非直交 AO 基底での展開係数」であることの代償**であり、係数を実際の内積（ブラケット）へ変換するために、ブラ側・ケット側の双方で $\mathbf S^x$ を要する。

| 基底 | $\langle\mu\vert\varphi_i^{h}\rangle$ | 射影項 |
|------|-----------------------------------|--------|
| 正規直交（ $\mathbf S^x=\mathbf I$） | $A^{x}_{\mu i}$ | $\mathbf A^x\mathbf B(\mathbf A^x)^{\dagger}$ （ $\mathbf S^x$ なし） |
| 非直交 AO（一般） | $(\mathbf S^x\mathbf A^x)_{\mu i}$ | $\mathbf S^x\mathbf A^x\mathbf B(\mathbf A^x)^{\dagger}\mathbf S^x$ |

### 補足：なぜ AO 基底は非直交なのか（直交化の可能性）

上で $\mathbf S^x$ が2個現れるのは「 AO 基底が非直交（ $\mathbf S^x\neq\mathbf I$）」だからである。ではなぜ非直交なのか、直交にはできないのか。

**(a) 非直交は原子中心基底の空間的重なりに由来する。** AO $\chi_\mu$ は原子核に中心を持つ関数（ GTO / STO ）であり、重なり行列の要素は

```math
S^{x}_{\mu\nu} = \int \chi_\mu^{*}(\mathbf r)\,\chi_\nu(\mathbf r)\,d\mathbf r
```

で与えられる。異なる原子上の2つの関数は空間的に広がって重なり合うため、一般に $S^{x}_{\mu\nu}\neq 0$（ $\mu\neq\nu$）となる。これは化学結合を生む軌道重なりそのものであり、原理的にゼロにはできない。同一原子上でも、角運動量が異なれば（ s と p など）球面調和関数の直交性により重なりは消えるが、同じ角運動量で動径部が異なる関数（ multiple-zeta など）は一般に重なりをもつ。したがって全体として $\mathbf S^x\neq\mathbf I$ である。

**(b) 直交化は可能だが、 AO の利点を失う。** 対称（ Löwdin ）直交化 $\tilde{\boldsymbol\chi}=\boldsymbol\chi\,\mathbf S^{-1/2}$ や正準直交化により正規直交 AO（ OAO ）を作れば $\mathbf S^x=\mathbf I$ となり、射影項は $\mathbf A^x\mathbf B(\mathbf A^x)^{\dagger}$ に簡約される。しかし直交化された関数 $\tilde\chi_\mu$ はもとの AO の線形結合（非局所的に広がった関数）であり、「特定の原子に属する基底関数」という素直な意味と移植性を失う。さらに $\mathbf S^{-1/2}$ の構成自体に $\mathbf S^x$ を要するため、 $\mathbf S^x$ から完全に逃れられるわけではない。

**(c) 慣例として非直交 AO のまま扱う。** このため標準的な HF / Roothaan–Hall では、非直交 AO 基底のまま $\mathbf S^x$ を陽に保持し（ $\mathbf{FC}=\mathbf{SC}\boldsymbol\varepsilon$）、原子中心関数の利点を生かす。本ノートの $\mathbf A^x$ も実際に計算で用いる非直交 AO 基底での展開係数であり、その帰結として射影項に $\mathbf S^x$ が2個現れる。なお第7節の $S$-正規直交化（式8）は、許容空間の縮約基底 $\mathbf T^x$ に対してこの直交化を施し、縮約 SCF を標準固有値問題に帰着させる操作にあたる。

---

射影項 $\mathbf P^x=\mathbf S^{x}\mathbf A^{x}\mathbf B(\mathbf A^{x})^{\dagger}\mathbf S^{x}$ を Fock 行列に加えて Roothaan–Hall を解くと、
$\{\varphi_i^{h}\}$ 方向の軌道エネルギーが $\sim B_i$ だけ持ち上がり、被占有空間から押し出される。
典型的には $B_i = 10^{6}\,E_{\mathrm h}$ 程度を用いる。

---

## 3. 有限レベルシフトの残差

有限 $B$ では除去は近似的でしかない。被占有 MO に残る $\{\varphi_i^{h}\}$ 射影成分は摂動論的に $O(B^{-1})$ で、
$B$ を上げると $\widetilde{\mathbf F}^x$ の条件数が悪化し、SCF 収束・数値安定性とトレードオフになる。
この残差を厳密にゼロにするのが次の矩形変換である。

---

## 4. 矩形変換による厳密 deflation（式2）

許容空間（ $\{\varphi_i^{h}\}$ と $S$-直交する部分空間）を張る矩形行列
$\mathbf T^x\in\mathbb{R}^{N\times(N-m)}$ を、

```math
(\mathbf A^{x})^{\dagger}\mathbf S^{x}\mathbf T^{x} = \mathbf 0 \qquad (2)
```

を満たすように構成する。すなわち $\mathbf T^x$ の列は、 $m\times N$ 行列 $(\mathbf A^{x})^{\dagger}\mathbf S^{x}$ の
**核 $\ker\big((\mathbf A^{x})^{\dagger}\mathbf S^{x}\big)$（次元 $N-m$）の基底**にとる。

MO 係数を $\mathbf C^x = \mathbf T^x \mathbf C_{\mathrm r}^x$ と縮約パラメータ化すると、任意の $\mathbf C_{\mathrm r}^x$ に対して

```math
(\mathbf A^{x})^{\dagger}\mathbf S^{x}\mathbf C^{x} = (\mathbf A^{x})^{\dagger}\mathbf S^{x}\mathbf T^x\mathbf C_{\mathrm r}^x = \mathbf 0
```

が恒等的に成立する。除去対象は最初から変分空間に含まれない。

---

## 5. 合同変換と射影項の消失（式3）

式 $(1')$ を許容空間へ合同変換すると、

```math
(\mathbf T^x)^{\dagger}\widetilde{\mathbf H}^{x}\mathbf T^x = (\mathbf T^x)^{\dagger}(\mathbf H^x+\mathbf V^x)\mathbf T^x + \big[(\mathbf A^x)^{\dagger}\mathbf S^x\mathbf T^x\big]^{\dagger}\mathbf B\,\big[(\mathbf A^x)^{\dagger}\mathbf S^x\mathbf T^x\big]
```

第2項は式(2)より厳密に $\mathbf 0$。 $\mathbf B$ の大きさに依らず消える。したがって

```math
\mathbf H_{\mathrm r}^{x} = (\mathbf T^x)^{\dagger}\left(\mathbf H^x+\mathbf V^x\right)\mathbf T^x \qquad (3)
```

レベルシフト項は完全に消去され、その役割は拘束条件(2)が引き継ぐ。すなわち
**「巨大 $\mathbf B$ を付与する」操作が「部分空間を陽に取り除く」操作に置き換わる**。

---

## 6. 縮約 SCF 方程式（式4–7）

縮約空間（次元 $N-m$）での一般化固有値問題：

```math
\mathbf F_{\mathrm r}^{x}\mathbf C_{\mathrm r}^{x} = \mathbf S_{\mathrm r}^{x}\mathbf C_{\mathrm r}^{x}\boldsymbol\varepsilon^{x} \qquad (4)
```

各行列は元空間の量を $\mathbf T^x$ で合同変換したもの：

```math
\mathbf S_{\mathrm r}^{x} = (\mathbf T^x)^{\dagger}\mathbf S^x\mathbf T^x \qquad (5)
```

```math
\mathbf F_{\mathrm r}^{x} = (\mathbf T^x)^{\dagger}\big[\mathbf H^x+\mathbf V^x+\mathbf G^x[\mathbf D^x]\big]\mathbf T^x \qquad (6)
```

```math
\mathbf D^x = \mathbf T^x\mathbf D_{\mathrm r}^{x}(\mathbf T^x)^{\dagger} \qquad (7)
```

$\mathbf G^x[\mathbf D^x]$ は二電子部分（Coulomb + exchange）。SCF は縮約空間で密度 $\mathbf D_{\mathrm r}^x$ を更新し、
$\mathbf D^x$ を介して $\mathbf G^x$ を組む。除去軌道は密度に一切寄与しない。

---

## 7. $S$-正規直交化（式8）

$\mathbf T^x$ を $S$-正規直交、すなわち

```math
(\mathbf T^x)^{\dagger}\mathbf S^x\mathbf T^x = \mathbf I
```

にとれば $\mathbf S_{\mathrm r}^x=\mathbf I$ となり、(4) は標準固有値問題に帰着する：

```math
\mathbf F_{\mathrm r}^{x}\mathbf C_{\mathrm r}^{x} = \mathbf C_{\mathrm r}^{x}\boldsymbol\varepsilon^{x} \qquad (8)
```

実装上は $\ker\big((\mathbf A^{x})^{\dagger}\mathbf S^{x}\big)$ の基底に対して $S$-計量での Löwdin 直交化
（あるいは正準直交化）を施して $\mathbf T^x$ を得る。

---

## 8. 次元削減

$sp^3$ 炭素を 1 本切断する場合、局所価電子空間では $4\to 3$。全体では
$\mathbf T^x:\ N\times(N-1)$、変分次元は $N\to N-1$。一般に $m$ 本除去で $N\to N-m$ となり、
Fock 構築後の対角化コストもこの縮約次元で評価される。

---

## 9. レベルシフト法との関係

矩形変換による厳密 deflation は、レベルシフトの $B\to\infty$ 極限と等価である：

```math
\text{矩形変換 (deflation)} \quad\Longleftrightarrow\quad B \to +\infty
```

| | レベルシフト射影 (式1) | 矩形変換 (式3) |
|---|---|---|
| 除去の機構 | $\{\varphi_i^h\}$ に $\sim B_i$ を加え被占有空間から排除 | $\ker((\mathbf A^x)^\dagger\mathbf S^x)$ へ変分空間を制限 |
| 拘束 | $\mathbf P^x=\mathbf S^x\mathbf A^x\mathbf B(\mathbf A^x)^\dagger\mathbf S^x$ を Fock に加算 | $(\mathbf A^x)^\dagger\mathbf S^x\mathbf T^x=\mathbf 0$ |
| 残留成分 | $O(B^{-1})$（近似的） | 機械精度でゼロ（厳密） |
| 次元 | $N$ のまま | $N-m$ に縮約 |
| 数値特性 | 大 $B$ で条件数悪化 | 条件数悪化なし・対角化も縮約次元 |
| 関係 | — | $B\to\infty$ 極限が矩形変換 |

---

## 付録：記号一覧

| 記号 | 意味 |
|------|------|
| $\mathbf H^x,\ \mathbf V^x$ | コアハミルトニアン、環境 ESP |
| $\mathbf S^x$ | AO 重なり行列 |
| $\mathbf A^x$ | 除去混成軌道 $\varphi_i^h$ の AO 係数（ $N\times m$、列ベクトル） |
| $\mathbf B=\mathrm{diag}(B_i)$ | レベルシフト |
| $\mathbf P^x=\mathbf S^x\mathbf A^x\mathbf B(\mathbf A^x)^\dagger\mathbf S^x$ | AO 表現の射影項 |
| $\mathbf T^x$ | 矩形変換行列（ $N\times(N-m)$、 $\ker((\mathbf A^x)^\dagger\mathbf S^x)$ の基底） |
| $\mathbf C^x,\ \mathbf C_{\mathrm r}^x$ | MO 係数（元空間 / 縮約空間） |
| $\mathbf F^x,\ \mathbf G^x[\mathbf D^x]$ | Fock 行列、二電子部分 |
| $\mathbf D^x,\ \mathbf D_{\mathrm r}^x$ | 密度行列（元 / 縮約） |
| $\boldsymbol\varepsilon^x$ | 軌道エネルギー |
| $N,\ m$ | AO 基底数、除去軌道本数 |
