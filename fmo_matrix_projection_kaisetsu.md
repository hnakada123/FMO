# FMO境界処理における不要混成軌道の除去：レベルシフト射影と矩形変換

> 想定読者：Hartree–Fock / Roothaan–Hall 法（ $\mathbf{FC}=\mathbf{SC}\boldsymbol\varepsilon$、密度行列、一般化固有値問題、射影演算子）を理解している方。
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
**核 $\ker\!\big((\mathbf A^{x})^{\dagger}\mathbf S^{x}\big)$（次元 $N-m$）の基底**にとる。

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

実装上は $\ker\!\big((\mathbf A^{x})^{\dagger}\mathbf S^{x}\big)$ の基底に対して $S$-計量での Löwdin 直交化
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
