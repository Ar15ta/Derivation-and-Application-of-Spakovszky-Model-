|  属性  |             [[法术表]]              |                                |
| :--: | :------------------------------: | :----------------------------: |
|  关联  | [[Spakovszky稳定性方法——模型Ⅱ——径向无叶空间]] | [[Spakovszky稳定性方法——模型Ⅱ——轴向管道]] |
| 链接解释 |              准备施加在这              |             也施加在这              |
***
# $Model\space Ⅱ$
	一头一尾。

| 符号                                    | 含义                  | 无量纲化                   |
| :------------------------------------ | :------------------ | :--------------------- |
| $\hat{x}$                             | 轴向坐标                | $x/R$                  |
| $s$                                   | 拉普拉斯变量              | $s = \sigma - j\omega$ |
| $n$                                   | 周向波数                | 整数 $\ge 0$             |
| $\hat{\overline{v}}_x, \hat{\overline{v}}_\theta$ | 平均轴向/周向速度           | $v/U$                  |
| $A_n, B_n, C_n$                       | 势模态（向上/向下游衰减）及涡模态系数 | —                      |
| $\mathbb{X}_{\text{sys},n}$           | 系统传递矩阵（3×3）         | —                      |
| $\mathbf{IC}, \mathbf{EC}$            | 进口、出口边界条件矩阵         | —                      |

---
# 1. 边界条件概述
在 $Model \space II$ 中，边界条件用于封闭由组件传递矩阵串联而成的系统方程。根据系统两端（上游进口、下游出口）的物理配置，主要有两类边界条件：

| 边界类型             | 物理情境                    | 数学约束                               |
| :--------------- | :---------------------- | :--------------------------------- |
| **无限长管道**        | 上游、下游管道足够长，远离压缩部件的扰动可忽略 | 上游：无入**射**势波、无向上游传播的涡；下游：无出**射**势波 |
| **下游容腔（Plenum）** | 下游有体积很大的容腔，压力扰动空间均匀     | 容腔入口处：周向非均匀压力扰动为零（$n>0$ 的谐波压力幅值为零） |

此外，在组件内部（如转子‑定子间隙）的**连续性/匹配条件**也属于局部边界条件，但在传递矩阵法中它们被隐式地包含在组件自身的传递矩阵定义中（如转子出口与间隙入口的匹配）。本文主要讨论**系统层面**的边界条件。

# 2. 无限长管道边界条件
## 2.1 物理假设

- 上游管道无限长：扰动源（压缩机）是唯一能量来源，没有来自上游无穷远处的扰动。
- 下游管道无限长：扰动能量必须向下游传播，不能从下游无穷远处反射回来。
- 涡量只能沿主流方向对流，因此上游区域不可能有涡（因为涡产生于叶片排，且不能逆流而上）。
### 2.2 数学约束
对于第 $n$ 阶周向谐波（$n>0$），流函数通解（见轴向管道部分）包含：
- 势模态 $A_n e^{n\hat{x}}$（向上游衰减）、$B_n e^{-n\hat{x}}$（向下游衰减）
- 涡模态 $C_n e^{-k_n\hat{x}}$（向下游对流，$k_n = \frac{s}{\hat{\overline{v}}_x} + j n \frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x}$）

**上游区域**（$x \to -\infty$）：
- 不允许有向下游衰减的势波（$B_n e^{-n\hat{x}}$ 在 $x \to -\infty$ 时发散），故 $B_n^{up}=0$。
- 不允许有涡（涡只能向下游传播），故 $C_n^{up}=0$。

故上游边界条件应为 $B_n^{up} = 0$。同理，下游区域 $x \to +\infty$，$e^{-n x}$ 衰减，$e^{n x}$ 增长，故下游边界条件为 $A_n^{dn} = 0$。此外，下游可以有涡（$C_n^{dn}$ 自由），上游无涡（$C_n^{up}=0$）。所以**常用边界条件**为：

$$
\boxed{
\begin{aligned}
\text{上游（$x \to -\infty$）}&:\quad B_n^{up} = 0,\quad C_n^{up} = 0 \\
\text{下游（$x \to +\infty$）}&:\quad A_n^{dn} = 0
\end{aligned}
}
$$

## 2.3 在传递矩阵中的约束形式
定义 **进口边界条件矩阵** $\mathbf{IC}$ 和 **出口边界条件矩阵** $\mathbf{EC}$，使得：
$$
\mathbf{IC} \cdot \mathbf{c}^{up} = \mathbf{0}, \qquad \mathbf{EC} \cdot \mathbf{c}^{dn} = \mathbf{0}
$$
对于无限长管道：

$$
\mathbf{IC} = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}, \quad
\mathbf{EC} = \begin{bmatrix} 1 & 0 & 0 \end{bmatrix}
$$

因为 $B_n^{up}$ 和 $C_n^{up}$ 为零由 $IC$ 强制，$A_n^{dn}=0$ 由 $EC$ 强制。

结合 $\mathbf{c}^{dn} = \mathbb{X}_{\text{sys},n} \mathbf{c}^{up}$，得到总齐次系统：

$$
\begin{bmatrix} \mathbf{EC} \cdot \mathbb{X}_{\text{sys},n} \\ \mathbf{IC} \end{bmatrix} \cdot \mathbf{c}^{up} = \mathbf{0}
$$

非零解条件：

$$
\boxed{
\det \left( \begin{bmatrix} \mathbf{EC} \cdot \mathbb{X}_{\text{sys},n} \\ \mathbf{IC} \end{bmatrix} \right) = 0
}
$$

这就是特征值方程。
# 3. 下游容腔边界条件（Downstream Plenum）
## 3.1 物理假设
- 下游有一个大容腔（如燃烧室、集气室），其尺寸远大于管道尺寸。
- 容腔内的压力在空间上视为**均匀**（即忽略波动），因此容腔入口处（管道出口）的**非定常压力扰动**必须满足：**周向非均匀部分为零**。
- 注意：$n=0$（轴对称）压力扰动可以存在，它对应容腔整体的压强变化（声学模式）。
## 3.2 数学约束

设管道出口位于 $x = L$，则对于 $n>0$：

$$
\tilde{P}_n(L, s) = 0
$$

对于 $n=0$，无此约束（$P_0$ 自由）。

## 3.3 在传递矩阵中的约束形式
将容腔视为一个特殊的边界元件，直接给出出口系数之间的关系。在论文中，对于有容腔的情况，出口边界条件 $\mathbf{EC}$ 改为：

$$
\mathbf{EC} = \begin{bmatrix} 
\left( -\frac{s}{n} - \hat{\overline{v}}_x - j\hat{\overline{v}}_\theta \right) e^{nL}, &
\left( \frac{s}{n} - \hat{\overline{v}}_x + j\hat{\overline{v}}_\theta \right) e^{-nL}, &
0
\end{bmatrix}
$$

其中 $\hat{\overline{v}}_x, \hat{\overline{v}}_\theta$ 是出口管道的平均流速。这一行向量的物理含义是：它正是通过轴向管道传递矩阵计算出的 $\tilde{P}_n(L)$ 的表达式（忽略涡模态对压力的贡献，因为涡模态不产生压力扰动）。因此设定 $\tilde{P}_n(L)=0$ 即 $\mathbf{EC} \cdot \mathbf{c}^{dn} = 0$。

同时，上游进口边界条件 $\mathbf{IC}$ 仍为 $[0\ 1\ 0;\ 0\ 0\ 1]$（假设上游为无限长管道，无入射波和涡）。于是总齐次方程为：

$$
\begin{bmatrix} \mathbf{EC} \cdot \mathbb{X}_{\text{sys},n} \\ \mathbf{IC} \end{bmatrix} \cdot \mathbf{c}^{up} = 0
$$

行列式为零给出特征值方程。
# 6. 特征值方程的统一形式
无论哪种边界条件，最终都归结为：

$$
\det\!\left( \begin{bmatrix}
\mathbf{EC} \cdot \mathbb{X}_{\text{sys},n}(s) \\
\mathbf{IC}
\end{bmatrix} \right) = 0
$$

解此方程得到的 $s$ 即为系统模态特征值（$\sigma_n$ 为增长率，$\omega_n$ 为旋转角频率）。Model II 采用两种数值方法求解：
# 5. 示例：单级轴流压气机无限长管道边界条件
考虑一个转子——间隙——静子系统，上游为无限长管道，下游为无限长管道。则：

- 上游：$B^{up}=0, C^{up}=0$
- 下游：$A^{dn}=0$

系统传递矩阵 $\mathbb{X}_{\text{sys},n} = \mathbb{T}_{\text{ax},n}(x_4)^{-1} \cdot \mathbb{B}_{\text{sta},n} \cdot \mathbb{B}_{\text{gap},n} \cdot \mathbb{B}_{\text{rot},n} \cdot \mathbb{T}_{\text{ax},n}(x_1)$

构造：

$$
\mathbf{IC} = \begin{bmatrix}0&1&0\\0&0&1\end{bmatrix}, \quad
\mathbf{EC} = \begin{bmatrix}1&0&0\end{bmatrix}
$$

齐次方程：

$$
\begin{bmatrix} \mathbf{EC} \cdot \mathbb{X}_{\text{sys},n} \\ \mathbf{IC} \end{bmatrix} \begin{bmatrix} A^{up} \\ B^{up} \\ C^{up} \end{bmatrix} = 0 
$$

