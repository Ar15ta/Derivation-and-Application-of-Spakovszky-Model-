
|  属性  |            [[法术表]]             |                                  |                                   |                                |                                |
| :--: | :----------------------------: | :------------------------------: | :-------------------------------: | :----------------------------: | :----------------------------: |
|  关联  |       [[环绕积分法求解非线性特征值]]        |  [[Spakovszky稳定性方法——模型Ⅱ——边界条件]]  |  [[Spakovszky稳定性方法——模型Ⅱ——轴向管道]]   | [[Spakovszky稳定性方法——模型Ⅱ——轴流转子]] | [[Spakovszky稳定性方法——模型Ⅱ——离心转子]] |
| 链接解释 |           求解器（环绕积分法）           |              边界条件组装              |             轴向管道传递矩阵              |            轴流转子执行盘             |            离心转子执行盘             |
|      | [[Spakovszky稳定性方法——模型Ⅱ——轴向静叶]] | [[Spakovszky稳定性方法——模型Ⅱ——径向无叶空间]] | [[Spakovszky稳定性方法——模型Ⅱ——径向带叶扩压器]] |    [[SpakoVszky稳定性方法——模型Ⅰ]]    |                                |
|      |            轴向静叶执行盘             |            径向无叶空间传递矩阵            |            径向有叶扩压器执行盘             |          是$ModelⅠ$的进化          |                                |
***
# $Model\space Ⅱ$ 集大成者

    本文档整合 Model II 的所有部件传递矩阵，给出从输入参数到系统特征值求解的完整流程。重点包括：部件输入/输出定义、矩阵串联方法、边界条件嵌入以及特征值问题的构造。

***
# $n>0$
## 1. 无量纲基准约定
- **特征速度** $U$：取系统中**转子**的特征速度。  
  - 若系统包含**离心转子**：$U = U_2 = \Omega R_2$（叶轮出口切线速度）。  
  - 若系统为纯轴流：$U = \Omega R$（平均半径处切线速度）。  
- **特征长度** $R$：与 $U$ 对应，取同一半径（离心时用 $R_2$，轴流时用平均半径）。  
- **特征密度** $\rho$：进口密度（常数）。  
- 所有无量纲量均基于 $(U,R,\rho)$ 定义。  
  - 坐标：$\hat x = x/R$，$\hat r = r/R$  
  - 时间：$\tau = tU/R$  
  - 速度：$\hat v = v/U$  
  - 压力：$\hat p = p/(\rho U^2)$  
  - 拉普拉斯变量：$s$ 对应 $\partial/\partial\tau \to s$

**注意**：上下游轴向管道内的平均速度 $\hat{\overline{v}}_x, \hat{\overline{v}}_\theta$ 也采用同一 $U$ 无量纲化。若管道内平均速度源于叶轮出口参数，则须预先转换。
## 2. 部件输入与输出参数

每个部件均以**入口状态向量**映射到**出口状态向量**。状态向量的分量按部件类型选取（轴向或径向速度）。
### 2.1 轴向管道（`Axial Duct`）
- **输入**：系数向量 $\mathbf{c}_n = [A_n, B_n, C_n]^{\mathsf T}$（定义在入口 $\hat x=0$ 处）。  
- **输出**：状态向量 $\mathbf{q}_n(\hat x) = [\delta\tilde V_x, \delta\tilde V_\theta, \delta\tilde P]^{\mathsf T}$ 在任意 $\hat x$ 处。  
- **传递矩阵**：$\mathbf{q}_n(\hat x) = \mathbb{T}_{\text{ax},n}(\hat x) \cdot \mathbf{c}_n$（无参考点形式，入口为 0）。  
- **所需参数**：管道长度 $L/R$，背景速度 $\hat{\overline{v}}_x, \hat{\overline{v}}_\theta$（常数）。
### 2.2 轴流转子（`Axial Rotor`）
- **输入**：入口状态向量 $\mathbf{q}_{\text{in}} = [\delta\tilde V_{x1}, \delta\tilde V_{\theta1}, \delta\tilde P_1]^{\mathsf T}$。  
- **输出**：出口状态向量 $\mathbf{q}_{\text{out}} = [\delta\tilde V_{x2}, \delta\tilde V_{\theta2}, \delta\tilde P_2]^{\mathsf T}$。  
- **传递矩阵**：$\mathbf{q}_{\text{out}} = \mathbb{B}_{\text{rot},n}(s) \cdot \mathbf{q}_{\text{in}}$。  
- **所需参数**：
  - 几何：$\tan\beta_2$（出口相对气流角）、$\lambda_{\text{rot}}$（惯性系数）、$\tau_R$（损失滞后）  
  - 背景流：$\hat{\overline{v}}_x, \hat{\overline{v}}_{\theta1}, \hat{\overline{v}}_{\theta2}, \tan\beta_1, \tan\alpha_1$  
  - 损失导数：$\left.\frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\right|_{\text{avg}}$（工作点值）
### 2.3 离心转子（`Centrifugal Impeller`）
- **输入**：入口状态向量 $\mathbf{q}_{\text{in}} = [\delta\tilde V_{x1}, \delta\tilde V_{\theta1}, \delta\tilde P_1]^{\mathsf T}$（轴向）。  
- **输出**：出口状态向量 $\mathbf{q}_{\text{out}} = [\delta\tilde V_{r2}, \delta\tilde V_{\theta2}, \delta\tilde P_2]^{\mathsf T}$（径向）。  
- **传递矩阵**：$\mathbf{q}_{\text{out}} = \mathbb{B}_{\text{imp},n}(s) \cdot \mathbf{q}_{\text{in}}$。  
- **所需参数**：
  - 几何：$AR_{\text{imp}}$（面积-密度比）、$R_1/R_2$、$\tan\beta_2$、$\tan\beta_1$、$\tan\alpha_1$、$\lambda_{\text{imp}}$、$\tau_{\text{imp}}$  
  - 背景流：$\hat{\overline{v}}_{x1}, \hat{\overline{v}}_{r2}, \hat{\overline{v}}_{\theta1}, \hat{\overline{v}}_{\theta2}$  
  - 损失导数：$\left.\frac{d\hat{l}_{\text{imp}}^{ss}}{d\tan\beta_1}\right|_{\text{avg}}$
### 2.4 轴向静叶（`Axial Stator`）
- **输入**：入口状态向量 $\mathbf{q}_{\text{in}} = [\delta\tilde V_{x1}, \delta\tilde V_{\theta1}, \delta\tilde P_1]^{\mathsf T}$。  
- **输出**：出口状态向量 $\mathbf{q}_{\text{out}} = [\delta\tilde V_{x2}, \delta\tilde V_{\theta2}, \delta\tilde P_2]^{\mathsf T}$。  
- **传递矩阵**：$\mathbf{q}_{\text{out}} = \mathbb{B}_{\text{sta},n}(s) \cdot \mathbf{q}_{\text{in}}$。  
- **所需参数**：$\tan\alpha_2$、$\lambda_{\text{sta}}$、$\tau_S$、$\hat{\overline{v}}_x, \hat{\overline{v}}_{\theta1}, \hat{\overline{v}}_{\theta2}, \tan\alpha_1$、损失导数 $\left.\frac{d\hat{l}_S^{ss}}{d\tan\alpha_1}\right|_{\text{avg}}$。
### 2.5 径向无叶空间（`Radial Vaneless Space`）
- **输入**：入口状态向量 $\mathbf{q}_{\text{in}} = [\delta\tilde V_{r1}, \delta\tilde V_{\theta1}, \delta\tilde P_1]^{\mathsf T}$。  
- **输出**：出口状态向量 $\mathbf{q}_{\text{out}} = [\delta\tilde V_{r2}, \delta\tilde V_{\theta2}, \delta\tilde P_2]^{\mathsf T}$。  
- **传递矩阵**：  
$$
  \mathbb{M}_{\text{vls},n} = \mathbb{T}_{\text{rad},n}(r_2) \cdot \mathbb{T}_{\text{rad},n}(r_1)^{-1}, \quad \mathbf{q}_{\text{out}} = \mathbb{M}_{\text{vls},n} \cdot \mathbf{q}_{\text{in}}
  $$  
  其中 $\mathbb{T}_{\text{rad},n}$ 将系数向量 $[D_n,E_n,F_n]^{\mathsf T}$ 映射到状态。  
- **所需参数**：内外半径 $r_1, r_2$、流量常数 $Q$、环量常数 $\Gamma$（由背景自由涡给出）。
### 2.6 径向有叶扩压器（`Radial Vaned Diffuser`）
- **输入**：入口状态向量 $\mathbf{q}_{\text{in}} = [\delta\tilde V_{r1}, \delta\tilde V_{\theta1}, \delta\tilde P_1]^{\mathsf T}$。  
- **输出**：出口状态向量 $\mathbf{q}_{\text{out}} = [\delta\tilde V_{r2}, \delta\tilde V_{\theta2}, \delta\tilde P_2]^{\mathsf T}$。  
- **传递矩阵**：$\mathbf{q}_{\text{out}} = \mathbb{B}_{\text{dif},n}(s) \cdot \mathbf{q}_{\text{in}}$。  
- **所需参数**：
  - 几何：$AR'_{\text{dif}}$（径向投影面积-密度比）、$\tan\alpha_2$、$\lambda_{\text{dif}}$、$\tau_{\text{dif}}$  
  - 背景流：$\hat{\overline{v}}_{r1}, \hat{\overline{v}}_{\theta1}, \hat{\overline{v}}_{r2}, \hat{\overline{v}}_{\theta2}, \tan\alpha_1$  
  - 损失导数：$\left.\frac{d\hat{l}_{\text{dif}}^{ss}}{d\tan\alpha_1}\right|_{\text{avg}}$
### 2.7 边界条件（`Boundary Conditions`）
- **上游无限长管道**：  
  $$
  \mathbf{IC} = \begin{bmatrix}0 & 1 & 0 \\ 0 & 0 & 1\end{bmatrix}, \quad \mathbf{IC} \cdot \mathbf{c}^{\text{up}} = \mathbf{0}
  $$  
  即 $B_n^{\text{up}} = 0,\ C_n^{\text{up}} = 0$。  
- **下游集气室（$Plenum$）**：  
  设出口管道长度为 $L$，背景速度 $\hat{\overline{v}}_x, \hat{\overline{v}}_\theta$，则  
$$
  \mathbf{EC} = \begin{bmatrix}
  \left( -\frac{s}{n} - \hat{\overline{v}}_x - j\hat{\overline{v}}_\theta \right) e^{nL}, &
  \left( \frac{s}{n} - \hat{\overline{v}}_x + j\hat{\overline{v}}_\theta \right) e^{-nL}, &
  0
  \end{bmatrix}
  $$
  条件 $\tilde P_n(L)=0$ 等价于 $\mathbf{EC} \cdot \mathbf{c}^{\text{dn}} = 0$。
## 3. 矩阵堆叠示例（离心压缩系统）

考虑如下串联顺序（从左到右为流动方向）：
$$
\text{[上游无限长管道]} \;\to\; \text{[离心转子]} \;\to\; \text{[径向无叶空间]} \;\to\; \text{[径向有叶扩压器]} \;\to\; \text{[下游轴向管道]} \;\to\; \text{[集气室]}
$$
### 3.1 状态分量映射（关键）
- 有叶扩压器出口的状态向量为 $[\delta\tilde V_r,\ \delta\tilde V_\theta,\ \delta\tilde P]^\mathsf{T}$。
- 下游轴向管道进口的状态向量应为 $[\delta\tilde V_x,\ \delta\tilde V_\theta,\ \delta\tilde P]^\mathsf{T}$。

**近似处理**：忽略 90° 弯管的动力学，直接将径向速度赋给轴向速度：
$$
\delta\tilde V_x^{\text{(duct,in)}} = \delta\tilde V_r^{\text{(diffuser,out)}}
$$

周向速度和压力保持不变。在代码实现中，需在扩压器传递矩阵之后、轴向管道传递矩阵之前插入一个 **置换矩阵** $P_{\text{map}}$：
$$
P_{\text{map}} = \begin{bmatrix}
0 & 1 & 0 \\
0 & 0 & 1 \\
1 & 0 & 0
\end{bmatrix}
$$
其作用是将 $[V_r, V_\theta, P]^\mathsf{T}$ 转换为 $[V_x, V_\theta, P]^\mathsf{T}$（第一行取原第二分量，第二行取原第三分量，第三行取原第一分量）。
### 3.2 各段传递关系（$n\ge1$）
1. **上游管道**（长度 $L_{\text{up}}$）：  
$$
   \mathbf{q}_{\text{rot,in}} = \mathbb{T}_{\text{ax},n}(L_{\text{up}}) \cdot \mathbf{c}^{\text{up}}
   $$

2. **离心转子**：  
   $$
   \mathbf{q}_{\text{rot,out}} = \mathbb{B}_{\text{imp},n}(s) \cdot \mathbf{q}_{\text{rot,in}}
   $$

3. **径向无叶空间**（进口半径 $r_{\text{vls,in}}$，出口半径 $r_{\text{vls,out}}$）：  
   $$
   \mathbf{q}_{\text{vls,out}} = \mathbb{T}_{\text{rad},n}(r_{\text{vls,out}}) \cdot \mathbb{T}_{\text{rad},n}(r_{\text{vls,in}})^{-1} \cdot \mathbf{q}_{\text{rot,out}}
   $$

4. **径向有叶扩压器**：  
   $$
   \mathbf{q}_{\text{dif,out}} = \mathbb{B}_{\text{dif},n}(s) \cdot \mathbf{q}_{\text{vls,out}}
   $$

5. **状态映射**（$V_r \to V_x$）：  
   $$
   \mathbf{q}_{\text{duct,in}} = P_{\text{map}} \cdot \mathbf{q}_{\text{dif,out}}
   $$

6. **下游轴向管道**（长度 $L_{\text{dn}}$）：  
   $$
   \mathbf{c}^{\text{dn}} = \mathbb{T}_{\text{ax},n}(L_{\text{dn}})^{-1} \cdot \mathbf{q}_{\text{duct,in}}
   $$

### 3.3 系统传递矩阵

串联后从上游系数 $\mathbf{c}^{\text{up}}$ 到下游系数 $\mathbf{c}^{\text{dn}}$ 的映射为：
$$
\boxed{
\mathbb{X}_{\text{sys},n}(s) = 
\mathbb{T}_{\text{ax},n}(L_{\text{dn}})^{-1} \cdot P_{\text{map}} \cdot
\mathbb{B}_{\text{dif},n} \cdot
\Bigl[ \mathbb{T}_{\text{rad},n}(r_{\text{vls,out}}) \cdot \mathbb{T}_{\text{rad},n}(r_{\text{vls,in}})^{-1} \Bigr] \cdot
\mathbb{B}_{\text{imp},n} \cdot
\mathbb{T}_{\text{ax},n}(L_{\text{up}})
}
$$
### 3.4 边界条件
- 上游无限长管道：$\mathbf{IC} \cdot \mathbf{c}^{\text{up}} = \mathbf{0}$，$\mathbf{IC} = \begin{bmatrix}0&1&0\\0&0&1\end{bmatrix}$。  
- 下游集气室（$n\ge1$）：出口管道末端压力扰动为零，即 $\tilde P_n(L_{\text{dn}})=0$。  
  该条件等价于 $\mathbf{EC} \cdot \mathbf{c}^{\text{dn}} = 0$，其中
$$
  \mathbf{EC} = \begin{bmatrix}
  \left( -\frac{s}{n} - \hat{\overline{v}}_x - j\hat{\overline{v}}_\theta \right) e^{nL_{\text{dn}}}, &
  \left( \frac{s}{n} - \hat{\overline{v}}_x + j\hat{\overline{v}}_\theta \right) e^{-nL_{\text{dn}}}, &
  0
  \end{bmatrix}
  $$

  物理意义：压力节点产生全反射，波可在管道内来回传播。
### 3.5 特征值问题
结合 $\mathbf{c}^{\text{dn}} = \mathbb{X}_{\text{sys},n} \cdot \mathbf{c}^{\text{up}}$，得齐次系统：

$$
\underbrace{
\begin{bmatrix}
\mathbf{EC} \cdot \mathbb{X}_{\text{sys},n}(s) \\
\mathbf{IC}
\end{bmatrix}
}_{\displaystyle \mathbb{Y}_{\text{sys},n}(s)}
\cdot \mathbf{c}^{\text{up}} = \mathbf{0}
$$
特征值方程：
$$
\boxed{
\det\!\bigl( \mathbb{Y}_{\text{sys},n}(s) \bigr) = 0, \qquad n = 1,2,\dots
}
$$

## 4. 特征值求解方法：围道积分法（Beyn 算法）

采用 **围道积分法** 一次性求解围道内所有特征值，避免逐点扫描。
### 4.1 算法步骤

1. **选择围道** $\Gamma$（通常为圆心 $c$、半径 $R$ 的圆），使内部包含待求特征值，且 $\mathbb{Y}_{\text{sys},n}(z)$ 在 $\Gamma$ 上可逆。

2. **离散围道计算矩矩阵**  
   取 $N$ 个等距点 $z_\ell = c + R e^{i\theta_\ell}$，$\theta_\ell = 2\pi\ell/N$，用梯形公式：
   $$
   S_p = \frac{1}{2\pi i}\oint_\Gamma z^p \,\mathbb{Y}_{\text{sys},n}(z)^{-1} dz
        \approx \frac{1}{N} \sum_{\ell=1}^{N} z_\ell^{\,p} \,\mathbb{Y}_{\text{sys},n}(z_\ell)^{-1}, \qquad p=0,1.
   $$

3. **奇异值分解**：$S_0 = U \Sigma V^H$。根据奇异值阈值确定有效秩 $r$，取 $U_r$ 为 $U$ 的前 $r$ 列，$V_r$ 为 $V$ 的前 $r$ 列。

4. **投影**：$A_0 = U_r^H S_0 V_r$，$A_1 = U_r^H S_1 V_r$。

5. **求解小型广义特征值问题**：$A_1 \boldsymbol{x} = \lambda A_0 \boldsymbol{x}$。特征值 $\lambda$ 即为围道内的原特征值 $s_j$。
### 4.2 注意事项
- 围道圆心 $c$ 应位于待求特征值分布中心，半径 $R$ 足够大以包络它们。
- 积分点数 $N$ 通常取 $32\sim 64$，梯形公式具有谱收敛性。
- 计算 $\mathbb{Y}_{\text{sys},n}(z_\ell)^{-1}$ 时需使用稳定的线性求解器（LU 分解）。
- 若围道内无特征值，$S_0$ 的数值秩为零。
## 5. 关于 $n=0$ 的说明
本版本仅处理 $n\ge1$ 的非轴对称模态。$n=0$（轴对称喘振模态）需采用不同的集总参数模型（如 Moore‑Greitzer），其边界条件也不同（集气室允许整体压力脉动）。相关实现留作占位符，后续补充。