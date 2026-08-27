|  属性  |            [[法术表]]             |
| :--: | :----------------------------: |
|  关联  | [[Spakovszky稳定性方法——模型Ⅱ——轴向管道]] |
| 链接解释 |            我来组成头部！             |
***
# $Model\space Ⅱ$
	二维、定常、不可压缩、无粘流动。
`在柱坐标系`$(r,\theta)$ `中，径向空间（如离心压气机叶轮与扩压器之间的无叶空间）内的平均背景流动。假设流场轴对称`（$\partial/\partial\theta = 0$）`且无旋`（$\nabla \times \mathbf{V}=0$）`，壁面效应忽略。

## 主要变量与符号

| 变量名 | 符号 | 单位 | 量纲 | 说明 |
| :-: | :-: | :-: | :-: | :-: |
| 径向速度 | $v_r$ | $\text{m/s}$ | $L T^{-1}$ | 沿半径方向的速度分量 |
| 周向速度 | $v_\theta$ | $\text{m/s}$ | $L T^{-1}$ | 垂直于半径的切向速度 |
| 半径 | $r$ | $\text{m}$ | $L$ | 径向坐标 |
| 密度 | $\rho$ | $\text{kg/m}^3$ | $M L^{-3}$ | 常数（不可压） |
| 环量 | $\Gamma$ | $\text{m}^2/\text{s}$ | $L^2 T^{-1}$ | 绕某闭合曲线的速度环量，此处为常数 |
| 流量常数 | $Q$ | $\text{m}^2/\text{s}$ | $L^2 T^{-1}$ | 单位高度上的体积流量 |

***
# 1. 建模部分
## 1.1 自由涡背景
### 1.1.1 理想流动的速度分布推导
在无叶空间内，流动为**二维、定常、不可压、无粘、轴对称**。柱坐标系 $(r,\theta)$ 中，速度场 $\mathbf{v} = (v_r(r), v_\theta(r), 0)$。
- **连续方程（不可压）**：
$$
  \nabla \cdot \mathbf{v} = \frac{1}{r}\frac{\partial (r v_r)}{\partial r} + \frac{1}{r}\frac{\partial v_\theta}{\partial \theta} + \frac{\partial v_z}{\partial z} = 0
$$
轴对称且无轴向流动（$v_z=0$，$\partial/\partial\theta=0$）简化为：
$$
  \frac{1}{r}\frac{d}{dr}(r v_r) = 0 \quad \Longrightarrow \quad r v_r = \text{常数} \tag{1}
$$
定义该常数为 $Q$（流量常数），量纲为 $\text{m}^2/\text{s}$：
$$
  \boxed{v_r(r) = \frac{Q}{r}} \tag{2}
$$
- **无旋条件**：这一假设其实可以由上游管道里推导过的涡量输运方程体现：这两者都是二维，无黏的。但是这里的假设要比上游管道里的那个均匀背景假设弱一些。

自由涡的核心特征是**涡量为零**（$\nabla \times \mathbf{v}=0$）。在二维极坐标中，只有 $z$ 方向涡量分量：
$$
  \omega_z = \frac{1}{r}\frac{\partial (r v_\theta)}{\partial r} - \frac{1}{r}\frac{\partial v_r}{\partial \theta} = 0
$$
由于轴对称（$\partial v_r/\partial \theta = 0$），无旋条件简化为：
$$
  \frac{1}{r}\frac{d}{dr}(r v_\theta) = 0 \quad \Longrightarrow \quad r v_\theta = \text{常数} \tag{3}
$$
定义该常数为 $\Gamma$（环量常数），量纲为 $\text{m}^2/\text{s}$：
$$
  \boxed{v_\theta(r) = \frac{\Gamma}{r}} \tag{4}
$$
注意：$\Gamma$ 正是绕半径为 $r$ 的圆周的速度环量 $\oint \mathbf{v}\cdot d\mathbf{l} = 2\pi r v_\theta = 2\pi \Gamma$，与 $r$ 无关。

这是一个无源无旋流！
### 1.1.2 动量方程检验
对于上述速度分布，需验证其是否满足无粘流动的欧拉方程（径向与周向分量）。不可压、定常、无体积力时，柱坐标下的欧拉方程为：
- **径向分量**：
$$
  v_r\frac{\partial v_r}{\partial r} + \frac{v_\theta}{r}\frac{\partial v_r}{\partial \theta} - \frac{v_\theta^2}{r} = -\frac{1}{\rho}\frac{\partial p}{\partial r}
$$
由于轴对称，$\partial v_r/\partial \theta = 0$，代入 $v_r = Q/r$，$v_\theta = \Gamma/r$：
$$
  \frac{Q}{r}\cdot \frac{d}{dr}\left(\frac{Q}{r}\right) - \frac{1}{r}\left(\frac{\Gamma}{r}\right)^2 = -\frac{1}{\rho}\frac{d p}{d r}
$$
计算导数：
$$
  \frac{d}{dr}\left(\frac{Q}{r}\right) = -\frac{Q}{r^2}
$$
因此：
$$
  \frac{Q}{r}\cdot\left(-\frac{Q}{r^2}\right) - \frac{\Gamma^2}{r^3} = -\frac{Q^2}{r^3} - \frac{\Gamma^2}{r^3} = -\frac{1}{\rho}\frac{d p}{d r}
$$
即：
$$
  \frac{d p}{d r} = \rho \frac{Q^2 + \Gamma^2}{r^3} \tag{5}
$$
- **周向分量**：
$$
  v_r\frac{\partial v_\theta}{\partial r} + \frac{v_\theta}{r}\frac{\partial v_\theta}{\partial \theta} + \frac{v_r v_\theta}{r} = -\frac{1}{\rho r}\frac{\partial p}{\partial \theta}
$$
轴对称下 $\partial p/\partial \theta = 0$，且 $\partial v_\theta/\partial \theta = 0$，代入得：
  $$
  v_r\frac{d v_\theta}{d r} + \frac{v_r v_\theta}{r} = 0
  $$
将 $v_r = Q/r$，$v_\theta = \Gamma/r$ 代入：
  $$
  \frac{Q}{r}\cdot\left(-\frac{\Gamma}{r^2}\right) + \frac{Q}{r}\cdot\frac{\Gamma}{r}\cdot\frac{1}{r} = -\frac{Q\Gamma}{r^3} + \frac{Q\Gamma}{r^3} = 0
  $$
自动满足。

因此，式 (2) 和 (4) 同时满足连续方程、无旋条件和欧拉方程的周向分量，并导出一个可积分的径向压力梯度（式(5)）。积分可得压力分布：
$$
p(r) = p_0 - \frac{1}{2}\rho \frac{Q^2+\Gamma^2}{r^2}
$$
其中 $p_0$ 为无穷远处的参考压力。压力随半径减小而升高（向心流动减速增压），符合扩压器特性。

### 1.1.3 $\Gamma$ 与 $Q$ 的物理意义

- **$Q$（流量常数）**：由质量守恒决定。考虑无叶空间的高度为 $b$（垂直于 $r\theta$ 平面），则通过半径为 $r$ 的圆柱面的体积流量为：
$$
  \dot q_{v} = (2\pi r b) \cdot v_r(r) = 2\pi b Q
$$
因此 $Q = \dot q_{v} / (2\pi b)$，即 **单位高度上的流量（$b=1$ 时的流量）**。对于给定的工作点，$Q$ 由上游叶轮出口的径向流量决定。

- **$\Gamma$（环量常数）**：由叶轮出口的切向速度分布决定。叶轮对流体做功，赋予流体周向动量。对于无叶空间入口（$r = r_2$，叶轮出口半径），平均周向速度 $v_{\theta 2}$ 可通过速度三角形和滑移因子计算。则 $\Gamma = r_2 v_{\theta 2}$。在无叶空间中，由于无粘且无外力矩，角动量守恒（即 $r v_\theta = \text{常数}$），因此 $\Gamma$ 沿径向不变。这也可以从周向欧拉方程直接推导：定常、无粘、轴对称下，周向动量方程为 $v_r \frac{d (r v_\theta)}{dr} = 0$，故 $r v_\theta = \text{常数}$。

自由涡的速度分布是一个对数螺旋线：
![[Spakovszky稳定性方法——模型Ⅱ——径向空间配图.png]]
## 1.2 得到小扰动
### 1.2.1 线性化
将总流场分解为背景（上划线）与小扰动（$\delta$）：
$$
v_r(r,\theta,t) = \overline{v}_r(r) + \delta v_r(r,\theta,t), \quad
v_\theta(r,\theta,t) = \overline{v}_\theta(r) + \delta v_\theta(r,\theta,t), \quad
p(r,\theta,t) = \overline{p}(r) + \delta p(r,\theta,t)
$$
涡量也分解为背景涡量 $\overline{\omega}_z$ 与扰动 $\delta\zeta$。由于背景无旋，$\overline{\omega}_z = 0$，故：
$$
\omega_z = \delta\zeta
$$
代入之前已经推导过的涡量输运方程：
$$
\frac{\partial \omega_z}{\partial t} + v_r \frac{\partial \omega_z}{\partial r} + \frac{v_\theta}{r} \frac{\partial \omega_z}{\partial \theta} = 0
$$

得到：
$$
\frac{\partial \delta\zeta}{\partial t}
+ \left( \overline{v}_r + \delta v_r \right) \frac{\partial \delta\zeta}{\partial r}
+ \frac{\overline{v}_\theta + \delta v_\theta}{r} \frac{\partial \delta\zeta}{\partial \theta} = 0
$$
忽略二阶小量（$\delta v_r \partial \delta\zeta/\partial r$ 等），得到线性化方程：
$$
\frac{\partial \delta\zeta}{\partial t} + \overline{v}_r \frac{\partial \delta\zeta}{\partial r} + \frac{\overline{v}_\theta}{r} \frac{\partial \delta\zeta}{\partial \theta} = 0 
$$
也可以写成：
$$
\frac{\partial \delta\zeta}{\partial t} + \frac{Q}{r}\frac{\partial \delta\zeta}{\partial r} + \frac{\Gamma}{r^2} \frac{\partial \delta\zeta}{\partial \theta} = 0 \tag{6}
$$
然后，就该用流函数表示涡量并线性化了。

这里我应该说，流函数就是一个超级换元工具。原本的封闭系统由不可压缩连续性$+$二维动量方程解两个速度一个压力，流函数引入后把两个速度压缩成一个变量（在连续性方程的解里面继续找动量方程的解），然后涡量方程（通过取旋度，压力梯度消去了）把压力也换走，整个系统用阶次换变量，变得好解了。

我们需要把涡量用流函数表达，首先给出这里的流函数定义：
$$
v_r = \frac{1}{r}\frac{\partial \psi}{\partial \theta}, \qquad v_\theta = -\frac{\partial \psi}{\partial r}
$$
将总流函数分解：$\psi = \overline{\psi}(r) + \delta\psi(r,\theta,t)$。背景流函数 $\overline{\psi}$ 可由背景速度积分得到（但不必要）。

由涡量定义，唯一还存活的分量是：
$$
\omega_z = \frac{1}{r}\frac{\partial (r v_\theta)}{\partial r} - \frac{1}{r}\frac{\partial v_r}{\partial \theta}
$$

将 $v_r, v_\theta$ 用 $\psi$ 表示，并线性化。对于扰动部分，由于背景无旋（$\overline{\omega}_z=0$），有：
$$
\delta\zeta = \frac{1}{r}\frac{\partial}{\partial r}\left( -r \frac{\partial \delta\psi}{\partial r} \right) - \frac{1}{r^2}\frac{\partial^2 \delta\psi}{\partial \theta^2}
$$
### 1.2.2 无量纲化
引入特征量，具体定义和上游管道中保持一致：
- 特征速度 $U$
- 特征半径 $R$
- 特征密度 $\rho$
- 特征时间 $R/U$

定义**无量纲变量**：
$$
\hat{r} = \frac{r}{R}, \quad \tau = \frac{t}{R/U} = \frac{tU}{R}, \quad \hat{\overline{v}}_r = \frac{\overline{v}_r}{U}, \quad \hat{\overline{v}}_\theta = \frac{\overline{v}_\theta}{U}, \quad \delta\hat{p} = \frac{\delta p}{\rho U^2}
$$
$$
\delta\hat{\psi} = \frac{\delta\psi}{U R}, \quad \delta\hat{\zeta} = \frac{\delta\zeta R}{U}, \quad \hat{Q} = \frac{Q}{U R}, \quad \hat{\Gamma} = \frac{\Gamma}{U R}
$$
背景自由涡的无量纲形式：
$$
\hat{\overline{v}}_r = \frac{\hat{Q}}{\hat{r}}, \qquad \hat{\overline{v}}_\theta = \frac{\hat{\Gamma}}{\hat{r}}
$$
将式(6)两边乘以 $R/U^2$：
$$
\frac{\partial \delta\hat{\zeta}}{\partial \tau} + \frac{\hat{Q}}{\hat{r}} \frac{\partial \delta\hat{\zeta}}{\partial \hat{r}} + \frac{\hat{\Gamma}}{\hat{r}^2} \frac{\partial \delta\hat{\zeta}}{\partial \hat{\theta}} = 0 \tag{7}
$$
这就是我们的无量纲线性化涡量方程。

同理，涡量和流函数的关系：
$$
-\delta\hat\zeta = \nabla^2 \delta\hat\psi, \quad \text{其中} \quad \nabla^2 = \frac{1}{\hat r}\frac{\partial}{\partial \hat r}\left( \hat r \frac{\partial}{\partial \hat r} \right) + \frac{1}{\hat r^2}\frac{\partial^2}{\partial \hat\theta^2} \tag{8}
$$

建模部分到此结束。
***
# 2. 求解部分
## 2.1时频变换
由式(7)，对时间 $\tau$ 作单边拉普拉斯变换，记 $s$ 为无量纲拉普拉斯变量，$\delta\hat{Z}(\hat{r},\hat{\theta},s) = \mathcal{L}\{\delta\hat{\zeta}(\hat{r},\hat{\theta},\tau)\}$。
$$
s \delta\hat{Z} + \frac{\hat{Q}}{\hat{r}} \frac{\partial \delta\hat{Z}}{\partial \hat{r}} + \frac{\hat{\Gamma}}{\hat{r}^2} \frac{\partial \delta\hat{Z}}{\partial \hat{\theta}} = 0 \tag{9}
$$
由于流动在周向具有周期性（$2\pi$ 周期），将 $\delta\hat{Z}$ 展开为傅里叶级数：
$$
\delta\hat{Z}(\hat{r},\hat{\theta},s) = \sum_{n=0}^{\infty} \delta\tilde{Z}_n(\hat{r},s) e^{j n \hat{\theta}} \tag{10}
$$
代入式(9)，利用各谐波正交性，对每个 $n$ 得到常微分方程：
$$
s \delta\tilde{Z}_n + \frac{\hat{Q}}{\hat{r}} \frac{d \delta\tilde{Z}_n}{d \hat{r}} + j n \frac{\hat{\Gamma}}{\hat{r}^2} \delta\tilde{Z}_n = 0, \quad n \ge 0 \tag{11}
$$
类似地，拉普拉斯变换后，无量纲流函数扰动 $\delta\hat{\Psi}(\hat{r},\hat{\theta},s) = \mathcal{L}\{\delta\hat{\psi}\}$ 与速度扰动的关系可仍沿用时域定义（拉普拉斯变换与空间导数可交换）：
$$
\delta\hat{V}_r = \frac{1}{\hat{r}} \frac{\partial \delta\hat{\Psi}}{\partial \hat{\theta}}, \qquad
\delta\hat{V}_\theta = -\frac{\partial \delta\hat{\Psi}}{\partial \hat{r}}
$$
其中 $\delta\hat{V}_r = \mathcal{L}\{\delta\hat{v}_r\}$，$\delta\hat{V}_\theta = \mathcal{L}\{\delta\hat{v}_\theta\}$。

涡量扰动 $\delta\hat{Z}$ 与 $\delta\hat{\Psi}$ 的关系（由式(8)拉普拉斯变换得到）为：
$$
\delta\hat{Z} = -\nabla^2 \delta\hat{\Psi}, \qquad
\nabla^2 = \frac{1}{\hat{r}} \frac{\partial}{\partial \hat{r}} \left( \hat{r} \frac{\partial}{\partial \hat{r}} \right) + \frac{1}{\hat{r}^2} \frac{\partial^2}{\partial \hat{\theta}^2}
$$
## 2.2 求解涡量扰动 $\delta\tilde{Z}_n$
对每个谐波 $n \ge 0$，方程 (11) 为一阶线性齐次 ODE：
$$
s \delta\tilde{Z}_n + \frac{\hat{Q}}{\hat{r}} \frac{d \delta\tilde{Z}_n}{d \hat{r}} + j n \frac{\hat{\Gamma}}{\hat{r}^2} \delta\tilde{Z}_n = 0
$$
整理为：
$$
\frac{d \delta\tilde{Z}_n}{d \hat{r}} + \left( \frac{s \hat{r}}{\hat{Q}} + j n \frac{\hat{\Gamma}}{\hat{Q} \hat{r}} \right) \delta\tilde{Z}_n = 0
$$
分离变量并积分（从某个参考半径 $\hat{r}_0$ 到 $\hat{r}$）：
$$
\int_{\delta\tilde{Z}_n(\hat{r}_0)}^{\delta\tilde{Z}_n(\hat{r})} \frac{d \delta\tilde{Z}_n}{\delta\tilde{Z}_n}
= -\int_{\hat{r}_0}^{\hat{r}} \left( \frac{s \xi}{\hat{Q}} + j n \frac{\hat{\Gamma}}{\hat{Q} \xi} \right) d\xi
$$
得：
$$
\ln \frac{\delta\tilde{Z}_n(\hat{r})}{\delta\tilde{Z}_n(\hat{r}_0)} = -\frac{s}{2\hat{Q}} (\hat{r}^2 - \hat{r}_0^2) - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \frac{\hat{r}}{\hat{r}_0}
$$
因此：
$$
\delta\tilde{Z}_n(\hat{r},s) = \delta\tilde{Z}_n(\hat{r}_0,s) \; \exp\!\left( -\frac{s}{2\hat{Q}} (\hat{r}^2 - \hat{r}_0^2) - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \frac{\hat{r}}{\hat{r}_0} \right)
$$
记 $\tilde{F}_n(s) = \delta\tilde{Z}_n(\hat{r}_0,s) \, e^{\frac{s}{2\hat{Q}} \hat{r}_0^2 + j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \hat{r}_0}$，则：

$$
\delta\tilde{Z}_n(\hat{r},s) = \tilde{F}_n(s) \; \exp\!\left( -\frac{s}{2\hat{Q}} \hat{r}^2 - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \hat{r} \right) \tag{12}
    $$
该式适用于所有 $n \ge 0$。对于 $n=0$，指数中 $j n$ 项消失，变成：
$$
\delta\tilde{Z}_0(\hat{r},s) = \tilde{F}_0(s) \; \exp\!\left( -\frac{s}{2\hat{Q}} \hat{r}^2 \right)
$$
于是：
$$
\boxed{\delta\hat{Z}(\hat{r},\hat{\theta},s) = \sum_{n=0}^{\infty}\tilde{F}_n(s) \; \exp\!\left( -\frac{s}{2\hat{Q}} \hat{r}^2 - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \hat{r} + j n \hat\theta\right)} 
$$
## 2.3 求解流函数泊松方程
我们从频域泊松方程出发：
$$
\nabla^2 \delta\hat{\Psi} = -\delta\hat{Z}, \qquad
\nabla^2 = \frac{1}{\hat{r}} \frac{\partial}{\partial \hat{r}} \left( \hat{r} \frac{\partial}{\partial \hat{r}} \right) + \frac{1}{\hat{r}^2} \frac{\partial^2}{\partial \hat{\theta}^2}
$$
其中 $\delta\hat{Z}(\hat{r},\hat{\theta},s)$ 由式 (12) 给出。将方程分解为齐次部分和非齐次部分，分别求解。
### 2.3.1 齐次解
$$
\frac{1}{\hat{r}}\frac{\partial}{\partial \hat{r}}\left( \hat{r} \frac{\partial \delta\hat{\Psi}_h}{\partial \hat{r}} \right) + \frac{1}{\hat{r}^2}\frac{\partial^2 \delta\hat{\Psi}_h}{\partial \hat{\theta}^2} = 0 
$$
由分离变量法，齐次解为：
$$
\begin{aligned}
\delta\hat{\Psi}_h(\hat{r},\hat{\theta}) = &\alpha_0 + \alpha_1 \ln \hat{r} + \alpha_2 \hat{\theta} + \alpha_3 \hat{\theta} \ln \hat{r} \\
&+ \sum_{n=1}^{\infty} \left( \tilde{D}_n \hat{r}^n + \tilde{E}_n \hat{r}^{-n} \right)
\left( \delta\tilde{A}_n \cos n\hat{\theta} + \delta\tilde{B}_n \sin n\hat{\theta} \right)
\end{aligned}
$$
这里常数项由于解的任意性可以弃掉，两个对数项会导致后面积分出来的压力多值，所以也弃掉。最后角度项虽然本身导致流函数多值性，但它贡献了 $n=0$ 的径向速度场，不能丢。
### 2.3.2 特解（常数变易法）
对于每个谐波 $n$，设 $\delta\hat{\Psi}_p = \sum_{n=0}^{\infty} \delta\tilde{\Psi}_{p,n}(\hat{r},s) e^{jn\hat{\theta}}$，代入得：
$$
\frac{1}{\hat{r}} \frac{d}{d\hat{r}}\left( \hat{r} \frac{d \delta\tilde{\Psi}_{p,n}}{d\hat{r}} \right) - \frac{n^2}{\hat{r}^2} \delta\tilde{\Psi}_{p,n} = -\tilde{F}_n(s) e^{-\frac{s}{2\hat{Q}} \hat{r}^2 - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \hat{r}}
$$

> [!NOTE]
> 常数变易法给出特解：
> $$
> y_p = y_1 \int \frac{-y_2 f}{W} d\hat{r} + y_2 \int \frac{y_1 f}{W} d\hat{r}
> $$
> 其中 $f = -\delta\tilde{Z}_n$，$W = y_1 y_2' - y_2 y_1'$。
- 对于 $n=0$：
$$
\delta\tilde{\Psi}_{p,0} = \tilde{F}_0(s) \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2} \cdot \xi\,(\ln \hat{r} - \ln \xi) \, d\xi$$
- 对于 $n>0$：
  $$
  \delta\tilde{\Psi}_{p,n} = \tilde{F}_n(s) \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}} \xi^2 - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \xi} \left( \hat{r}^n \xi^{-n+1} - \hat{r}^{-n} \xi^{n+1} \right) d\xi
  $$
记积分结果为 $\tilde{\mathcal{R}}_n(\hat{r},s)$，则 $\delta\tilde{\Psi}_{p,n} = \tilde{F}_n(s) \tilde{\mathcal{R}}_n(\hat{r},s)$。
### 2.3.3 叠加齐次解与特解
完整流函数（频域）：
$$
\delta\tilde{\Psi}_n(\hat{r},s) = \delta\tilde{\Psi}_{h,n}(\hat{r},s) + \tilde{F}_n(s) \tilde{\mathcal{R}}_n(\hat{r},s)
$$
- 当 $n>0$：$\delta\tilde{\Psi}_{h,n} = \tilde{D}_n(s) \hat{r}^n + \tilde{E}_n(s) \hat{r}^{-n}$；
- 当 $n=0$：$\delta\tilde{\Psi}_{h,0} = \tilde{D}_0(s)\hat\theta$。
因此：
$$
\delta\tilde{\Psi}_n(\hat{r},s) = 
\begin{cases}
\tilde{D}_0(s)\hat \theta + \tilde{F}_0(s) \delta\tilde{\mathcal{R}}_0(\hat{r},s), & n=0 \\
\tilde{D}_n(s) \hat{r}^n + \tilde{E}_n(s) \hat{r}^{-n} + \tilde{F}_n(s) \tilde{\mathcal{R}}_n(\hat{r},s), & n>0
\end{cases}
$$
最后有[^1]：
$$
\begin{aligned}
\delta\hat{\Psi}(\hat{r},\hat{\theta},s) = &\ \tilde{D}_0(s)\,\hat{\theta} \\
&+ \sum_{n=1}^{\infty} \Bigl[ \tilde{D}_n(s)\,\hat{r}^{n} + \tilde{E}_n(s)\,\hat{r}^{-n} \Bigr] e^{j n \hat{\theta}} \\
&+ \tilde{F}_0(s) \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2} \cdot \xi\,(\ln \hat{r} - \ln \xi) \, d\xi\\
&+ \sum_{n=1}^{\infty} \tilde{F}_n(s) \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}} \xi^2 - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \xi} \left( \hat{r}^n \xi^{-n+1} - \hat{r}^{-n} \xi^{n+1} \right) d\xi \; e^{j n \hat{\theta}}
\end{aligned}\tag{13}
$$
## 2.4 由流函数扰动导出速度扰动
扰动速度与扰动流函数的关系为：
$$
\delta\hat{V}_r = \frac{1}{\hat{r}}\frac{\partial (\delta\hat{\Psi})}{\partial \hat{\theta}}, \qquad
\delta\hat{V}_\theta = -\frac{\partial (\delta\hat{\Psi})}{\partial \hat{r}}
$$
频域中 $\delta\hat{\Psi}(\hat{r},\hat{\theta},s) = \mathcal{L}\{\delta\hat{\psi}\}$：
$$
\begin{aligned}
\delta\hat{\Psi}(\hat{r},\hat{\theta},s) = &\ \tilde{D}_0(s)\,\hat{\theta} \\
&+ \sum_{n=1}^{\infty} \Bigl[ \tilde{D}_n(s)\,\hat{r}^{n} + \tilde{E}_n(s)\,\hat{r}^{-n} \Bigr] e^{j n \hat{\theta}} \\
&+ \tilde{F}_0(s)\,\tilde{\mathcal{I}}_0(\hat{r},s) \\
&+ \sum_{n=1}^{\infty} \tilde{F}_n(s)\,\tilde{\mathcal{R}}_n(\hat{r},s)\; e^{j n \hat{\theta}}
\end{aligned}
$$
其中辅助函数：
$$
\tilde{\mathcal{I}}_0(\hat{r},s) = \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2} \cdot \xi\,(\ln \hat{r} - \ln \xi) \, d\xi,
$$
$$
\tilde{\mathcal{R}}_n(\hat{r},s) = \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2 - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \xi} \left( \hat{r}^{n} \xi^{-n+1} - \hat{r}^{-n} \xi^{n+1} \right) d\xi.
$$
### 2.4.1 径向速度扰动 $\delta\hat{V}_r$
由 $\delta\hat{V}_r = \frac{1}{\hat{r}}\frac{\partial \delta\hat{\Psi}}{\partial \hat{\theta}}$，先计算 $\partial\delta\hat{\Psi}/\partial\hat{\theta}$。
- **第一项** $\tilde{D}_0(s)\,\hat{\theta}$：对 $\hat{\theta}$ 导数为 $\tilde{D}_0(s)$。
- **第二项**（$n\ge1$）：$\bigl[\tilde{D}_n \hat{r}^{n}+\tilde{E}_n \hat{r}^{-n}\bigr] e^{jn\hat{\theta}}$ 的导数为 $j n \bigl[\tilde{D}_n \hat{r}^{n}+\tilde{E}_n \hat{r}^{-n}\bigr] e^{jn\hat{\theta}}$。
- **第三项** $\tilde{F}_0 \delta\tilde{\mathcal{I}}_0$：与 $\hat{\theta}$ 无关，导数为 0。
- **第四项**（$n\ge1$）：$\tilde{F}_n \tilde{\mathcal{R}}_n e^{jn\hat{\theta}}$ 的导数为 $j n \tilde{F}_n \tilde{\mathcal{R}}_n e^{jn\hat{\theta}}$。

于是：
$$
\frac{\partial \delta\hat{\Psi}}{\partial \hat{\theta}} = \tilde{D}_0(s) + \sum_{n=1}^{\infty} j n \Bigl[ \tilde{D}_n(s) \hat{r}^{n} + \tilde{E}_n(s) \hat{r}^{-n} + \tilde{F}_n(s) \tilde{\mathcal{R}}_n(\hat{r},s) \Bigr] e^{j n \hat{\theta}}
$$

除以 $\hat{r}$ 即得径向速度扰动：
$$
\boxed{
\delta\hat{V}_r(\hat{r},\hat{\theta},s) = \frac{\tilde{D}_0(s)}{\hat{r}} + \sum_{n=1}^{\infty} \frac{j n}{\hat{r}} \Bigl[ \tilde{D}_n(s) \hat{r}^{n} + \tilde{E}_n(s) \hat{r}^{-n} + \tilde{F}_n(s) \tilde{\mathcal{R}}_n(\hat{r},s) \Bigr] e^{j n \hat{\theta}}
} \tag{14}
$$
### 2.4.2 周向速度扰动 $\delta\hat{V}_\theta$
由 $\delta\hat{V}_\theta = -\frac{\partial \delta\hat{\Psi}}{\partial \hat{r}}$，计算 $\partial\delta\hat{\Psi}/\partial \hat{r}$ 逐项求导。
**(a) 第一项** $\tilde{D}_0(s)\,\hat{\theta}$：对 $\hat{r}$ 导数为 0。
**(b) 第二项（$n\ge1$）**：
$$
\frac{\partial}{\partial \hat{r}}\Bigl[ \bigl(\tilde{D}_n \hat{r}^{n}+\tilde{E}_n \hat{r}^{-n}\bigr) e^{jn\hat{\theta}} \Bigr] = \bigl( n \tilde{D}_n \hat{r}^{n-1} - n \tilde{E}_n \hat{r}^{-n-1} \bigr) e^{jn\hat{\theta}}.
$$
**(c) 第三项（$n=0$ 特解）**：
$$
\frac{\partial}{\partial \hat{r}}\bigl[ \tilde{F}_0 \tilde{\mathcal{I}}_0 \bigr] = \tilde{F}_0 \frac{\partial \tilde{\mathcal{I}}_0}{\partial \hat{r}}.
$$
由莱布尼茨法则：
$$
\frac{\partial \delta\tilde{\mathcal{I}}_0}{\partial \hat{r}} = e^{-\frac{s}{2\hat{Q}}\hat{r}^2} \hat{r} (\ln \hat{r} - \ln \hat{r}) + \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2} \xi \cdot \frac{1}{\hat{r}} d\xi = \frac{1}{\hat{r}} \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2} \xi d\xi \equiv \frac{\tilde{J}_0(\hat{r},s)}{\hat{r}},
$$
其中定义
$$
\tilde{J}_0(\hat{r},s) = \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2} \xi d\xi = \frac{\hat{Q}}{s}\left( e^{-\frac{s}{2\hat{Q}}\hat{r}_0^2} - e^{-\frac{s}{2\hat{Q}}\hat{r}^2} \right) \quad (s \neq 0).
$$
**(d) 第四项（$n\ge1$ 特解）**：
$$
\frac{\partial}{\partial \hat{r}}\bigl[ \tilde{F}_n \tilde{\mathcal{R}}_n e^{jn\hat{\theta}} \bigr] = \tilde{F}_n \frac{\partial \tilde{\mathcal{R}}_n}{\partial \hat{r}} e^{jn\hat{\theta}}.
$$
$\tilde{\mathcal{R}}_n$ 的定义中，被积函数显含 $\hat{r}$ 于括号内。求导时，上限代入的边界项相消（因为 $\hat{r}^n \hat{r}^{-n+1} - \hat{r}^{-n} \hat{r}^{n+1}=0$），故只需对括号内求导：
$$
\frac{\partial \tilde{\mathcal{R}}_n}{\partial \hat{r}} = \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2 - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \xi} \left( n \hat{r}^{n-1} \xi^{-n+1} + n \hat{r}^{-n-1} \xi^{n+1} \right) d\xi = n \tilde{\mathcal{S}}_n(\hat{r},s),
$$
其中
$$
\tilde{\mathcal{S}}_n(\hat{r},s) = \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2 - j n \frac{\hat{\Gamma}}{\hat{Q}} \ln \xi} \left( \hat{r}^{n-1} \xi^{-n+1} + \hat{r}^{-n-1} \xi^{n+1} \right) d\xi.
$$
现在将以上导数组合成 $-\partial\delta\hat{\Psi}/\partial\hat{r}$：
$$
\delta\hat{V}_\theta = - \left[ \sum_{n=1}^{\infty} (n \tilde{D}_n \hat{r}^{n-1} - n \tilde{E}_n \hat{r}^{-n-1}) e^{jn\hat{\theta}} + \tilde{F}_0 \frac{\delta\tilde{J}_0}{\hat{r}} + \sum_{n=1}^{\infty} \tilde{F}_n \, n \delta\tilde{\mathcal{S}}_n e^{jn\hat{\theta}} \right].
$$
整理得：
$$
\boxed{
\delta\hat{V}_\theta(\hat{r},\hat{\theta},s) = -\frac{\tilde{F}_0(s) \tilde{J}_0(\hat{r},s)}{\hat{r}} - \sum_{n=1}^{\infty} n \Bigl[ \tilde{D}_n(s) \hat{r}^{n-1} - \tilde{E}_n(s) \hat{r}^{-n-1} + \tilde{F}_n(s) \tilde{\mathcal{S}}_n(\hat{r},s) \Bigr] e^{j n \hat{\theta}}.
}\tag{15}
$$
## 2.5 压力扰动推导
### 2.5.1 轴对称扰动（$n=0$）
从线性化的径向动量方程出发，推导压力扰动 $\delta\tilde{P}_0$ 的表达式。

**步骤 1：线性化径向动量方程**
柱坐标系下，不可压无粘流动的径向动量方程为：
$$
\frac{\partial v_r}{\partial t} + v_r \frac{\partial v_r}{\partial r} + \frac{v_\theta}{r}\frac{\partial v_r}{\partial \theta} - \frac{v_\theta^2}{r} = -\frac{1}{\rho}\frac{\partial p}{\partial r}
$$
将流动分解为背景（上划线）与小扰动（$\delta$）：
$$
v_r = \overline{v}_r(r) + \delta v_r(r,\theta,t),\quad
v_\theta = \overline{v}_\theta(r) + \delta v_\theta(r,\theta,t),\quad
p = \overline{p}(r) + \delta p(r,\theta,t)
$$
代入方程，减去背景部分（满足定常欧拉方程），忽略二阶小量，得到线性化方程：
$$
\frac{\partial \delta v_r}{\partial t}
+ \overline{v}_r \frac{\partial \delta v_r}{\partial r}
+ \frac{\overline{v}_\theta}{r}\frac{\partial \delta v_r}{\partial \theta}
+ \delta v_r \frac{d\overline{v}_r}{dr}
- \frac{2\overline{v}_\theta}{r}\delta v_\theta
= -\frac{1}{\rho}\frac{\partial \delta p}{\partial r}
\tag{16}
$$
**步骤 2：无量纲化与频域变换**
引入同样的无量纲变量
$$
\hat r = \frac{r}{R},\quad \tau = \frac{tU}{R},\quad
\delta\hat{v}_r = \frac{\delta v_r}{U},\quad
\delta\hat{v}_\theta = \frac{\delta v_\theta}{U},\quad
\delta\hat{p} = \frac{\delta p}{\rho U^2},\quad
\hat Q = \frac{Q}{UR},\quad
\hat\Gamma = \frac{\Gamma}{UR}
$$
背景自由涡的无量纲形式：
$$
\hat{\overline{v}}_r = \frac{\hat Q}{\hat r},\qquad
\hat{\overline{v}}_\theta = \frac{\hat\Gamma}{\hat r},\qquad
\frac{d\hat{\overline{v}}_r}{d\hat r} = -\frac{\hat Q}{\hat r^2}
$$
对方程 (16) 进行无量纲化（密度直接被吃掉成1），并对时间作拉普拉斯变换（$\partial/\partial\tau \to s$），同时考虑轴对称扰动（$\partial/\partial\theta = 0$），得：
$$
s\,\delta\hat{V}_r
+ \hat{\overline{v}}_r \frac{d\delta\hat{V}_r}{d\hat r}
+ \delta\hat{V}_r \frac{d\hat{\overline{v}}_r}{d\hat r}
- \frac{2\hat{\overline{v}}_\theta}{\hat r}\,\delta\hat{V}_\theta
= -\frac{d\delta\hat{p}}{d\hat r}
\tag{17}
$$
**步骤 3：代入 $n=0$ 的速度扰动**
由 2.4 节，$n=0$ 时：
$$
\delta\hat{V}_r = \frac{\tilde{D}_0(s)}{\hat r},\qquad
\delta\hat{V}_\theta = -\frac{\tilde{F}_0(s)\,\delta\tilde{J}_0(\hat r,s)}{\hat r}
$$
其中
$$
\delta\tilde{J}_0(\hat r,s) = \int_{\hat r_0}^{\hat r} e^{-\frac{s}{2\hat Q}\xi^2}\,\xi\,d\xi
$$
将上述表达式代入式 (17)，逐项计算：
- $s\,\delta\hat{V}_r = \dfrac{s\tilde{D}_0}{\hat r}$
- $\hat{\overline{v}}_r \dfrac{d\delta\hat{V}_r}{d\hat r} = \dfrac{\hat Q}{\hat r}\cdot\left(-\dfrac{\tilde{D}_0}{\hat r^2}\right) = -\dfrac{\hat Q\tilde{D}_0}{\hat r^3}$
- $\delta\hat{V}_r \dfrac{d\hat{\overline{v}}_r}{d\hat r} = \dfrac{\tilde{D}_0}{\hat r}\cdot\left(-\dfrac{\hat Q}{\hat r^2}\right) = -\dfrac{\hat Q\tilde{D}_0}{\hat r^3}$
- $-\dfrac{2\hat{\overline{v}}_\theta}{\hat r}\delta\hat{V}_\theta = -\dfrac{2(\hat\Gamma/\hat r)}{\hat r}\cdot\left(-\dfrac{\tilde{F}_0 \delta\tilde{J}_0}{\hat r}\right) = +\dfrac{2\hat\Gamma\tilde{F}_0 \delta\tilde{J}_0}{\hat r^3}$
合并前三项，得到：
$$
\frac{s\tilde{D}_0}{\hat r} - \frac{2\hat Q\tilde{D}_0}{\hat r^3} + \frac{2\hat\Gamma\tilde{F}_0 \delta\tilde{J}_0}{\hat r^3}
= -\frac{d\delta\hat{p}}{d\hat r}
$$
移项：
$$
\frac{d\delta\hat{p}}{d\hat r}
= -\frac{s\tilde{D}_0}{\hat r} + \frac{2\hat Q\tilde{D}_0}{\hat r^3} - \frac{2\hat\Gamma\tilde{F}_0 \delta\tilde{J}_0}{\hat r^3}
\tag{18}
$$
**步骤 4：将 $\delta\tilde{J}_0$ 写为显式指数形式**
由涡量解已知：
$$
\delta\tilde{J}_0(\hat r,s) = \frac{\hat Q}{s}\left(e^{-\frac{s}{2\hat Q}\hat r_0^2} - e^{-\frac{s}{2\hat Q}\hat r^2}\right)
$$
代入式 (17)：
$$
\frac{d\delta\hat{p}}{d\hat r}
= -\frac{s\tilde{D}_0}{\hat r}
+ \frac{2\hat Q\tilde{D}_0}{\hat r^3}
- \frac{2\hat\Gamma\tilde{F}_0}{\hat r^3}\cdot\frac{\hat Q}{s}\left(e^{-\frac{s}{2\hat Q}\hat r_0^2} - e^{-\frac{s}{2\hat Q}\hat r^2}\right)
$$
展开：
$$
\frac{d\delta\hat{p}}{d\hat r}
= -\frac{s\tilde{D}_0}{\hat r}
+ \frac{2\hat Q\tilde{D}_0}{\hat r^3}
- \frac{2\hat\Gamma\hat Q}{s}\tilde{F}_0\frac{e^{-\frac{s}{2\hat Q}\hat r_0^2}}{\hat r^3}
+ \frac{2\hat\Gamma\hat Q}{s}\tilde{F}_0\frac{e^{-\frac{s}{2\hat Q}\hat r^2}}{\hat r^3}
\tag{18'}
$$

其中第三项为与 $\hat r_0$ 有关的常数乘以了 $1/\hat r^3$，与第二项形式相同。将其合并进第二项，重新定义系数：
$$
\tilde{D}_0' = \tilde{D}_0 - \frac{\hat\Gamma}{s} e^{-\frac{s}{2\hat Q}\hat r_0^2} \tilde{F}_0
$$
则压力梯度简化为：
$$
\frac{d\delta\hat{p}}{d\hat r}
= \left( \frac{2\hat Q}{\hat r^3} - \frac{s}{\hat r} \right) \tilde{D}_0'(s)
+ \frac{2\hat\Gamma}{\hat r^3 s \hat Q} e^{-\frac{s}{2\hat Q}\hat r^2} \tilde{F}_0(s)
$$
**步骤5：得到积分形式**
对 $\hat r$ 从 $\hat r_0$ 到 $\hat r$ 积分，并将积分下限贡献和 $\delta\hat{p}(\hat r_0)$ 吸收进常数 $C_0(s)$，得：
$$
\delta\hat{p}(\hat r,s) = C_0(s) + \int_{\hat r_0}^{\hat r} \left[ \left( \frac{2\hat Q}{\xi^3} - \frac{s}{\xi} \right) \tilde{D}_0'(s) + \frac{2\hat\Gamma}{\xi^3 s \hat Q} e^{-\frac{s}{2\hat Q}\xi^2} \tilde{F}_0(s) \right] d\xi
$$
在传输矩阵中，常数 $C_0(s)$ 可吸收进 $\tilde{D}_0'$ 的重定义，因此只保留与 $\hat r$ 相关的积分部分。重新标记 $\tilde{D}_0' \to \tilde{D}_0$，并恢复压力符号为 $\delta\tilde{P}_0$，得：
$$
\boxed{
\delta\tilde{P}_0(\hat r,s) = \int_{\hat r_0}^{\hat r} \left[ \left( \frac{2\hat Q}{\xi^3} - \frac{s}{\xi} \right) \tilde{D}_0(s) + \frac{2\hat\Gamma}{\xi^3 s \hat Q} e^{-\frac{s}{2\hat Q}\xi^2} \tilde{F}_0(s) \right] d\xi
}
\tag{19}
$$

### 2.5.2 压力扰动推导（$n>0$）—— 基于周向动量方程与流函数
对于 $n>0$ 的周向谐波扰动，压力扰动 $\delta\tilde{P}_n$ 可以通过线性化的**周向动量方程**直接导出，且结果可用频域流函数 $\delta\hat{\Psi}$ 的谐波分量 $\delta\tilde{\Psi}_n$ 及其导数表示，无需显式求速度。
#### 线性化周向动量方程
柱坐标系下，不可压无粘流动的周向动量方程为：
$$
\frac{\partial v_\theta}{\partial t} + v_r \frac{\partial v_\theta}{\partial r} + \frac{v_\theta}{r} \frac{\partial v_\theta}{\partial \theta} + \frac{v_r v_\theta}{r} = -\frac{1}{\rho r} \frac{\partial p}{\partial \theta}
$$
将流动分解为背景（上划线）与小扰动（$\delta$），代入并减去背景方程（背景定常、轴对称，$\partial\overline{p}/\partial\theta=0$），忽略二阶小量，得到线性化方程：
$$
\frac{\partial \delta v_\theta}{\partial t}
+ \overline{v}_r \frac{\partial \delta v_\theta}{\partial r}
+ \frac{\overline{v}_\theta}{r} \frac{\partial \delta v_\theta}{\partial \theta}
+ \delta v_r \frac{d\overline{v}_\theta}{dr}
+ \frac{\overline{v}_\theta}{r} \delta v_r
+ \frac{\overline{v}_r}{r} \delta v_\theta
= -\frac{1}{\rho r} \frac{\partial \delta p}{\partial \theta}
\tag{20}
$$
#### 无量纲化与拉普拉斯变换

回顾一下，记频域中的无量纲扰动量为：
$$
\delta\hat{V}_r(\hat r,\theta,s) = \mathcal{L}\{\delta\hat{v}_r\},\quad
\delta\hat{V}_\theta(\hat r,\theta,s) = \mathcal{L}\{\delta\hat{v}_\theta\},\quad
\delta\hat{P}(\hat r,\theta,s) = \mathcal{L}\{\delta\hat{p}\}
$$
同时，频域流函数扰动为 $\delta\hat{\Psi}(\hat r,\theta,s) = \mathcal{L}\{\delta\hat{\psi}\}$，且满足：
$$
\delta\hat{V}_r = \frac{1}{\hat r}\frac{\partial \delta\hat{\Psi}}{\partial \theta},\qquad
\delta\hat{V}_\theta = -\frac{\partial \delta\hat{\Psi}}{\partial \hat r}
$$
无量纲化后的周向动量方程（频域）为：
$$
s\,\delta\hat{V}_\theta
+ \hat{\overline{v}}_r \frac{\partial \delta\hat{V}_\theta}{\partial \hat r}
+ \frac{\hat{\overline{v}}_\theta}{\hat r} \frac{\partial \delta\hat{V}_\theta}{\partial \theta}
+ \delta\hat{V}_r \frac{d\hat{\overline{v}}_\theta}{d\hat r}
+ \frac{\hat{\overline{v}}_\theta}{\hat r} \delta\hat{V}_r
+ \frac{\hat{\overline{v}}_r}{\hat r} \delta\hat{V}_\theta
= -\frac{1}{\hat r} \frac{\partial \delta\hat{P}}{\partial \theta}
\tag{21}
$$
#### 代入谐波形式

对于第 $n>0$ 次空间谐波，设频域流函数扰动为：
$$
\delta\hat{\Psi}(\hat r,\theta,s) = \delta\tilde{\Psi}_n(\hat r,s) e^{jn\theta}
$$
其中波浪号 $\tilde{}$ 表示空间傅里叶系数，且 $\delta\tilde{\Psi}_n$ 本身仍是频域量（依赖 $s$）。则速度扰动和压力扰动的谐波分量为：
$$
\delta\hat{V}_r = \frac{1}{\hat r}\frac{\partial \delta\hat{\Psi}}{\partial \theta}
= \frac{jn}{\hat r} \delta\tilde{\Psi}_n e^{jn\theta}, \qquad
\delta\hat{V}_\theta = -\frac{\partial \delta\hat{\Psi}}{\partial \hat r}
= -\frac{d\delta\tilde{\Psi}_n}{d\hat r} e^{jn\theta}, \qquad
\delta\hat{P} = \delta\tilde{P}_n(\hat r,s) e^{jn\theta}
$$
代入方程 (21)，并利用背景自由涡 $\hat{\overline{v}}_r = \hat Q/\hat r$，$\hat{\overline{v}}_\theta = \hat\Gamma/\hat r$，$\frac{d\hat{\overline{v}}_\theta}{d\hat r} = -\hat\Gamma/\hat r^2$。各项计算如下（已约去公因子 $e^{jn\theta}$）：
- $s\,\delta\hat{V}_\theta \to s\left(-\frac{d\delta\tilde{\Psi}_n}{d\hat r}\right)$
- $\hat{\overline{v}}_r \frac{\partial \delta\hat{V}_\theta}{\partial \hat r} \to \frac{\hat Q}{\hat r} \left(-\frac{d^2\delta\tilde{\Psi}_n}{d\hat r^2}\right)$
- $\frac{\hat{\overline{v}}_\theta}{\hat r} \frac{\partial \delta\hat{V}_\theta}{\partial \theta} \to \frac{\hat\Gamma/\hat r}{\hat r} \cdot (jn) \left(-\frac{d\delta\tilde{\Psi}_n}{d\hat r}\right) = -\frac{jn\hat\Gamma}{\hat r^2} \frac{d\delta\tilde{\Psi}_n}{d\hat r}$
- $\delta\hat{V}_r \frac{d\hat{\overline{v}}_\theta}{d\hat r} \to \frac{jn}{\hat r}\delta\tilde{\Psi}_n \cdot \left(-\frac{\hat\Gamma}{\hat r^2}\right) = -\frac{jn\hat\Gamma}{\hat r^3} \delta\tilde{\Psi}_n$
- $\frac{\hat{\overline{v}}_\theta}{\hat r} \delta\hat{V}_r \to \frac{\hat\Gamma/\hat r}{\hat r} \cdot \frac{jn}{\hat r}\delta\tilde{\Psi}_n = \frac{jn\hat\Gamma}{\hat r^3} \delta\tilde{\Psi}_n$
- $\frac{\hat{\overline{v}}_r}{\hat r} \delta\hat{V}_\theta \to \frac{\hat Q/\hat r}{\hat r} \left(-\frac{d\delta\tilde{\Psi}_n}{d\hat r}\right) = -\frac{\hat Q}{\hat r^2} \frac{d\delta\tilde{\Psi}_n}{d\hat r}$

注意第四项与第五项互为相反数，**相互抵消**（物理上源于自由涡的角动量守恒）。将所有左边项合并，得到：
$$
- s \frac{d\delta\tilde{\Psi}_n}{d\hat r}
- \frac{\hat Q}{\hat r} \frac{d^2\delta\tilde{\Psi}_n}{d\hat r^2}
- \frac{jn\hat\Gamma}{\hat r^2} \frac{d\delta\tilde{\Psi}_n}{d\hat r}
- \frac{\hat Q}{\hat r^2} \frac{d\delta\tilde{\Psi}_n}{d\hat r}
= -\frac{jn}{\hat r} \delta\tilde{P}_n
$$
两边乘以 $-1$：
$$
s \frac{d\delta\tilde{\Psi}_n}{d\hat r}
+ \frac{\hat Q}{\hat r} \frac{d^2\delta\tilde{\Psi}_n}{d\hat r^2}
+ \left( \frac{jn\hat\Gamma}{\hat r^2} + \frac{\hat Q}{\hat r^2} \right) \frac{d\delta\tilde{\Psi}_n}{d\hat r}
= \frac{jn}{\hat r} \delta\tilde{P}_n
$$
#### 解出 $\delta\tilde{P}_n$

将左边提取公因子，并乘以 $\hat r/(jn)$，得：
$$
\delta\tilde{P}_n(\hat r,s) = \frac{\hat r}{jn} \left[ \frac{\hat Q}{\hat r} \frac{d^2\delta\tilde{\Psi}_n}{d\hat r^2}
+ \left( s + \frac{jn\hat\Gamma}{\hat r^2} + \frac{\hat Q}{\hat r^2} \right) \frac{d\delta\tilde{\Psi}_n}{d\hat r} \right]
$$
即：
$$
\boxed{
\delta\tilde{P}_n(\hat r,s) = \frac{1}{jn} \left( \hat Q \frac{d^2\delta\tilde{\Psi}_n}{d\hat r^2}
+ \left( s\hat r + \frac{jn\hat\Gamma}{\hat r} + \frac{\hat Q}{\hat r} \right) \frac{d\delta\tilde{\Psi}_n}{d\hat r} \right)
}
\tag{21}
$$
其中 $\delta\tilde{\Psi}_n(\hat r,s)$ 是频域流函数扰动的第 $n$ 次空间谐波分量。
## 2.6 传递矩阵 $\mathbb{T}_{\text{rad},n}$ 的推导
径向空间背景为自由涡：$\hat{\overline{v}}_r = \hat Q / \hat r$，$\hat{\overline{v}}_\theta = \hat\Gamma / \hat r$。无量纲流函数 $\delta\hat{\Psi}$ 满足泊松方程，其通解由齐次势模态和旋转模态叠加而成（2.3节）。速度与流函数的关系为：
$$
\delta\hat{V}_r = \frac{1}{\hat r}\frac{\partial \delta\hat{\Psi}}{\partial \theta},\qquad
\delta\hat{V}_\theta = -\frac{\partial \delta\hat{\Psi}}{\partial \hat r}
$$
压力扰动由动量方程导出：$n=0$ 时用径向动量方程积分；$n>0$ 时用周向动量方程。
将第 $n$ 次谐波的流函数通解写为：
$$
\delta\tilde{\Psi}_n(\hat r,s) = \tilde{D}_n \hat r^{\,n} + \tilde{E}_n \hat r^{-n} + \tilde{F}_n \tilde{\mathcal{R}}_n(\hat r,s)
$$
其中 $\tilde{\mathcal{R}}_n$ 为特解积分。代入速度及压力表达式，整理系数即得传递矩阵（可以自己回代验证）。
$$
\begin{bmatrix}
\delta \hat V_r \\[2pt]
\delta \hat V_\theta \\[2pt]
\delta \hat P
\end{bmatrix}
(\hat r,\hat\theta,s) 
= 
\mathbb{T}_{\text{rad},0} \cdot
\begin{bmatrix}
\tilde D_0(s) \\[2pt]
\tilde F_0(s)
\end{bmatrix}
\;+\;
\sum_{n=1}^{\infty} \mathbb{T}_{\text{rad},n} \cdot e^{jn\hat\theta}\cdot
\begin{bmatrix}
\tilde D_n(s) \\[2pt]
\tilde E_n(s) \\[2pt]
\tilde F_n(s)
\end{bmatrix}
$$
### 一、$n = 0$（轴对称）
状态向量 $[\delta\tilde{V}_{r0},\ \delta\tilde{V}_{\theta0},\ \delta\tilde{P}_0]^{\mathsf{T}}$，系数向量 $[\tilde{D}_0,\ \tilde{F}_0]^{\mathsf{T}}$（$\tilde{E}_0$ 不产生物理贡献）：
$$
\mathbb{T}_{\text{rad},0} =
\begin{bmatrix}
\dfrac{1}{\hat r} & 0 \\\\
0 & \dfrac{\hat Q}{s\hat r}e^{-\frac{s}{2\hat Q}\hat r^2} \\\\
-\dfrac{\hat Q}{\hat r^2} - s\ln\hat r & \dfrac{2\hat\Gamma}{s\hat Q}\tilde{\mathcal R}_0(\hat r,s)
\end{bmatrix}
$$
其中 $\tilde{\mathcal R}_0(\hat r,s) = \displaystyle\int_{\hat r_0}^{\hat r} e^{-\frac{s}{2\hat Q}\xi^2}\,\xi^{-3}d\xi$，常数吸收进势流系数。
### 二、$n > 0$
状态向量 $[\delta\tilde{V}_{rn},\ \delta\tilde{V}_{\theta n},\ \delta\tilde{P}_n]^{\mathsf{T}}$，系数向量 $[\tilde{D}_n,\ \tilde{E}_n,\ \tilde{F}_n]^{\mathsf{T}}$：
$$
\mathbb{T}_{\text{rad},n} =
\begin{bmatrix}
jn\hat r^{\,n-1} & jn\hat r^{-n-1} & \dfrac{jn}{\hat r}\tilde{\mathcal{R}}_n \\
-n\hat r^{\,n-1} & n\hat r^{-n-1} & -n\tilde{\mathcal{S}}_n \\
\mathbb{T}_{31} & \mathbb{T}_{32} & \mathbb{T}_{33}
\end{bmatrix}

$$
其中
$$
\begin{aligned}
\mathbb{T}_{31} &= -js\hat r^{\,n} - jn(\hat Q + j\hat\Gamma)\hat r^{\,n-2},\\
\mathbb{T}_{32} &= js\hat r^{-n} - jn(\hat Q - j\hat\Gamma)\hat r^{-n-2},\\
\mathbb{T}_{33} &= \frac{1}{jn}\left( \hat Q \frac{d^2\tilde{\mathcal{R}}_n}{d\hat r^2}
+ \left( s\hat r + \frac{jn\hat\Gamma}{\hat r} + \frac{\hat Q}{\hat r} \right) \frac{d\tilde{\mathcal{R}}_n}{d\hat r} \right)
\end{aligned}
$$
辅助函数 $$\tilde{\mathcal{S}}_n = \int_{\hat r_0}^{\hat r} e^{-\frac{s}{2\hat Q}\xi^2 - jn\frac{\hat\Gamma}{\hat Q}\ln\xi}
\left( \hat r^{\,n-1}\xi^{-n+1} + \hat r^{-n-1}\xi^{\,n+1} \right) d\xi$$

[^1]: 当 $n\to 0$ 时，积分核中的括号内表达式：
	$$
	\hat{r}^n \xi^{-n+1} - \hat{r}^{-n} \xi^{n+1} \;\to\; \xi - \xi = 0
	$$
	出现 $0/0$ 型不定式（若保留分母 $2n$）。为此考虑极限：
	
	令 $A(n) = \hat{r}^n \xi^{-n+1} = \xi\,(\hat{r}/\xi)^n$，$B(n) = \hat{r}^{-n} \xi^{n+1} = \xi\,(\xi/\hat{r})^n$。
	
	则
	$$
	\frac{A(n)-B(n)}{2n} = \xi \cdot \frac{(\hat{r}/\xi)^n - (\xi/\hat{r})^n}{2n}
	$$
	利用 $(\hat{r}/\xi)^n = e^{n\ln(\hat{r}/\xi)} \approx 1 + n\ln(\hat{r}/\xi) + O(n^2)$，$(\xi/\hat{r})^n = e^{n\ln(\xi/\hat{r})} \approx 1 - n\ln(\hat{r}/\xi) + O(n^2)$。相减得 $2n\ln(\hat{r}/\xi) + O(n^3)$，除以 $2n$ 后极限为 $\ln(\hat{r}/\xi) = \ln\hat{r} - \ln\xi$。
	
	因此：
	$$
	\lim_{n\to0} \frac{\hat{r}^n \xi^{-n+1} - \hat{r}^{-n} \xi^{n+1}}{2n} = \xi\,(\ln\hat{r} - \ln\xi)
	$$
	代入 $n>0$ 特解表达式（保留分母 $2n$），并令 $\tilde{F}_0(s) = \lim_{n\to0} \tilde{F}_n(s)$，得到：
	$$
	\delta\hat{\Psi}_{p,0} = \tilde{F}_0(s) \int_{\hat{r}_0}^{\hat{r}} e^{-\frac{s}{2\hat{Q}}\xi^2} \cdot \xi\,(\ln\hat{r} - \ln\xi) \, d\xi
	$$