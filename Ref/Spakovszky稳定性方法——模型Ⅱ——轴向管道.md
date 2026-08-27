|  属性  |   [[法术表]]    |          |           |           |            |           |
| :--: | :----------: | :------: | :-------: | :-------: | :--------: | :-------: |
|  关联  | [[涡量、旋度恒等式]] | [[欧拉方程]] | [[ODE求解]] | [[傅里叶变换]] | [[拉普拉斯变换]] | [[PDE求解]] |
| 链接解释 |    两个都用上了    |   基础方程   |    解！     |  用来消时间导数  |  用来消周向导数   |    解！     |
***
# $Model\space Ⅱ$
	Model II 的轴向和径向传递矩阵构建了一个灵活，简约的稳定性预测艺术品。这里是它的开始，最上游的组件。

|    变量名    |                        符号                        |           单位           |        量纲         |                        说明                         |
| :-------: | :----------------------------------------------: | :--------------------: | :---------------: | :-----------------------------------------------: |
|   速度矢量    |           $\mathbf{v}$ 或 $\mathbf{V}$            |  $\text{m}/\text{s}$   |    $L T^{-1}$     |          二维速度场，分量为 $(V_x, V_\theta, 0)$           |
|   涡量矢量    | $\boldsymbol{\omega} = \nabla \times \mathbf{V}$ |    $\text{s}^{-1}$     |     $T^{-1}$      |              速度的旋度，二维流动中仅有 $z$ 分量非零               |
| 涡量 $z$ 分量 |                    $\omega_z$                    |    $\text{s}^{-1}$     |     $T^{-1}$      |                  二维流动中唯一的非零涡量分量                   |
|    静压     |                       $p$                        |      $\text{Pa}$       | $M L^{-1} T^{-2}$ |                       流体静压                        |
|    总压     |                      $p_t$                       |      $\text{Pa}$       | $M L^{-1} T^{-2}$ | 不可压流中 $p_t = p + \frac{1}{2}\rho\|\mathbf{V}\|^2$ |
|    密度     |                      $\rho$                      | $\text{kg}/\text{m}^3$ |    $M L^{-3}$     |                    假设为常数（不可压缩）                    |
|    势函数    |                      $\Phi$                      | $\text{m}^2/\text{s}$  |   $L^2 T^{-1}$    |         无旋流的速度势，$\mathbf{V} = \nabla\Phi$         |
|  无因次势函数   |                   $\hat{\Phi}$                   |           –            |         1         |            $\hat{\Phi} = \Phi / (UR)$             |
|   特征速度    |                       $U$                        |  $\text{m}/\text{s}$   |    $L T^{-1}$     |                     转子中径切线速度                      |
|   特征长度    |                       $R$                        |       $\text{m}$       |        $L$        |                       平均半径                        |
|   轴向坐标    |                       $x$                        |     $\text{m}$ 或 –     |      $L$ 或 1      |                  有量纲或已由 $R$ 无因次化                  |
|   流量系数    |                 $\phi = V_x / U$                 |           –            |         1         |                      无因次轴向速度                      |
|   无因次时间   |                  $\tau = tU/R$                   |           –            |         1         |                   以转子转动为特征的时间尺度                   |

***
# 主体
# 1. 建模部分
## 1.1 从欧拉方程到涡量方程
无粘不可压流动的动量方程（欧拉方程）：
$$ \frac{\partial \mathbf{v}}{\partial t} + (\mathbf{v} \cdot \nabla) \mathbf{v} = -\frac{1}{\rho} \nabla p  $$
利用涡量$\nabla$恒等式：
$$ (\mathbf{v} \cdot \nabla) \mathbf{v} = \nabla \left( \frac{1}{2} |\mathbf{v}|^2 \right) - \mathbf{v} \times (\nabla \times \mathbf{v}) $$
代入欧拉方程：
$$ \frac{\partial \mathbf{v}}{\partial t} + \nabla \left( \frac{1}{2} |\mathbf{v}|^2 \right) - \mathbf{v} \times \boldsymbol{\omega} = -\frac{1}{\rho} \nabla p $$
引入总压 $p_t = p + \frac{1}{2}\rho |\mathbf{v}|^2$，两边乘 $\rho$ 并整理：
$$ \rho \frac{\partial \mathbf{v}}{\partial t} - \rho(\mathbf{v} \times \boldsymbol{\omega}) = -\nabla p_t $$
除以 $\rho$ 的形式：
$$ \frac{\partial \mathbf{v}}{\partial t} - \mathbf{v} \times \boldsymbol{\omega} = -\frac{1}{\rho}\nabla p_t $$
对上述两边取旋度 ($\nabla \times$)：
$$ \nabla \times \left( \frac{\partial \mathbf{v}}{\partial t} \right) - \nabla \times (\mathbf{v} \times \boldsymbol{\omega}) = \nabla \times \left( -\frac{1}{\rho}\nabla p_t \right) $$
逐项分析：
*   **第一项**：$\nabla \times (\partial \mathbf{v} / \partial t) = \partial (\nabla \times \mathbf{v}) / \partial t = \partial \boldsymbol{\omega} / \partial t$
*   **总压项**：标量函数梯度的旋度恒为零，即 $\nabla \times (\nabla f) \equiv 0$。因此 $\nabla \times (-\frac{1}{\rho}\nabla p_t) = 0$。
*   **第二项**：利用旋度$\nabla$恒等式替换：
$$ \nabla \times (\mathbf{v} \times \boldsymbol{\omega}) = (\boldsymbol{\omega} \cdot \nabla) \mathbf{v} - (\mathbf{v} \cdot \nabla) \boldsymbol{\omega} + \mathbf{v} (\nabla \cdot \boldsymbol{\omega}) - \boldsymbol{\omega} (\nabla \cdot \mathbf{v}) $$
## 1.2 假设简化
在轴-周向的二维流动中：
*   **涡线拉伸项 $(\boldsymbol{\omega} \cdot \nabla) \mathbf{v} = 0$**：因为涡量只有*垂直于二维平面的*分量，即 $\boldsymbol{\omega} = [0, 0, \omega_z]$，而速度只有$x$和$\theta$两个分量，所以这一整项为0了。
*   **旋度的散度恒为零**：$\nabla \cdot \boldsymbol{\omega} = 0$。
*   **不可压连续性**：$\nabla \cdot \mathbf{v} = 0$。

因此，旋度$\nabla$恒等式右边四项中的第一、三、四项全部为零，只剩下 $-(\mathbf{v} \cdot \nabla) \boldsymbol{\omega}$ 一项。

代入上述结果，得到：
$$ \frac{\partial \boldsymbol{\omega}}{\partial t} + (\mathbf{v} \cdot \nabla) \boldsymbol{\omega} = 0 $$
左边正是 $\boldsymbol{\omega}$ 的物质导数 $D\boldsymbol{\omega} / Dt$。因为只有 $z$ 分量非零，最终得到涡量输运方程：
$$ \boxed{\frac{D\omega_z}{Dt} = 0} $$
## 1.3 线性化
### 第一步：流函数定义
为了自动满足不可压连续方程（$\frac{\partial v_x}{\partial x} + \frac{1}{R}\frac{\partial v_\theta}{\partial \theta} = 0$），定义流函数 $\psi(x,\theta,t)$：
$$
\boxed{v_x = \frac{1}{R}\frac{\partial \psi}{\partial \theta}, \qquad v_\theta = -\frac{\partial \psi}{\partial x}} 
$$
**曲率项说明**：这里 $1/R$ (R为常数）的出现是因为对 $\theta$ 的物理长度微分是 $R d\theta$，所以 $v_x = \frac{1}{R}\frac{\partial \psi}{\partial \theta}$ 保证了 $v_x$ 的量纲正确（m/s）。

连续性方程变成：
$$
\frac{\partial v_x}{\partial x} + \frac{1}{R}\frac{\partial v_\theta}{\partial \theta}
= \frac{\partial}{\partial x}\left( \frac{1}{R}\frac{\partial \psi}{\partial \theta} \right) + \frac{1}{R}\frac{\partial}{\partial \theta}\left( -\frac{\partial \psi}{\partial x} \right) = 0
$$
### 第二步：重新表示涡量输运
在柱坐标 $(x, r, \theta)$ 中，对于二维流动，涡量的唯一非零分量为：
$$
\omega_z = \frac{\partial v_\theta}{\partial x} - \frac{1}{R}\frac{\partial v_x}{\partial \theta}
$$
将刚才的流函数代入：
$$
\omega_z = \frac{\partial}{\partial x}\left( -\frac{\partial \psi}{\partial x} \right) - \frac{1}{R}\frac{\partial}{\partial \theta}\left( \frac{1}{R}\frac{\partial \psi}{\partial \theta} \right)
$$
$$
\omega_z = -\frac{\partial^2 \psi}{\partial x^2} - \frac{1}{R^2}\frac{\partial^2 \psi}{\partial \theta^2}
$$
定义拉普拉斯算子（在 $x,\theta$ 空间，保留曲率）为：
$$
\nabla^2 \equiv \frac{\partial^2}{\partial x^2} + \frac{1}{R^2}\frac{\partial^2}{\partial \theta^2}
$$
则：
$$
\boxed{\omega_z = -\nabla^2 \psi}
$$
### 第三步：涡量输运方程的线性化
由**涡量输运方程（涡量守恒）**：
$$
\frac{D\omega_z}{Dt} = 0 \quad\Longrightarrow\quad
\frac{\partial \omega_z}{\partial t} + v_x \frac{\partial \omega_z}{\partial x} + \frac{v_\theta}{R}\frac{\partial \omega_z}{\partial \theta} = 0 
$$
注意：最后一项中，对弧长的导数是 $\frac{v_\theta}{R}\frac{\partial}{\partial \theta}$，因为 $v_\theta = R\dot{\theta}$，而 $\frac{\partial}{\partial s} = \frac{1}{R}\frac{\partial}{\partial \theta}$。

将各变量分解为 **定常均匀背景流** 与 **小扰动**：
$$
\begin{aligned}
v_x &= \overline{v}_x + \delta v_x, &\quad v_\theta &= \overline{v}_\theta + \delta v_\theta,\\
\omega_z &= \overline{\omega}_z + \delta\zeta, &\quad \psi &= \overline{\psi} + \delta\psi
\end{aligned}
$$
涡量输运方程就变成：
$$
\frac{\partial (\overline{\omega}_z + \delta\zeta)}{\partial t}
+ (\overline{v}_x + \delta v_x)\frac{\partial (\overline{\omega}_z + \delta\zeta)}{\partial x}
+ \frac{\overline{v}_\theta + \delta v_\theta}{R}\frac{\partial (\overline{\omega}_z + \delta\zeta)}{\partial \theta} = 0
$$
**逐项展开**：
- 时间项：$\frac{\partial \overline{\omega}_z}{\partial t}=0$，剩余 $\frac{\partial \delta\zeta}{\partial t}$。
- 第一对流项第一部分：$(\overline{v}_x + \delta v_x)\frac{\partial \overline{\omega}_z}{\partial x}$。由于背景均匀，$\frac{\partial \overline{\omega}_z}{\partial x}=0$，此项为 0。
- 第一对流项第二部分：$(\overline{v}_x + \delta v_x)\frac{\partial \delta\zeta}{\partial x} = \overline{v}_x\frac{\partial \delta\zeta}{\partial x} + \delta v_x\frac{\partial \delta\zeta}{\partial x}$。
- 第二对流项（含曲率）：$\frac{\overline{v}_\theta + \delta v_\theta}{R}\frac{\partial \overline{\omega}_z}{\partial \theta} = 0$（因为 $\frac{\partial \overline{\omega}_z}{\partial \theta}=0$）。
- 第二对流项第二部分：$\frac{\overline{v}_\theta + \delta v_\theta}{R}\frac{\partial \delta\zeta}{\partial \theta} = \frac{\overline{v}_\theta}{R}\frac{\partial \delta\zeta}{\partial \theta} + \frac{\delta v_\theta}{R}\frac{\partial \delta\zeta}{\partial \theta}$。

再忽略二阶小量（$\delta v_x\frac{\partial \delta\zeta}{\partial x}$ 和 $\frac{\delta v_\theta}{R}\frac{\partial \delta\zeta}{\partial \theta}$），得到线性化方程：
$$
\frac{\partial \delta\zeta}{\partial t} + \overline{v}_x\frac{\partial \delta\zeta}{\partial x} + \frac{\overline{v}_\theta}{R}\frac{\partial \delta\zeta}{\partial \theta} = 0 
$$
和$Model\space Ⅰ$一样**无量纲化**：取参考速度 $U$（转子中径切向速度），参考长度 $R$（平均半径）。

定义无量纲量：
$$
\tau = \frac{tU}{R},\quad \hat{x} = \frac{x}{R},\quad \hat{\theta} = \theta,\quad
\hat{v}_x = \frac{v_x}{U},\quad \hat{v}_\theta = \frac{v_\theta}{U},\quad \delta\hat{\zeta} = \frac{\delta\zeta R}{U},\quad \delta\hat{\psi} = \frac{\delta\psi}{UR}
$$
背景速度无量纲化后记为 $\hat{\overline{v}}_x, \hat{\overline{v}}_\theta$。将上式无量纲化后有公因子$\frac{U^2}{R^2}$，约掉：
$$
\frac{\partial \delta\hat{\zeta}}{\partial \tau} + \hat{\overline{v}}_x\frac{\partial \delta\hat{\zeta}}{\partial \hat{x}} + \hat{\overline{v}}_\theta\frac{\partial \delta\hat{\zeta}}{\partial \hat{\theta}} = 0
$$
### 第四步：流函数的线性化
做类似分解：
$$
\omega_z = -\nabla^2\psi \;\Longrightarrow\; \overline{\omega}_z + \delta\zeta = -\nabla^2(\overline{\psi} + \delta\psi)
$$
背景流满足 $\overline{\omega}_z = -\nabla^2\overline{\psi}$。将两式相减：
$$
\delta\zeta = -\nabla^2\delta\psi
$$
其中拉普拉斯算子保持原样（包含曲率项）：
$$
\nabla^2\delta\psi = \frac{\partial^2 \delta\psi}{\partial x^2} + \frac{1}{R^2}\frac{\partial^2 \delta\psi}{\partial \theta^2}
$$
无量纲化后：
$$
\delta\hat\zeta = -\nabla^2 \delta\hat\psi,\quad \nabla^2 = \frac{\partial^2}{\partial \hat x^2}+\frac{\partial^2}{\partial \hat \theta^2}
$$

我必须补充说明的是，因为这里二维直管段大轮毂比无量纲化后可以看成圆环面展开成一个平面，本质上是笛卡尔坐标系，所以才有这种形式。另外，这里的无量纲化取了管道的中径，和其他组件耦合时要注意参考值。

现在我们就得到了两个控制方程，并且经过巧妙的变量替换以后，正好两个未知数，是可解的（平均流速度可以由稳态方程得到）。

到这里，建模部分就完结了，我应该给个总结：引入流函数和涡量的核心思想是将二维无粘不可压流动的原始变量（速度、压力）重新配对：流函数自动满足连续性，涡量则通过取旋度消去压力梯度，将欧拉方程转化为涡量守恒方程（$D\boldsymbol{\omega}/Dt=0$）。其最大好处在于：压力完全被消除，非线性对流项简化为涡量的物质导数，线性化后成为一个纯粹的一阶双曲输运方程，可以解析求解（如特征线法）；同时流函数与涡量之间由泊松方程联系，实现了双曲型（涡量输运）与椭圆型（流函数重构）的分离，极大简化了稳定性分析中的模态求解。

这一方法依赖的假设包括：二维流动（无径向变化）、无粘、不可压缩、大轮毂比（轴向管道半径视为常数），且在线性化时要求扰动幅值远小于背景流。其主要限制是：无法直接处理粘性（涡量扩散）、三维效应（涡线拉伸）以及可压缩性；此外，当流动存在强非线性或间断（如激波）时，涡量守恒的简单形式不再成立。
***
# 2. 求解部分
# 2.1 涡量方程部分
首先是涡量方程（无量纲时域）：
$$
\frac{\partial \delta\hat{\zeta}}{\partial \tau} + \hat{\overline{v}}_x\frac{\partial \delta\hat{\zeta}}{\partial \hat{x}} + \hat{\overline{v}}_\theta\frac{\partial \delta\hat{\zeta}}{\partial \hat{\theta}} = 0 \tag{A}
$$
$$
-\delta\hat\zeta = \frac{\partial^2\delta\hat\psi}{\partial \hat x^2}+\frac{\partial^2\delta\hat\psi}{\partial \hat \theta^2} \qquad\Longrightarrow\qquad \delta\hat\zeta = -\nabla^2 \delta\hat\psi,\quad \nabla^2 = \frac{\partial^2}{\partial \hat x^2} + \frac{\partial^2}{\partial \hat \theta^2} \tag{B}
$$
### 2.1.1 变换（时域 → 频域）
用傅里叶和拉普拉斯变换处理是非常重要的步骤，这里不仅让涡量输运方程求解方便，还让后面流函数解的形式非常清晰！

先对时间 $\tau$ 做单边拉普拉斯变换，记
$$
\mathcal{L}\{ f(\tau) \} = F(s) = \int_0^\infty e^{-s\tau} f(\tau)\,d\tau
$$
并假设初始条件为零（$\tau=0$ 时扰动为零）。变换后 $\frac{\partial}{\partial\tau} \to s$，且变换与空间导数可交换。

定义像函数（频域无量纲扰动）：
$$
\delta\hat{Z}(\hat x,\hat\theta,s) = \mathcal{L}\{\delta\hat\zeta(\hat x,\hat\theta,\tau)\},\qquad
\delta\hat{\Psi}(\hat x,\hat\theta,s) = \mathcal{L}\{\delta\hat\psi(\hat x,\hat\theta,\tau)\}
$$
对方程 (A) 两边做拉普拉斯变换（基流只是数值，不变）：
$$
s \delta\hat{Z} + \hat{\overline{v}}_x \frac{\partial \delta\hat{Z}}{\partial \hat x} + \hat{\overline{v}}_\theta \frac{\partial \delta\hat{Z}}{\partial \hat \theta} = 0 \tag{C}
$$
方程 (B) 不含时间导数，变换后形式不变：
$$
\delta\hat{Z} = -\nabla^2 \delta\hat{\Psi} \tag{D}
$$
将涡量扰动按圆周方向展开为傅里叶级数（物理量在周向必定是周期的）：
$$
\delta\hat{Z}(\hat x, \hat \theta, s) = \sum_{n=0}^{\infty} \delta\tilde{Z}_n(\hat x, s) \cdot e^{j n \hat \theta}
$$
其中：
- $n$ 是周向波数（整数，$n=0$ 表示轴对称扰动）
- $\delta\tilde{Z}_n(\hat x, s)$ 是第 $n$ 阶谐波的复振幅（沿轴向变化）

将级数展开代入方程 (C)：
$$
s \sum_{n=0}^{\infty} \delta\tilde{Z}_n e^{j n \hat\theta}
+ \hat{\overline{v}}_x \sum_{n=0}^{\infty} \frac{d \delta\tilde{Z}_n}{d \hat x} e^{j n \hat\theta}
+ \hat{\overline{v}}_\theta \sum_{n=0}^{\infty} (j n) \delta\tilde{Z}_n e^{j n \hat\theta}
= 0
$$
由于各谐波 $e^{j n \hat\theta}$ 线性独立，可以对每个 $n$ 分别写出：
$$
\boxed{
s \delta\tilde{Z}_n + \hat{\overline{v}}_x \frac{d \delta\tilde{Z}_n}{d \hat x} + j n \hat{\overline{v}}_\theta \delta\tilde{Z}_n = 0
}
$$
## 2.1.2 求解涡量输运方程ODE
现在的方程是一阶线性齐次 ODE（$\hat x$ 为自变量）。其通解为：
$$
\delta\tilde{Z}_n(\hat x, s) = \tilde{C}_n(s) \cdot \exp\left( -\frac{s + j n \hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} \hat x \right)
$$
其中 $\tilde{C}_n(s)$ 是由边界条件确定的复常数。

将上式代回级数，得到涡量扰动的频域-空间表达式：
$$
\delta\hat{Z}(\hat x, \hat \theta, s) = \sum_{n=0}^{\infty} \tilde{C}_n(s) \cdot e^{j n \hat\theta} \cdot \exp\left( -\frac{s + j n \hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} \hat x \right)
$$
物理意义：涡量扰动被平均流向下游输运，且不同周向波数的扰动以不同的复波数传播。
## 2.2 流函数部分
涡量扰动与流函数扰动的关系是泊松方程：
$$
\delta\hat{Z}(\hat{x}, \hat{\theta}, s) = -\nabla^2 \delta\hat{\Psi}(\hat{x}, \hat{\theta}, s), \qquad
\nabla^2 = \frac{\partial^2}{\partial \hat{x}^2} + \frac{\partial^2}{\partial \hat{\theta}^2} \tag{D}
$$

涡量扰动刚才才解出来：
$$
\delta\hat{Z}(\hat{x}, \hat{\theta}, s) = \sum_{n=0}^{\infty} \tilde{C}_n(s) e^{-k_n \hat{x}} e^{j n \hat{\theta}}, \quad
k_n = \frac{s + j n \hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} \tag{E}
$$
### 2.2.1 求解流函数PDE
将 $\delta\hat{Z}$ 代入 (D) 得：
$$
\frac{\partial^2 \delta\hat{\Psi}}{\partial \hat{x}^2} + \frac{\partial^2 \delta\hat{\Psi}}{\partial \hat{\theta}^2} = - \sum_{n=0}^{\infty} \tilde{C}_n(s) e^{-k_n \hat{x}} e^{j n \hat{\theta}} \tag{F}
$$

这是一个线性非齐次偏微分方程。其通解 = 齐次解 + 特解。
#### 齐次方程：
$$
\frac{\partial^2 \delta\hat{\Psi}_h}{\partial \hat{x}^2} + \frac{\partial^2 \delta\hat{\Psi}_h}{\partial \hat{\theta}^2} = 0
$$
这是二维拉普拉斯方程。在直角坐标 $(\hat{x}, \hat{\theta})$ 中，其通解可分离变量，一般解为（见关联文件中，笛卡尔坐标系下的PDE求解部分）：
$$
\delta\hat{\Psi}_h(\hat{x}, \hat{\theta}, s) = \sum_{n=0}^{\infty} \left[ \tilde{A}_n(s) e^{n \hat{x}} + \tilde{B}_n(s) e^{-n \hat{x}} \right] \cos(n \hat{\theta} + \phi_n)
$$
但为了与指数形式的傅里叶级数匹配，通常采用复数形式：
$$
\delta\hat{\Psi}_h = \sum_{n=-\infty}^{\infty} \left( \tilde P_n e^{n \hat{x}} + \tilde Q_n e^{-n \hat{x}} \right) e^{j n \hat{\theta}}
$$
然而，当 $n=0$ 时，$e^{0\cdot \hat{x}} = 1$，此时拉普拉斯方程的解还包含 **$\hat{\theta}$ 和 $\hat{x}\hat{\theta}$ 的线性项（特征值为0时的项）**，因为：
$$
\frac{\partial^2}{\partial \hat{x}^2} (a \hat{\theta}) = 0,\quad
\frac{\partial^2}{\partial \hat{\theta}^2} (a \hat{\theta}) = 0
$$
$$
\frac{\partial^2}{\partial \hat{x}^2} (b \hat{x}\hat{\theta}) = 0,\quad
\frac{\partial^2}{\partial \hat{\theta}^2} (b \hat{x}\hat{\theta}) = 0
$$
所以完整齐次解应包含：
$$
\delta\hat{\Psi}_h = \alpha_0(s) \hat{x} + \alpha_1(s) \hat{\theta} + \alpha_2(s) \hat{x}\hat{\theta} + \alpha_3(s)
+ \sum_{n=1}^{\infty} \left( \tilde{A}_n(s) e^{n \hat{x}} + \tilde{B}_n(s) e^{-n \hat{x}} \right) e^{j n \hat{\theta}} + \text{共轭项}
$$
其中 $\alpha_0, \alpha_1, \alpha_2, \alpha_3$ 是 $s$ 的任意函数。
#### 特解
对于非齐次项中的每一谐波 $n$，设特解形式为：
$$
\delta\hat{\Psi}_{p,n} = \tilde{D}_n(s) e^{-k_n \hat{x}} e^{j n \hat{\theta}}
$$
代入 (F) 左边：
$$
\frac{\partial^2}{\partial \hat{x}^2} (\tilde{D}_n e^{-k_n \hat{x}} e^{j n \hat{\theta}}) = \tilde{D}_n k_n^2 e^{-k_n \hat{x}} e^{j n \hat{\theta}}
$$
$$
\frac{\partial^2}{\partial \hat{\theta}^2} (\tilde{D}_n e^{-k_n \hat{x}} e^{j n \hat{\theta}}) = \tilde{D}_n (-n^2) e^{-k_n \hat{x}} e^{j n \hat{\theta}}
$$
总和：
$$
(k_n^2 - n^2) \tilde{D}_n e^{-k_n \hat{x}} e^{j n \hat{\theta}} = - \tilde{C}_n e^{-k_n \hat{x}} e^{j n \hat{\theta}}
$$
所以：
$$
\tilde{D}_n = - \frac{\tilde{C}_n}{k_n^2 - n^2}, \quad \text{当 } k_n^2 \neq n^2
$$
对于 $n=0$：
$$
\tilde{D}_0 = - \frac{\tilde{C}_0}{k_0^2}, \quad k_0 = s / \hat{\overline{v}}_x
$$
因此特解总和为：
$$
\delta\hat{\Psi}_p = - \sum_{n=0}^{\infty} \frac{\tilde{C}_n}{k_n^2 - n^2} e^{-k_n \hat{x}} e^{j n \hat{\theta}}
$$
（这里对 $n=0$，分母 $k_0^2 - 0^2$，成立）
### 2.2.2 通解
合并齐次解与特解：
$$
\begin{aligned}
\delta\hat{\Psi}(\hat{x}, \hat{\theta}, s) = &\ \alpha_0(s) \hat{x} + \alpha_1(s) \hat{\theta} + \alpha_2(s) \hat{x}\hat{\theta} + \alpha_3(s) \\
&+ \sum_{n=1}^{\infty} \left( \tilde{A}_n(s) e^{n \hat{x}} + \tilde{B}_n(s) e^{-n \hat{x}} \right) e^{j n \hat{\theta}} \\
&- \sum_{n=0}^{\infty} \frac{\tilde{C}_n(s)}{k_n^2 - n^2} e^{-k_n \hat{x}} e^{j n \hat{\theta}}
\end{aligned} \tag{G}
$$

其中当 $n=0$ 时，$e^{j0\hat{\theta}}=1$，且 $\tilde{C}_0$ 项单独处理。
### 2.2.3 舍弃非物理项
Spakovszky博士论文中注释指出：$\alpha_0$ 和 $\alpha_2$ 项会导致从两个动量方程导出的压力扰动相互矛盾，故舍弃。$\alpha_3$ 是常数，可保留但不影响速度。$\alpha_1 \hat{\theta}$ 项是物理的（产生均匀轴向速度），保留。

因此简化：
$$
\delta\hat{\Psi}(\hat{x}, \hat{\theta}, s) = \alpha_1(s) \hat{\theta} + \alpha_3(s) + \sum_{n=1}^{\infty} \left( \tilde{A}_n(s) e^{n \hat{x}} + \tilde{B}_n(s) e^{-n \hat{x}} \right) e^{j n \hat{\theta}} - \sum_{n=0}^{\infty} \frac{\tilde{C}_n(s)}{k_n^2 - n^2} e^{-k_n \hat{x}} e^{j n \hat{\theta}}
$$

同时，注意到对于 $n=0$，$\tilde{A}_0 e^{0\cdot \hat{x}} = \tilde{A}_0$，可与 $\alpha_3$ 合并。论文将 $\alpha_1(s) \hat{\theta}$ 单独写出，并与 $n=0$ 的模态并列。所以我们写成：

$$
\begin{aligned}
\delta\hat{\Psi}(\hat{x}, \hat{\theta}, s) = &\ \alpha_1(s) \hat{\theta} \\
&+ \sum_{n=1}^{\infty} \tilde{A}_n(s) e^{n \hat{x}} e^{j n \hat{\theta}} + \sum_{n=1}^{\infty} \tilde{B}_n(s) e^{-n \hat{x}} e^{j n \hat{\theta}} \\
&- \sum_{n=1}^{\infty} \frac{\tilde{C}_n(s)}{k_n^2 - n^2} e^{-k_n \hat{x}} e^{j n \hat{\theta}} \\
&- \frac{\tilde{C}_0(s)}{k_0^2} e^{-k_0 \hat{x}} + \tilde{A}_0(s) \quad (\text{常数项})
\end{aligned} \tag{H}
$$

### 2.2.4 系数重定义
重新定义常数：
- 对于 $n>0$：
$$
  A_n(s) \equiv j n \tilde{A}_n(s), \quad B_n(s) \equiv j n \tilde{B}_n(s), \quad C_n(s) \equiv - j n \frac{\tilde{C}_n(s)}{k_n^2 - n^2}
  $$
- 对于 $n=0$：
$$
  C_0(s) \equiv - \frac{\tilde{C}_0(s)}{k_0^2} = - \tilde{C}_0(s) \frac{\hat{\overline{v}}_x^2}{s^2}, \quad A_0(s) \equiv \alpha_1(s)
$$并忽略常数 $\tilde{A}_0(s)$。

代入 (H) 得到：
$$
\begin{aligned}
\delta\hat{\Psi}(\hat{x}, \hat{\theta}, s) = &\ A_0(s) \hat{\theta} \\
&+ \sum_{n=1}^{\infty} \frac{A_n(s)}{j n} e^{n \hat{x} + j n \hat{\theta}} 
   + \sum_{n=1}^{\infty} \frac{B_n(s)}{j n} e^{-n \hat{x} + j n \hat{\theta}} \\
&+ C_0(s) e^{-k_0 \hat{x}} 
   + \sum_{n=1}^{\infty} \frac{C_n(s)}{j n} e^{-k_n \hat{x} + j n \hat{\theta}}
\end{aligned} \tag{I}
$$
其中 $k_0 = s / \hat{\overline{v}}_x$，$k_n = \frac{s + j n \hat{\overline{v}}_\theta}{\hat{\overline{v}}_x}$。
### 2.2.5 最终形式
写完整就是：
$$
\boxed{
\begin{aligned}
\delta\hat{\Psi}(\hat{x}, \hat{\theta}, s) = &\ A_0(s)\,\hat{\theta} \;+\; \sum_{n=1}^{\infty} \frac{A_n(s)}{j n} e^{n \hat{x} + j n \hat{\theta}} \\
&+ \sum_{n=1}^{\infty} \frac{B_n(s)}{j n} e^{-n \hat{x} + j n \hat{\theta}} \\
&+ C_0(s) e^{-\frac{s}{\hat{\overline{v}}_x} \hat{x}} \;+\; \sum_{n=1}^{\infty} \frac{C_n(s)}{j n} e^{-\left( \frac{s}{\hat{\overline{v}}_x} + j n \frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} \right) \hat{x} + j n \hat{\theta}}
\end{aligned}
} 
$$
### 2.2.6 物理解释
- $A_0(s)\hat{\theta}$：由 $n=0$ 齐次解中的 $\alpha_1 \hat{\theta}$ 项，对应 **均匀轴向速度扰动** $\delta\hat{V}_x = A_0(s)$。
- 含 $A_n(s)$ 的项：向上游衰减的势模态。
- 含 $B_n(s)$ 的项：向下游衰减的势模态。
- 含 $C_0(s)$ 的项：轴对称的涡模态，向下游对流。
- 含 $C_n(s)$ 的项：高阶涡模态，同时向下游传播并周向旋转。
## 2.3 从流函数解得到速度扰动
流函数解（式 (I)）为：

$$
\begin{aligned}
\delta\hat{\Psi}(\hat{x}, \hat{\theta}, s) = &\ A_0(s)\,\hat{\theta} \;+\; \sum_{n=1}^{\infty} \frac{A_n(s)}{j n} e^{n \hat{x} + j n \hat{\theta}} \\
&+ \sum_{n=1}^{\infty} \frac{B_n(s)}{j n} e^{-n \hat{x} + j n \hat{\theta}} \\
&+ C_0(s) e^{-\frac{s}{\hat{\overline{v}}_x} \hat{x}} \;+\; \sum_{n=1}^{\infty} \frac{C_n(s)}{j n} e^{-\left( \frac{s}{\hat{\overline{v}}_x} + j n \frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} \right) \hat{x} + j n \hat{\theta}}
\end{aligned}
$$

其中 $A_n(s), B_n(s), C_n(s)$ 为待定系数，$\hat{\overline{v}}_x, \hat{\overline{v}}_\theta$ 为无量纲背景速度常数。

根据流函数定义（无量纲化后）：
$$
\delta\hat{V}_x = \frac{\partial \delta\hat{\Psi}}{\partial \hat{\theta}}, \qquad
\delta\hat{V}_\theta = -\frac{\partial \delta\hat{\Psi}}{\partial \hat{x}}
$$
### 2.3.1 轴向速度扰动 $\delta\hat{V}_x$
对 $\hat\theta$ 求偏导：
- 第一项：$\displaystyle \frac{\partial}{\partial\hat\theta}\bigl( A_0(s)\,\hat\theta \bigr) = A_0(s)$
- 第二项（$n\ge1$）：$\displaystyle \frac{\partial}{\partial\hat\theta}\left( \frac{A_n(s)}{j n} e^{n\hat{x} + j n\hat\theta} \right) = A_n(s) e^{n\hat{x} + j n\hat\theta}$
- 第三项（$n\ge1$）：$\displaystyle \frac{\partial}{\partial\hat\theta}\left( \frac{B_n(s)}{j n} e^{-n\hat{x} + j n\hat\theta} \right) = B_n(s) e^{-n\hat{x} + j n\hat\theta}$
- 第四项：$\displaystyle \frac{\partial}{\partial\hat\theta}\left( C_0(s) e^{-\frac{s}{\hat{\overline{v}}_x}\hat{x}} \right) = 0$
- 第五项（$n\ge1$）：$\displaystyle \frac{\partial}{\partial\hat\theta}\left( \frac{C_n(s)}{j n} e^{-k_n\hat{x} + j n\hat\theta} \right) = C_n(s) e^{-k_n\hat{x} + j n\hat\theta}$，其中 $k_n = \frac{s}{\hat{\overline{v}}_x} + j n \frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x}$

合并得：

$$
\boxed{
\delta\hat{V}_x(\hat{x},\hat\theta,s) = A_0(s) \;+\; \sum_{n=1}^{\infty} \Bigl[ A_n(s) e^{n\hat{x}} + B_n(s) e^{-n\hat{x}} + C_n(s) e^{-k_n\hat{x}} \Bigr] e^{j n\hat\theta}
} \tag{J}
$$
### 2.3.2 周向速度扰动 $\delta\hat{V}_\theta$
先计算 $-\partial/\partial\hat{x}$ 作用于各项：
- 第一项：$\displaystyle -\frac{\partial}{\partial\hat{x}}\bigl( A_0(s)\,\hat\theta \bigr) = 0$
- 第二项（$n\ge1$）：$\displaystyle -\frac{\partial}{\partial\hat{x}}\left( \frac{A_n(s)}{j n} e^{n\hat{x} + j n\hat\theta} \right) = j A_n(s) e^{n\hat{x} + j n\hat\theta}$
- 第三项（$n\ge1$）：$\displaystyle -\frac{\partial}{\partial\hat{x}}\left( \frac{B_n(s)}{j n} e^{-n\hat{x} + j n\hat\theta} \right) = -j B_n(s) e^{-n\hat{x} + j n\hat\theta}$
- 第四项：$\displaystyle -\frac{\partial}{\partial\hat{x}}\left( C_0(s) e^{-\frac{s}{\hat{\overline{v}}_x}\hat{x}} \right) = \frac{s}{\hat{\overline{v}}_x} C_0(s) e^{-\frac{s}{\hat{\overline{v}}_x}\hat{x}}$
- 第五项（$n\ge1$）：$\displaystyle -\frac{\partial}{\partial\hat{x}}\left( \frac{C_n(s)}{j n} e^{-k_n\hat{x} + j n\hat\theta} \right) = \frac{k_n}{j n} C_n(s) e^{-k_n\hat{x} + j n\hat\theta} = \left( \frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} - j\frac{s}{n\hat{\overline{v}}_x} \right) C_n(s) e^{-k_n\hat{x} + j n\hat\theta}$

因此：

$$
\boxed{
\begin{aligned}
\delta\hat{V}_\theta(\hat{x},\hat\theta,s) = &\ j\sum_{n=1}^{\infty} A_n(s) e^{n\hat{x} + j n\hat\theta} \;-\; j\sum_{n=1}^{\infty} B_n(s) e^{-n\hat{x} + j n\hat\theta} \\
&+ \frac{s}{\hat{\overline{v}}_x} C_0(s) e^{-\frac{s}{\hat{\overline{v}}_x}\hat{x}} \;+\; \sum_{n=1}^{\infty} \left( \frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} - j\frac{s}{n\hat{\overline{v}}_x} \right) C_n(s) e^{-k_n\hat{x} + j n\hat\theta}
\end{aligned}
} \tag{K}
$$
## 2.4 压力扰动 $\delta\hat{P}$
无量纲线性化欧拉方程的 $x$ 分量（频域）为：
$$
s \delta\hat{V}_x + \hat{\overline{v}}_x \frac{\partial \delta\hat{V}_x}{\partial \hat{x}} + \hat{\overline{v}}_\theta \frac{\partial \delta\hat{V}_x}{\partial \hat{\theta}} = -\frac{\partial \delta\hat{P}}{\partial \hat{x}} \tag{L}
$$
将 $\delta\hat{V}_x$ 的表达式 (J) 代入。由于方程线性且各谐波独立，分别处理 $n=0$ 和 $n \ge 1$ 的 Fourier 系数。采用**从参考点 $\hat{x}_0$ 到 $\hat{x}$ 的定积分**。
### 2.4.1 $n = 0$：轴对称扰动
对于 $n=0$，$\delta\hat{V}_x^{(0)} = A_0(s)$（与 $\hat{x}$ 无关），且 $\partial/\partial\hat{\theta}=0$。代入 (L)：
$$
s A_0(s) = -\frac{d \delta\tilde{P}_0}{d \hat{x}}
$$
从 $\hat{x}_0$ 到 $\hat{x}$ 积分：
$$
\delta\tilde{P}_0(\hat{x}) - \delta\tilde{P}_0(\hat{x}_0) = -s A_0(s) (\hat{x} - \hat{x}_0)
$$

通常取参考点 $\hat{x}_0 = 0$，并记 $\delta\hat{P}_0(0)$ 为已知值，则：

$$
\boxed{ \delta\tilde{P}_0(\hat{x}) = \delta\tilde{P}_0(0) - s A_0(s) \,\hat{x} } \tag{M}
$$
### 2.4.2 $n \ge 1$：非轴对称扰动
对于固定的 $n$，$\delta\hat{V}_x^{(n)}(\hat{x}) = A_n e^{n\hat{x}} + B_n e^{-n\hat{x}} + C_n e^{-k_n\hat{x}}$，其中 $k_n = \dfrac{s}{\hat{\overline{v}}_x} + j n \dfrac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x}$。代入 (L) 的 $n$ 次谐波分量（$\partial/\partial\hat{\theta} \to j n$）：
$$
s \delta\hat{V}_x^{(n)} + \hat{\overline{v}}_x \frac{d \delta\hat{V}_x^{(n)}}{d \hat{x}} + j n \hat{\overline{v}}_\theta \delta\hat{V}_x^{(n)} = -\frac{d \delta\tilde{P}_n}{d \hat{x}} \tag{N}
$$
逐项计算左端：
- 对 $A_n e^{n\hat{x}}$：贡献 $\lambda_A A_n e^{n\hat{x}}$，其中 $\lambda_A = s + n(\hat{\overline{v}}_x + j \hat{\overline{v}}_\theta)$
- 对 $B_n e^{-n\hat{x}}$：贡献 $\lambda_B B_n e^{-n\hat{x}}$，其中 $\lambda_B = s - n\hat{\overline{v}}_x + j n \hat{\overline{v}}_\theta$
- 对 $C_n e^{-k_n\hat{x}}$：贡献为零

因此：
$$
\frac{d \delta\tilde{P}_n}{d \hat{x}} = -\left( \lambda_A A_n e^{n\hat{x}} + \lambda_B B_n e^{-n\hat{x}} \right) \tag{O}
$$
从参考点 $\hat{x}_0$ 到 $\hat{x}$ 积分：
$$
\delta\tilde{P}_n(\hat{x}) - \delta\tilde{P}_n(\hat{x}_0) = -\int_{\hat{x}_0}^{\hat{x}} \left( \lambda_A A_n e^{n\xi} + \lambda_B B_n e^{-n\xi} \right) d\xi
$$
计算积分：
$$
\int_{\hat{x}_0}^{\hat{x}} \lambda_A A_n e^{n\xi} d\xi = \frac{\lambda_A A_n}{n} \left( e^{n\hat{x}} - e^{n\hat{x}_0} \right)
$$
$$
\int_{\hat{x}_0}^{\hat{x}} \lambda_B B_n e^{-n\xi} d\xi = -\frac{\lambda_B B_n}{n} \left( e^{-n\hat{x}} - e^{-n\hat{x}_0} \right)
$$
因此：
$$
\delta\tilde{P}_n(\hat{x}) = \delta\tilde{P}_n(\hat{x}_0) - \frac{\lambda_A A_n}{n} \left( e^{n\hat{x}} - e^{n\hat{x}_0} \right) + \frac{\lambda_B B_n}{n} \left( e^{-n\hat{x}} - e^{-n\hat{x}_0} \right) \tag{P}
$$

通常取参考点 $\hat{x}_0 = 0$，则：
$$
\boxed{ \delta\tilde{P}_n(\hat{x}) = \delta\tilde{P}_n(0) - \frac{\lambda_A A_n}{n} \left( e^{n\hat{x}} - 1 \right) + \frac{\lambda_B B_n}{n} \left( e^{-n\hat{x}} - 1 \right) } \tag{P'}
$$
### 2.4.3 完整压力扰动表达式
综合 $n=0$ 和 $n \ge 1$ 的结果：
$$
\boxed{
\begin{aligned}
\delta\hat{P}(\hat{x},\hat{\theta},s) = &\ \bigl[ \delta\tilde{P}_0(0) - s A_0(s) \,\hat{x} \bigr] \\
&+ \sum_{n=1}^{\infty} \Biggl[ \delta\tilde{P}_n(0) - \frac{\lambda_A A_n}{n} \bigl( e^{n\hat{x}} - 1 \bigr) + \frac{\lambda_B B_n}{n} \bigl( e^{-n\hat{x}} - 1 \bigr) \Biggr] e^{j n \hat{\theta}}
\end{aligned}
}
$$
## 2.5 传递矩阵
### 2.5.1 构成
**回顾之前的定义：**
$$
k_n = \frac{s}{\hat{\overline{v}}_x} + j n \frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x},\qquad
\lambda_A = s + n(\hat{\overline{v}}_x + j\hat{\overline{v}}_\theta),\qquad
\lambda_B = s - n\hat{\overline{v}}_x + j n\hat{\overline{v}}_\theta
$$

**状态向量与系数向量的关系：**
$$
\begin{bmatrix}
\delta\hat{V}_x \\
\delta\hat{V}_\theta \\
\delta\hat{P}
\end{bmatrix}(\hat{x},\hat{\theta},s)
=
\mathbb{T}_0(\hat{x},s)
\begin{bmatrix}
A_0(s) \\ C_0(s)
\end{bmatrix}
+
\sum_{n=1}^{\infty}
\mathbb{T}_n(\hat{x},s)
\begin{bmatrix}
A_n(s) \\ B_n(s) \\ C_n(s)
\end{bmatrix}
e^{j n \hat{\theta}}
+
\begin{bmatrix}
0 \\ 0 \\ \delta\hat{P}_0(0) + \sum_{n=1}^{\infty} \delta\hat{P}_n(0)\,e^{j n \hat{\theta}}
\end{bmatrix}
$$

**矩阵 $\mathbb{T}_0$（轴对称，$n=0$）**
$$
\mathbb{T}_0(\hat{x},s) =
\begin{bmatrix}
1 &0 \\\\
0 & \displaystyle\frac{s}{\hat{\overline{v}}_x}\,e^{-k_0 \hat{x}}\\ \\
-s \hat{x} & 0\\
\end{bmatrix},
\qquad k_0 = \frac{s}{\hat{\overline{v}}_x}
$$
这里$\frac{s}{\hat{\overline v}_x}$可以被$C_0$以重定义的方式吸收，于是有：
$$
\mathbb{T}_0(\hat{x},s) =
\begin{bmatrix}
1 &0 \\\\
0 &e^{-k_0 \hat{x}}\\ \\
-s \hat{x} & 0\\
\end{bmatrix},
\qquad k_0 = \frac{s}{\hat{\overline{v}}_x}
$$

**矩阵 $\mathbb{T}_n$（非轴对称，$n\ge 1$）**
$$
\mathbb{T}_n(\hat{x},s) =
\begin{bmatrix}
e^{n \hat{x}} & e^{-n \hat{x}} & e^{-k_n \hat{x}}\\
j e^{n \hat{x}} & -j e^{-n \hat{x}} & \left(\displaystyle\frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} - j\frac{s}{n\hat{\overline{v}}_x}\right) e^{-k_n \hat{x}} 
\\-\displaystyle\frac{\lambda_A}{n}(e^{n \hat{x}}-1) & \displaystyle\frac{\lambda_B}{n}(e^{-n \hat{x}}-1) & 0
\end{bmatrix}
$$
> 积分常数 $\delta\hat{P}_n(0)$ 由边界条件确定，在构造两点间传递矩阵时会自动消去。

这里，也可以通过设定压力参考点$x_0$压力为0，使矩阵非常简洁，只是要注意，这种形式下自由项是不需要加的：
$$
\mathbb{T}_n(\hat{x},s) =
\begin{bmatrix}
e^{n \hat{x}} & e^{-n \hat{x}} & e^{-k_n \hat{x}}\\
j e^{n \hat{x}} & -j e^{-n \hat{x}} & \left(\displaystyle\frac{\hat{\overline{v}}_\theta}{\hat{\overline{v}}_x} - j\frac{s}{n\hat{\overline{v}}_x}\right) e^{-k_n \hat{x}} 
\\-\displaystyle\frac{\lambda_A}{n}e^{n \hat{x}} & \displaystyle\frac{\lambda_B}{n}e^{-n \hat{x}} & 0
\end{bmatrix}
$$
这是简化为了和论文对上，便于使用。