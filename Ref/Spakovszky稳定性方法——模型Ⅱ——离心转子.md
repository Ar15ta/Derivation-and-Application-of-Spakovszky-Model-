
|  属性  |            [[法术表]]             |              |             |                          |
| :--: | :----------------------------: | ------------ | ----------- | ------------------------ |
|  关联  | [[Spakovszky稳定性方法——模型Ⅱ——轴流转子]] | [[涡量、旋度恒等式]] | [[叶轮机械总参数]] | [[Spakovszky稳定性方法——模型Ⅰ]] |
| 链接解释 |       差不多，都是用匹配条件建立传递矩阵        | 又用           | 出现其中的定义     | 复用非定常总压降，但是速度变成相对速度      |
***
# $Model\space Ⅱ$
    激动的盘。

| 变量名                 | 符号                     | 单位          | 量纲               | 说明                                       |
| ------------------- | ---------------------- | ----------- | ---------------- | ---------------------------------------- |
| 绝对速度                | $\mathbf{v}$           | $\text{m/s}$ | $L T^{-1}$       | 静止坐标系下的速度                              |
| 相对速度                | $\mathbf{w}$           | $\text{m/s}$ | $L T^{-1}$       | 随叶轮旋转的坐标系下的速度                         |
| 叶轮角速度矢量             | $\boldsymbol{\Omega}$  | $\text{rad/s}$ | $T^{-1}$         | 大小 $\Omega$，方向沿旋转轴                     |
| 位置矢量                | $\mathbf{r}$           | $\text{m}$     | $L$              | 从旋转中心指向质点                               |
| 绝对涡量                | $\boldsymbol{\omega}_a = \nabla \times \mathbf{v}$ | $\text{s}^{-1}$ | $T^{-1}$         | 绝对坐标系下的涡量                              |
| 相对涡量                | $\boldsymbol{\omega}_r = \nabla \times \mathbf{w}$ | $\text{s}^{-1}$ | $T^{-1}$         | 旋转坐标系下的涡量                              |
| 相对涡量轴向分量           | $\omega_x$             | $\text{s}^{-1}$ | $T^{-1}$         | 旋转轴方向的分量（这里 $x$ 轴为轴向）                  |
| 叶轮出口半径             | $R_2$                  | $\text{m}$     | $L$              | 叶轮出口平均半径                               |
| 叶轮出口圆周速度           | $U_2 = \Omega R_2$     | $\text{m/s}$ | $L T^{-1}$       | 叶片出口处的切向速度（理想无滑移）                      |
| 叶轮出口通道宽度           | $b$                    | $\text{m}$     | $L$              | 垂直于叶片表面的通道宽度                           |
| 叶片总数                | $N$                    | –              | 1                | 主叶片 + 分叶片（若有）                          |
| 出口叶片金属角（后弯角）      | $\chi_2$               | $\text{rad}$   | 1                | 叶片出口切线与圆周切线之间的夹角，通常 $0 \le \chi_2 \le 90^\circ$ |
| 滑移速度（大小）           | $v_s$                  | $\text{m/s}$ | $L T^{-1}$       | $U_2 - v_{\theta 2}$，方向与旋转相反              |
| 滑移因子                | $\sigma_W$             | –              | 1                | $\sigma_W = v_{\theta 2} / U_2$           |

***
# 主体
![[Spakovszky稳定性方法——模型Ⅱ——离心转子配图.png]]

# 1.旋转总压守恒
假设旋转总压（Rotary Stagnation Pressure）
$$
p^* = p_t - \rho \Omega R v_\theta
$$
沿流线守恒。$p_t = p + \frac{1}{2}\rho v^2$ 为绝对总压，$\Omega$ 为叶轮角速度，$R$ 为半径，$v_\theta$ 为绝对切向速度。

在叶轮进口（STA1）和出口（STA2）处分别写出：
$$
p_{t1} - \rho \Omega R_1 v_{\theta1} = p_{t2} - \rho \Omega R_2 v_{\theta2}
$$
移项：
$$
p_{t2} - p_{t1} = \rho \Omega (R_2 v_{\theta2} - R_1 v_{\theta1})
$$
两边除以 $\rho U_2^2$，其中 $U_2 = \Omega R_2$：
$$
\frac{p_{t2} - p_{t1}}{\rho U_2^2} = \frac{\Omega (R_2 v_{\theta2} - R_1 v_{\theta1})}{(\Omega R_2)^2}
= \frac{v_{\theta2}}{\Omega R_2} - \frac{R_1}{R_2}\cdot\frac{v_{\theta1}}{\Omega R_2}
$$
注意到 $\frac{v_{\theta2}}{\Omega R_2} = \frac{v_{\theta2}}{U_2}$，正是出口处的无量纲切向速度，记作 $\hat v_{\theta2}$；同样 $\frac{ v_{\theta1}}{\Omega R_2} = \frac{v_{\theta1}}{U_2}$ 也是无量纲量，记作 $\hat v_{\theta1}$。因此：
$$
\boxed{\hat \psi_I = \hat v_{\theta 2} - \frac{R_1}{R_2} \hat v_{\theta 1}}\tag{1}
$$
# 2. 非定常总压降
在离心叶轮中，叶片通道长且弯曲，流体微团沿流线运动时存在明显的非定常加速度。这部分加速度对应的惯性力需要由叶片做功来提供（或反之，阻碍流动变化），表现为总压的额外变化。该修正用于连接 **稳态等熵总压升**（由欧拉方程给出）与 **实际总压升**（计入损失和非定常惯性）。

在 **相对坐标系** 中，沿一条流线 $s$（从进口 $s1$ 到出口 $s2$）积分非定常加速度，得到的总压降为：
$$
\Delta p_{t,\text{unsteady}} = -\int_{s1}^{s2} \left( \frac{\partial}{\partial t}\bigg|_{\text{rel}} \right) W\, ds
$$
其中：
- $W$ —— 相对速度的大小（有量纲，$\text{m/s}$）；
- $\tau = t \, U_2 / R_2$ —— 无量纲时间，以叶轮出口轮缘速度 $U_2 = \Omega R_2$ 和出口半径 $R_2$ 为参考；
- $\frac{\partial}{\partial\tau}\big|_{\text{rel}}$ —— 在旋转坐标系中观察的当地时间导数。
# 3 坐标变换：绝对坐标系中的表达
在小扰动分析中，扰动通常是在**绝对坐标系**中描述的。旋转坐标系与绝对坐标系的时间导数关系为：
$$
\frac{\partial}{\partial\tau}\bigg|_{\text{rel}} = \frac{\partial}{\partial t}\bigg|_{\text{abs}} + \Omega \frac{\partial}{\partial\theta}\bigg|_{\text{abs}}
$$
其中 $\theta$ 为周向坐标，$\frac{\partial}{\partial\theta}$ 项来源于旋转角速度 $\Omega$ 带来的效应。

如果两边都用$U_2$和$R_2$无量纲化，就有：
$$
\frac{\partial}{\partial\tau}\bigg|_{\text{rel}} = \frac{\partial}{\partial\tau}\bigg|_{\text{abs}} +\frac{\partial}{\partial\hat\theta}\bigg|_{\text{abs}}
$$
因此，非定常惯性总压降可写作：
$$
\boxed{
\Delta P_{t,\text{unsteady}} = -\int_{\hat s_1}^{\hat s_2} \left( \frac{\partial}{\partial\tau} + \frac{\partial}{\partial\hat\theta} \right) \hat W \, d\hat s \tag{2}
}
$$

在离心叶轮中，我们希望将分布积分简化为一个集总表达式，仅依赖于 **出口径向速度 $V_{r2}$**（无量纲）。这需要利用 **质量守恒（连续性方程）** 沿流管建立 $W(s)$ 与 $V_{r2}$ 的关系。
# 4 流管视角
基本假设：
- 流动沿流管一维化（截面上参数均匀）
- 流管截面积 $A(s)$ 定义为 **垂直于当地流线** 的横截面积
- 流管壁面为流线，无流体穿越
- 允许密度变化，但沿流线 $s$ 的密度分布 $\rho(s)$ 视为已知（从平均流计算或状态方程得到）

通过流管任意截面的 **质量流量** 为常数：
$$
\dot m = \rho(s) \, W(s) \, A(s) = \text{常数} \tag{3}
$$
其中：
- $W(s)$ —— 相对速度的大小（方向沿流线，因而垂直于 $A(s)$）
- $A(s)$ —— 垂直于流线的流管横截面积
- $\rho(s)$ —— 当地密度

在叶轮 **出口**（下标 2），我们通常采用 **径向截面**（法向沿径向）来度量流量。径向截面的质量为：
$$
\dot m = \rho_2 \, V_{r2} \, A_{2,\text{radial}} \tag{4}
$$
其中：
- $V_{r2}$ —— 出口径向速度
- $A_{2,\text{radial}}$ —— 出口处垂直于径向的环形面积（例如 $2\pi R_2 b_2$）

由于质量守恒，式 (1) 与 (2) 相等：

$$
\rho(s) \, W(s) \, A(s) = \rho_2 \, v_{r2} \, A_{2,\text{radial}} \tag{5}
$$
解得：
$$
W(s) = \frac{\rho_2 \, A_{2,\text{radial}}}{\rho(s) \, A(s)} \, v_{r2} \tag{6}
$$
定义 **面积‑密度比** 沿流线的分布，跟前面的系列不一样，这里$s$不是拉普拉斯变量，而是相对流线位置：
$$
AR(s) \equiv \frac{\rho_2 \, A_{2,\text{radial}}}{\rho(s) \, A(s)} \tag{7}
$$
顺带地把无量纲表达式给出来：
$$
\boxed{\hat W(s) = AR(\hat s) \cdot \hat v_{r2}} \tag {8}
$$
代回式(5)：
$$
\Delta P_{t,\text{unsteady}} = -\int_{\hat s_1}^{\hat s_2} \left( \frac{\partial}{\partial\tau} + \frac{\partial}{\partial\hat \theta} \right) \bigl( AR(\hat s) \, \hat v_{r2}(\tau,\hat \theta) \bigr) \, d\hat s
$$

由于偏导数只作用于 $\hat v_{r2}$，可以提到积分号外：
$$
\Delta P_{t,\text{unsteady}} = -\left( \frac{\partial}{\partial\tau} + \frac{\partial}{\partial\theta} \right) \hat v_{r2} \; \cdot \; \int_{\hat s_1}^{\hat s_2} AR(\hat s) \, \hat ds
$$

为了搞定这个积分，我们令：

$$
\lambda_{\text{imp}} = \int_{\hat s_1}^{\hat s_2} AR(\hat s) \, d\hat s \tag{9}
$$
设：
$$
\rho(s)A(s) = \rho_1 A_{1} + k\,(s - s_1), \qquad k = \frac{\rho_2 A_{2,\text{radial}} - \rho_1 A_{1}}{s_{\text{imp}}}
$$
其中，平均的流道长度$\hat s_{\text{imp}} = \hat s_2 - \hat s_1$。

将 $\rho(s)A(s)$ 代入 $AR(s)$ 表达式：
$$
AR(s) = \frac{\rho_2 A_{2,\text{radial}}}{\rho_1 A_{1} + k\,(s - s_1)} \tag{10}
$$
使用换元法，令 $\xi = s - s_1$，则 $\xi \in [0, s_{\text{imp}}]$，分母记作 $a + k\xi$，其中$a = \rho_1 A_{1}$：


于是：
$$
\lambda_{\text{imp}} = \int_0^{s_{\text{imp}}} \frac{\rho_2 A_{2,\text{radial}}}{a + k\xi} \, d\xi
$$
积分基本公式 $\int \frac{d\xi}{a+k\xi} = \frac{1}{k}\ln(a+k\xi)$，故：
$$
\lambda_{\text{imp}} = \rho_2 A_{2,\text{radial}} \cdot \frac{1}{k} \left[ \ln(a+k\xi) \right]_0^{s_{\text{imp}}}
$$
代入 $\xi = s_{\text{imp}}$：$a + b s_{\text{imp}} = \rho_2 A_{2,\text{radial}}$；$\xi = 0$：$a = \rho_1 A_{1}$。因此：
$$
\lambda_{\text{imp}} = \frac{\rho_2 A_{2,\text{radial}}}{k} \ln\left( \frac{\rho_2 A_{2,\text{radial}}}{\rho_1 A_{1}} \right)
$$
将 $k = \frac{\rho_2 A_{2,\text{radial}} - \rho_1 A_{1}}{s_{\text{imp}}}$ 代入：

$$
\lambda_{\text{imp}} = \rho_2 A_{2,\text{radial}} \cdot \frac{s_{\text{imp}}}{\rho_2 A_{2,\text{radial}} - \rho_1 A_{1}} \cdot \ln\left( \frac{\rho_2 A_{2,\text{radial}}}{\rho_1 A_{1}} \right)\tag{11}
$$

定义 **总面积‑密度比**（进口到出口）：
$$
AR_{\text{imp}} = \frac{\rho_2 A_{2,\text{radial}}}{\rho_1 A_{1}}
$$
代入上式：
$$
\lambda_{\text{imp}} = s_{\text{imp}} \cdot \frac{AR_{\text{imp}} \ln(AR_{\text{imp}})}{AR_{\text{imp}} - 1}\tag{12}
$$
所以:
$$
\Delta P_{t,\text{unsteady}} = -\lambda_{\text{imp}}\left( \frac{\partial \hat v_{r2}}{\partial\tau} + \frac{\partial \hat v_{r2}}{\partial\theta} \right) 
$$
上式中的面积计算方式如下：
$$
A_1 =h_1\cdot\frac{2\pi R_1}{N_{main}}\qquad{\space}\qquad{\space}
A_2 =h_2\cdot\left(\frac{2\pi R_2}{N_{main}+N_{splitter}}-t_{imp}\right)
$$
这里$h$是叶高、$N$是叶片数量、$t_{imp}$是叶片尾缘的厚度。这很明显的近似，入口处咱就没考虑叶片厚度。开心可以把前缘厚度也算上的。另外，这里直接用的环形面积，其实也可以采用流道真实面积，多乘一个$cos\chi$就好。
# 5. 匹配条件
由连续性得到：
$$
\hat v_{x1}=AR_{\text{imp}}\hat v_{r2}\tag{13}
$$
由滑移因子得到：
$$
\hat v_{\theta2}=\sigma_W+\hat v_{r2}\cdot\tan\chi_2=1+\hat v_{r2}\cdot tan\beta_2\tag{14}
$$
这里的1是因为圆周速度$R_2\Omega$被归一化了(这里的归一化标准也和前面几份文档的不一样，使用时务必统一），还有注意金属角和气流角不同。

结合式(4)和式(17)，考虑上我们的损失和非定常需要的额外总压降，我们有：
$$
\hat p_{t2}-\hat p_{t1}=1+\hat v_{r2}\cdot\tan\beta_2-\frac{R1}{R2}\hat v_{x1}\cdot \tan\alpha_1 -\hat l_{imp} -\lambda_{imp}\left( \frac{\partial}{\partial\tau} + \frac{\partial}{\partial\theta} \right)\hat v_{r2}\tag{15}
$$

和轴流转子一样，这里的损失可以表示成：
$$
\left[1+\tau_{imp }\cdot\left(\frac{\partial}{\partial \tau}+\frac{\partial}{\partial\hat\theta}\right)\right]\hat l_{imp}=l_{imp}^{ss}(tan\beta_1)\tag{16}
$$
类似地，这里的一阶迟滞估算由下式给出：
$$
\tau_{imp}=\tau_u\cdot\frac{2s_{imp}}{\hat{W}_1+\hat{W}_2}\tag{17}
$$
# 6. 小扰动
## 6.1 传递式展开
展开式（18）的总压：
$$
\hat p_{2}+\frac{1}{2}\hat v_{\theta 2}^2+\frac{1}{2}\hat v_{r 2}^2
-\hat p_{1}-\frac{1}{2}\hat v_{\theta 1}^2-\frac{1}{2}\hat v_{x 1}^2
=
1+\hat v_{r2}\cdot\tan\beta_2-\frac{R1}{R2}\hat v_{x1}\cdot \tan\alpha_1 -\hat l_{imp} -\lambda_{imp}\left( \frac{\partial}{\partial\tau} + \frac{\partial}{\partial\theta} \right)\hat v_{r2}
$$
# 6.2 线性化
进行线性化：
$$
\hat p_{2}+\hat {\overline v}_{\theta 2}\delta\hat v_{\theta 2}+\hat {\overline v}_{r 2}\delta\hat v_{r 2}
-\hat p_{1}-\hat {\overline v}_{\theta 1}\delta\hat v_{\theta 1}-\hat {\overline v}_{x 1}\delta\hat v_{x 1}
=
\delta\hat v_{r2}\cdot\tan\beta_2-\frac{R1}{R2}\delta\hat v_{x1}\cdot \tan\alpha_1 -\delta\hat l_{imp} -\lambda_{imp}\left( \frac{\partial}{\partial\tau} + \frac{\partial}{\partial\theta} \right)\delta\hat v_{r2}
$$
然后，我们只希望保留$\hat v_{x1}$和$\hat v_{r2}$作为传递关系，所以利用前面的速度三角形关系做点变换：
$$
\hat p_{2}+\frac{\hat {\overline v}_{\theta 2}\tan{\beta_2}}{AR_{imp}}\delta\hat v_{x 1}+\frac{\hat {\overline v}_{r 2}}{AR_{imp}}\delta\hat v_{x 1}
-\hat p_{1}-\hat {\overline v}_{\theta 1}\delta\hat v_{\theta 1}-\hat {\overline v}_{x 1}\delta\hat v_{x 1}
=
\frac{\tan\beta_2}{AR_{imp}}\delta\hat v_{x1}-\frac{R1}{R2}\delta\hat v_{x1}\cdot \tan\alpha_1 -\delta\hat l_{imp} -\frac{\lambda_{imp}}{AR_{imp}}\left( \frac{\partial}{\partial\tau} + \frac{\partial}{\partial\theta} \right)\delta\hat v_{x1}
$$
这样我们就完全消除了得到$\hat p_2$式中所需的所有“激盘出口参数”，让我们取拉普拉斯变换，在周向进行傅里叶展开，并把视角聚焦到第$n$阶谐波上，同样$e^{jn\theta}$暂时约掉，重构物理场时再拿回来，有：
$$
\delta\tilde P_{2}+\frac{\hat {\overline v}_{\theta 2}\tan{\beta_2}}{AR_{imp}}\delta\tilde V_{x 1}+\frac{\hat {\overline v}_{r 2}}{AR_{imp}}\delta\tilde V_{x 1}
-\delta\tilde P_{1}-\hat {\overline v}_{\theta 1}\delta\tilde V_{\theta 1}-\hat {\overline v}_{x 1}\delta\tilde V_{x 1}
=
\frac{\tan\beta_2}{AR_{imp}}\delta\tilde V_{x1}-\frac{R1}{R2}\delta\tilde V_{x1}\cdot \tan\alpha_1 -\delta\tilde L_{imp} -\frac{\lambda_{imp}}{AR_{imp}}\left( s+jn\right)\delta\tilde V_{x1}
$$
损失项可以直接用轴流转子的解：
$$
\delta\tilde{L}_R =  \frac{d\hat{l}_{imp}^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}}\cdot \frac{   \delta\tilde{V}_{\theta1} - \tan\beta_1 \delta\tilde{V}_{x1}} {\hat{\overline v}_{x1}[1 + \tau_{imp}(s+jn)]}
$$
这里和轴流转子不一样的地方在于，参考的解变成（链式法则）：
$$
\frac{d\hat{l}_{imp}^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}}
=
\frac{d\hat{l}_{imp}^{ss}}{d\hat{\overline v}_{r2}}\cdot\frac{R_1}{R_2}\frac{1}{(\tan{\alpha_1}-\tan{\beta_1})^2AR_{imp}}
$$
这样，整个式子就变成：
$$
\delta\tilde P_{2}+\frac{\hat {\overline v}_{\theta 2}\tan{\beta_2}}{AR_{imp}}\delta\tilde V_{x 1}+\frac{\hat {\overline v}_{r 2}}{AR_{imp}}\delta\tilde V_{x 1}
-\delta\tilde P_{1}-\hat {\overline v}_{\theta 1}\delta\tilde V_{\theta 1}-\hat {\overline v}_{x 1}\delta\tilde V_{x 1}
=
\frac{\tan\beta_2}{AR_{imp}}\delta\tilde V_{x1}-\frac{R1}{R2}\delta\tilde V_{x1}\cdot \tan\alpha_1 - \frac{d\hat{l}_{imp}^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}}\cdot \frac{   \delta\tilde{V}_{\theta1} - \tan\beta_1 \delta\tilde{V}_{x1}} {\hat{\overline v}_{x1}[1 + \tau_{imp}(s+jn)]} -\frac{\lambda_{imp}}{AR_{imp}}\left( s+jn\right)\delta\tilde V_{x1}\tag{21}
$$
# 7. 传递矩阵
输入向量（进口 STA 1）：  
$$
\begin{bmatrix} \delta\hat V_{r2} \\ \delta\hat V_{\theta2} \\ \delta\hat P_2 \end{bmatrix}(\hat\theta,s)
=
\sum_{n=0}^\infty\mathbb{B}_{\text{imp},n}\cdot e^{jn\theta}\cdot\begin{bmatrix} \delta\tilde V_{x1} \\ \delta\tilde V_{\theta1} \\ \delta\tilde P_1 \end{bmatrix}(s)
$$

传递矩阵：
$$
\boxed{
\mathbb{B}_{\text{imp},n} =
\begin{bmatrix}
\dfrac{1}{AR_{\text{imp}}} & 0 & 0 \\
\dfrac{\tan\beta_2}{AR_{\text{imp}}} & 0 & 0 \\
B_{31} & B_{32} & 1
\end{bmatrix}
}
$$
其中
$$
\begin{aligned}
B_{31} &= \frac{\tan\beta_2-\lambda_{\text{imp}}(s+jn)-\hat{\overline v}_{r2} - \hat{\overline v}_{\theta2}\tan\beta_2}{AR_{\text{imp}}}
        - \frac{R_1}{R_2}\tan\alpha_1
        - \hat{\overline v}_{x1}
        + K_{\text{imp}}\tan\beta_1,\\
B_{32} &= \hat{\overline v}_{\theta1} - K_{\text{imp}},\\
K_{\text{imp}} &= \frac{1}{\hat{\overline v}_{x1}}\,
            \frac{d\hat l_{\text{imp}}^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}}
                \frac{1}{1 + \tau_{\text{imp}}(s+jn)}
\end{aligned}
$$