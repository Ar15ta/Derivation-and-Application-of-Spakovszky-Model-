|  属性  | [[法术表]]  |                  |               |
| :--: | :------: | :--------------: | :-----------: |
|  关联  | [[欧拉方程]] | [[连续性方程——比体积形式]] |  [[欧拉涡轮方程]]   |
| 链接解释 |    基础    |        约束        | 理论的稳态压升特性与此有关 |
***
#  $Model\space Ⅰ$

`该模型将整个压气机视为一个集中的半激盘，是以往线性动力系统分析的正统延续。`

|  变量名   |        符号        |  单位   |        量纲         |                                      说明                                       |
| :----: | :--------------: | :---: | :---------------: | :---------------------------------------------------------------------------: |
|  出口静压  |      $p_2$       |  Pa   | $M L^{-1} T^{-2}$ |                              压气机出口（最后一个叶片排下游）静压                               |
|  进口总压  |     $p_{t1}$     |  Pa   | $M L^{-1} T^{-2}$ |                               压气机进口（第一个叶片排上游）总压                               |
|   密度   |      $\rho$      | kg/m³ |    $M L^{-3}$     |                                 气流密度（假设不可压缩）                                  |
| 叶片切向速度 |       $U$        |  m/s  |    $L T^{-1}$     |                                  转子叶片中径处的线速度                                  |
| 总静压升系数 |   $\psi^{ts}$    |   –   |         1         |                     稳态总‑静压升特性，通常由实验或计算给出，是流量系数 $\phi$ 的函数                     |
|  流量系数  | $\phi = v_x / U$ |   –   |         1         |                                气流轴向速度与叶片切线速度之比                                |
| 转子惯性因子 |    $\lambda$     |   –   |         1         |  所有转子叶片排的流体惯性之和，$\lambda = \sum_{\text{rotors}} \frac{c_x}{R \cos^2\gamma}$   |
| 总惯性因子  |      $\mu$       |   –   |         1         | 所有叶片排（转子+静子）的流体惯性之和，$\mu = \sum_{\text{all rows}} \frac{c_x}{R \cos^2\gamma}$ |
|  轴向弦长  |      $c_x$       |   m   |        $L$        |                                   叶片排的轴向弦长                                    |
|  平均半径  |       $R$        |   m   |        $L$        |                                     叶片排中径                                     |
|  安装角   |     $\beta$      |  rad  |         –         |                                 叶片排的安装角（无因次）                                  |
| 无因次时间  | $\tau = t U / R$ |   –   |         1         |                                 以转子转动为特征的时间尺度                                 |
|  周向角度  |     $\theta$     |  rad  |         –         |                                     周向坐标                                      |

***
# 前提
## 几何与运动学假设
| 假设         | 说明                                                      |
| :--------- | :------------------------------------------------------ |
| **高轮毂比**   | 假设轮毂比足够高，径向流动效应可忽略，流动可视为二维 $(x, \theta)$                |
| **恒定环形面积** | 上下游管道截面积不变，无突扩或突缩                                       |
| **无轴向间隙**  | Model I 假设各叶片排紧密耦联，无叶排间轴向间隙。Spakovszky Model II 中此假设被放宽 |
| **无限长管道**  | 上下游管道足够长，非均匀压力场在远端完全衰减，不与远场边界干涉                         |
## 流体性质假设
| 假设          | 说明                            |
| :---------- | :---------------------------- |
| **不可压缩流动**  | 马赫数足够低，密度 $\rho$ 为常数，声学效应可忽略  |
| **无粘（管道中）** | 管道流动中忽略粘性；粘性效应仅在叶片排内以损失系数形式体现 |
| **忽略热传导**   | 管道中无热传导，气流总温不变                |
| **忽略雷诺数效应** | 假设损失特性等对雷诺数不敏感                |
## 叶片排建模假设
| 假设          | 说明                                                  |
| :---------- | :-------------------------------------------------- |
| **半激盘**     | 叶片排被压缩为无厚度的作用盘，进出口匹配由代数关系描述                         |
| **长波长**     | 扰动周向波长 $\gg$ 叶片栅距。由 Longley(1994) 的论证，直至第5～6阶谐波仍可满足 |
| **准稳态叶片响应** | 低约化频率下，叶片排的气流转折和损失可瞬时响应局部流动条件。此假设可通过一阶迟滞修正          |
| **完美导叶**    | 进口导叶出口角恒定，不受进口流动条件影响                                |
| **高稠度出口叶片** | 末排叶片出口角恒定，出口旋流不受扰动影响                                |
| **损失可加性**   | 在设计反动度下，级总损失可分解为转子损失和静子损失的线性叠加                      |
## 扰动假设
| 假设        | 说明                                                                                     |
| :-------- | :------------------------------------------------------------------------------------- |
| **小扰动**   | 扰动幅度足够小，略去二阶及以上非线性项                                                                    |
| **模态解形式** | 扰动可表示为 $\delta \phi = \Re\{a e^{in\theta + \sigma \tau + i\omega \tau}\}$，各谐波在线性阶段独立演化 |
| **上游无旋**  | 上游流场无涡量来源，允许使用势函数描述                                                                    |
| **下游有旋**  | 下游流场包含从叶片排脱落的涡量，必须用欧拉方程描述                                                              |
## 系统组件假设
| 序号  | 假设          | 说明                           |
| :-- | :---------- | :--------------------------- |
| S1  | **集总集气室**   | 集气室容积大，内部压力空间均匀，仅为时间的函数      |
| S2  | **等熵集气室过程** | 集气室内压力变化是等熵的                 |
| S3  | **平方律节气门**  | 节气门压降与当地速度平方成正比，$k_T$ 为节气门系数 |
# 1. 压气机

## 1.1 单叶片排的非定常惯性效应
![[Spakovszky稳定性方法——模型Ⅰ配图.png]]
考虑一个叶片排，其叶片通道可近似为一条长为 $c_x / \cos \beta$ 的平行管道，流体的平均相对速度为 $v_x / \cos \beta$。当叶片排经历非均匀（这里特指转子旋转坐标系带来的周向非均匀）、非定常流动时，通道内流体加速所需的额外压差可表示为：
$$
\left( p_{\text{in}} - p_{\text{out}} \right)_{\text{inertia}} = \rho \cdot (\text{通道长度}) \cdot \frac{\partial}{\partial t} (\text{相对速度})
$$
如果没能理解，对于一维分析，取一个控制体作为微元，此时通道长度变成微元长度。式两边同乘以通流面积，就有==额外的惯性压力差$*$面积 = 控制体质量$*$非定常效应引起的额外相对加速度==，是不是就很明显了呢！

代入几何和速度关系，得到：
$$
\left( \frac{p_{\text{in}} - p_{\text{out}}}{\rho U^2} \right)_{\text{inertia}} = \Lambda \cdot \left. \frac{\partial}{\partial \tau} \left( \frac{v_x}{U} \right) \right|_{\text{rel to blade}}
$$
其中 $\Lambda = \frac{c_x}{R \cos^2\gamma}$ 为单排惯性参数，$\tau = t U / R$ 为无因次时间。

对于**转子**，相对坐标系中的时间导数变换到绝对坐标系需加上周向对流项（周向非均匀就在这）：
$$
\left. \frac{\partial}{\partial \tau} \right|_{\text{rel}} = \left. \left(\frac{\partial}{\partial \tau} + \frac{\partial}{\partial \theta}\right)\right|_{\text{abs}}
$$
对于**静子**，相对坐标系即为绝对坐标系，时间导数保持 $\partial / \partial \tau$。
## 1.2 叶片排稳态性能的描述
稳态时，叶片排的总静压升由其气动特性决定，用 $\psi^{ts}(\phi)$ 表示（为简化，这里假设进口无预旋或采用完美导叶，使特性仅依赖于 $\phi$）。稳态总-静压升定义为：
$$
\psi^{ts}(\phi) = \frac{p_2 - p_{t1}}{\rho U^2} \bigg|_{\text{steady}}
$$
## 1.3 全压气机的合成
若忽略叶片排之间的轴向间隙（即假设所有叶片排紧密耦联），则各排进出口的轴向速度相同，扰动在周向均匀传播。将所有叶片排的稳态压力升与惯性压差代数相加，得到压气机瞬时总静压升方程：
$$
\frac{p_2 - p_{t1}}{\rho U^2} = \psi^{ts}(\phi)-\; \underbrace{\sum_{\text{rotors}} \Lambda_{\text{row}} \left( \frac{\partial \phi}{\partial \tau} + \frac{\partial \phi}{\partial \theta} \right) \;-\; \sum_{\text{stators}} \Lambda_{\text{row}} \frac{\partial \phi}{\partial \tau}}_{\text{inertia terms}}
$$
将同类项合并，定义：
$$
\lambda = \sum_{\text{rotors}} \Lambda_{\text{row}}, \qquad 
\mu = \sum_{\text{rotors}} \Lambda_{\text{row}} + \sum_{\text{stators}} \Lambda_{\text{row}}
$$
最终得到：
$$
\boxed{\frac{p_2 - p_{t1}}{\rho U^2} = \psi^{\mathrm{ts}}(\phi) - \lambda \frac{\partial\phi}{\partial\theta} - \mu \frac{\partial\phi}{\partial\tau}}
$$
该式明确表明：压气机出口静压与进口总压之差（即实际总静压升）由两部分构成——稳态特性提供的压力升，以及因非定常流动引起的惯性压降。后者取决于扰动流量系数的周向梯度（与转子惯性 $\lambda$ 相关）和时间梯度（与总惯性 $\mu$ 相关）。这一方程是压气机动态稳定性分析的核心，将失速先兆波的传播与增长同压气机特性及其几何参数直接联系了起来。
## 1.4 稳态特性的进一步分解：理想增压与损失
 $\psi^{ts}(\phi)$ 并非一个纯粹的经验黑箱。在Model I中，它被拆解为理想等熵部分与粘性损失两部分，以便单独捕捉非定常响应。
$$
\psi^{\mathrm{ts}}(\phi) = \psi_{I}^{\mathrm{ts}}(\phi) - L_R(\phi) - L_S(\phi)
$$
- $\psi_{I}^{\mathrm{ts}}(\phi)$ —— **理想总静压升系数**，对应气流在叶片排中等熵转折时所能达到的最大压力升。
- $L_R(\phi), L_S(\phi)$ —— 转子、静子的**总压损失系数**，分别刻画了转子与静子中的真实粘性耗散。

这一分离的物理动机在于：理想压力升源于气流转折（欧拉功），其变化几乎是瞬时的；而损失主要源于边界层分离等粘性过程，对流动变化存在明显的时间滞后。将 $L_R、L_S$ 单独写出，为后续引入损失时间滞后提供了直接接口。讨论系统稳定性时，损失特性斜率 $\partial L / \partial \phi$ 将直接决定系统的阻尼。
## 1.5 一阶迟滞
单纯的准稳态模型假设损失能瞬间响应流动变化，但这在高频扰动下过于理想化。因此，引入了一阶迟滞。其核心思想是：**叶片排的总压损失变化，并非即刻完成，而是以指数形式弛豫到一个新的稳态值**。

该过程用一个一阶常微分方程来精确描述：
$$
\left( 1 + \tau_{\text{loss}} \cdot \frac{\partial}{\partial \tau} \right) L = L^{\text{ss}}
$$
其中：
- $L$ 是当前的真实损失值。
- $L^{\text{ss}}$ 是基于当前流场条件查表得到的稳态目标损失。
- $\tau_{\text{loss}}$ 是损失响应的时间常数，可以由基流和特征尺度估计，代表了弛豫时间。

此外，记反力度为$\sigma$，可以得到：
$$
\left( 1 + \tau_{\text{ROTORloss}} \cdot (\frac{\partial \phi}{\partial \tau} + \frac{\partial \phi}{\partial \theta}  )\right)L_R(\phi) 
=
\sigma(\psi_{I}^{\mathrm{ts}}(\phi) - \psi^{\mathrm{ts}}(\phi))
$$
$$
\left( 1 + \tau_{\text{STATORloss}} \cdot \frac{\partial \phi}{\partial \tau}  \right)L_S(\phi) 
=
(1-\sigma)\left(\psi_{I}^{\mathrm{ts}}(\phi) - \psi^{\mathrm{ts}}(\phi)\right)
$$

这个方程的含义是：稳态损失 $L^{\text{ss}}$ 是"目标"，而真实的损失 $L$ 以 $\tau_{\text{loss}}$ 为时间尺度，**逐渐迫近**这个目标。这极大地增强了模型对非定常效应的捕捉能力。到这里，惯性和黏性导致的迟滞两个非定常部分就都建模完毕了

# 2. 管道流场与集气室动力学
压气机自身的特性必须与上下游管道中的流场相匹配。换句话说，**压气机产生的压力差，必须等于管道流场中对应两点之间的实际压力差**。因此，需要建立上下游管道以及集气室的流动模型。

核心假设：
- 管道中的流动是**不可压、无粘、二维**的（轴向 $x$ 和周向 $\theta$）。
- 上下游管道足够长，压气机产生的非均匀压力场在远前方和远后方完全衰减，不与远端边界发生干涉。
## 2.1 上游管道：无旋势流
上游来流均匀且定常，管道中没有涡量来源（因为前面假设了叶片排以外没有黏性）。根据开尔文定理，整个上游区域的流动是**无旋的**。

无旋流动的速度可以表示为一个势函数 $\Phi$ 的梯度：
$$ v_x = \frac{\partial \Phi}{\partial x}, \quad v_\theta = \frac{1}{R} \frac{\partial \Phi}{\partial \theta} $$
代入不可压连续性方程 $\nabla \cdot \mathbf{v} = 0$：
$$ \frac{\partial v_x}{\partial x} + \frac{1}{R} \frac{\partial v_\theta}{\partial \theta} = 0 $$
即得到拉普拉斯方程：
$$ \frac{\partial^2 \Phi}{\partial x^2} + \frac{1}{R^2} \frac{\partial^2 \Phi}{\partial \theta^2} = 0   $$
**物理作用**：这是上游流场的控制方程。结合"扰动在上游无限远处衰减至零"的边界条件，可求解出上游管道中任意位置的速度和压力分布。对于从压气机向上游传播的势流扰动，该方程决定了它们以指数形式衰减（$e^{n x/R}$）的规律。这意味着，高阶谐波（$n$ 越大）的扰动衰减越快，因此上游传感器只能检测到低阶模态。
## 2.2 下游管道：有旋欧拉方程
下游情况与上游截然不同。气流经过压气机叶片排后，由于叶片载荷沿周向不均匀，会从叶片尾缘**脱落涡量**。因此，下游管道中的流动是**有旋的**，不能再用简单的势函数来描述。

需要直接使用描述无粘流动的**二维欧拉方程**：

**轴向动量方程**：

$$ \frac{\partial v_x}{\partial t} + v_x \frac{\partial v_x}{\partial x} + \frac{v_\theta}{R} \frac{\partial v_x}{\partial \theta} = -\frac{1}{\rho} \frac{\partial p}{\partial x} $$

**周向动量方程**：

$$ \frac{\partial v_\theta}{\partial t} + v_x \frac{\partial v_\theta}{\partial x} + \frac{v_\theta}{R} \frac{\partial v_\theta}{\partial \theta} = -\frac{1}{\rho} \frac{\partial p}{\partial \theta} $$

**连续性方程**：

$$ \frac{\partial v_x}{\partial x} + \frac{1}{R} \frac{\partial v_\theta}{\partial \theta} = 0 $$

**物理作用**：这是下游流场的完整控制方程组。与上游不同，这里的解包含两类扰动：
1.  **势流扰动**：同样以指数形式向下游衰减（$e^{-n x/R}$），在远下游消失。
2.  **涡量扰动**：由压气机脱落，随主流向下游对流，**不衰减**。

这两类扰动的叠加，使下游管道中的"压力-速度"关系远比上游复杂。涡量扰动的存在，是下游流场能够携带能量、影响系统稳定性的关键原因。
## 2.3 集气室：集总参数动力学
压气机下游连接着**集气室**（plenum），其几何尺寸远大于管道，内部流动可近似为**空间均匀**——压力 $p_p$ 处处相同，仅是时间的函数。

集气室如同一个"气动弹簧"：当流入和流出的质量流量不平衡时，其内部压力就会发生变化。假设过程是**等熵**的：

$$ \frac{dp_p}{dt} = \frac{\gamma p_a}{V_P} \left( A_C v_{x3} - A_T v_{x4} \right) \quad $$
# 3. 界面耦合匹配条件
要将各组件（上游管道、压气机、下游管道、集气室）的控制方程组合成一个封闭的整机模型，必须在组件间的界面处施加物理匹配条件。在 $Model I$ 中，压气机被集成为一个整体半激盘，因此只需在其**进口站 1** 和**出口站 2** 与管道流场之间建立匹配关系。
## 3.1 压气机进口匹配（上游管道 → 站 1）
上游管道中的流动是**无旋势流**，由拉普拉斯方程 (2.5) 描述。在站 1 处，管道出口的流动变量必须与压气机进口变量相匹配。
### 3.1.1 速度连续（质量守恒）
轴向速度连续：
$$ v_{x1}^{\text{duct}} = v_{x1}^{\text{comp}} \equiv v_{x1} $$
### 3.1.2 总压匹配
上游管道中的非定常伯努利方程给出扰动总压与扰动速度势的关系。对于 $n$ 阶谐波，衰减的势流解给出[^1]：
$$ \frac{\delta p_{t1}}{\rho U^2} = -\frac{1}{|n|} \frac{\partial}{\partial \tau} (\delta \phi) $$
该关系表明：上游管道中的压力扰动与流量扰动的时间变化率成正比，且与谐波数成反比。高阶谐波在上游衰减更快，因而对上游流场的影响更弱。
### 3.1.3 进口旋流角
在“完美导叶”的简化假设下，进口绝对旋流角 $\alpha_1$ 为常数，其扰动为零：
$$ \alpha_1 = \text{const} \quad \text{或} \quad \delta v_{\theta 1} = 0 $$
若无导叶，则为轴向进气 $\alpha_1 = 0$。
## 3.2 压气机出口匹配（站 2 → 下游管道）
在站 2 处，压气机出口变量必须与下游管道进口变量匹配。
### 3.2.1 速度连续（质量守恒）
轴向速度连续：
$$ v_{x2}^{\text{comp}} = v_{x2}^{\text{duct}} \equiv v_{x2} $$
### 3.2.2 出口旋流角
在“高稠度出口叶片”的简化假设下，出口绝对旋流角 $\alpha_2$ 为常数，其扰动为零：
$$ \alpha_2 = \text{const} \quad \text{或} \quad \delta \alpha_2 = 0 $$
此假设意味着出口切向速度扰动完全由轴向速度扰动决定：$\delta v_{\theta 2} = \delta v_{x2} \tan \alpha_2$。
### 3.2.3 静压匹配
压气机出口静压等于下游管道进口静压：
$$ p_2^{\text{comp}} = p_2^{\text{duct}} \equiv p_2 $$
在下游管道中，有旋欧拉方程的解给出静压扰动与流量扰动的关系。对于高稠度出口假设，该关系为：
$$ \frac{\delta p_2}{\rho U^2} = \frac{1}{|n|} \frac{\partial}{\partial \tau} (\delta \phi) $$
注意，与上游总压关系相比，下游静压关系**符号相反**。此外，下游还存在不衰减的涡量扰动，其幅值由叶片排脱落涡量决定。
### 3.2.4 脱落涡量
压气机叶片排脱落的总涡量由叶片载荷的周向分布决定。涡量强度由环量差确定：
$$ \Gamma_{\text{shed}} = 2\pi R (v_{\theta 2} - v_{\theta 1}) $$
在线化形式中，脱落涡量扰动与流量扰动直接相关。
## 3.3 集气室界面匹配
### 3.3.1 下游管道出口 → 集气室进口（站 3）
下游管道出口静压等于集气室压力（集气室内空间均匀）：
$$ p_3 = p_p $$
在远下游，势流扰动已完全衰减，此处只有涡量扰动起作用。该压力条件为下游欧拉方程提供了必要的出口边界条件。
### 3.3.2 集气室出口 → 节气门（站 4）
集气室出口速度 $v_{x4}$ 通过节气门特性与集气室压力 $p_p$ 及环境背压 $p_6$ 关联：
$$ p_p - p_6 = \frac{1}{2} k_T v_{x4}^2 $$
质量连续（不可压）给出：
$$ A_C v_{x3} = A_T v_{x4} $$
其中 $A_C$、$A_T$ 分别为压气机下游管道和节气门管道的截面积。

# 4. 完整求解流程
从一台具体的压气机出发，到最终判断其稳定性，需要经历三个步骤：**获取特性曲线**、**求解稳态工作点**、**小扰动特征值分析**。
## 第一步：获取总静压升特性曲线 $ψ^{ts}(φ)$
特性曲线是模型的**核心输入**。获取它的方法有三种，精度逐级递增：
### 1、基于速度三角形的理想特性（零损失）
即使没有任何实验或CFD数据，仅凭几何设计参数，就可以画出压气机的"骨架"——速度三角形，并推导出理想（等熵）总静压升特性。

**转子出口速度三角形（考虑滑移）**：
$$ V_{\theta 2} = \sigma U_2 - V_{x2} \tan \beta_{2b} $$
其中 $\sigma_{Slice}$ 为滑移因子（Wiesner公式：$\sigma_{Slice} = 1 - \sqrt{\cos \beta_{2b}} / N^{0.7}$），$N$ 为叶片数，$\beta_{2b}$ 为出口金属角。

**欧拉功**：
$$ \Delta h_t = U_2 V_{\theta 2} - U_1 V_{\theta 1} $$
**无因次理想总静压升**（假设进口无预旋 $V_{\theta 1} = 0$）：
$$ \psi_I^{ts}(\phi) = \frac{\Delta h_t - \frac{1}{2}V_2^2}{U_2^2} $$
代入速度三角形关系后，$\psi_I^{ts}$ 是 $\phi$ 的**线性函数**，且斜率为负：
$$ \boxed{\psi_I^{ts}(\phi) = A - B \phi} $$
其中 $A$、$B$ 是由叶片角和滑移因子决定的常数。纯理想压气机天生是**气动稳定**的（特性线向右下倾斜）。是粘性损失改变了这一斜率。
### 2、引入经验损失模型：从理想走向真实
真实的压气机存在粘性损失。在无CFD和实验的条件下，可用经典经验关联式

**总压损失系数 $\omega$ 的Koch & Smith模型**：
总压损失系数定义为 $\omega = \Delta p_t / (0.5 \rho V_1^2)$。对于设计工况附近：
$$ \omega \approx \omega_{\min} \cdot \left[ 1 + K \cdot (i - i_{\text{opt}})^2 \right] $$
其中 $i = \beta_1 - \beta_{1b}$ 是攻角，$i_{\text{opt}}$ 是最小损失攻角，$K$ 是经验常数，可从公开文献获取。

**总压损失的无因次形式**：
$$ L_R(\phi) = \frac{\Delta p_{t,\text{rotor}}}{\rho U^2} = \omega \cdot \frac{1}{2} \phi^2 (1 + \tan^2 \beta_1) $$
**真实总静压升特性**：
$$ \psi^{ts}(\phi) = \psi_I^{ts}(\phi) - L_R(\phi) - L_S(\phi) $$
在高流量时，攻角合理，损失小，$\psi^{ts}$ 接近 $\psi_I^{ts}$，斜率为负。
当流量减小时，攻角增大，损失迅速增加，$L_R + L_S$ 的斜率超过 $\psi_I^{ts}$ 的下降速度，导致 **$\psi^{ts}$ 曲线向下弯曲，斜率从负变正**。

**斜率由负变正的转折点，就是模型预测的失速边界**：
$$ \frac{d\psi^{ts}}{d\phi} = 0 \quad \text{（失速临界点）} $$
### 3、利用实验或CFD数据（最精确）
若条件允许，直接从实验或CFD提取特性曲线及其导数，精度最高。
## 第二步：稳态工作点求解（耦合迭代）
在获得 $\psi^{ts}(\phi)$ 曲线后，需要在给定的系统条件下求出稳态工作流量 $\overline{\phi}$ 和对应的稳态压力 $\overline{p}_p$。

**已知量**：
- 上游大气总压 $p_0$，下游环境背压 $p_6$
- 节气门系数 $k_T$
- 特性曲线 $\psi^{ts}(\phi)$、损失特性 $L_R(\phi)$、$L_S(\phi)$
- 惯性参数 $\lambda$、$\mu$

**稳态假设**：所有时间导数 $\partial/\partial\tau = 0$，所有周向导数 $\partial/\partial\theta = 0$。

**控制方程**：
$$ \begin{cases} \psi^{ts}(\phi) = \dfrac{p_2 - p_{t1}}{\rho U^2} \quad \text{（压气机总静压升）} \\ p_{t1} = p_0 - \dfrac{1}{2}\rho v_x^2 \quad \text{（上游伯努利）} \\ p_2 = p_p \quad \text{（下游管道无损失，出口静压等于集气室压力）} \\ p_p - p_6 = \dfrac{1}{2} k_T v_x^2 \quad \text{（节气门特性）} \end{cases} $$

**迭代求解流程**：
1. **猜测**：给定初始流量系数 $\phi^{(0)}$
2. **上游计算**：$p_{t1} = p_0 - 0.5 \rho (\phi U)^2$
3. **压气机计算**：$p_2 = p_{t1} + \rho U^2 \cdot \psi^{ts}(\phi)$
4. **下游计算**：$p_p = p_2$
5. **节气门反算**：由 $p_p - p_6 = 0.5 k_T (\phi_T U)^2$，解出 $\phi_T$
6. **收敛判断**：若 $|\phi - \phi_T| < \epsilon$，则收敛，稳态工作点为 $\overline{\phi} = \phi$；否则，更新 $\phi$ 并返回步骤2（例如采用二分法或牛顿法）

**输出**：稳态工作点流量 $\overline{\phi}$ 和集气室压力 $\overline{p}_p$。

**重要**：在此稳态工作点处，求出所有特性曲线的**导数值**——$d\psi^{ts}/d\phi$、$\partial L_R/\partial\phi$、$\partial L_S/\partial\phi$。这些导数值将直接输入下一步的特征值计算。
## 第三步：小扰动稳定性分析（特征值求解）
在稳态工作点附近，假设存在小扰动，分析其随时间增长或衰减，以及周向传播的特性。

**扰动假设**：
$$ \phi(\theta, \tau) = \overline{\phi} + \delta\phi(\theta, \tau), \quad \delta\phi = \Re\{ a e^{in\theta + s\tau} \} $$
其中 $s = \sigma + i\omega$ 为复特征值。
### 1——线化压气机方程
线化形式为：
$$ \frac{\delta p_2}{\rho U^2} - \frac{\delta p_{t1}}{\rho U^2} = \frac{d\psi^{ts}}{d\phi} \delta\phi - \lambda \frac{\partial}{\partial\theta}(\delta\phi) - \mu \frac{\partial}{\partial\tau}(\delta\phi) $$
### 2——代入上下游阻抗关系
**上游**（衰减势流解，见上游总压匹配推导）：
$$ \frac{\delta p_{t1}}{\rho U^2} = -\frac{1}{|n|} \frac{\partial}{\partial\tau}(\delta\phi) $$
**下游**（有旋流，高稠度出口假设）：
$$ \frac{\delta p_2}{\rho U^2} = \frac{1}{|n|} \frac{\partial}{\partial\tau}(\delta\phi) $$
### 3——模态代入与特征值求解
将 $\delta\phi = a e^{in\theta + s\tau}$ 代入，注意 $\partial/\partial\tau \to s$，$\partial/\partial\theta \to in$：
$$ \frac{1}{|n|} s - \left( -\frac{1}{|n|} s \right) = \frac{d\psi^{ts}}{d\phi} - \lambda(in) - \mu s $$
化简：
$$ \left( \mu + \frac{2}{|n|} \right) s + i\lambda n = \frac{d\psi^{ts}}{d\phi} $$
分离实部（$\sigma$）和虚部（$\omega$）：
$$ \boxed{\sigma_n = \frac{\frac{d\psi^{ts}}{d\phi}}{\mu + \frac{2}{|n|}}} \qquad \boxed{\omega_n = \frac{\lambda n}{\mu + \frac{2}{|n|}}} $$
### 4——物理解读
| 结果 | 含义 |
| :--- | :--- |
| $\sigma_n < 0$ | 第 $n$ 阶谐波衰减，流动稳定 |
| $\sigma_n > 0$ | 第 $n$ 阶谐波增长，失速先兆 |
| $\sigma_n = 0$ | 临界稳定，$d\psi^{ts}/d\phi = 0$ 为失速边界 |
| $\omega_n > 0$ | 扰动沿转子旋转方向传播 |
| $\omega_n \propto n$ | 高阶谐波旋转更快 |
### 5——考虑损失迟滞（更一般的形式）
若加入一阶损失迟滞效应，特征值方程会稍有不同的形式，其核心效果是使不同谐波的稳定性边界发生分离，高阶模态往往更稳定。

[^1]: ## 势流分解：拉普拉斯方程与谐波解
	柱坐标，二维无旋流动的速度可以表示为一个势函数 $\Phi$ 的梯度：
	$$ v_x = \frac{\partial \Phi}{\partial x}, \quad v_\theta = \frac{1}{R} \frac{\partial \Phi}{\partial \theta} $$
	代入不可压连续性方程 $\nabla \cdot \mathbf{v} = 0$，得到拉普拉斯方程：
	$$ \frac{\partial^2 \Phi}{\partial x^2} + \frac{1}{R^2} \frac{\partial^2 \Phi}{\partial \theta^2} = 0 \quad  $$
	考虑小扰动，将流动分解为稳态背景流加上小扰动：
	$$ \Phi(x, \theta, \tau) = \overline{\Phi}(x) + \delta \Phi(x, \theta, \tau) $$
	背景流是轴向均匀的（$\overline{v}_x = \text{const}$），$\overline{\Phi} = \overline{v}_x x$。扰动势函数 $\delta \Phi$ 同样满足拉普拉斯方程。
	
	由于拉普拉斯方程是线性的，且周向是周期性的，可将扰动分解为周向傅里叶谐波。对于第 $n$ 阶谐波，设解的形式为：
	$$ \delta \Phi_n(x, \theta, \tau) = \tilde{\Phi}_n(x) \cdot e^{in\theta + s\tau} $$
	其中 $s = \sigma + i\omega$ 是复特征值。
	
	将此形式代入拉普拉斯方程：
	$$ \frac{d^2 \tilde{\Phi}_n}{dx^2} + \frac{(in)^2}{R^2} \tilde{\Phi}_n = 0 \;\Longrightarrow\; \frac{d^2 \tilde{\Phi}_n}{dx^2} - \frac{n^2}{R^2} \tilde{\Phi}_n = 0 $$
	通解为：
	$$ \tilde{\Phi}_n(x) = A_n e^{n x / R} + B_n e^{-n x / R} $$
	## 边界条件：衰减解的确定
	物理约束——扰动在远上游必须消失（$x \to -\infty$ 时，$\delta \Phi \to 0$）：
	- 当 $x \to -\infty$ 时，$e^{n x / R} \to 0$（衰减）
	- 当 $x \to -\infty$ 时，$e^{-n x / R} \to \infty$（发散）
	因此，**必须令 $B_n = 0$**，只保留指数衰减的 $A_n$ 项。上游势函数解为：
	$$ \boxed{\delta \Phi_n(x, \theta, \tau) = A_n \cdot e^{n x / R} \cdot e^{in\theta + s\tau}} $$
	这就是**衰减势流解**——扰动在向远上游（$x \to -\infty$）传播时，振幅按指数 $e^{n x/R}$ 衰减。这也解释了为什么高阶谐波（$n$ 大）在上游衰减极快。
	## 总压匹配关系
	对有了上面的解，对于无旋流动，还需要线性化的非定常伯努利方程（无视重力）——全流场适用。对于扰动形式，非定常的伯努利方程等价于：
	$$ \frac{\partial}{\partial t} (\delta \Phi) + \frac{\delta p_t}{\rho U^2} = 0 $$
	或：
	$$ \frac{\delta p_t}{\rho U^2} = -\frac{\partial}{\partial t} (\delta \Phi) $$
	现在将已求得的势函数解代入（ $\delta \Phi_n = A_n e^{n x / R} e^{in\theta + s\tau}$）。从流量系数的扰动形式出发，发现流量系数和势函数之间的关系：
	$$ \delta \phi = \frac{\delta v_x}{U} = \frac{1}{U} \frac{\partial}{\partial x}(\delta \Phi) = \frac{n}{R U} \delta \Phi \;\Longrightarrow\; \delta \Phi = \frac{UR}{n} \delta \phi $$
	代入上面的伯努利方程，并注意抽象出的长度尺度经过无因次化后（$x$ 坐标已由 $R$ 无因次化），得到最终总压关系(常数被分母的无因次时间吸收）：
	$$ \boxed{\frac{\delta p_{t1}}{\rho U^2} = -\frac{1}{|n|} \frac{\partial}{\partial \tau} (\delta \phi)} $$
	其中 $|n|$ 确保该关系对任意符号的谐波均成立（衰减速率仅取决于 $|n|$，与旋转方向无关）。
