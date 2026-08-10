
|  属性  |            [[法术表]]             |
| :--: | :----------------------------: |
|  关联  | [[Spakovszky稳定性方法——模型Ⅱ——轴流转子]] |
| 链接解释 |              一脉相承              |
***
# $Model\space Ⅱ$
    适用范围：轴向压气机（或风扇）的静子叶片排，作为半执行盘（semi‑actuator disk），连接上游和下游轴向管道。模型假设叶片稠度足够，可将叶片排视为无厚度的周向均匀装置，仅改变流动的周向速度（去旋）和总压。适用于不可压、小扰动、频域稳定性分析。

| 变量名 | 符号 | 单位 | 量纲 | 说明 |
| :-: | :-: | :-: | :-: | :-: |
| 轴向速度 | $v_x$ | $\text{m/s}$ | $L T^{-1}$ | 绝对坐标系下轴向速度分量 |
| 周向速度 | $v_\theta$ | $\text{m/s}$ | $L T^{-1}$ | 绝对坐标系下周向速度分量 |
| 静压 | $p$ | $\text{Pa}$ | $M L^{-1} T^{-2}$ | 当地静压 |
| 总压 | $p_t$ | $\text{Pa}$ | $M L^{-1} T^{-2}$ | $p_t = p + \frac12 \rho (v_x^2+v_\theta^2)$ |
| 密度 | $\rho$ | $\text{kg/m}^3$ | $M L^{-3}$ | 常数（不可压） |
| 静子总压损失（有量纲） | $l_S$ | $\text{Pa}$ | $M L^{-1} T^{-2}$ | 通过静子的总压降 |
| 无量纲总压损失 | $\hat l_S$ | – | – | $\hat l_S = l_S / (\rho U^2)$ |
| 稳态损失特性 | $l_S^{ss}(\cdot)$ | $\text{Pa}$ | $M L^{-1} T^{-2}$ | 稳态下损失与进口参数的关系 |
| 静子惯性系数 | $\lambda_{\text{sta}}$ | – | – | $\lambda_{\text{sta}} = c_x / (R \cos^2\gamma)$ |
| 损失时间滞后（有量纲） | $\tau_S$ | $\text{s}$ | $T$ | 损失响应的延迟时间 |
| 无量纲损失滞后 | $\tau_S$ | – | – | $\tau_S = \tau_u \cdot \dfrac{c_x}{R\cos\gamma\,\overline v_x}$ |
| 经验滞后系数 | $\tau_u$ | – | – | 通常取 $1.0$～$1.5$ |
| 轴向弦长 | $c_x$ | $\text{m}$ | $L$ | 叶片轴向投影长度 |
| 平均半径 | $R$ | $\text{m}$ | $L$ | 特征长度（转子中径） |
| 叶片安装角 | $\gamma$ | $\text{rad}$ | – | 叶片弦线与周向（或轴向）的夹角，几何常数 |
| 特征速度 | $U$ | $\text{m/s}$ | $L T^{-1}$ | 转子中径切线速度，无量纲化基准 |
| 无量纲轴向速度 | $\hat v_x = v_x/U$ | – | – | |
| 无量纲周向速度 | $\hat v_\theta = v_\theta/U$ | – | – | |
| 无量纲静压 | $\hat p = p/(\rho U^2)$ | – | – | |
| 无量纲时间 | $\tau = tU/R$ | – | – | |
| 拉普拉斯变量（无量纲） | $s$ | – | – | $\partial/\partial\tau \to s$ |
| 周向波数 | $n$ | – | – | 整数，$n=0$ 轴对称 |
| 绝对气流角 | $\alpha_1, \alpha_2$ | $\text{rad}$ | – | $\tan\alpha = v_\theta/v_x$ |
| 相对气流角（参考） | $\beta_1$ | $\text{rad}$ | – | $\tan\beta_1 = (U - v_{\theta1})/v_{x1}$ |
| 背景流（上划线） | $\overline{(\cdot)}$ | – | – | 定常、均匀平均量 |
| 小扰动（时域） | $\delta(\cdot)$ | – | – | |
| 频域小扰动 | $\delta\hat V_x, \delta\hat V_\theta, \delta\hat P$ | – | – | |
| 第 $n$ 谐波系数 | $\delta\tilde V_{x,n}, \delta\tilde V_{\theta,n}, \delta\tilde P_n$ | – | – | $\delta\hat V_x = \sum_n \delta\tilde V_{x,n} e^{jn\theta}$ |
***
# 主体
这个很简单，让我们来看图说话！
![[Spakovszky稳定性方法——模型Ⅱ——轴向静叶配图.png]]
# 1. 匹配条件
咱还没无量纲，来看速度匹配和压力关系：
$$
v_{x1}=v_{x2}\tag{1}
$$
$$
v_{\theta 2}=v_{x 2}\cdot\tan{\alpha_2}\tag{2}
$$
静叶没有扭速关系，也没有旋转带来的非定常性（周向导数项）：
$$
p_{t2}-p_{t1}=-l_S-\lambda_{sta}\frac{\partial v_{x1}}{\partial t}\tag{3}
$$
不过，还是要考虑迟滞的：
$$
\left(1+\tau_S\cdot\frac{\partial}{\partial t}\right)l_S=l_S^{ss}(v_{x1},v_{\theta 1})\tag{4}
$$
$$
\tau_S=\tau_u\cdot\frac{c_x}{R\cos\gamma\overline v_x}
$$
# 2. 无量纲化、线性化......
无量纲化参照量和其他文档相同即可，如果没有，对于轴流设备，就取转子中径的速度参数即可。
$$
\hat v_{x1}=\hat v_{x2}\tag{5}
$$
$$
\hat v_{\theta 2}=\hat v_{x 2}\cdot\tan{\alpha_2}\tag{6}
$$
$$
\hat p_{t2}-\hat p_{t1}=-\hat l_S-\lambda_{sta}\frac{\partial \hat v_{x1}}{\partial \tau}\tag{7}
$$
展开总压，由于匹配关系，轴向速度会被约掉：
$$
\hat p_{2}+\frac{1}{2}v_{\theta 2}^2-\hat p_{1}-\frac{1}{2}v_{\theta 1}^2=-\hat l_S-\lambda_{sta}\frac{\partial \hat v_{x1}}{\partial \tau}\tag{8}
$$
线性化：
$$
\delta\hat p_{2}+\hat{\overline v}_{\theta2}\delta \hat v_{\theta 2}-\delta\hat p_{1}-\hat{\overline v}_{\theta1}\delta \hat v_{\theta1}=-\delta\hat l_S-\lambda_{sta}\frac{\partial \delta\hat v_{x1}}{\partial \tau}\tag{8}
$$
拉普拉斯变换，周向傅里叶展开，对于$n$阶谐波有：
$$
\delta\tilde P_{2}+\hat{\overline v}_{\theta2}\delta\tilde V_{\theta 2}-\delta\tilde P_{1}-\hat{\overline v}_{\theta1}\delta\tilde V_{\theta1}=-\delta\tilde L_S-\lambda_{sta}s \delta\tilde V_{x1}\tag{9}
$$
这个损失的扰动，继续抄转子作业，由于迟滞公式里没有对$\hat\theta$的导数了，把分母里的$jn$给删掉。除此之外，咱没有相对速度，只有绝对速度角了：
$$
\boxed{
\delta\tilde{L}_S =  \frac{d\hat{l}_S^{ss}}{d\tan\alpha_1}\bigg|_{\text{avg}}\cdot \frac{    \delta\tilde{V}_{\theta1}-\tan\alpha_1 \delta\tilde{V}_{x1}} {\hat{\overline v}_{x1}[1 + s\cdot\tau_S]}}
\tag{10}
$$
其中$\dfrac{d\hat{l}_R^{ss}}{d\tan\alpha_1}\bigg|_{\text{avg}}$ 在工作点取值，同样有$M-G$模型可以参考，这两者在速度三角形里互余，所以抄过来还得再变个号：
$$
\frac{d\hat{l}_S^{ss}}{d\tan\alpha_1}\bigg|_{\text{avg}}=\frac{d\hat{l}_S^{ss}}{d\hat v_{x1}}\bigg|_{\text{avg}}\cdot\frac{-1}{(\tan{\alpha_1}-\tan{\beta_1})^2}
$$
# 3. 传递矩阵
有
$$
\mathbb{B}_{\text{sta},n}(s)=
\begin{bmatrix}
1 &&& 
0 &&& 
0 &\\[4pt]
\tan\alpha_2 &&& 
0 &&& 
0 &\\[12pt]
\displaystyle
-\hat{\overline v}_{\theta2}\tan\alpha_2 + \left.\frac{d\hat l_S^{ss}}{d\tan\alpha_1}\right|_{\text{avg}}\frac{\tan\alpha_1}{\hat{\overline v}_{x1}(1+s\tau_S)} \;-\; s\lambda_{\text{sta}}
&&&
\displaystyle
\hat{\overline v}_{\theta1} - \left.\frac{d\hat l_S^{ss}}{d\tan\alpha_1}\right|_{\text{avg}}\frac{1}{\hat{\overline v}_{x1}(1+s\tau_S)}
&&& 1
\end{bmatrix}
$$
$$
\begin{bmatrix}
\delta\hat V_{x2,n} \\[2pt]
\delta\hat V_{\theta 2,n} \\[2pt]
\delta\hat P_{2,n}
\end{bmatrix}(\hat\theta,s)
=
\mathbb{B}_{\text{sta},n}(s)\cdot e^{jn\hat\theta}
\cdot
\begin{bmatrix}
\delta\tilde V_{x1,n} \\[2pt]
\delta\tilde V_{\theta 1,n} \\[2pt]
\delta\tilde P_{1,n}
\end{bmatrix}(s)
,\qquad
\forall n \ge 0
$$
