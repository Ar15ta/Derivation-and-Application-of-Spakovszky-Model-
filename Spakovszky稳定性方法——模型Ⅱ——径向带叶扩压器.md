
|  属性  |            [[法术表]]             |                                |
| :--: | :----------------------------: | :----------------------------: |
|  关联  | [[Spakovszky稳定性方法——模型Ⅱ——离心转子]] | [[Spakovszky稳定性方法——模型Ⅱ——轴流转子]] |
| 链接解释 |           用了其中的匹配关系            |              用了解               |
***
# $Model\space Ⅱ$
	`适用范围：离心压缩机中的径向有叶扩压器（vaned diffuser），作为半执行盘（semi‑actuator disk），连接上游径向空间（如叶轮出口或径向间隙）和下游轴向管道/集气室。模型假设叶片稠度足够，可将叶片排视为无厚度的周向均匀装置，改变流动的径向速度（通过面积‑密度比）和周向速度（去旋），并计入总压损失。适用于不可压、小扰动、频域稳定性分析。`

|     变量名      |                            符号                             |       单位        |        量纲         | 说明                                                           |
| :----------: | :-------------------------------------------------------: | :-------------: | :---------------: | :----------------------------------------------------------- |
|      密度      |                          $\rho$                           | $\text{kg/m}^3$ |    $M L^{-3}$     | 常数（不可压）                                                      |
|      静压      |                            $p$                            |   $\text{Pa}$   | $M L^{-1} T^{-2}$ | 当地静压                                                         |
|      总压      |                           $p_t$                           |   $\text{Pa}$   | $M L^{-1} T^{-2}$ | $p_t = p + \frac{1}{2}\rho (v_r^2+v_\theta^2)$               |
|     径向速度     |                           $v_r$                           |  $\text{m/s}$   |    $L T^{-1}$     | 绝对坐标系下径向速度分量                                                 |
|     周向速度     |                        $v_\theta$                         |  $\text{m/s}$   |    $L T^{-1}$     | 绝对坐标系下周向速度分量                                                 |
| 扩压器总压损失（有量纲） |                     $l_{\text{dif}}$                      |   $\text{Pa}$   | $M L^{-1} T^{-2}$ | 通过扩压器的总压降                                                    |
|    稳态损失特性    |               $l_{\text{dif}}^{ss}(\cdot)$                |   $\text{Pa}$   | $M L^{-1} T^{-2}$ | 稳态下损失与进口参数的关系                                                |
|  惯性系数（无量纲）   |                  $\lambda_{\text{dif}}$                   |        –        |         –         | 由式(2)定义，积分使用通流面积 $A$                                         |
| 损失时间滞后（有量纲）  |                    $\tau_{\text{dif}}$                    |   $\text{s}$    |        $T$        | 损失响应的延迟时间                                                    |
|   无量纲损失滞后    |                    $\tau_{\text{dif}}$                    |        –        |         –         | $\tau_{\text{dif}} = \tau_u \cdot 2s_{\text{dif}}/(V_1+V_2)$ |
|    经验滞后系数    |                         $\tau_u$                          |        –        |         –         | 通常取 $1.0$                                                    |
|    径向投影面积    |                           $A'$                            |  $\text{m}^2$   |       $L^2$       | 不乘 $\cos\chi$，用于连续性匹配                                        |
|     通流面积     |                            $A$                            |  $\text{m}^2$   |       $L^2$       | $A = A' \cos\chi$，用于惯性系数积分                                   |
|  径向投影面积‑密度比  |                    $AR'_{\text{dif}}$                     |        –        |         –         | $AR'_{\text{dif}} = \frac{\rho_2 A_2'}{\rho_1 A_1'}$         |
|   通流面积‑密度比   |                     $AR_{\text{dif}}$                     |        –        |         –         | $AR_{\text{dif}} = \frac{\rho_2 A_2}{\rho_1 A_1}$            |
|    叶片金属角     |                     $\chi_1, \chi_2$                      |  $\text{rad}$   |         –         | 叶片切线与径向的夹角（在 $r$‑$\theta$ 平面内）                               |
|    绝对气流角     |                   $\alpha_1, \alpha_2$                    |  $\text{rad}$   |         –         | $\tan\alpha = v_\theta / v_r$                                |
|     叶片高度     |                            $h$                            |   $\text{m}$    |        $L$        | 径向扩压器叶片高度（轴向尺寸）                                              |
|      半径      |                        $R_1, R_2$                         |   $\text{m}$    |        $L$        | 扩压器进口、出口半径                                                   |
|     叶片数      |                     $N_{\text{dif}}$                      |        –        |         –         | 扩压器叶片总数                                                      |
|     叶片厚度     |                     $t_{\text{dif}}$                      |   $\text{m}$    |        $L$        | 叶片尾缘厚度（周向扣除）                                                 |
|     流道长度     |                     $s_{\text{dif}}$                      |   $\text{m}$    |        $L$        | 沿叶片通道的平均长度                                                   |
|     特征速度     |                            $U$                            |  $\text{m/s}$   |    $L T^{-1}$     | 叶轮出口切线速度（无量纲化基准）                                             |
|     特征长度     |                            $R$                            |   $\text{m}$    |        $L$        | 平均半径（通常取 $R_2$）                                              |
|   无量纲径向速度    |                    $\hat v_r = v_r/U$                     |        –        |         –         |                                                              |
|   无量纲周向速度    |               $\hat v_\theta = v_\theta/U$                |        –        |         –         |                                                              |
|    无量纲静压     |                  $\hat p = p/(\rho U^2)$                  |        –        |         –         |                                                              |
|    无量纲总压     |                $\hat p_t = p_t/(\rho U^2)$                |        –        |         –         |                                                              |
|    无量纲时间     |                       $\tau = tU/R$                       |        –        |         –         |                                                              |
| 拉普拉斯变量（无量纲）  |                            $s$                            |        –        |         –         | $\partial/\partial\tau \to s$                                |
|     周向波数     |                            $n$                            |        –        |         –         | 整数，$n=0$ 轴对称                                                 |
|   平均量（上划线）   |                   $\overline{(\cdot)}$                    |        –        |         –         | 定常、均匀背景流                                                     |
|   小扰动（时域）    |                      $\delta(\cdot)$                      |        –        |         –         |                                                              |
|    频域小扰动     | $\delta\tilde V_r, \delta\tilde V_\theta, \delta\tilde P$ |        –        |         –         |                                                              |
***
# 主体
![[Spakovszky稳定性方法——模型Ⅱ——径向静叶配图.png]]
# 1. 匹配条件
定义：$AR_\text{dif}=\frac{\rho_2 A_2}{\rho_1 A_1}$和$AR'_\text{dif}=\frac{\rho_2 A_2'}{\rho_1 A_1'}$。同样，由于非定常性而需要的额外总压降可以表示为：
$$
\triangle p_{t,unsteady}=-\int_{s1}^{s2}\frac{\rho_2 A_2}{\rho(s) A(s)}ds\cdot\frac{\partial v_{r2}}{\partial t}=-\lambda_\text{dif}\cdot\frac{\partial v_{r2}}{\partial t}\tag{1}
$$
计算惯性系数时，不用径向面积，而是扩压器中的通流面积，公式抄过来，但多了金属角(写到这里，其实可以看到，这惯性系数把出口的圆环形面积和当地流道的面积做了个近似，希望影响不大）：
$$
A_1 =h_1\cdot\frac{2\pi R_1}{N_\text{dif}}\cos\chi_1\qquad{\space}\qquad{\space}
A_2 =h_2\cdot\left(\frac{2\pi R_2}{N_\text{dif}}-t_\text{dif}\right)\cos{\chi_2}
$$

在密度面积乘积线性变化的假设下，直接把离心叶轮的积分解拿过来用：
$$
\lambda_{\text{dif}} = s_{\text{dif}} \cdot \frac{AR_{\text{dif}} \ln(AR_{\text{dif}})}{AR_{\text{dif}} - 1}\tag{2}
$$
由径向质量守恒：
$$
   v_{r1} = AR'_{\text{dif}} \cdot v_{r2}, \quad AR'_{\text{dif}} = \frac{\rho_2 A_2'}{\rho_1 A_1'} \tag{3}
$$

$$
A_1' =h_1\cdot\frac{2\pi R_1}{N_\text{dif}}\qquad{\space}\qquad{\space}
A_2' =h_2\cdot\left(\frac{2\pi R_2}{N_\text{dif}}-t_\text{dif}\right)
$$
周向速度匹配：
$$
v_{\theta 2} = v_{r2} \cdot \tan\alpha_2 \tag{4}
$$
总压变化（稳态损失 + 非定常惯性）：
$$
p_{t2} - p_{t1} = -l_{\text{dif}} - \lambda_{\text{dif}} \frac{\partial v_{r2}}{\partial t} \tag{5}
$$
其中 $l_{\text{dif}}$ 为稳态总压损失，$\lambda_{\text{dif}}$ 由式 (2) 定义，其积分使用通流面积 $A$（含 $\cos\chi$）。

损失滞后（一阶滞后）：   
$$
   \left(1 + \tau_{\text{dif}} \frac{\partial}{\partial t}\right) l_{\text{dif}} = l_{\text{dif}}^{ss}(v_{r1}, v_{\theta 1}) \tag{6}
$$
其中$\tau_{\text{dif}}$ 为损失响应时间，近似为扩压器内平均流通时间（也就是说可以用基流参数来算，是个常数）：
$$
   \tau_{\text{dif}} = \tau_u \cdot \frac{2s_{\text{dif}}}{\boldsymbol{v}_1 + \boldsymbol{v}_2}, \quad \boldsymbol{v} = \sqrt{v_r^2 + v_\theta^2} \tag{7}
$$
# 2. 无量纲化与线性化
## 2.1 无量纲基准  
- 特征速度 $U$（叶轮出口切线速度）  
- 特征长度 $R$（平均半径，通常取叶轮出口半径 $R_2$）  
- 密度 $\rho$（参考密度，通常取进口密度）  

定义无量纲量：
$$
\hat v = \frac{v}{U},\quad \hat p = \frac{p}{\rho U^2},\quad \hat p_t = \frac{p_t}{\rho U^2},\quad \hat L = \frac{L}{\rho U^2},\quad \tau = \frac{tU}{R},\quad \hat s_\text{dif} = \frac{s_\text{dif}}{R}.
$$
## 2.2 变换
对无量纲匹配关系中，时间 $\tau$ 做拉普拉斯变换 $\partial/\partial\tau \to s$，并对周向 $\theta$ 做傅里叶展开 $e^{jn\theta}$。记频域小扰动量为 $\delta\tilde V_{r1,n}, \delta\tilde V_{\theta1,n}, \delta\tilde P_{1,n}$ 等，对第$n$阶谐波有。
$$
\delta\tilde P_{t2} - \delta\tilde P_{t1} = -\delta\tilde L_{\text{dif}} - \lambda_{\text{dif}} \frac{\partial \delta\tilde V_{r2}}{\partial \tau} \tag{8}
$$
将总压表示为静压与动压之和（不可压）：
$$
\delta\tilde P_{2} + \hat{\overline v}_{r2}\delta\tilde V_{r2} + \hat{\overline v}_{\theta2}\delta\tilde V_{\theta2} 
- \left( \delta\tilde P_{1} + \hat{\overline v}_{r1}\delta\tilde V_{r1} + \hat{\overline v}_{\theta1}\delta\tilde V_{\theta1} \right)
= -\delta\tilde L_{\text{dif}} - \lambda_{\text{dif}} \cdot s \, \delta\tilde V_{r2}
$$
代入速度匹配关系：
$$
\delta\tilde P_{2} + \hat{\overline v}_{r2}\frac{\delta\tilde V_{r1}}{AR'_{\text{dif}}} + \hat{\overline v}_{\theta2}\frac{\delta\tilde V_{r1}}{AR'_{\text{dif}}} \cdot \tan\alpha_2  
- \left( \delta\tilde P_{1} + \hat{\overline v}_{r1}\delta\tilde V_{r1} + \hat{\overline v}_{\theta1}\delta\tilde V_{\theta1} \right)
= -\delta\tilde L_{\text{dif}} - \lambda_{\text{dif}} \cdot s \, \frac{\delta\tilde V_{r1}}{AR'_{\text{dif}}}
$$
以及从轴向静叶拿来的非定常损失解：
$$
\boxed{
\delta\tilde{L}_\text{dif} =  \frac{d\hat{l}_\text{dif}^{ss}}{d\tan\alpha_1}\bigg|_{\text{avg}}\cdot \frac{    \delta\tilde{V}_{\theta1}-\tan\alpha_1 \delta\tilde{V}_{r1}} {\hat{\overline v}_{x1}[1 + s\cdot\tau_S]}}
\tag{9}
$$
最后有：
$$
\delta\tilde P_{2} = \delta\tilde P_{1} 
+ \left[ \hat{\overline v}_{\theta1} - \frac{1}{\hat{\overline v}_{r1}(1+s\tau_{\text{dif}})} \frac{\partial \hat l_{\text{dif}}^{ss}}{\partial \tan\alpha_1} \right] \delta\tilde V_{\theta1}
+ \left[ \hat{\overline v}_{r1} - \frac{\hat{\overline v}_{r2} + \hat{\overline v}_{\theta2}\tan\alpha_2}{AR'_{\text{dif}}} - \frac{\lambda_{\text{dif}} \cdot s}{AR'_{\text{dif}}} + \frac{\tan\alpha_1}{\hat{\overline v}_{r1}(1+s\tau_{\text{dif}})} \frac{\partial \hat l_{\text{dif}}^{ss}}{\partial \tan\alpha_1} \right] \delta\tilde V_{r1}\tag{10}
$$
# 3. 传递矩阵
有：
$$
\begin{bmatrix}
\delta\hat V_{r2,n} \\[2pt]
\delta\hat V_{\theta 2,n} \\[2pt]
\delta\hat P_{2,n}
\end{bmatrix}(\hat\theta,s)
=
\mathbb{B}_{\text{dif},n}(s)\cdot e^{jn\hat\theta}
\cdot
\begin{bmatrix}
\delta\tilde V_{r1,n} \\[2pt]
\delta\tilde V_{\theta 1,n} \\[2pt]
\delta\tilde P_{1,n}
\end{bmatrix}(s)
$$
其中：
$$
\mathbb{B}_{\text{dif},n}(s) =
\begin{bmatrix}
\displaystyle \frac{1}{AR'_{\text{dif}}} & 0 & 0 &\\[10pt]
\displaystyle \frac{\tan\alpha_2}{AR'_{\text{dif}}} & 0 & 0 &\\[12pt]
\displaystyle
\hat{\overline v}_{r1} - \frac{\hat{\overline v}_{r2} + \hat{\overline v}_{\theta2}\tan\alpha_2}{AR'_{\text{dif}}} - \frac{\lambda_{\text{dif}} \cdot s}{AR'_{\text{dif}}} + \frac{\tan\alpha_1}{\hat{\overline v}_{r1}(1+s\tau_{\text{dif}})} \frac{\partial \hat l_{\text{dif}}^{ss}}{\partial \tan\alpha_1}
&
\displaystyle
\hat{\overline v}_{\theta1} - \frac{1}{\hat{\overline v}_{r1}(1+s\tau_{\text{dif}})} \frac{\partial \hat l_{\text{dif}}^{ss}}{\partial \tan\alpha_1}
& 1
\end{bmatrix}
$$
