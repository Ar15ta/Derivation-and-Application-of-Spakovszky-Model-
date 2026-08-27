
|  属性  |            [[法术表]]             |            |                          |
| :--: | :----------------------------: | ---------- | ------------------------ |
|  关联  | [[Spakovszky稳定性方法——模型Ⅱ——轴向管道]] | [[欧拉涡轮方程]] | [[Spakovszky稳定性方法——模型Ⅰ]] |
| 链接解释 |    转子执行盘的上游连接轴向管道，下游连接间隙或静子    | 直接用        | 惯性导致的额外总压，以及损失在这里有提及     |

***
# $Model\space Ⅱ$
	适用范围：轴向压气机（或风扇）的转子叶片排，作为半执行盘（semi‑actuator disk），连接上游和下游轴向管道。假设叶片稠度足够，可将叶片排视为无厚度的周向均匀装置，仅改变流动的角动量和总压。

|          变量名          |                               符号                               |       单位        |       量纲        |                                      说明                                       |
| :-------------------: | :------------------------------------------------------------: | :-------------: | :-------------: | :---------------------------------------------------------------------------: |
|         轴向速度          |                             $v_x$                              |  $\text{m/s}$   |    $lT^{-1}$    |                                    绝对坐标系下                                     |
|         周向速度          |                           $v_\theta$                           |  $\text{m/s}$   |    $lT^{-1}$    |                                    绝对坐标系下                                     |
|          静压           |                              $p$                               |   $\text{Pa}$   | $Ml^{-1}T^{-2}$ |                                                                               |
|          密度           |                             $\rho$                             | $\text{kg/m}^3$ |    $Ml^{-3}$    |                                    常数（不可压）                                    |
|        转子轮缘速度         |                              $U$                               |  $\text{m/s}$   |    $lT^{-1}$    |                                  中径处切线速度，常数                                   |
|         平均半径          |                              $R$                               |   $\text{m}$    |       $l$       |                                     特征长度                                      |
|        无量纲轴向坐标        |                        $\hat{x} = x/R$                         |       $-$       |       $1$       |                                                                               |
|         无量纲时间         |                         $\tau = tU/R$                          |       $-$       |       $1$       |                                                                               |
|        拉普拉斯变量         |                              $s$                               |       $-$       |       $1$       |                       对应 $\partial/\partial\tau \to s$                        |
|        无量纲轴向速度        |                      $\hat{v}_x = v_x/U$                       |       $-$       |       $1$       |                                                                               |
|        无量纲周向速度        |                 $\hat{v}_\theta = v_\theta/U$                  |       $-$       |       $1$       |                                                                               |
|         无量纲静压         |                    $\hat{p} = p/(\rho U^2)$                    |       $-$       |       $1$       |                                                                               |
|       稳态背景（无量纲）       |       $\hat{\overline{v}}_x,\ \hat{\overline{v}}_\theta$       |       $-$       |       $1$       |                                      常数                                       |
|      小扰动（时域无量纲）       |    $\delta\hat{v}_x,\ \delta\hat{v}_\theta,\ \delta\hat{p}$    |       $-$       |       $1$       |                                                                               |
|        频域无量纲扰动        |    $\delta\hat{V}_x,\ \delta\hat{V}_\theta,\ \delta\hat{P}$    |       $-$       |       $1$       |                                    拉普拉斯变换后                                    |
|   第 $n$ 谐波傅里叶系数（频域）   | $\delta\tilde{V}_x,\ \delta\tilde{V}_\theta,\ \delta\tilde{P}$ |       $-$       |       $1$       |            $\delta\hat{V}_x = \sum \delta\tilde{V}_x e^{jn\theta}$            |
|     转子进口相对气流角（稳态）     |                           $\beta_1$                            |  $\text{rad}$   |       $1$       | $\tan\beta_1 = \frac{1 - \hat{\overline{v}}_{\theta1}}{\hat{\overline{v}}_x}$ |
|       转子出口相对气流角       |                           $\beta_2$                            |  $\text{rad}$   |       $1$       |                                     给定常数                                      |
|        转子惯性系数         |                     $\lambda_{\text{rot}}$                     |       $-$       |       $1$       |              $\lambda_{\text{rot}} = \frac{c_x/R}{\cos^2\gamma}$              |
|     损失时间滞后常数（无量纲）     |                            $\tau_R$                            |       $-$       |       $1$       |   通常 $\tau_R = \tau_u \cdot \frac{c_x}{R\cos\gamma\,\hat{\overline{v}}_x}$    |
|       总压损失（无量纲）       |                  $\hat{l}_R = l_R/(\rho U^2)$                  |       $-$       |       $1$       |                                                                               |
|        稳态损失特性         |                 $\hat{l}_R^{ss}(\tan\beta_1)$                  |       $-$       |       $1$       |                                                                               |
| 损失对 $\tan\beta_1$ 的导数 |          $\hat{l}_R' = d\hat{l}_R^{ss}/d\tan\beta_1$           |       $-$       |       $1$       |                                    在工作点取值                                     |

***
# 主体
![[Spakovszky稳定性方法——模型Ⅱ——轴流转子配图.png]]

# 1. 基本假设与几何
- 转子叶片排简化为**半执行盘**：轴向和周向尺度远小于波长，因此可视为无厚度，位于 $\hat{x}=0$。
- 流动无粘、不可压、等熵（除了总压损失）。
- 叶片出口气流角 $\beta_2$ 与进口条件无关（即**无攻角扰动**，忽略非定常偏差）。
- 叶片通道内的流体惯性用集总参数 $\lambda_{\text{rot}}$ 模拟，它表示转子叶片排对轴向速度扰动的非定常压力响应。
- 总压损失 $\hat{l}_R$ 是进口相对气流角 $\beta_1$ 的函数，且对非定常扰动有一阶滞后（时间常数 $\tau_R$）。
# 2. 时域匹配条件（无量纲）
在转子盘前后（STA 1 进口，STA 2 出口），连续性和出口流动角给出：
$$
\begin{aligned}
\hat{v}_{x2} &= \hat{v}_{x1} \\
\hat{v}_{\theta 2} &= 1 + \hat{v}_{x2}\cdot\tan\beta_2   \\
\end{aligned}
$$
角度定义中，分子分母都已除以叶轮中径线速度$U$：
$$
\begin{aligned}
\tan\beta_2&=\frac{\hat{v}_{\theta 2} -1}{\hat{v}_{x2}}  \\
\tan\alpha_1 &= \frac{\hat{v}_{\theta 1}}{\hat{v}_{x1}}
\end{aligned}
$$
考虑损失和迟滞的总压升方程（Euler 涡轮方程——两边除以了无量纲参数$\rho U^2$，考虑了无量纲形式的损失和惯性项，这两个都在$Model\space Ⅰ$里有）：
$$
\hat{p}_{t2} - \hat{p}_{t1} = \hat{v}_{\theta 2} - \hat{v}_{\theta 1} - \hat{l}_R - \lambda_{\text{rot}}\left(\frac{\partial \hat{v}_{x1}}{\partial \tau}+\frac{\partial \hat{v}_{x1}}{\partial \hat\theta}\right) \tag{1}
$$
也就是：
$$
\hat{p}_{t2} - \hat{p}_{t1} = 1 + (\tan\beta_2 - \tan\alpha_1)\hat{v}_{x1} - \hat{l}_R - \lambda_{\text{rot}}\left(\frac{\partial \hat{v}_{x1}}{\partial \tau}+\frac{\partial \hat{v}_{x1}}{\partial \hat\theta}\right) \tag{2}
$$
这俩工况依赖的角度是“已知的变量”。

损失的一阶滞后模型：

$$
\left[1+\tau_\text{R}\cdot\left(\frac{\partial \hat{v}_{x1}}{\partial \tau}+\frac{\partial \hat{v}_{x1}}{\partial \hat\theta}\right)\right]\hat{l}_R = \hat{l}_R^{ss}\bigl(\hat v_{x1},\hat v_{\theta})=\hat{l}_R^{ss}\bigl(\tan\beta_1) \tag{3}
$$
然后，这里迟滞常数的定义为：
$$
\tau_R=\tau_u\cdot\frac{c_x}{R \cos{\gamma}\cdot\overline{v}_x}
$$
其中，$\tau_u$的范围为$1.0$~$1.5$这样子，这个数的来源是其他文献的实验。
# 3. 转换
## 3.1 总压项与惯性项
对式(1)线性化：
$$
\hat{\overline{p}}_{2}+\delta\hat{p}_{2}+\frac{1}{2}(\hat{\overline{v}}_{\theta 2}+\delta{\hat{v}}_{\theta 2}+\hat{\overline{v}}_{x 2}+\delta{\hat{v}}_{x2})^2 -\hat{\overline{p}}_{1}- \delta\hat{p}_{1}-\frac{1}{2}(\hat{\overline{v}}_{\theta 1}+\delta{\hat{v}}_{\theta 1}+\hat{\overline{v}}_{x 1}+\delta{\hat{v}}_{x 1})^2  
= 
 1 + (\tan\beta_2 - \tan\alpha_1)(\delta\hat{v}_{x1}+\hat{\overline v}_{x1})- \delta\hat{{l}}_R - \lambda_{\text{rot}}\left[\frac{\partial (\delta\hat{v}_{x1}+\hat{\overline{v}}_{x1})}{\partial \tau}+\frac{\partial (\delta\hat{v}_{x1}+\hat{\overline{v}}_{x1})}{\partial \hat\theta}\right] \tag{4}
$$
减去基流方程有：
$$
\delta\hat{p}_{2}+\hat{\overline{v}}_{\theta 2}\delta{\hat{v}}_{\theta 2}- \delta\hat{p}_{1}-\hat{\overline{v}}_{\theta 1}\delta{\hat{v}}_{\theta 1}  
= 
 (\tan\beta_2 - \tan\alpha_1)\delta\hat{v}_{x1} - \delta\hat{{l}}_R - \lambda_{\text{rot}}\left[\frac{\partial \delta\hat{v}_{x1}}{\partial \tau}+\frac{\partial \delta\hat{v}_{x1}}{\partial \hat\theta}\right] \tag{4}
$$
两边取拉普拉斯变换：
$$
\delta\hat{P}_{2}+\hat{\overline{v}}_{\theta 2}\delta\hat{V}_{\theta 2}- \delta\hat{P}_{1}-\hat{\overline{v}}_{\theta 1}\delta\hat{V}_{\theta 1}  
= 
 (\tan\beta_2 - \tan\alpha_1)\delta\hat{v}_{x1} - \delta\hat{\mathcal{l}}_R - \lambda_{\text{rot}}\left[s \delta\hat{v}_{x1}+\frac{\partial \delta\hat{v}_{x1}}{\partial \hat\theta}\right] \tag{5}
$$
周向傅里叶变换展开，并取$n$阶谐波,公因子$e^{jn\theta}$消去：
$$
\delta\tilde{P}_{2}+\hat{\overline{v}}_{\theta 2}\delta\tilde{V}_{\theta 2}- \delta\tilde{P}_{1}-\hat{\overline{v}}_{\theta 1}\delta\tilde{V}_{\theta 1}  
= 
 (\tan\beta_2 - \tan\alpha_1)\delta\tilde{v}_{x1} - \delta\tilde{{L}}_R - \lambda_{\text{rot}}\left(s +jn\right)\delta\tilde{v}_{x1}\tag{6}
$$
先放在这，让我们去看式(3)经过这番操作后变成什么了。
## 3.2 损失项
### 3.2.1 速度表示形式
左边：

$$
\left[1 + \tau_R\left(\frac{\partial}{\partial\tau} + \frac{\partial}{\partial\hat\theta}\right)\right] (\hat{\overline{l}}_R + \delta \hat l_R)
= \hat{\overline{l}}_R + \left[1 + \tau_R\left(\frac{\partial}{\partial\tau} + \frac{\partial}{\partial\hat\theta}\right)\right] \delta \hat l_R
$$
（因为 $\overline{l}_R$ 为常数，其导数为零）

右边：
$$
\hat{l}_R^{ss}(\hat{\overline{v}}_{x1} + \delta \hat v_{x1},\; \hat{\overline{v}}_{\theta 1} + \delta \hat v_{\theta1})
$$
可以将其分解为“背景值 + 扰动增量”：
$$
\hat{l}_R^{ss}(\hat{\overline{v}}_x + \delta\hat v_{x1},\; \hat{\overline{v}}_\theta + \delta\hat v_{\theta1})
= \underbrace{\hat{l}_R^{ss}(\hat{\overline{v}}_{x1}, \hat{\overline{v}}_{\theta 1})}_{\text{背景值}}
+ \underbrace{\Big[ \hat{l}_R^{ss}(\hat{\overline{v}}_{x1} + \delta \hat v_{x1},\; \hat{\overline{v}}_{\theta 1} + \delta \hat v_{\theta1}) - \hat{l}_R^{ss}(\hat{\overline{v}}_x, \hat{\overline{v}}_\theta) \Big]}_{\text{扰动增量} \equiv \delta\hat{l}_R^{ss}}
$$
把这俩代入式(3)得到：
$$
\hat{\overline{l}}_R + \left[1 + \tau_R(\partial_\tau + \partial_\theta)\right] \delta l_R
= \underbrace{\hat{l}_R^{ss}(\hat{\overline{v}}_{x1}, \hat{\overline{v}}_{\theta 1})}_{\hat{\overline{l}}_R} + \delta\hat{l}_R^{ss}
$$
约去 $\hat{\overline{l}}_R$：

$$
\left[1 + \tau_R\left(\frac{\partial}{\partial\tau} + \frac{\partial}{\partial\hat\theta}\right)\right] \delta l_R = \delta\hat{l}_R^{ss}\tag{7}

$$

$\delta\hat{l}_R^{ss}$ 是 $\delta\hat v_{x1}, \delta\hat v_{\theta1}$ 的非线性函数。为得到线性方程，保留其一阶近似（线性主部），这里需要注意的是，我们冻结的是偏导数部分，而不是乘上的扰动量：
$$
\delta\hat{l}_R^{ss} \approx \frac{\partial \hat{l}_R^{ss}}{\partial \hat v_{x1}}\bigg|_{\text{avg}} \delta\hat v_{x1}
+ \frac{\partial \hat{l}_R^{ss}}{\delta\hat v_{\theta1}}\bigg|_{\text{avg}} \delta\hat v_{\theta1}
\tag{8}
$$
将 (8) 代入 (7)，得到线性化扰动方程：
$$
\left[1 + \tau_R(\partial_\tau+\partial_\theta)\right] \delta\hat{l}_R
= \frac{\partial \hat{l}_R^{ss}}{\partial \hat v_{x1}}\bigg|_{\text{avg}}\delta\hat v_{x1}
+ \frac{\partial \hat{l}_R^{ss}}{\partial\hat v_{\theta1}}\bigg|_{\text{avg}} \delta\hat v_{\theta1}
\tag{9}
$$
对时间 $\tau$ 作拉普拉斯变换：$\partial/\partial\tau \to s$。  

对周向 $\hat\theta$ 作傅里叶级数展开，取第 $n$ 阶谐波：$\partial/\partial\hat\theta \to jn$，且记
$$
\delta\hat v_{x1} = \delta\tilde{V}_{x1} e^{jn\hat\theta},\quad
\delta\hat v_{\theta1} = \delta\tilde{V}_{\theta1} e^{jn\hat\theta},\quad
\delta\hat{l}_R = \delta\tilde{L}_R e^{jn\hat\theta}
$$

代入 (9)，把公因子$e^{jn\theta}$约掉：
$$
\left[1 + \tau_R(s + jn)\right] \delta\tilde{L}_R
= \frac{\partial \hat{l}_R^{ss}}{\partial \hat v_{x 1}}\bigg|_{\text{avg}} \delta\tilde{V}_{x1}
+ \frac{\partial \hat{l}_R^{ss}}{\partial \hat v_{\theta 1}}\bigg|_{\text{avg}} \delta\tilde{V}_{\theta1}
\tag{10}
$$

解得：
$$
\boxed{
\delta\tilde{L}_R = \frac{ \frac{\partial \hat{l}_R^{ss}}{\partial\hat v_{x 1}}\bigg|_{\text{avg}}\delta \tilde{V}_{x1}
+ \frac{\partial \hat{l}_R^{ss}}{\partial \hat v_{\theta 1}}\bigg|_{\text{avg}}\delta \tilde{V}_{\theta1}}{1 + \tau_R(s+ jn)}}\tag{11}
$$
### 3.2.2 角度表示形式
我们想把这个冻结的导数给换掉，由无量纲相对气流角定义：
$$
\tan\beta_1 = \frac{\hat v_{\theta1} - 1}{\hat v_{x1}}
$$
根据链式法则：
$$
\frac{\partial \hat{l}_R^{ss}}{\partial \hat v_{x1}}\bigg|_{\text{avg}} 
= \frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot \frac{\partial \tan\beta_1}{\partial \hat {\overline v}_{x1}}
= \frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot \left( -\frac{\hat{\overline v}_{\theta1} - 1}{\hat {\overline v}_{x1}^2} \right)
= -\,\frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot \frac{\tan\beta_1}{\hat{\overline v}_{x1}}
$$
$$
\frac{\partial \hat{l}_R^{ss}}{\partial \hat v_{\theta1}} \bigg|_{\text{avg}}= \frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot \frac{\partial \tan\beta_1}{\partial \hat{
\overline v}_{\theta1}}
= \frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot \frac{1}{\hat{\overline v}_{x1}}
$$

将上述偏导数代入 (11)：
$$
\delta\tilde{L}_R = \frac{ \frac{d\hat{l}_R^{ss}}{d\tan\beta_1} \bigg|_{\text{avg}}\left( -\frac{\tan\beta_1}{\hat{\overline v}_{x1}} \delta\tilde{V}_{x1} + \frac{1}{\hat{\overline v}_{x1}} \delta\tilde{V}_{\theta1} \right) }{1 + \tau_R(s+jn)}
$$

整理得：
$$
\boxed{
\delta\tilde{L}_R =  \frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}}\cdot \frac{   \delta\tilde{V}_{\theta1} - \tan\beta_1 \delta\tilde{V}_{x1}} {\hat{\overline v}_{x1}[1 + \tau_R(s+jn)]}}
\tag{12}
$$
其中$\dfrac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}}$ 在工作点取值，有$M-G$模型可以参考：
$$
\frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}}=\frac{d\hat{l}_R^{ss}}{d\hat v_{x1}}\bigg|_{\text{avg}}\cdot\frac{1}{(\tan{\alpha_1}-\tan{\beta_1})^2}
$$

## 4. 轴流转子转子的第 $n$ 阶谐波匹配条件
将损失表达式 (12) 代入压力平衡方程 (6)：
$$
\delta\tilde{P}_2 + \hat{\overline{v}}_{\theta2}\,\delta\tilde{V}_{\theta2}
- \delta\tilde{P}_1 - \hat{\overline{v}}_{\theta1}\,\delta\tilde{V}_{\theta1}
=  [\tan\beta_2 - \tan\alpha_1- \lambda_{\text{rot}}(s+jn)]\delta\tilde{V}_{x1} -\frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}}\cdot \frac{   \delta\tilde{V}_{\theta1} - \tan\beta_1 \delta\tilde{V}_{x1}} {\hat{\overline v}_{x1}[1 + \tau_R(s+jn)]}
$$
# 5. 压力
将第 4 节得到的匹配方程与速度匹配条件结合，可写出转子执行盘在频域、第 $n$ 阶谐波下的传递矩阵形式
已知速度匹配条件：
$$
\delta\tilde{V}_{x2} = \delta\tilde{V}_{x1}, \qquad
\delta\tilde{V}_{\theta2} = \tan\beta_2 \cdot \delta\tilde{V}_{x1}
$$


其中
$$
\delta\tilde{L}_R = 
\frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot
\frac{ \delta\tilde{V}_{\theta1} - \tan\beta_1 \delta\tilde{V}_{x1} }
{ \hat{\overline{v}}_{x1}\bigl[1 + \tau_R(s+jn)\bigr] }
\tag{12}
$$

将 $\delta\tilde{V}_{\theta2} = \tan\beta_2 \delta\tilde{V}_{x1}$ 代入上式左边，并移项解出 $\delta\tilde{P}_2$：
$$
\delta\tilde{P}_2 = \delta\tilde{P}_1 + \hat{\overline{v}}_{\theta1}\delta\tilde{V}_{\theta1}
- \hat{\overline{v}}_{\theta2}\tan\beta_2 \delta\tilde{V}_{x1}
+ \left( \tan\beta_2 - \tan\alpha_1 - \lambda_{\text{rot}}(s+jn) \right) \delta\tilde{V}_{x1}
- \delta\tilde{L}_R
$$

将 $\delta\tilde{L}_R$ 代入，合并 $\delta\tilde{V}_{x1}$ 和 $\delta\tilde{V}_{\theta1}$ 的系数：
$$
\delta\tilde{P}_2 = \delta\tilde{P}_1
+ \left( \hat{\overline{v}}_{\theta1} - \frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot \frac{1}{\hat{\overline{v}}_{x1}[1+\tau_R(s+jn)]} \right) \delta\tilde{V}_{\theta1}
+ \Biggl[ \tan\beta_2 - \tan\alpha_1 - \lambda_{\text{rot}}(s+jn) - \hat{\overline{v}}_{\theta2}\tan\beta_2
+ \frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot \frac{\tan\beta_1}{\hat{\overline{v}}_{x1}[1+\tau_R(s+jn)]} \Biggr] \delta\tilde{V}_{x1}
$$
# 5. 传递矩阵
记
$$
K = \frac{d\hat{l}_R^{ss}}{d\tan\beta_1}\bigg|_{\text{avg}} \cdot \frac{1}{\hat{\overline{v}}_{x1}[1+\tau_R(s+jn)]}
$$
则：
$$
B_{31} = \tan\beta_2 - \tan\alpha_1 - \lambda_{\text{rot}}(s+jn) - \hat{\overline{v}}_{\theta2}\tan\beta_2 + K \tan\beta_1
$$
$$
B_{32} = \hat{\overline{v}}_{\theta1} - K
$$

因此，转子执行盘的传递矩阵为：
$$
\mathbb{B}_n
=
\begin{bmatrix}
1 & 0 & 0 \\[2pt]
\tan\beta_2 & 0 & 0 \\[2pt]
B_{31} & B_{32} & 1
\end{bmatrix}


$$
$$
\boxed{
\begin{bmatrix}
\delta\hat{V}_{x2} \\[2pt]
\delta\hat{V}_{\theta2} \\[2pt]
\delta\hat{P}_2
\end{bmatrix}
(\hat\theta,s)
=
\sum_{n=0}^{\infty}\mathbb{B}_n\cdot e^{jn\hat\theta}\cdot
\begin{bmatrix}
\delta\tilde{V}_{x1} \\[2pt]
\delta\tilde{V}_{\theta1} \\[2pt]
\delta\tilde{P}_1
\end{bmatrix}
(s)}
$$