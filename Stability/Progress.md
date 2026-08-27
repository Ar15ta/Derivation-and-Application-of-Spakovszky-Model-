# 项目进展：NCEPU 离心压气机稳定性模型

基于 Spakovszky CC3 执行盘理论的旋转失速/喘振模型。数据链路：
CFX 后处理 → `Geo_and_Base_Flow/NCEPU_CFX.csv` → `NCEPU_CFX_Loader`
（几何 + CoolProp 物性 → 无量纲参数字典）→ `Prediction/NCEPU/Compute_NCEPU_Margin`
（逐工况特征值求解，存结果 CSV）→ `Post_Process/Plot_*` 绘图。

## 一、基础约定（已完成）

- **站位与两层命名**：STA2=叶轮进口、STA3=叶轮出口、STA4=扩压器出口。几何/CSV/loader 用站位下标（`R_3`、`b_3`、`c_m2`、`radial_sta3` 等，"物理意义在前、站位在后"）；Matrix 组件接口用部件号 1=入口/2=出口（`Vx_bar1_imp`、`beta2_imp_deg`）。Matrix 层保持不动，由 loader 做两层翻译。
- **速度符号**（Cumpsty/Greitzer）：c=绝对速度，w=相对速度，u=轮缘速度；分量 m=子午、u=周向、r=径向、x=轴向；p/pt=静/总压，T/Tt=静/总温，a=静声速。
- **切向翻译**：模型统一 `u3=|ω|·R_3>0`、切向正方向=旋转方向、后掠 χ<0；loader 取 `c_u_sign=sign(rpm)`，`c_u(model)=c_u_sign·c_u(CFX)`，u3 取绝对值。rpm 与 c_u 同号（如 -32500 rpm 配负 c_u）可正确处理；`z_axis` 仅记录坐标系，不参与切向翻译。
- **出口角闭合**：α/β 与相对速度 w 全部由分量反算（w_θ=c_u−u(r)，β=atan2(w_θ,c_m)），CSV 不填角度列；金属角 χ + Wiesner 滑移 σ_W=1−√(cos χ)/(Z+Z_split)^0.7 仅在 c_u3 缺失时降级使用。
- **params.txt**：数值解析器只接受纯数字（不支持算式）；关键值 `l_comp=0.0409`、`t_imp=0.0018`、`tau_u=1.0`、`fluid=CO2`。

## 二、CFX 数据流（loader）

- **CSV 24 列**：每站位 7 个量（`c_m,c_u,rho,p,T,pt,Tt`）+ `rpm,m` + `dL_dtanbeta1`。角度 α/β、相对速度 w、π_tt、a、dh_t 全部由代码派生，不填列；eta 仅打印用已全链路删除。`write_csv_template` 只出干净表头。
- **平均方式（最新结论）**：CFX 截面导出量（含 c_m/c_u/rho 与 p/T/pt/Tt）统一用 **areaAveraged（面积平均）**，实测全部面积平均更准确；c_u 填原始带符号值，loader 按 rpm 符号翻转。
- **真实气体物性**：工质由 params.txt `fluid=CO2` 指定（loader 用 `read_fluid` 单独解析字符串行）。静声速 `a=PropsSI('A',T,p)`、总焓 `h=PropsSI('H',Tt,pt)`、`dh_t=h3−h2`；CoolProp 缺失/失败返回 0，不中断（sCO₂ 不做 γ=1.4 回退）。
- **几何/惯性派生**：`r_2rms=√((r_2s²+r_2h²)/2)`；`A3=b_3·(2πR_3−Z_total·t_imp)` 含尾缘堵塞；面积密度比 `AR_imp=ρ3A3/(ρ2A2)`；`s_imp=l_comp/R_3`；`λ_imp=s·AR·ln(AR)/(AR−1)`；`τ_imp=τu·2s_imp/(Ŵ1+Ŵ2)`，Ŵ=w/u3（τu 是经验系数，默认 1.0）。
- **损失导数 dL/dtanβ1**：优先级 CSV 直填 > 自动差分 > 置零。自动差分用 Spakovszky 间接法 μ=Δh_t/u3²、ψ=(pt3−pt2)/(ρ2 u3²)、L=μ−ψ。邻点取同转速高/低流量侧各自**最近点（无流量容差）**：两侧都有做中心差分；端点只有单侧时退化为单侧差分（最小流量点用高流量侧迎风前向，最大流量点用低流量侧后向）。

## 三、特征值求解与绘图

- **扫描顺序**（Compute_NCEPU_Margin）：同转速内按流量**降序**（大流量稳定侧 → 小流量喘振侧），保证热启动从稳定侧连续追踪、早停在越过 σ=0 时才触发（原升序会让最小流量点首点即早停）。
- **结果 CSV**：`Prediction/NCEPU/NCEPU_Eigenvalues_<rpm>.csv`，表头 `rpm,m,pi_tt,n,sigma,omega`；关键参数（dL/AR/lambda）缺失的工况自动 SKIP。
- **绘图**：
  - 特征值图 `Plot_NCEPU_Stability_Margain.py`：`--csv` 可缺省（自动取 `Prediction/NCEPU` 下最新结果）；设计流量 `M_DESIGN=13.5 kg/s`（可用 `--m-design` 覆盖）；边界标注取稳定侧 σ 最大点（失稳前）与不稳定侧 σ 最小点（失稳点），二者分居 σ=0 两侧，流量标注精确到 4 位小数。
  - 能量分布图 `Plot_NCEPU_Distribution.py`：同样把 `M_DESIGN` 更新为 13.5，消费 acoustics（CoolProp 提供的 a/rho/p），接口不变。

## 四、改动文件清单

- `Geo_and_Base_Flow/NCEPU_params.txt`：`t_imp=0.0018`、`tau_u=1.0`、`fluid=CO2`、`l_comp=0.0409`。
- `Geo_and_Base_Flow/NCEPU_CFX_Loader.py`：24 列表头；全部量 areaAveraged 约定；CoolProp 查声速/焓；角度/w/π_tt/dh_t 派生；几何/惯性/尾缘堵塞；自动差分改最近邻（无容差、端点单侧）；删除 eta。
- `Geo_and_Base_Flow/NCEPU_CFX.csv`：已填真实 -32500 rpm 工况。
- `Prediction/NCEPU/Compute_NCEPU_Margin.py`：流量降序扫描；结果表头删 eta。
- `Post_Process/Plot_NCEPU_Stability_Margain.py`：--csv 缺省自动查找；M_DESIGN=13.5；边界标注分居两侧、4 位小数；清理冗余（m_design 默认值、合并 track 单遍绘制等）。
- `Post_Process/Plot_NCEPU_Distribution.py`：M_DESIGN=13.5。
- `Model/Model_NCEPU.py`、`Matrix/*`：**未改动**（保持组件 1/2 接口）。

## 五、验证（真实数据，-32500 rpm）

- 切向翻译正常：rpm=−32500、c_u 同号为负，翻转后 u3>0、β3 后掠为负，π_tt≈1.54。
- 自动差分对 3 个流量点（8.696/8.644/6.603 kg/s）全部生效（端点单侧、中间中心），不再整表 SKIP。
- 扫描方向为大流量→小流量；8.696 点 σ 全负（稳定），8.644 点 σ 转正（n=1~5 失稳），早停在该点触发——喘振边界落在 m/13.5 ≈ 0.64（8.696→8.644 之间）。

## 六、待办 / 备注

- 边界附近补更密流量点（当前 8.644 与 8.696 仅差 0.05 kg/s，且 6.603 单侧差分跨距较大），提高边界与差分精度。
- `rho4` 建议直接给出（缺省回退 rho3）；`pt4/Tt4` 已提取但暂未被消费，留作扩压器总压恢复扩展。
- `b4`、plenum 体积、进口金属角当前模型不需要（plenum 边界直接加在扩压器出口）。
- 直接跑 loader 文件会因 `from Matrix...` 报 `No module named 'Matrix'`；应从 Stability 根目录用 `python -m Geo_and_Base_Flow.NCEPU_CFX_Loader`（Compute/Plot 脚本已带 sys.path 引导，可直接运行）。
- 特征值图 PNG 默认存到 `Spa_RGP/` 根目录（三层 dirname）；如需落 `Stability/` 下再调路径。
