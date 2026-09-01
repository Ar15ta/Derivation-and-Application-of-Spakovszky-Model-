"""Plot_Eckardt_O_Distribution.py - Eckardt O 模态能量空间分布绘图

从已有特征值 CSV 读取特征值, 对指定工况的各 n 模态计算其在物理空间的
状态向量分布, 然后绘制 Energy Function (声学扰动能量密度) 沿物理空间
的分布曲线, 复现 Spakovszky 论文 Figure 5-18 风格。

──────────────────────────────────────────────────────────────────
物理空间分段 (归一化坐标, R2=0.2m, 全程未改动物理模型):
  - 上游管道: x ∈ (-∞, 0]   (无限长管道, x=0 = 叶轮进口 STA2, 负值)
  - 叶轮:     x = 0           (执行盘, 零厚度, 能量跳变)
  - 扩压器:   s ∈ [0, 0.69]   (径向弧长 = r - r3, 正值, r3=1.0, r4=1.69)
──────────────────────────────────────────────────────────────────

Energy Function (声学扰动能量流密度, 论文 Eq.5-23):
  δε_n(x, s) = P̃_n(x, s) · Ṽ_n(x, s)
  = (压力扰动) × (流向速度扰动)
  物理意义: 声学功率通量密度 (acoustic intensity), 取实部保留能量方向
  归一化: 除以本模态上游入口 δε_inlet, 即 δε/δε_inlet

──────────────────────────────────────────────────────────────────
关于 x 轴范围 [-0.15, +0.15] 的实现机制 (重要: 不改动物理空间与边界条件):

  物理模型仍是 "无限长上游管道 + 零厚度叶轮 + 0.69 弧长扩压器 + 气室".
  x 轴 [-0.15, +0.15] 的显示范围由两层独立机制实现, 均为后处理性质:

  1) 计算采样层 (Model_Eckardt_O.compute_modal_shape 的 x_upstream_max 参数):
     x_up = np.linspace(-x_upstream_max, 0.0, n_points)
     x_upstream_max=0.15 仅决定在哪些物理点上 "读出" 模态形状 (类似探针位置),
     SVD 提取、特征方程 det(Y)=0、状态传递矩阵 M@coeffs 全部不变,
     上游仍按无限长管道建模.

  2) 视觉显示层 (matplotlib 视图裁剪):
     ax.set_xlim(-0.15, 0.15)
     扩压器实际计算到 x=+0.69, set_xlim 只显示前 ~12% 弧长 (r ∈ [1.0, ~1.12]).
     这只是 matplotlib 的视图裁剪, 不影响底层数组与物理结果.

  边界条件完全未变:
    - 上游: extract_modal_coeffs(Y, force_upstream_bc=True) → B_n=C_n=0
      (强制消除上游反射波 = 无限长管道无反射条件)
    - 下游: EC_downstream_plenum 在扩压器出口直接施加气室边界 (无下游管道)
──────────────────────────────────────────────────────────────────

用法:
  python Plot_Eckardt_O_Distribution.py --rpm 14000 --n-max 3
  python Plot_Eckardt_O_Distribution.py --rpm 14000 --m 5.0 --n 1
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import argparse
import sys
import os

# 将 Stability 根目录加入 sys.path
STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if STABILITY_ROOT not in sys.path:
    sys.path.insert(0, STABILITY_ROOT)

from Model.Model_Eckardt_O import compute_modal_shape
from Matrix.rotor import load_params
from Prediction.Eckardt_O.Compute_Eckardt_O_Stability_Margin import (
    load_results, load_yaml, find_record
)
from Geo_and_Base_Flow.Eckardt_O_Geo_Calculator import compute_base_flow

# 参数文件路径
PARAMS_FILE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Geo_and_Base_Flow", "Eckardt_O_params.txt"
)

# 均值线包路径
MEANLINE_PKG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "MeanLIne_Design", "radcomp-main"
)
YAML_DATA = os.path.join(MEANLINE_PKG, "data", "known_compressors.yml")

# 理想气体比热比 (空气)
GAMMA = 1.4

# 设计点参考质量流量 (14000 RPM, 归一化用)
M_DESIGN = 5.31   # kg/s


# ════════════════════════════════════════════════════════════════
# 基流声学参数计算
# ════════════════════════════════════════════════════════════════
def compute_station_acoustics(geo_data):
    """从 geo_data 计算各站位静态声速 c [m/s] 和密度 ρ [kg/m³]

    各站位数据可用性:
      sta1: 有 a0(总声速), M, P        → c = a0/√(1+0.2M²),  ρ = γP/c²
      sta2: 有 P, ρ                     → c = √(γP/ρ)
      sta3: 有 P, ρ                     → c = √(γP/ρ)
      sta4: 有 c(速度), M, P            → c_sound = V/M,      ρ = γP/c²
    """
    U2 = geo_data['U2']

    def _sound_speed(sta):
        if sta == 'sta1':
            a0 = geo_data['sta1_a0']
            M = geo_data['sta1_M']
            return a0 / np.sqrt(1.0 + 0.2 * M * M)
        elif sta in ('sta2', 'sta3'):
            P = geo_data[f'{sta}_P']
            rho = geo_data[f'{sta}_rho']
            return np.sqrt(GAMMA * P / rho)
        elif sta == 'sta4':
            c_vel = geo_data['sta4_c']      # 气流速度
            M = geo_data['sta4_M']
            return c_vel / M if M > 1e-6 else 340.0
        raise ValueError(sta)

    def _density(sta):
        if sta == 'sta1':
            P = geo_data['sta1_P']
            c = _sound_speed('sta1')
            return GAMMA * P / (c * c)
        elif sta in ('sta2', 'sta3'):
            return geo_data[f'{sta}_rho']
        elif sta == 'sta4':
            P = geo_data['sta4_P']
            c = _sound_speed('sta4')
            return GAMMA * P / (c * c)
        raise ValueError(sta)

    acoustics = {
        'U2': U2,
        'rho': {s: _density(s) for s in ('sta1', 'sta2', 'sta3', 'sta4')},
        'c':   {s: _sound_speed(s) for s in ('sta1', 'sta2', 'sta3', 'sta4')},
    }
    # 参考马赫数 M_ref = U2/c
    acoustics['M_ref'] = {s: U2 / acoustics['c'][s]
                          for s in ('sta1', 'sta2', 'sta3', 'sta4')}
    return acoustics


# ════════════════════════════════════════════════════════════════
# 能量函数 (论文定义: δε_n = P̃_n · Ṽ_n)
# ════════════════════════════════════════════════════════════════
def energy_function(q_array, rho, U2, c):
    """计算声学扰动能量函数 δε_n = P̃_n · Ṽ_n (论文 Eq.5-23 定义)

    物理意义:
      声学扰动功率通量密度 (acoustic intensity), 即单位面积上声能流的方向 & 大小.
      与能量密度 (1/2·ρ·U²·|δV|²) 不同, 能量函数刻画的是 "能量流向", 因此
      必须保留符号 (实部), 而非取幅值. 论文 Figure 5-18 即采用此定义.

    定义式:
      δε_n(x, s) = P̃_n(x, s) · Ṽ_n(x, s)
                 = (复数压力扰动) × (复数流向速度扰动)
      取实部 → 功率通量的物理方向 (正=向下游, 负=向上游)

    输入:
      q_array : (N, 3) complex, 归一化状态向量 [δV_stream, δV_θ, δP]
                (速度归一化 /U2, 压力归一化 /(ρ·U2²))
                - 上游管道: q[0]=δVx (轴向, 流向)
                - 扩压器:   q[0]=δVr (径向, 流向)
                - 叶轮出口: q[0]=δVr (径向, 流向)
      rho     : 局部密度 [kg/m³]   (本函数未直接使用, 保留接口便于扩展)
      U2      : 叶轮出口轮缘速度 [m/s]  (本函数未直接使用, 归一化常数在调用处消去)
      c       : 局部静态声速 [m/s]  (本函数未直接使用, 保留接口便于扩展)

    返回:
      δε : (N,) float, 能量函数的实部 (保留正负号 → 反映能量流向)

    注: 归一化常数 (ρ·U2³) 在 plot_energy_multi_n 的 scale 处统一消去,
        故本函数直接用归一化状态量相乘即可.
    """
    # 论文定义: δε = P̃ · Ṽ = δP · δV_stream
    # q_array[:, 2] = δP (归一化压力扰动)
    # q_array[:, 0] = δV_stream (归一化流向速度扰动)
    delta_eps = q_array[:, 2] * q_array[:, 0]
    # 取实部: 复数乘积的实部 = 功率通量在时间平均下的方向分量
    return delta_eps.real


# ════════════════════════════════════════════════════════════════
# 模态能量分布计算
# ════════════════════════════════════════════════════════════════
def compute_energy_distribution(shape_result, acoustics, eps_imp=0.03):
    """把 compute_modal_shape 的结果转换为能量分布 (按物理空间三段切分)

    物理空间分段 (归一化坐标, R2=0.2m):
      段0 Upstream : x ∈ [-0.2, 0]   (仅显示窗口, 物理模型仍为无限长管道)
      段1 Impeller : x = 0            (零厚度执行盘, 入口/出口状态跳变)
      段2 Diffuser : x ∈ [0, 0.69]   (s = r - r3, r3=1.0, r4=1.69)

    各段对应站位声学参数:
      Upstream  → sta1 (诱导轮入口 ≈ 上游管道远端基流)
      Impeller  → 入口用 sta2, 出口用 sta3 (执行盘跨越叶轮, 进出口不同)
      Diffuser  → sta3 (入口) 线性插值到 sta4 (出口)

    输入:
      shape_result : compute_modal_shape 的返回 dict (含 'segments' 三段)
      acoustics    : compute_station_acoustics 的返回 dict (含各站位 ρ, c, M_ref)
      eps_imp      : 叶轮显示半宽 (保留参数, 当前实现已改为 x=0 垂直阶跃, 不再用)

    返回 dict:
      'x'        : list of np.ndarray, 各段物理坐标 (上游负值, 扩压器正值)
      'E'        : list of np.ndarray, 各段能量函数值 (实部, 保留符号)
      'bounds'   : list of (x_start, x_end, seg_name), 各段坐标范围
      'seg_meta' : list of dict, 各段物理参数 (rho, c, M_ref)
    """
    U2 = acoustics['U2']
    segments = shape_result['segments']
    # segments[0]=Upstream, [1]=Impeller(出口), [2]=Diffuser

    x_list, E_list, bounds, seg_meta = [], [], [], []

    # ── 上游管道: x ∈ [-L, 0], 用 sta1 基流 (诱导轮入口 ≈ 上游管道) ──
    up = segments[0]
    x_up = up['coord']                      # 负值 [-L, 0]
    q_up = up['q']
    rho_up = acoustics['rho']['sta1']
    c_up = acoustics['c']['sta1']
    E_up = energy_function(q_up, rho_up, U2, c_up)
    x_list.append(x_up)
    E_list.append(E_up)
    bounds.append((x_up[0], x_up[-1], 'Upstream'))
    seg_meta.append({'rho': rho_up, 'c': c_up, 'M_ref': U2 / c_up})

    # ── 叶轮: x=0 处能量跳变 (入口状态 → 出口状态) ──
    # 入口状态 = 上游末端; 出口状态 = segments[1].q[0]
    q_imp_in = segments[0]['q'][-1]          # 叶轮入口 [δVx, δVθ, δP]
    q_imp_out = segments[1]['q'][0]          # 叶轮出口 [δVr, δVθ, δP]
    rho2, c2 = acoustics['rho']['sta2'], acoustics['c']['sta2']
    rho3, c3 = acoustics['rho']['sta3'], acoustics['c']['sta3']

    E_imp_in = energy_function(q_imp_in.reshape(1, 3), rho2, U2, c2)[0]
    E_imp_out = energy_function(q_imp_out.reshape(1, 3), rho3, U2, c3)[0]

    # 叶轮: 零厚度执行盘, 入口/出口均在 x=0 (垂直阶跃)
    x_list.append(np.array([0.0, 0.0]))
    E_list.append(np.array([E_imp_in, E_imp_out]))
    bounds.append((0.0, 0.0, 'Impeller'))
    seg_meta.append({'rho_in': rho2, 'c_in': c2, 'rho_out': rho3, 'c_out': c3})

    # ── 扩压器: s = r - r3 ∈ [0, r4-r3], 密度声速线性插值 ──
    diff = segments[2]
    r_diff = diff['coord']                   # [1.0, 1.69]
    q_diff = diff['q']
    s_diff = r_diff - r_diff[0]              # 弧长 [0, 0.69]
    x_diff = s_diff                           # 从 x=0 开始 (紧接叶轮出口)

    rho3, rho4 = acoustics['rho']['sta3'], acoustics['rho']['sta4']
    c3, c4 = acoustics['c']['sta3'], acoustics['c']['sta4']
    t = (r_diff - r_diff[0]) / (r_diff[-1] - r_diff[0])
    rho_diff = rho3 + (rho4 - rho3) * t
    c_diff = c3 + (c4 - c3) * t

    E_diff = np.array([
        energy_function(q_diff[i:i + 1], rho_diff[i], U2, c_diff[i])[0]
        for i in range(len(q_diff))
    ])
    x_list.append(x_diff)
    E_list.append(E_diff)
    bounds.append((x_diff[0], x_diff[-1], 'Diffuser'))
    seg_meta.append({'rho': rho_diff, 'c': c_diff})

    return {'x': x_list, 'E': E_list, 'bounds': bounds, 'seg_meta': seg_meta}


# ════════════════════════════════════════════════════════════════
# 绘图: 单 n 能量分布
# ════════════════════════════════════════════════════════════════
def plot_energy_single_mode(shape_result, acoustics, m_flow, rpm,
                             save_path=None, show=True, normalize=True):
    """绘制单个模态的 Energy Function 沿物理空间分布"""
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'

    dist = compute_energy_distribution(shape_result, acoustics)
    x_list, E_list, bounds = dist['x'], dist['E'], dist['bounds']

    s_star = shape_result['eig']
    n_val = shape_result['n']

    fig, ax = plt.subplots(figsize=(11, 6))

    # 归一化: 除以上游入口值 (使入口处 E/E_in=1, 对数坐标下起点=0 dB)
    E_inlet = E_list[0][0]   # 上游管道最远端 (= 入口边界条件处)
    scale = E_inlet if (normalize and E_inlet > 0) else 1.0

    colors = {'Upstream': '#1f77b4', 'Impeller': '#d62728', 'Diffuser': '#2ca02c'}
    labels = {'Upstream': 'Upstream duct',
              'Impeller': 'Impeller (actuator disk)',
              'Diffuser': 'Vaneless diffuser'}

    # 画各段
    for x_seg, E_seg, (_, _, name) in zip(x_list, E_list, bounds):
        c = colors.get(name, 'gray')
        E_plot = E_seg / scale
        if name == 'Impeller':
            # 叶轮: 用阶梯线连接入口和出口 (显示跳变)
            x_step = np.array([x_seg[0], x_seg[1], x_seg[1]])
            E_step = np.array([E_plot[0], E_plot[0], E_plot[1]])
            ax.plot(x_step, E_step, '-', color=c, linewidth=2.5, label=labels[name])
            ax.scatter(x_seg, E_plot, color=c, s=90, zorder=5,
                       edgecolors='black', linewidths=0.8)
            # 阶跃百分比标注
            step_pct = abs(E_plot[1] - E_plot[0]) / E_plot[0] * 100
            step_dir = '↑' if E_plot[1] > E_plot[0] else '↓'
            x_mid_step = x_seg[1] + 0.06
            y_mid_step = np.sqrt(E_plot[0] * E_plot[1])
            ax.annotate(f'Impeller step: {step_dir} {step_pct:.2f}%',
                        xy=(x_seg[1], E_plot[1]),
                        xytext=(x_mid_step, y_mid_step),
                        fontsize=9, color=c, fontweight='bold',
                        arrowprops=dict(arrowstyle='->', color=c, lw=1.2),
                        va='center')
        else:
            ax.plot(x_seg, E_plot, '-', color=c, linewidth=2.5, label=labels[name])

    # 叶轮分界竖线
    ax.axvline(x=0, color='gray', linestyle=':', linewidth=1.3, alpha=0.7)

    # 段背景着色
    ylim = ax.get_ylim()
    for (s, e, name) in bounds:
        c = colors.get(name, 'gray')
        if name == 'Impeller':
            ax.axvspan(s, e, color=c, alpha=0.10)
        else:
            ax.axvspan(s, e, color=c, alpha=0.06)

    # 段名标注
    for (s, e, name) in bounds:
        mid = 0.5 * (s + e) if e > s else s
        ax.annotate(labels.get(name, name),
                    xy=(mid, 0.96), xycoords=('data', 'axes fraction'),
                    ha='center', va='top', fontsize=10, style='italic',
                    color=colors.get(name, 'gray'))

    # 轴标签
    ax.set_xlabel(r'Streamwise coordinate  $x/R_2$  (upstream)  /  '
                  r'$s/R_2$  (diffuser)', fontsize=13)
    ylabel = r'Energy function  $E(x)/E_{\mathrm{inlet}}$  (log scale)'
    ax.set_ylabel(ylabel, fontsize=12)

    # 标题
    sigma = s_star.real
    omega = -s_star.imag
    stability = 'UNSTABLE' if sigma > 0 else 'stable'
    title = (f'Energy function distribution  |  {int(rpm)} rpm  '
             r'$\dot{m}$' + f'={m_flow:.3f} kg/s  |  n={n_val}  |  '
             r'$\sigma$' + f'={sigma:+.5f}, '
             r'$\omega$' + f'={omega:+.5f}  ({stability})')
    ax.set_title(title, fontsize=11.5, pad=12)

    # 刻度
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.tick_params(axis='both', which='minor', length=4)

    # 对数 y 轴
    ax.set_yscale('log')

    ax.legend(loc='upper right', fontsize=10, framealpha=1.0,
              facecolor='white', edgecolor='black')
    ax.grid(True, which='both', alpha=0.25)
    ax.grid(True, which='minor', alpha=0.12)

    # 底部标注各站位声速/密度
    info = (f"sta1: ρ={acoustics['rho']['sta1']:.3f} kg/m³, "
            f"c={acoustics['c']['sta1']:.1f} m/s  |  "
            f"sta2: ρ={acoustics['rho']['sta2']:.3f}, "
            f"c={acoustics['c']['sta2']:.1f}  |  "
            f"sta3: ρ={acoustics['rho']['sta3']:.3f}, "
            f"c={acoustics['c']['sta3']:.1f}  |  "
            f"sta4: ρ={acoustics['rho']['sta4']:.3f}, "
            f"c={acoustics['c']['sta4']:.1f}")
    fig.text(0.5, 0.01, info, ha='center', fontsize=8.5, color='gray')

    fig.subplots_adjust(top=0.90, bottom=0.16, left=0.10, right=0.95)
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f'[saved] {save_path}')
    if show:
        plt.show()
    else:
        plt.close(fig)


# ════════════════════════════════════════════════════════════════
# 绘图: 多 n 能量分布对比
# ════════════════════════════════════════════════════════════════
# n=1~3 线型方案 (论文风格, 纯线型区分 n)
N_LINESTYLES = {1: '-', 2: '--', 3: '-.'}

# 段颜色: 上游蓝, 叶轮红, 扩压器绿
SEG_COLORS = {'Upstream': '#1f77b4', 'Impeller': '#d62728', 'Diffuser': '#2ca02c'}


def plot_energy_multi_n(all_shapes, acoustics, m_flow, rpm,
                         save_path=None, show=True):
    """绘制多个 n 模态的 Energy Function 叠加对比 (n=1~3)

    各段用颜色区分 (蓝=上游, 红=叶轮, 绿=扩压器), 各 n 用线型区分.
    """
    from matplotlib.lines import Line2D

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'

    # 计算各模态能量分布
    all_dists = [compute_energy_distribution(s, acoustics) for s in all_shapes]

    fig, ax = plt.subplots(figsize=(11, 6.5))

    for shape, dist in zip(all_shapes, all_dists):
        n_val = shape['n']
        if n_val not in N_LINESTYLES:
            continue   # 只画 n=1~3
        ls = N_LINESTYLES[n_val]

        x_list, E_list, bounds = dist['x'], dist['E'], dist['bounds']

        # 归一化: 除以本模态入口值 (保持符号)
        E_inlet = E_list[0][0]
        scale = E_inlet if E_inlet != 0 else 1.0

        # 分别绘制各段 (避免拼接导致的交叉)
        for x_seg, E_seg, (_, _, name) in zip(x_list, E_list, bounds):
            E_plot = E_seg / scale
            color = SEG_COLORS.get(name, 'gray')

            if name == 'Impeller':
                # 叶轮: 垂直阶跃 at x=0
                ax.plot([0.0, 0.0], [E_plot[0], E_plot[1]],
                        ls, color=color, linewidth=1.5)
            else:
                ax.plot(x_seg, E_plot, ls, color=color, linewidth=1.2)

    # 叶轮分界线
    ax.axvline(x=0, color='gray', linestyle=':', linewidth=1.0, alpha=0.6)

    # 零基准线
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

    # 段标注 (顶部)
    seg_labels = {'Upstream': 'Inlet duct', 'Impeller': 'Impeller', 'Diffuser': 'Vaneless Diffuser'}
    ax.annotate('Inlet duct', xy=(-0.1, 0.97), xycoords=('data', 'axes fraction'),
                ha='center', va='top', fontsize=10, style='italic', color=SEG_COLORS['Upstream'])
    ax.annotate('Impeller', xy=(0.0, 0.97), xycoords=('data', 'axes fraction'),
                ha='center', va='top', fontsize=10, style='italic', color=SEG_COLORS['Impeller'])
    ax.annotate('Vaneless Diffuser', xy=(0.1, 0.97), xycoords=('data', 'axes fraction'),
                ha='center', va='top', fontsize=10, style='italic', color=SEG_COLORS['Diffuser'])

    # 叶轮处归一化半径标注 (黑色, 带箭头, 在 x 轴刻度与标题之间)
    # y=-0.095: 位于刻度标签下方, xlabel (labelpad=22) 上方, 避免重叠
    # arrowstyle='<-': 箭头尖在 xytext(文本处), 尾部在 xy(中间=叶轮),
    #                  即箭头从中间指向两边 (标注两侧归一化半径)
    r_norm = acoustics.get('r_norm', {'sta2': 0.52, 'sta3': 1.0, 'sta4': 1.69})
    ax.annotate(r'$r_2/R_2 = ' + f'{r_norm["sta2"]:.2f}$',
                xy=(0.0, -0.095), xycoords=('data', 'axes fraction'),
                xytext=(-0.12, -0.095), textcoords=('data', 'axes fraction'),
                ha='right', va='center', fontsize=9.5, color='black',
                arrowprops=dict(arrowstyle='<-', color='black', lw=1.0))
    ax.annotate(r'$r_3/R_2 = ' + f'{r_norm["sta3"]:.2f}$',
                xy=(0.0, -0.095), xycoords=('data', 'axes fraction'),
                xytext=(0.12, -0.095), textcoords=('data', 'axes fraction'),
                ha='left', va='center', fontsize=9.5, color='black',
                arrowprops=dict(arrowstyle='<-', color='black', lw=1.0))

    # labelpad=22: 把 xlabel 下推到 r 标注下方, 避免与 r 标注重叠
    ax.set_xlabel(r'Streamwise coordinate  $x/R_2$', fontsize=13, labelpad=22)
    ax.set_ylabel(r'Energy function  $\delta \varepsilon_n / \delta \varepsilon_{\mathrm{inlet}}$',
                  fontsize=13)

    # 限制 x 轴显示范围 [-0.15, +0.15]  (仅视觉裁剪, 不改动物理空间/边界条件)
    # ----------------------------------------------------------------
    # 物理模型仍为: 无限长上游管道 + 零厚度叶轮(x=0) + 扩压器(s∈[0,0.69])
    # 这里 set_xlim 只裁剪显示窗口:
    #   - 上游: 显示 x∈[-0.15, 0] (实际已由 compute_modal_shape 的 x_upstream_max=0.15
    #           控制采样点, 这里只是与之对齐)
    #   - 扩压器: 物理算到 x=+0.69, 这里只显示前 ~12% (r∈[1.0, ~1.12])
    # 不影响特征值、SVD 提取、边界条件 (上游 B_n=C_n=0, 下游气室).
    # ----------------------------------------------------------------
    ax.set_xlim(-0.15, 0.15)

    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(MultipleLocator(0.01))
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.tick_params(axis='both', which='minor', length=4)

    # 自定义图例: 线型区分 n (黑色代理线)
    legend_elements = [Line2D([0], [0], color='black', linestyle=N_LINESTYLES[n],
                              linewidth=1.2, label=f'n={n}')
                       for n in sorted(N_LINESTYLES.keys())]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=11,
              framealpha=1.0, facecolor='white', edgecolor='black')
    ax.grid(True, which='both', alpha=0.25)
    ax.grid(True, which='minor', alpha=0.12)

    # 用 subplots_adjust 替代 tight_layout, 保证底部有足够空间容纳
    # 刻度 + r 标注 + xlabel 三层 (tight_layout 会把 xlabel 拉到 r 标注位置导致重叠)
    fig.subplots_adjust(bottom=0.20, top=0.95, left=0.10, right=0.97)
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f'[saved] {save_path}')
    if show:
        plt.show()
    else:
        plt.close(fig)


# ════════════════════════════════════════════════════════════════
# 辅助
# ════════════════════════════════════════════════════════════════
def find_closest_result(all_results, target_m):
    """在 all_results 中找最接近 target_m 的工况"""
    best = None
    best_diff = float('inf')
    for res in all_results:
        diff = abs(res['m'] - target_m)
        if diff < best_diff:
            best_diff = diff
            best = res
    return best, best_diff


def select_by_growth_rate(all_results, rank=2):
    """选取增长率第 rank 大的工况点 (rank=1 最大, rank=2 第二大)

    按工况 (m) 分组, 每个工况取跨所有 n 的最大 σ, 按 σ 降序选第 rank 个.
    返回 (res, max_sigma) 或 None.
    """
    m_to_res = {}
    m_to_max_sigma = {}
    for res in all_results:
        m = res['m']
        max_sig = max((eig.real for eigs in res['eigenvalues'].values()
                       for eig in eigs), default=-999.0)
        m_to_res[m] = res
        m_to_max_sigma[m] = max_sig
    sorted_m = sorted(m_to_max_sigma.items(), key=lambda x: x[1], reverse=True)
    if len(sorted_m) >= rank:
        m, sig = sorted_m[rank - 1]
        return m_to_res[m], sig
    return None


def compute_acoustics_and_params(rpm, m_flow, record, base_params):
    """计算指定工况的声学参数和模型参数

    返回 (acoustics, params)
    """
    geo_data = compute_base_flow(record, rpm, m_flow)
    from Prediction.Eckardt_O.Compute_Eckardt_O_Stability_Margin import (
        update_params_from_geo
    )
    params = update_params_from_geo(dict(base_params), geo_data)
    acoustics = compute_station_acoustics(geo_data)
    return acoustics, params


# ════════════════════════════════════════════════════════════════
# 主入口
# ════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Eckardt O 模态能量空间分布绘图 (Energy Function)'
    )
    parser.add_argument('--rpm', type=float, default=None,
                        help='指定单转速 [rpm] (默认遍历全部 5 个转速)')
    parser.add_argument('--n-max', type=int, default=3,
                        help='多 n 对比时的最大 n (默认 3)')
    parser.add_argument('--no-show', action='store_true',
                        help='不弹窗显示, 仅保存')
    args = parser.parse_args()

    # 输出目录: e:\Data\Simulation\Spa\
    save_dir = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))

    # 预加载 YAML record 和基础参数 (所有转速共用)
    records = load_yaml(YAML_DATA)
    record = find_record(records, "Eckardt O")
    base_params = load_params(PARAMS_FILE)

    # 转速列表
    all_rpms = [10000, 12000, 14000, 16000, 18000]
    target_rpms = [int(args.rpm)] if args.rpm else all_rpms

    n_range = list(range(1, args.n_max + 1))
    total_plotted = 0

    for rpm in target_rpms:
        csv_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "Prediction", "Eckardt_O",
            f"Eckardt_O_Eigenvalues_{rpm}rpm.csv"
        )
        if not os.path.exists(csv_path):
            print(f"[skip] CSV not found: {csv_path}")
            continue

        print(f"\n{'='*60}")
        print(f"  {rpm} rpm")
        print(f"{'='*60}")
        all_results, _ = load_results(csv_path)

        # 选取增长率最大的工况点 (rank=1)
        selected = select_by_growth_rate(all_results, rank=1)
        if selected is None:
            print(f"  [skip] 工况数不足, 无法选取最大增长率工况")
            continue
        target_res, max_sigma = selected
        m_flow = target_res['m']
        print(f"  选取最大增长率工况: m={m_flow:.3f} kg/s, "
              f"PR={target_res['PR']:.4f}, max σ={max_sigma:+.5f}")

        # 计算基流参数
        acoustics, params = compute_acoustics_and_params(
            rpm, m_flow, record, base_params)
        # 添加归一化半径用于标注
        acoustics['r_norm'] = {
            'sta2': params['STA2_radial'],
            'sta3': params['STA3_radial'],
            'sta4': params['STA4_radial'],
        }

        # 计算 n=1~3 模态形状
        all_shapes = []
        for n_val in n_range:
            eigs = target_res['eigenvalues'].get(n_val, [])
            if not eigs:
                continue
            s_star = eigs[0]
            # x_upstream_max=0.15: 仅决定上游管道采样的物理坐标网格 x∈[-0.15, 0]
            # (np.linspace(-0.15, 0, 60)), 类似在无限长管道里只在 [-0.15,0] 段
            # 放置 60 个探针. SVD/特征方程/边界条件全部不变, 上游仍按无限长
            # 管道建模 (force_upstream_bc=True → B_n=C_n=0).
            shape = compute_modal_shape(s_star, n_val, params,
                                        n_points=60, x_upstream_max=0.15)
            all_shapes.append(shape)

        if not all_shapes:
            print(f"    无可用模态, 跳过")
            continue

        save_path = os.path.join(
            save_dir,
            f"Eckardt_O_EnergyFunction_{rpm}rpm_m{m_flow/M_DESIGN:.3f}.png")
        plot_energy_multi_n(all_shapes, acoustics, m_flow, rpm,
                            save_path=save_path, show=False)
        total_plotted += 1

    print(f"\nDone. 共生成 {total_plotted} 张图.")
