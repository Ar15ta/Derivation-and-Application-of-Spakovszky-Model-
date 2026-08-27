"""Plot_NCEPU_Distribution.py - NCEPU 模态能量空间分布绘图

从 Compute_NCEPU_Margin.py 输出的特征值 CSV 读取结果, 对指定工况的各 n
模态计算其在物理空间的状态向量分布, 然后绘制 Energy Function (声学扰动
能量密度) 沿物理空间的分布曲线, 复现 Spakovszky 论文 Figure 5-18 风格。

──────────────────────────────────────────────────────────────────
物理空间分段 (归一化坐标, R_3 = 叶轮出口半径):
  - 上游管道: x ∈ (-∞, 0]   (无限长管道, x=0 = 叶轮进口 STA2, 负值)
  - 叶轮:     x = 0           (执行盘, 零厚度, 能量跳变)
  - 扩压器:   s ∈ [0, s4]     (径向弧长 = r - r3, 正值, r3=1.0)
──────────────────────────────────────────────────────────────────

Energy Function (声学扰动能量流密度, 论文 Eq.5-23):
  δε_n(x, s) = P̃_n(x, s) · Ṽ_n(x, s)
  取实部保留能量方向; 归一化: 除以本模态上游入口 δε_inlet.

用法:
  python Plot_NCEPU_Distribution.py --eig ../Prediction/NCEPU/NCEPU_Eigenvalues_14000rpm.csv
  python Plot_NCEPU_Distribution.py --eig result.csv --m 5.0 --n-max 3
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

from Model.Model_NCEPU import compute_modal_shape
from Geo_and_Base_Flow.NCEPU_CFX_Loader import (
    load_all_operating_points, DEFAULT_CSV_FILE, DEFAULT_PARAMS_FILE
)
from Prediction.NCEPU.Compute_NCEPU_Margin import load_results

# 设计点参考质量流量 [kg/s] (设计转速下, 归一化/文件名用, 可由命令行覆盖)
M_DESIGN = 13.5


# ════════════════════════════════════════════════════════════════
# 能量函数 (论文定义: δε_n = P̃_n · Ṽ_n)
# ════════════════════════════════════════════════════════════════
def energy_function(q_array):
    """计算声学扰动能量函数 δε_n = P̃_n · Ṽ_n (论文 Eq.5-23)

    输入:
      q_array : (N, 3) complex, 归一化状态向量 [δV_stream, δV_θ, δP]
    返回:
      δε : (N,) float, 能量函数的实部 (保留正负号 → 反映能量流向)
    """
    return (q_array[:, 2] * q_array[:, 0]).real


# ════════════════════════════════════════════════════════════════
# 模态能量分布计算
# ════════════════════════════════════════════════════════════════
def compute_energy_distribution(shape_result, acoustics):
    """把 compute_modal_shape 的结果转换为能量分布 (按物理空间三段切分)

    返回 dict: x, E, bounds, seg_meta
    """
    u3 = acoustics['u3']
    segments = shape_result['segments']

    # 上游基流声速/密度: 用 STA2 (叶轮进口) 近似
    rho_up = acoustics['rho']['sta2']
    a_up = acoustics['a']['sta2']
    rho2, a2 = acoustics['rho']['sta2'], acoustics['a']['sta2']
    rho3, a3 = acoustics['rho']['sta3'], acoustics['a']['sta3']
    rho4, a4 = acoustics['rho']['sta4'], acoustics['a']['sta4']

    x_list, E_list, bounds, seg_meta = [], [], [], []

    # ── 上游管道 ──
    up = segments[0]
    x_up = up['coord']
    q_up = up['q']
    E_up = energy_function(q_up)
    x_list.append(x_up)
    E_list.append(E_up)
    bounds.append((x_up[0], x_up[-1], 'Upstream'))
    seg_meta.append({'rho': rho_up, 'a': a_up, 'M_ref': u3 / a_up if a_up > 0 else 0.0})

    # ── 叶轮: x=0 处能量跳变 ──
    q_imp_in = segments[0]['q'][-1]
    q_imp_out = segments[1]['q'][0]
    E_imp_in = energy_function(q_imp_in.reshape(1, 3))[0]
    E_imp_out = energy_function(q_imp_out.reshape(1, 3))[0]

    x_list.append(np.array([0.0, 0.0]))
    E_list.append(np.array([E_imp_in, E_imp_out]))
    bounds.append((0.0, 0.0, 'Impeller'))
    seg_meta.append({'rho_in': rho2, 'a_in': a2, 'rho_out': rho3, 'a_out': a3})

    # ── 扩压器: 密度声速线性插值 ──
    diff = segments[2]
    r_diff = diff['coord']
    q_diff = diff['q']
    s_diff = r_diff - r_diff[0]
    x_diff = s_diff

    t = (r_diff - r_diff[0]) / (r_diff[-1] - r_diff[0])
    rho_diff = rho3 + (rho4 - rho3) * t
    a_diff = a3 + (a4 - a3) * t

    E_diff = np.array([energy_function(q_diff[i:i + 1])[0]
                       for i in range(len(q_diff))])
    x_list.append(x_diff)
    E_list.append(E_diff)
    bounds.append((x_diff[0], x_diff[-1], 'Diffuser'))
    seg_meta.append({'rho': rho_diff, 'a': a_diff})

    return {'x': x_list, 'E': E_list, 'bounds': bounds, 'seg_meta': seg_meta}


# ════════════════════════════════════════════════════════════════
# 绘图: 多 n 能量分布对比
# ════════════════════════════════════════════════════════════════
N_LINESTYLES = {1: '-', 2: '--', 3: '-.'}
SEG_COLORS = {'Upstream': '#1f77b4', 'Impeller': '#d62728', 'Diffuser': '#2ca02c'}


def plot_energy_multi_n(all_shapes, acoustics, m_flow, rpm, r_norm,
                         save_path=None, show=True):
    """绘制多个 n 模态的 Energy Function 叠加对比 (n=1~3)"""
    from matplotlib.lines import Line2D

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'

    all_dists = [compute_energy_distribution(s, acoustics) for s in all_shapes]

    fig, ax = plt.subplots(figsize=(11, 6.5))

    for shape, dist in zip(all_shapes, all_dists):
        n_val = shape['n']
        if n_val not in N_LINESTYLES:
            continue
        ls = N_LINESTYLES[n_val]

        x_list, E_list, bounds = dist['x'], dist['E'], dist['bounds']
        E_inlet = E_list[0][0]
        scale = E_inlet if E_inlet != 0 else 1.0

        for x_seg, E_seg, (_, _, name) in zip(x_list, E_list, bounds):
            E_plot = E_seg / scale
            color = SEG_COLORS.get(name, 'gray')
            if name == 'Impeller':
                ax.plot([0.0, 0.0], [E_plot[0], E_plot[1]],
                        ls, color=color, linewidth=1.5)
            else:
                ax.plot(x_seg, E_plot, ls, color=color, linewidth=1.2)

    ax.axvline(x=0, color='gray', linestyle=':', linewidth=1.0, alpha=0.6)
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)

    # 段标注 (顶部)
    ax.annotate('Inlet duct', xy=(-0.1, 0.97), xycoords=('data', 'axes fraction'),
                ha='center', va='top', fontsize=10, style='italic',
                color=SEG_COLORS['Upstream'])
    ax.annotate('Impeller', xy=(0.0, 0.97), xycoords=('data', 'axes fraction'),
                ha='center', va='top', fontsize=10, style='italic',
                color=SEG_COLORS['Impeller'])
    ax.annotate('Vaneless Diffuser', xy=(0.1, 0.97),
                xycoords=('data', 'axes fraction'),
                ha='center', va='top', fontsize=10, style='italic',
                color=SEG_COLORS['Diffuser'])

    # 归一化半径标注
    r2 = r_norm.get('sta2', 0.52)
    r3 = r_norm.get('sta3', 1.0)
    ax.annotate(r'$r_2/R_2 = ' + f'{r2:.2f}$',
                xy=(0.0, -0.095), xycoords=('data', 'axes fraction'),
                xytext=(-0.12, -0.095), textcoords=('data', 'axes fraction'),
                ha='right', va='center', fontsize=9.5, color='black',
                arrowprops=dict(arrowstyle='<-', color='black', lw=1.0))
    ax.annotate(r'$r_3/R_2 = ' + f'{r3:.2f}$',
                xy=(0.0, -0.095), xycoords=('data', 'axes fraction'),
                xytext=(0.12, -0.095), textcoords=('data', 'axes fraction'),
                ha='left', va='center', fontsize=9.5, color='black',
                arrowprops=dict(arrowstyle='<-', color='black', lw=1.0))

    ax.set_xlabel(r'Streamwise coordinate  $x/R_2$', fontsize=13, labelpad=22)
    ax.set_ylabel(r'Energy function  $\delta \varepsilon_n / \delta \varepsilon_{\mathrm{inlet}}$',
                  fontsize=13)

    title = (f'Energy function distribution  |  {int(rpm)} rpm  '
             r'$\dot{m}$' + f'={m_flow:.3f} kg/s')
    ax.set_title(title, fontsize=11.5, pad=12)

    ax.set_xlim(-0.15, 0.15)
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(MultipleLocator(0.01))
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.tick_params(axis='both', which='minor', length=4)

    legend_elements = [Line2D([0], [0], color='black', linestyle=N_LINESTYLES[n],
                              linewidth=1.2, label=f'n={n}')
                       for n in sorted(N_LINESTYLES.keys())]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=11,
              framealpha=1.0, facecolor='white', edgecolor='black')
    ax.grid(True, which='both', alpha=0.25)
    ax.grid(True, which='minor', alpha=0.12)

    # 底部声学参数标注
    info = (f"sta2: ρ={acoustics['rho']['sta2']:.3f} kg/m³, "
            f"a={acoustics['a']['sta2']:.1f} m/s  |  "
            f"sta3: ρ={acoustics['rho']['sta3']:.3f}, "
            f"a={acoustics['a']['sta3']:.1f}  |  "
            f"sta4: ρ={acoustics['rho']['sta4']:.3f}, "
            f"a={acoustics['a']['sta4']:.1f}")
    fig.text(0.5, 0.01, info, ha='center', fontsize=8.5, color='gray')

    fig.subplots_adjust(bottom=0.20, top=0.92, left=0.10, right=0.97)
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
def find_operating_point(all_op, rpm_target, m_target, m_tol=0.01):
    """在 CFX 工况列表中找匹配 (rpm, m) 的工况"""
    for op in all_op:
        if abs(op['rpm'] - rpm_target) < 1.0 and abs(op['m'] - m_target) < m_tol:
            return op
    # 退而求其次: 同转速下最接近的流量
    same_rpm = [op for op in all_op if abs(op['rpm'] - rpm_target) < 1.0]
    if same_rpm:
        return min(same_rpm, key=lambda op: abs(op['m'] - m_target))
    return None


def select_by_growth_rate(all_results, rank=1):
    """选取增长率第 rank 大的工况点"""
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


# ════════════════════════════════════════════════════════════════
# 主入口
# ════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='NCEPU 模态能量空间分布绘图 (Energy Function)'
    )
    parser.add_argument('--eig', type=str, required=True,
                        help='特征值结果 CSV 路径 (Compute_NCEPU_Margin 输出)')
    parser.add_argument('--csv', type=str, default=DEFAULT_CSV_FILE,
                        help='CFX 气动参数 CSV 路径')
    parser.add_argument('--params', type=str, default=DEFAULT_PARAMS_FILE,
                        help='NCEPU_params.txt 路径')
    parser.add_argument('--m', type=float, default=None,
                        help='指定质量流量 [kg/s] (默认选增长率最大工况)')
    parser.add_argument('--rpm', type=float, default=None,
                        help='指定转速 [rpm] (默认从特征值 CSV 推断)')
    parser.add_argument('--n-max', type=int, default=3,
                        help='多 n 对比时的最大 n (默认 3)')
    parser.add_argument('--m-design', type=float, default=M_DESIGN,
                        help='设计点流量 [kg/s], 用于文件名归一化')
    parser.add_argument('--no-show', action='store_true',
                        help='不弹窗显示, 仅保存')
    args = parser.parse_args()

    if not os.path.exists(args.eig):
        print(f"Eigenvalue CSV not found: {args.eig}")
        sys.exit(1)
    if not os.path.exists(args.csv):
        print(f"CFX CSV not found: {args.csv}")
        sys.exit(1)

    # 加载 CFX 全部工况 (含 params 和声学量)
    print(f"Loading CFX CSV: {args.csv}")
    all_op = load_all_operating_points(args.csv, args.params)

    # 加载特征值结果
    print(f"Loading eigenvalues: {args.eig}")
    all_results = load_results(args.eig)

    rpms_in_eig = sorted(set(r['rpm'] for r in all_results))
    rpm_target = args.rpm if args.rpm is not None else rpms_in_eig[0]
    same_rpm = [r for r in all_results if abs(r['rpm'] - rpm_target) < 1.0]
    if not same_rpm:
        print(f"特征值 CSV 中没有转速 {rpm_target} 的数据")
        sys.exit(1)

    # 选定工况
    if args.m is not None:
        target_res = min(same_rpm, key=lambda r: abs(r['m'] - args.m))
    else:
        sel = select_by_growth_rate(same_rpm, rank=1)
        if sel is None:
            print("无法选取工况")
            sys.exit(1)
        target_res, max_sigma = sel

    m_flow = target_res['m']
    print(f"  工况: rpm={rpm_target:.0f} m={m_flow:.3f} kg/s, "
          f"π_tt={target_res['pi_tt']:.4f}")

    # 找到对应的 CFX 工况 params (含声学量)
    op = find_operating_point(all_op, rpm_target, m_flow)
    if op is None:
        print(f"CFX CSV 中没有匹配 rpm={rpm_target}, m={m_flow} 的工况")
        sys.exit(1)
    params = op['params']
    acoustics = params['acoustics']

    r_norm = {
        'sta2': params['radial_sta2'],
        'sta3': params['radial_sta3'],
        'sta4': params['radial_sta4'],
    }

    # 计算 n=1~n_max 模态形状
    n_range = list(range(1, args.n_max + 1))
    all_shapes = []
    for n_val in n_range:
        eigs = target_res['eigenvalues'].get(n_val, [])
        if not eigs:
            continue
        s_star = eigs[0]
        print(f"  n={n_val}: s* = {s_star.real:+.5f} {s_star.imag:+.5f}j")
        shape = compute_modal_shape(s_star, n_val, params,
                                    n_points=60, x_upstream_max=0.15)
        all_shapes.append(shape)

    if not all_shapes:
        print("无可用模态, 退出")
        sys.exit(1)

    save_dir = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    save_path = os.path.join(
        save_dir,
        f"NCEPU_EnergyFunction_{int(rpm_target)}rpm_m{m_flow/args.m_design:.3f}.png")
    plot_energy_multi_n(all_shapes, acoustics, m_flow, rpm_target, r_norm,
                        save_path=save_path, show=not args.no_show)
