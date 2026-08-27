"""Plot_NCEPU_Stability_Margin.py - NCEPU 特征值绘图

从 Compute_NCEPU_Margin.py 输出的 CSV 读取特征值结果, 绘制模态迁移图
(特征值随流量变化轨迹) 和稳定性裕度图。

用法:
  python Plot_NCEPU_Stability_Margain.py --csv ../Prediction/NCEPU/NCEPU_Eigenvalues_14000rpm.csv
  python Plot_NCEPU_Stability_Margain.py --csv result.csv --m-design 5.0
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import glob

# 将 Stability 根目录加入 sys.path
STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if STABILITY_ROOT not in sys.path:
    sys.path.insert(0, STABILITY_ROOT)

from Prediction.NCEPU.Compute_NCEPU_Margin import load_results, track_modes


# 设计点质量流量 [kg/s] (设计转速下), 用于流量归一化; 可由 --m-design 覆盖
M_DESIGN = 13.5


# ──────────────────────────────
# 绘图
# ──────────────────────────────
def plot_eigenvalue_map(all_results, rpm, n_range=None, paper_convention=True,
                        save_path=None, m_design=M_DESIGN, show=True):
    """绘制多个 n 的特征值分布图 (多工况叠加 + 模态追踪虚线)

    参数
    ----
    all_results : list of dict
        load_results 返回的特征值结果
    rpm : float
        当前转速 (仅用于标题/标注)
    n_range : list[int] or None
        要绘制的周向波数列表, 默认 [1, 2, 3, 4, 5]
    paper_convention : bool
        True 时 ω = -Im(s) (论文约定)
    m_design : float
        设计点流量 [kg/s], 用于颜色归一化
    """
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'

    if n_range is None:
        n_range = [1, 2, 3, 4, 5]

    fig, ax = plt.subplots(figsize=(9, 7))

    markers_mode = ['o', 's', '^', 'D', 'v', 'p']

    m_all_norm = [r['m'] / m_design for r in all_results]
    cmap = plt.get_cmap('viridis')
    if len(m_all_norm) > 1:
        norm = plt.Normalize(vmin=min(m_all_norm), vmax=max(m_all_norm))
    else:
        norm = plt.Normalize(vmin=m_all_norm[0] - 0.1, vmax=m_all_norm[0] + 0.1)

    all_points = []
    for n_idx, n_val in enumerate(n_range):
        mode_tracks = track_modes(all_results, n_val, max_modes=1)
        mkr = markers_mode[n_idx % len(markers_mode)]

        for track in mode_tracks:
            sigmas, omegas = [], []
            for m_val, eig in track:
                if eig is None:
                    sigmas.append(np.nan)
                    omegas.append(np.nan)
                    continue
                sigma = eig.real
                omega = -eig.imag if paper_convention else eig.imag
                sigmas.append(sigma)
                omegas.append(omega)
                all_points.append((m_val, sigma, omega))
                ax.scatter(sigma, omega, s=130,
                           color=cmap(norm(m_val / m_design)), marker=mkr,
                           edgecolors='black', linewidths=0.8, zorder=3)

            ax.plot(sigmas, omegas, '--', color='black', linewidth=1.5, alpha=0.6)

    # 稳定性边界
    ax.axvline(x=0.0, color='red', linestyle='--', linewidth=1.8, zorder=2)
    ax.axhline(y=0.0, color='gray', linestyle=':', linewidth=1.0, zorder=1)

    # 不稳定区域着色
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.fill_betweenx([ylim[0] - 100, ylim[1] + 100],
                     0, xlim[1] + 100,
                     color='red', alpha=0.05, zorder=0)
    ax.set_xlim(xlim[0] - abs(xlim[0]) * 0.08,
                xlim[1] + abs(xlim[1]) * 0.08)
    ax.set_ylim(ylim[0] - abs(ylim[0]) * 0.08,
                ylim[1] + abs(ylim[1]) * 0.08)

    # 失稳边界点标记: 稳定侧取 σ<0 中最大 σ (失稳前最后一个稳定点),
    # 不稳定侧取 σ>0 中最小 σ (首个失稳点), 二者必分居 σ=0 两侧
    ylim_final = ax.get_ylim()
    y_span = ylim_final[1] - ylim_final[0]
    stable_pts = [p for p in all_points if p[1] < 0.0]
    unstable_pts = [p for p in all_points if p[1] > 0.0]
    if stable_pts and unstable_pts:
        left_pt = max(stable_pts, key=lambda p: p[1])    # 最贴近 0 的稳定点
        right_pt = min(unstable_pts, key=lambda p: p[1])  # 最贴近 0 的失稳点
        for pt, place in [(left_pt, 'above'), (right_pt, 'below')]:
            m_pt, sigma_pt, omega_pt = pt
            m_norm_pt = m_pt / m_design
            if place == 'above':
                text_y = omega_pt + 0.06 * y_span
                va = 'bottom'
            else:
                text_y = omega_pt - 0.06 * y_span
                va = 'top'
            ax.annotate(rf'$\dot{{m}}$={m_norm_pt:.4f}',
                        xy=(sigma_pt, omega_pt),
                        xytext=(sigma_pt, text_y),
                        fontsize=11, fontweight='bold',
                        ha='center', va=va,
                        bbox=dict(facecolor='white', edgecolor='black',
                                  boxstyle='round,pad=0.3'),
                        arrowprops=dict(arrowstyle='->', color='black',
                                        linewidth=1.5, linestyle='-',
                                        shrinkA=2, shrinkB=4))

    omega_label = (r'$\omega$(Rotating Speed)'
                   if paper_convention else r'Im($s$)')
    ax.set_xlabel(r'$\sigma$ (Growth Rate)', fontsize=14)
    ax.set_ylabel(omega_label, fontsize=14)
    ax.set_title(f'NCEPU Eigenvalue Map  -  {int(rpm)} rpm', fontsize=14)

    from matplotlib.ticker import MultipleLocator

    x_step = 0.1
    y_step = 0.5
    ax.xaxis.set_major_locator(MultipleLocator(x_step))
    ax.xaxis.set_minor_locator(MultipleLocator(x_step / 5))
    ax.yaxis.set_major_locator(MultipleLocator(y_step))
    ax.yaxis.set_minor_locator(MultipleLocator(y_step / 5))
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', length=4)

    legend_handles = []
    for n_idx, n_val in enumerate(n_range):
        mkr = markers_mode[n_idx % len(markers_mode)]
        legend_handles.append(
            plt.Line2D([0], [0], marker=mkr, color='w',
                       markerfacecolor='white', markersize=14,
                       markeredgecolor='black', markeredgewidth=0.8,
                       label=f'n={n_val}')
        )
    ax.legend(handles=legend_handles, loc='lower left', fontsize=14,
              facecolor='white', edgecolor='black', framealpha=1.0)

    ax.grid(True, alpha=0.3)

    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r'$\dot{m}/\dot{m}_{\mathrm{design}}$', fontsize=14,
                   fontfamily='Times New Roman')
    cbar.ax.tick_params(labelsize=12)
    for label in cbar.ax.get_yticklabels():
        label.set_fontfamily('Times New Roman')

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f'[saved] {save_path}')
    if show:
        plt.show()
    else:
        plt.close(fig)


# ──────────────────────────────
# 主入口
# ──────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='NCEPU 特征值绘图 (从预计算 CSV 读取)'
    )
    parser.add_argument('--csv', type=str, default=None,
                        help='特征值结果 CSV 文件路径 (Compute_NCEPU_Margin 输出); '
                             '缺省时自动取 Prediction/NCEPU 下最新的 NCEPU_Eigenvalues_*.csv')
    parser.add_argument('--n-max', type=int, default=4,
                        help='最大周向波数 n (绘制 n=1..n_max)')
    parser.add_argument('--paper/--no-paper', dest='paper', default=True,
                        action='store_true',
                        help='使用论文约定 ω=-Im(s)')
    parser.add_argument('--save', type=str, default=None,
                        help='保存图片路径')
    parser.add_argument('--m-design', type=float, default=M_DESIGN,
                        help=f'设计点流量 [kg/s] (默认 {M_DESIGN})')
    parser.add_argument('--no-show', action='store_true',
                        help='不弹窗显示, 仅保存')
    args = parser.parse_args()

    csv_path = args.csv
    if csv_path is None:
        eig_dir = os.path.join(STABILITY_ROOT, 'Prediction', 'NCEPU')
        candidates = glob.glob(os.path.join(eig_dir, 'NCEPU_Eigenvalues_*.csv'))
        if not candidates:
            print(f"未找到特征值结果 CSV, 请先运行 Compute_NCEPU_Margin.py "
                  f"(查找目录: {eig_dir})")
            sys.exit(1)
        csv_path = max(candidates, key=os.path.getmtime)

    if not os.path.exists(csv_path):
        print(f"CSV not found: {csv_path}")
        sys.exit(1)

    print(f"Loading: {csv_path}")
    all_results = load_results(csv_path)
    if not all_results:
        print("CSV 中没有有效结果")
        sys.exit(1)

    rpms = sorted(set(r['rpm'] for r in all_results))
    rpm = rpms[0]
    n_range = list(range(1, args.n_max + 1))

    # 归一化设计流量 (默认设计点 13.5 kg/s, 可用 --m-design 覆盖)
    m_design = args.m_design

    save_path = args.save
    if save_path is None:
        project_root = os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))))
        save_path = os.path.join(
            project_root,
            f"NCEPU_Eigenvalues_n1-{args.n_max}_{int(rpm)}rpm.png")

    plot_eigenvalue_map(all_results, rpm, n_range=n_range,
                        paper_convention=args.paper,
                        save_path=save_path, m_design=m_design,
                        show=not args.no_show)
