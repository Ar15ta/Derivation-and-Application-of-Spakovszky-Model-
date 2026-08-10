"""Plot_Eckardt_O_Eigenvalues.py - Eckardt O 特征值绘图

从实验 CSV 读转速和流量范围，逐转速从中间流量向下扫直到找到不稳定点，
绘制模态迁移图。

用法:
  python Plot_Eckardt_O_Eigenvalues.py                          (遍历全部转速)
  python Plot_Eckardt_O_Eigenvalues.py --rpm 14000              (单转速)
  python Plot_Eckardt_O_Eigenvalues.py --npz path/to/file.csv   (从已有 CSV 绘制)
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import csv

# 将 Stability 根目录加入 sys.path
STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if STABILITY_ROOT not in sys.path:
    sys.path.insert(0, STABILITY_ROOT)

from Prediction.Eckardt_O.Compute_Eckardt_O_Stability_Margin import (
    load_results, track_modes, sweep_eigenvalues,
    save_results, load_yaml, find_record
)

# 均值线包路径
MEANLINE_PKG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "MeanLIne_Design", "RadCompressor_MeanLine", "radcomp-main"
)
YAML_DATA = os.path.join(MEANLINE_PKG, "data", "known_compressors.yml")

# 实验曲线路径
EXP_CSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "MeanLIne_Design", "RadCompressor_MeanLine", "Validation", "Eckardt_Exp_Curve.csv"
)


def parse_exp_curve(csv_path):
    """从 Eckardt_Exp_Curve.csv 解析各转速的 (m, PR) 数据

    返回 {rpm: [(m, pr), ...]}
    """
    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        lines = list(csv.reader(f))

    if len(lines) < 3:
        return {}

    # 第一行: "18000 RPM,,,16000 RPM,,,..."
    header = lines[0]
    rpm_cols = {}  # rpm -> start_col_index
    for i, cell in enumerate(header):
        cell = cell.strip()
        if cell.endswith('RPM'):
            rpm = int(cell.split()[0])
            if rpm not in rpm_cols:
                rpm_cols[rpm] = i

    # 数据行从第3行开始 (第2行是 Mass_Flow,Pressure_Ratio 子标题)
    data = {}
    for rpm in sorted(rpm_cols.keys()):
        data[rpm] = []

    for line in lines[2:]:
        for rpm, col_start in rpm_cols.items():
            m_str = line[col_start].strip() if col_start < len(line) else ''
            pr_str = line[col_start + 1].strip() if col_start + 1 < len(line) else ''
            if m_str and pr_str:
                try:
                    m_val = float(m_str)
                    pr_val = float(pr_str)
                    data[rpm].append((m_val, pr_val))
                except ValueError:
                    continue

    return data


def compute_sweep_below_mid(record, rpm, m_start, m_design,
                            step_norm=0.05, n_max=5, max_steps=15):
    """从 m_start 开始按归一化步长向下扫，直到出现不稳定根或步数耗尽

    返回 (all_results, mid_idx, boundary_msg)
    """
    all_results = []
    m_vals = [m_start]
    for i in range(max_steps):
        m_next = m_vals[-1] - step_norm * m_design
        if m_next <= 0:
            break
        m_vals.append(m_next)

    results = sweep_eigenvalues(record, rpm, m_vals, n_max=n_max)
    all_results.extend(results)

    # 检查是否出现不稳定根
    boundary_msg = None
    sorted_by_m = sorted(results, key=lambda r: r['m'], reverse=True)
    for res in sorted_by_m:
        for n_val, eigs in res['eigenvalues'].items():
            for eig in eigs:
                if eig.real > 0:
                    boundary_msg = (
                        f"{int(rpm)} rpm: 不稳定 m = {res['m']:.3f} kg/s "
                        f"(m/m_d = {res['m']/m_design:.3f}), "
                        f"n={n_val}, σ={eig.real:+.5f}"
                    )
                    break
            if boundary_msg:
                break
        if boundary_msg:
            break

    return all_results, boundary_msg


# ──────────────────────────────
# 绘图
# ──────────────────────────────
def plot_eigenvalue_map(all_results, rpm, n_range=None, paper_convention=True,
                        save_path=None, m_design=5.31, show=True):
    """绘制多个 n 的特征值分布图 (多工况叠加 + 模态追踪虚线)

    n_range: 要绘制的周向波数列表, 默认 [1, 2, 3, 4, 5]
    每个 n 用不同 marker 区分, 颜色仍按归一化流量映射。
    """
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'

    if n_range is None:
        n_range = [1, 2, 3, 4, 5]

    fig, ax = plt.subplots(figsize=(9, 7))

    markers_mode = ['o', 's', '^', 'D', 'v', 'p']

    m_all = [r['m'] for r in all_results]
    m_all_norm = [m / m_design for m in m_all]
    cmap = plt.get_cmap('viridis')
    if len(m_all) > 1:
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
                else:
                    sigmas.append(eig.real)
                    omegas.append(-eig.imag if paper_convention else eig.imag)
                    all_points.append((m_val, eig.real, -eig.imag if paper_convention else eig.imag))

            ax.plot(sigmas, omegas, '--', color='black', linewidth=1.5, alpha=0.6)

            for m_val, eig in track:
                if eig is None:
                    continue
                sigma = eig.real
                omega = -eig.imag if paper_convention else eig.imag
                m_norm = m_val / m_design
                ax.scatter(sigma, omega,
                           s=130, color=cmap(norm(m_norm)), marker=mkr,
                           edgecolors='black', linewidths=0.8, zorder=3)

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

    # 失稳边界点标记: 取|σ|最接近0的两个点, 按σ排序, 左侧点标签置上方, 右侧点标签置下方
    ylim_final = ax.get_ylim()
    y_span = ylim_final[1] - ylim_final[0]
    if len(all_points) >= 2:
        sorted_by_abs_sigma = sorted(all_points, key=lambda p: abs(p[1]))
        two_pts = sorted_by_abs_sigma[:2]
        two_pts_by_sigma = sorted(two_pts, key=lambda p: p[1])
        left_pt, right_pt = two_pts_by_sigma[0], two_pts_by_sigma[1]
        for pt, place in [(left_pt, 'above'), (right_pt, 'below')]:
            m_pt, sigma_pt, omega_pt = pt
            m_norm_pt = m_pt / m_design
            if place == 'above':
                text_y = omega_pt + 0.06 * y_span
                va = 'bottom'
            else:
                text_y = omega_pt - 0.06 * y_span
                va = 'top'
            ax.annotate(rf'$\dot{{m}}$={m_norm_pt:.2f}',
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

    from matplotlib.ticker import MultipleLocator

    # 固定刻度间隔: x轴 0.1, y轴 0.5
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
    leg = ax.legend(handles=legend_handles, loc='lower left', fontsize=14,
                    facecolor='white', edgecolor='black', framealpha=1.0)

    ax.grid(True, alpha=0.3)

    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label('Normalized Mass Flow', fontsize=14, fontfamily='Times New Roman')
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
        description='Eckardt O 特征值绘图 (从实验 CSV 读转速范围, 自动扫流量)'
    )
    parser.add_argument('--npz', type=str, default=None,
                        help='已有计算结果 CSV 文件路径 (跳过计算直接绘图)')
    parser.add_argument('--rpm', type=float, default=None,
                        help='指定单转速 (默认遍历全部转速)')
    parser.add_argument('--n-max', type=int, default=5,
                        help='最大周向波数 n (绘制 n=1..n_max)')
    parser.add_argument('--paper/--no-paper', dest='paper', default=True,
                        action='store_true',
                        help='使用论文约定 ω=-Im(s)')
    parser.add_argument('--save', type=str, default=None,
                        help='保存图片路径 (单文件模式)')
    parser.add_argument('--m-design', type=float, default=5.31,
                        help='设计点流量 (kg/s)')
    args = parser.parse_args()

    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # ── 模式 1: 从已有 CSV 直接绘图 ──
    if args.npz is not None:
        print(f"Loading: {args.npz}")
        all_results, rpm = load_results(args.npz)
        n_range = list(range(1, args.n_max + 1))
        save_path = args.save or os.path.join(
            project_root,
            f"Eckardt_O_Eigenvalues_n1-{args.n_max}_{int(rpm)}rpm.png")
        plot_eigenvalue_map(all_results, rpm, n_range=n_range,
                            paper_convention=args.paper,
                            save_path=save_path, m_design=args.m_design)
        sys.exit(0)

    # ── 模式 2: 从实验 CSV 读转速, 自动计算 ──
    if not os.path.exists(EXP_CSV):
        print(f"Experiment CSV not found: {EXP_CSV}")
        sys.exit(1)

    exp_data = parse_exp_curve(EXP_CSV)
    if not exp_data:
        print("Failed to parse experiment CSV")
        sys.exit(1)

    # 过滤转速
    target_rpms = [args.rpm] if args.rpm else sorted(exp_data.keys())
    target_rpms = [r for r in target_rpms if r in exp_data]
    if not target_rpms:
        print(f"No matching RPM data found")
        sys.exit(1)

    print(f"Target RPMs: {target_rpms}")
    print(f"Design mass flow (m_design): {args.m_design} kg/s")
    print(f"Step: 0.05 × m_design = {0.05 * args.m_design:.4f} kg/s\n")

    # 加载 YAML
    records = load_yaml(YAML_DATA)
    record = find_record(records, "Eckardt O")

    for rpm in target_rpms:
        points = exp_data[rpm]
        m_flows = [p[0] for p in points]
        m_mid = (min(m_flows) + max(m_flows)) / 2
        m_mid_norm = m_mid / args.m_design

        # 限制扫描起点: 以 14000rpm 为基准 (阈值 1.0), 按转速比例调整
        # 高转速失稳点流量更高 → 阈值拉高; 低转速 → 阈值降低
        start_threshold = rpm / 14000.0
        if m_mid_norm > start_threshold:
            m_start = start_threshold * args.m_design
            start_note = (f" (中点 {m_mid_norm:.3f} > {start_threshold:.3f}, "
                          f"起点限制为 {start_threshold:.3f})")
        else:
            m_start = m_mid
            start_note = ""

        print(f"\n{'='*60}")
        print(f"  {int(rpm)} rpm: 实验流量 {min(m_flows):.2f} ~ {max(m_flows):.2f} kg/s, "
              f"起点 m = {m_start:.3f} (m/m_d = {m_start/args.m_design:.3f}){start_note}")
        print(f"{'='*60}")

        all_results, boundary_msg = compute_sweep_below_mid(
            record, rpm, m_start, args.m_design,
            step_norm=0.05, n_max=5, max_steps=15)

        if boundary_msg:
            print(f"  [!] {boundary_msg}")

        # 保存 CSV
        pred_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "Prediction", "Eckardt_O"
        )
        os.makedirs(pred_dir, exist_ok=True)
        csv_path = os.path.join(pred_dir, f"Eckardt_O_Eigenvalues_{int(rpm)}rpm.csv")
        save_results(all_results, rpm, csv_path)

        # 绘图
        n_range = list(range(1, args.n_max + 1))
        save_path = os.path.join(
            project_root,
            f"Eckardt_O_Eigenvalues_n1-{args.n_max}_{int(rpm)}rpm.png")
        plot_eigenvalue_map(all_results, rpm, n_range=n_range,
                            paper_convention=args.paper,
                            save_path=save_path, m_design=args.m_design,
                            show=False)

    print("\nDone.")
