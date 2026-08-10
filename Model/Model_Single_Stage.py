# Model_Stage.py - 单级压缩系统模型：上游无限管道 + 转子 + 轴向间隙 + 静子 + 下游无限管道
import os
import sys

# 确保 Stability 根目录在 sys.path 中, 使 Matrix/ 等包可导入
_STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _STABILITY_ROOT not in sys.path:
    sys.path.insert(0, _STABILITY_ROOT)

import numpy as np

# 导入组件模块（每个模块自己读取参数）
from Matrix.rotor import B_rot_n_num, compute_rotor_params, load_params
from Matrix.stator import B_sta_n_num, compute_stator_params
from Matrix.axial import Tn_ax_coeff_num, axial_duct_state_transfer
from Matrix.boundary import IC_upstream_infinite, EC_downstream_infinite

# 默认配置文件（轴流单级）
DEFAULT_PARAMS_FILE = 'system_params.txt'


def compute_stage_params(params):
    """统一计算转子与静子派生参数，确保串联接口速度自洽。

    返回的 params 同时包含：
      - Vtheta_bar1_rot, Vtheta_bar2_rot  （转子进/出口周向速度）
      - Vtheta_bar1_sta, Vtheta_bar2_sta  （静子进/出口周向速度）
    其中 Vtheta_bar2_rot 应与 Vtheta_bar1_sta 一致（间隙中周向速度连续）。
    """
    params = compute_rotor_params(params)
    params = compute_stator_params(params)
    return params


def check_stage_interface(params, tol=1e-6, raise_on_fail=False):
    """检查转子出口与静子进口在间隙中的周向速度是否一致。

    串联系统要求间隙内平均流连续：
        Vθ_bar2_rot == Vθ_bar1_sta
    若两者不一致，说明 beta2_deg 与 alpha3_deg 在参数文件中不自洽，
    转子出口流出条件与静子进口接收条件代表不同物理工作点。
    """
    Vth_rot_out = params['Vtheta_bar2_rot']
    Vth_sta_in = params['Vtheta_bar1_sta']
    residual = Vth_rot_out - Vth_sta_in

    if abs(residual) > tol:
        msg = (
            "\n[!] 转子-静子间隙速度不连续\n"
            f"    转子出口 Vθ_bar2_rot = {Vth_rot_out:.6f}\n"
            f"    静子进口 Vθ_bar1_sta = {Vth_sta_in:.6f}\n"
            f"    残差 = {residual:+.3e}\n"
            "    请检查 system_params.txt 中 beta2_deg 与 alpha3_deg 是否自洽。\n"
        )
        if raise_on_fail:
            raise ValueError(msg)
        else:
            print(msg)
        return False
    else:
        print(f"[OK] 间隙速度连续 (Vθ_gap = {Vth_rot_out:.6f})")
        return True


def build_system_matrix(s_val, n_val=1, params=None, params_file=None):
    """
    构建完整单级系统传递矩阵

    串联顺序（从上游到下游，坐标约定）：
        x ∈ (-∞, 0)：上游无限长管道
        x = 0          ：转子执行盘
        x ∈ [0, L_gap] ：轴向间隙管道
        x = L_gap      ：静子执行盘
        x ∈ (L_gap, +∞)：下游无限长管道

    系统矩阵公式（参考 边界条件文档 §5）：
        X_sys = T_dn(L_gap; x0=L_gap)^{-1} · B_sta · G_gap(L_gap) · B_rot · T_up(0; x0=0)

    其中：
      - T_up        : 上游管道系数 → 转子进口状态 （系数空间 → 状态空间）
      - B_rot       : 转子进口状态 → 转子出口状态 （状态 → 状态，执行盘）
      - G_gap       : 转子出口状态 → 静子进口状态 （状态 → 状态，长度 L_gap）
      - B_sta       : 静子进口状态 → 静子出口状态 （状态 → 状态，执行盘）
      - inv_T_dn    : 静子出口状态 → 下游管道系数 （状态空间 → 系数空间）

    关键数值考虑：下游管道的参考点 x0 取为静子出口 x = L_gap，
    这样 T_dn(L_gap; L_gap) 中 exp(±n·Δx) = 1，避免高阶 n 下
    exp(n·L_gap) 带来的矩阵病态。
    间隙段 G_gap 使用 axial_duct_state_transfer (状态↔状态)，
    避免在状态、系数空间之间反复转换带来的数值误差。

    参数：
    --------
    s_val : complex
        拉普拉斯变量值
    n_val : int
        周向波数
    params : dict or None
        预计算的参数字典, 为 None 时自动从文件加载
    params_file : str or None
        参数文件路径，为 None 时使用默认文件

    返回：
    --------
    X_sys : np.ndarray (3×3)
        完整系统传递矩阵（上游系数 c^up → 下游系数 c^dn）
    """
    if params is None:
        params = compute_stage_params(load_params(params_file or DEFAULT_PARAMS_FILE))

    Vx = params['Vx']
    L_gap = params['L_gap']
    mu_up = params.get('mu_up', 1.0)
    mu_dn = params.get('mu_dn', 1.0)

    # 单级设定：
    #   上游远场无预旋 → Vθ_up = 0
    #   间隙中周向速度 = 转子出口周向速度 = 静子进口周向速度
    #   下游远场周向速度 = 静子出口周向速度
    Vtheta_up = 0.0
    Vtheta_gap = params['Vtheta_bar2_rot']
    Vtheta_dn = params['Vtheta_bar2_sta']

    # 转子执行盘传递矩阵（状态 → 状态）
    B_rot = B_rot_n_num(s=s_val, n=n_val, params=params)

    # 轴向间隙管道传递矩阵（状态 → 状态，从间隙入口 x=0 到间隙出口 x=L_gap）
    G_gap = axial_duct_state_transfer(x=L_gap, s=s_val, n=n_val, Vx_bar=Vx, Vtheta_bar=Vtheta_gap, mu=mu_up)
    # 静子执行盘传递矩阵（状态 → 状态）
    B_sta = B_sta_n_num(s=s_val, n=n_val, params=params)

    # 上游管道：转子进口 x=0⁻ 处的状态 = T_up(0;0) · c^up
    T_up = Tn_ax_coeff_num(x=0, s=s_val, n=n_val,
                           Vx_bar=Vx, Vtheta_bar=Vtheta_up, x0=0, mu=mu_up)

    # 下游管道：静子出口位于 x = L_gap，下游管道参考点也取 x0 = L_gap
    # 这样 T_dn(L_gap; L_gap) 中 exp(±n·Δx) = exp(0) = 1，避免病态 exp(n·L_gap) 因子
    # 反解：c^dn = T_dn(L_gap; L_gap)^{-1} · state(L_gap⁺)
    T_dn_at_exit = Tn_ax_coeff_num(x=L_gap, s=s_val, n=n_val,
                                   Vx_bar=Vx, Vtheta_bar=Vtheta_dn,
                                   x0=L_gap, mu=mu_dn)
    inv_T_dn = np.linalg.inv(T_dn_at_exit)

    # 系统传递矩阵（从上游到下游依次相乘）
    #   x = -∞ → x = 0  ： T_up         上游管道系数 → 转子进口状态
    #   x = 0  → x = 0  ： B_rot        转子执行盘
    #   x = 0  → x = L_gap : G_gap        间隙管道状态传递
    #   x = L_gap → x = L_gap : B_sta     静子执行盘
    #   x = L_gap → x = +∞ ： inv_T_dn   静子出口状态 → 下游管道系数
    X_sys = inv_T_dn @ B_sta @ G_gap @ B_rot @ T_up

    return X_sys


def build_characteristic_matrix(s_val, n_val=1, params=None, params_file=None):
    """构建特征矩阵 Y_sys = [EC·X_sys; IC]

    根据边界条件文档，特征值方程为：
        det(Y_sys) = 0，其中 Y_sys = [EC·X_sys; IC]

    上游无限长管道：IC = [[0,1,0],[0,0,1]]  (B_n=0, C_n=0)
    下游无限长管道：EC = [[1,0,0]]            (A_n=0)

    参数：
    --------
    s_val : complex
        拉普拉斯变量值
    n_val : int
        周向波数
    params : dict or None
        预计算的参数字典, 为 None 时自动从文件加载
    params_file : str or None
        参数文件路径，为 None 时使用默认文件
    """
    X_sys = build_system_matrix(s_val, n_val, params, params_file)

    IC = IC_upstream_infinite(n_val)  # 2×3
    EC = EC_downstream_infinite(n_val)  # 1×3

    Y_sys = np.vstack([
        EC @ X_sys,
        IC
    ])

    return Y_sys


def characteristic_equation(s_val, n_val=1, params=None, params_file=None):
    """计算特征方程的值：det(Y_sys) = det([EC·X_sys; IC])

    参数：
    --------
    s_val : complex
        拉普拉斯变量值
    n_val : int
        周向波数
    params : dict or None
        预计算的参数字典, 为 None 时自动从文件加载
    params_file : str or None
        参数文件路径，为 None 时使用默认文件
    """
    Y_sys = build_characteristic_matrix(s_val, n_val, params, params_file)
    return np.linalg.det(Y_sys)


def make_system_matrix_function(n_val=1, params_file=None, params=None):
    """为 solver 提供矩阵接口：预加载参数,避免每次求值重复读文件

    参数：
    --------
    n_val : int
        周向波数
    params_file : str or None
        参数文件路径，为 None 时使用默认文件
    params : dict or None
        预计算的参数字典, 为 None 时自动从文件加载
    """
    if params is None:
        params = compute_stage_params(load_params(params_file or DEFAULT_PARAMS_FILE))
    def func(s):
        return build_characteristic_matrix(s, n_val, params)
    return func


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from Eigen_Value_Solver.Eigen_Hunter import hybrid_hunt

    # 加载参数并计算转子+静子派生参数（使用默认配置文件）
    params = compute_stage_params(load_params(DEFAULT_PARAMS_FILE))

    # 检查转子-静子间隙的速度连续性
    check_stage_interface(params, raise_on_fail=False)

    print()
    print("  (｡♥‿♥｡)  单级 (Rotor + Gap + Stator) 特征值 (Eigen_Hunter)")
    L_gap_print = params['L_gap']
    if np.isinf(L_gap_print):
        print(f"  轴向间隙长度  L_gap = inf (将对比 L_gap=0 与 L_gap=2.4)")
    else:
        print(f"  轴向间隙长度  L_gap = {L_gap_print:.4f}")
    print()
    # 两个执行盘组件（转子 + 静子）→ 每个 n 期望有2个主特征值
    header = f"  {'n':>3s}  {'#':>2s}  {'σ':>10s}  {'ω':>10s}  {'s':>24s}  {'稳?':>6s}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    def run_sweep(params_input, label=""):
        """对给定参数运行 n=1..6 特征值扫描"""
        res = {}
        for n_val in range(1, 7):
            Y_func = make_system_matrix_function(n_val, params=params_input)
            eigs = hybrid_hunt(Y_func,
                               sigma_range=(-6.0, 2.0),
                               omega_range=(-6.0, 6.0),
                               verbose=False, true_root_tol=1e-10)
            eigs_sorted = sorted(eigs, key=lambda x: x.real, reverse=True)[:2]
            res[n_val] = eigs_sorted
            if label:
                if eigs_sorted:
                    for k, eig in enumerate(eigs_sorted, start=1):
                        sigma = eig.real
                        omega = -eig.imag
                        s_str = f"{sigma:+.5f}-{omega:+.5f}j" if omega >= 0 else f"{sigma:+.5f}+{-omega:+.5f}j"
                        stability = "不稳定" if sigma > 0 else "稳定"
                        n_show = f"{n_val:3d}" if k == 1 else "   "
                        print(f"  [{label}] {n_show}  {k:>2d}  {sigma:+10.5f}  {omega:+10.5f}  {s_str:>24s}  {stability:>6s}")
                    if len(eigs_sorted) < 2:
                        print(f"  [{label}] {'':>3s}  {2:>2d}  {'--':>10s}  {'--':>10s}  {'区域内未找到第2个':>24s}  {'?':>6s}")
                else:
                    print(f"  [{label}] {n_val:3d}  {'-':>2s}  {'--':>10s}  {'--':>10s}  {'未找到特征值':>24s}  {'?':>6s}")
        return res

    L_gap_val = params.get('L_gap', 0.0)

    # 若 L_gap=inf，则运行 L_gap=0 和 L_gap=2.4 两个工况并叠加绘图
    multi_case = np.isinf(L_gap_val)
    if multi_case:
        params_0 = params.copy()
        params_0['L_gap'] = 0.0
        print(f"\n  --- L_gap = 0 ---\n")
        results = run_sweep(params_0, label="L_gap=0")

        print(f"\n  --- L_gap = 2.4 ---\n")
        params_24 = params.copy()
        params_24['L_gap'] = 2.4
        results_24 = run_sweep(params_24, label="L_gap=2.4")
    else:
        results = run_sweep(params, label=f"L_gap={L_gap_val}")

    # ── 绘制特征值分布图 ──
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams['legend.fontsize'] = 14

    fig, ax = plt.subplots(figsize=(9, 7))
    n_list = sorted(results.keys())

    def _dist(p, q):
        return np.sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2)

    def plot_case(ax, results_dict, markers, facecolor, edgecolors,
                  track_colors, label_prefix):
        """在 ax 上绘制一组特征值结果，返回 legend handles"""
        legend_handles = []
        track_a, track_b = [], []

        for n_val in n_list:
            eigs = results_dict.get(n_val, [])
            mkr = markers[(n_val - 1) % len(markers)]

            if not eigs:
                track_a.append((np.nan, np.nan))
                track_b.append((np.nan, np.nan))
                continue

            pts = [(eig.real, -eig.imag) for eig in eigs]
            for sigma, omega in pts:
                ax.scatter(sigma, omega,
                           s=130, facecolor=facecolor, marker=mkr,
                           edgecolors=edgecolors, linewidths=1.5,
                           zorder=3)

            if len(track_a) == 0:
                track_a.append(pts[0])
                track_b.append(pts[1] if len(pts) > 1 else (np.nan, np.nan))
            else:
                la, lb = track_a[-1], track_b[-1]
                if len(pts) == 1:
                    da = _dist(pts[0], la) if not np.isnan(la[0]) else np.inf
                    db = _dist(pts[0], lb) if not np.isnan(lb[0]) else np.inf
                    if da <= db:
                        track_a.append(pts[0]); track_b.append((np.nan, np.nan))
                    else:
                        track_a.append((np.nan, np.nan)); track_b.append(pts[0])
                else:
                    d_aa = _dist(pts[0], la) if not np.isnan(la[0]) else np.inf
                    d_bb = _dist(pts[1], lb) if not np.isnan(lb[0]) else np.inf
                    d_ab = _dist(pts[1], la) if not np.isnan(la[0]) else np.inf
                    d_ba = _dist(pts[0], lb) if not np.isnan(lb[0]) else np.inf
                    if d_aa + d_bb <= d_ab + d_ba:
                        track_a.append(pts[0]); track_b.append(pts[1])
                    else:
                        track_a.append(pts[1]); track_b.append(pts[0])

            legend_handles.append(
                (plt.Line2D([0], [1], linestyle='--', color='gray',
                            marker=mkr, markerfacecolor=facecolor,
                            markeredgecolor=edgecolors, markeredgewidth=1.2,
                            linewidth=1.2, markersize=8),
                 f'{label_prefix} n={n_val}'))

        for track, color in zip([track_a, track_b], track_colors):
            sigmas = [p[0] for p in track]
            omegas = [p[1] for p in track]
            ax.plot(sigmas, omegas, '--', color=color,
                    linewidth=1.2, alpha=0.5, zorder=2)

        return legend_handles

    markers_all = ['o', 's', '^', 'D', 'v', 'p']

    # 当前 L_gap: 蓝色实心标记
    label_first = 'L_gap=0' if multi_case else f'L_gap={L_gap_val}'
    legend_handles = plot_case(ax, results, markers_all,
                               facecolor='lightblue', edgecolors='blue',
                               track_colors=['blue', 'cyan'],
                               label_prefix=label_first)

    if multi_case:
        # L_gap = 2.4: 黑色空心标记
        legend_handles += plot_case(ax, results_24, markers_all,
                                    facecolor='white', edgecolors='black',
                                    track_colors=['black', 'gray'],
                                    label_prefix='L_gap=2.4')

    ax.axvline(x=0.0, color='red', linestyle='--', linewidth=1.8, zorder=2)
    ax.axhline(y=0.0, color='gray', linestyle=':', linewidth=1.0, zorder=1)

    handles, labels = zip(*legend_handles)
    ax.legend(handles, labels, loc='upper left', fontsize=14,
              handlelength=2.5)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.fill_betweenx([ylim[0] - 100, ylim[1] + 100],
                     0, xlim[1] + 100,
                     color='red', alpha=0.05, zorder=0)
    ax.set_xlim(xlim[0] - abs(xlim[0]) * 0.08, xlim[1] + abs(xlim[1]) * 0.08)
    ax.set_ylim(ylim[0] - abs(ylim[0]) * 0.08, ylim[1] + abs(ylim[1]) * 0.08)

    title = 'Single-Stage Compression System Eigenvalues'
    if multi_case:
        title += ' (L_gap→inf and L_gap→0)'
    ax.set_xlabel(r'$\sigma$ (Growth Rate)', fontsize=14)
    ax.set_ylabel(r'$\omega$ (Rotating Speed)', fontsize=14)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    fname = 'Eigenvalues_Stage_compare.png' if multi_case else 'Eigenvalues_Stage.png'
    save_path = os.path.join(_STABILITY_ROOT, '..', fname)
    fig.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f'[saved] {os.path.abspath(save_path)}')
    plt.show()
