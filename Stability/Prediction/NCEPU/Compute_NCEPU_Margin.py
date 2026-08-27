"""Compute_NCEPU_Margin.py - NCEPU 多工况特征值计算 (基于 CFX 提取数据)

从 CFX 后处理导出的气动参数 CSV 读取各工况数据, 逐工况求解特征值,
将结果保存为 CSV 文件供绘图脚本使用。

气动参数 CSV 由 Geo_and_Base_Flow/NCEPU_CFX_Loader.py 加载,
几何参数在 Geo_and_Base_Flow/NCEPU_params.txt 中设置。

用法:
  python Compute_NCEPU_Margin.py --csv NCEPU_CFX.csv
  python Compute_NCEPU_Margin.py --csv NCEPU_CFX.csv --rpm 14000
  python Compute_NCEPU_Margin.py --csv NCEPU_CFX.csv --nmax 5 --no-early-stop
"""

import numpy as np
import argparse
import sys
import os

# 将 Stability 根目录加入 sys.path
STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if STABILITY_ROOT not in sys.path:
    sys.path.insert(0, STABILITY_ROOT)

from Geo_and_Base_Flow.NCEPU_CFX_Loader import (
    load_all_operating_points, DEFAULT_CSV_FILE, DEFAULT_PARAMS_FILE
)
from Model.Model_NCEPU import build_cc3_characteristic_matrix
from Eigen_Value_Solver.Eigen_Hunter import hybrid_hunt


# ──────────────────────────────
# 多工况特征值求解
# ──────────────────────────────
def sweep_eigenvalues_from_cfx(csv_path, params_file, rpm_filter=None,
                               n_max=5, sigma_range=(-6.0, 2.0),
                               omega_range=(-6.0, 6.0),
                               max_modes=2, true_root_tol=1e-10,
                               early_stop=True):
    """从 CFX CSV 逐工况计算特征值

    参数
    ----
    csv_path : str
        CFX 气动参数 CSV 路径
    params_file : str
        NCEPU_params.txt 路径
    rpm_filter : float or None
        若指定, 仅计算该转速的工况
    early_stop : bool
        若为 True, 一旦检测到不稳定特征值 (σ > 0) 即停止后续流量点计算

    返回
    ----
    all_results : list of dict
    """
    all_op = load_all_operating_points(csv_path, params_file)
    if rpm_filter is not None:
        all_op = [op for op in all_op if abs(op['rpm'] - rpm_filter) < 1.0]

    if not all_op:
        print("没有匹配的工况点")
        return []

    # 按 (rpm, m) 排序: 同转速内流量降序, 从大流量(堵塞/稳定侧)往
    # 小流量(喘振侧)扫, 使热启动从稳定侧开始、早停在越过 σ=0 时才触发
    all_op.sort(key=lambda x: (x['rpm'], -x['m']))

    all_results = []
    n_range = range(1, n_max + 1)
    _prev_primary = {}

    for i, op in enumerate(all_op):
        m_flow = op['m']
        rpm = op['rpm']
        params = op['params']
        extra = op['extra']

        print(f"  [{i+1}/{len(all_op)}] rpm={rpm:.0f} m={m_flow:.3f} kg/s ...",
              end=" ", flush=True)

        # 检查关键参数有效性
        if (params.get('dL_dtanbeta1_imp', 0.0) == 0.0
                or params.get('AR_imp', 0.0) <= 0.0
                or params.get('lambda_imp', 0.0) <= 0.0):
            print("SKIP (关键参数未设置: dL/AR/lambda)")
            continue

        eigs_n = {}
        unstable_detected = False
        for n_val in n_range:
            def Y_func(s, _n=n_val, _p=params):
                return build_cc3_characteristic_matrix(s, _n, _p)

            found = hybrid_hunt(Y_func,
                                sigma_range=sigma_range,
                                omega_range=omega_range,
                                verbose=False,
                                true_root_tol=true_root_tol)
            found_sorted = sorted(found, key=lambda x: x.real, reverse=True)

            # 热启动: 优先选择离上一流量点参考根最近的根
            ref = _prev_primary.get((rpm, n_val))
            if ref is not None and found_sorted:
                found_sorted.sort(key=lambda e: abs(e - ref))
            eigs_sorted = found_sorted[:max_modes]
            if eigs_sorted:
                _prev_primary[(rpm, n_val)] = eigs_sorted[0]
            eigs_n[n_val] = eigs_sorted
            if any(eig.real > 0 for eig in eigs_sorted):
                unstable_detected = True

        status = "UNSTABLE" if unstable_detected else "stable"
        pi_tt_str = f"π_tt={extra['pi_tt']:.3f}" if extra['pi_tt'] else "π_tt=N/A"
        print(f"OK ({status}, {pi_tt_str})")

        all_results.append({
            'rpm': rpm,
            'm': m_flow,
            'eigenvalues': eigs_n,
            'pi_tt': extra['pi_tt'] if extra['pi_tt'] else 0.0,
        })

        if early_stop and unstable_detected:
            remaining = len(all_op) - (i + 1)
            print(f"  [early stop] 检测到不稳定模态 (σ > 0), "
                  f"跳过剩余 {remaining} 个工况点")
            break

    return all_results


# ──────────────────────────────
# 保存 / 加载结果
# ──────────────────────────────
def save_results(all_results, filepath):
    """将计算结果保存为 CSV 文件"""
    with open(filepath, 'w', encoding='utf-8-sig') as f:
        f.write("rpm,m,pi_tt,n,sigma,omega\n")
        for res in all_results:
            m = res['m']
            rpm = res['rpm']
            pi_tt = res['pi_tt']
            for n_val in sorted(res['eigenvalues'].keys()):
                for eig in res['eigenvalues'][n_val]:
                    sigma = eig.real
                    omega = eig.imag
                    f.write(f"{rpm:.0f},{m:.6f},{pi_tt:.6f},"
                            f"{n_val},{sigma:.10f},{omega:.10f}\n")
    print(f"[saved] {filepath}")


def load_results(filepath):
    """从 CSV 文件加载计算结果, 返回 all_results 列表"""
    import csv
    all_results = []
    current_key = None
    current_res = None

    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        for row in reader:
            if not row:
                continue
            rpm = float(row[0])
            m = float(row[1])
            pi_tt = float(row[2])
            n_val = int(row[3])
            sigma = float(row[4])
            omega = float(row[5])

            key = (rpm, m)
            if key != current_key:
                if current_res is not None:
                    all_results.append(current_res)
                current_key = key
                current_res = {
                    'rpm': rpm, 'm': m,
                    'eigenvalues': {},
                    'pi_tt': pi_tt,
                }

            if n_val not in current_res['eigenvalues']:
                current_res['eigenvalues'][n_val] = []
            current_res['eigenvalues'][n_val].append(complex(sigma, omega))

    if current_res is not None:
        all_results.append(current_res)

    return all_results


# ──────────────────────────────
# 模态追踪
# ──────────────────────────────
def track_modes(all_results, n_val, max_modes=2):
    """对给定 n, 追踪特征值随流量的变化轨迹"""
    if not all_results:
        return []

    per_point = []
    for res in all_results:
        eigs = res['eigenvalues'].get(n_val, [])
        sorted_eigs = sorted(eigs, key=lambda x: x.real,
                             reverse=True)[:max_modes]
        per_point.append((res['m'], sorted_eigs))

    if not per_point:
        return []

    n_modes = len(per_point[0][1])
    tracks = [[(per_point[0][0], per_point[0][1][k])]
              for k in range(n_modes)]

    for i in range(1, len(per_point)):
        m_curr, eigs_curr = per_point[i]
        if not eigs_curr:
            for track in tracks:
                track.append((m_curr, None))
            continue

        used = set()
        track_matched = set()
        dists = []
        for ti, track in enumerate(tracks):
            if not track or track[-1][1] is None:
                continue
            last_eig = track[-1][1]
            for ei, eig in enumerate(eigs_curr):
                dists.append((abs(eig - last_eig), ti, ei))
        dists.sort()

        for dist, ti, ei in dists:
            if ti in track_matched or ei in used:
                continue
            tracks[ti].append((m_curr, eigs_curr[ei]))
            track_matched.add(ti)
            used.add(ei)

        for ti in range(len(tracks)):
            if ti not in track_matched:
                tracks[ti].append((m_curr, None))

        for ei, eig in enumerate(eigs_curr):
            if ei not in used:
                new_track = [(per_point[j][0], None) for j in range(i)]
                new_track.append((m_curr, eig))
                tracks.append(new_track)

    return tracks


# ──────────────────────────────
# 表格输出
# ──────────────────────────────
def print_eigenvalue_table(all_results):
    """打印各工况特征值表格"""
    for res in all_results:
        m = res['m']
        rpm = res['rpm']
        print(f"\n  rpm={rpm:.0f}  m = {m:.3f} kg/s  "
              f"π_tt = {res['pi_tt']:.4f}")
        print(f"  {'n':>3s}  {'#':>2s}  {'σ':>10s}  "
              f"{'ω':>10s}  {'稳?':>6s}")
        print("  " + "-" * 40)

        for n_val in sorted(res['eigenvalues'].keys()):
            eigs = res['eigenvalues'][n_val]
            for k, eig in enumerate(eigs, start=1):
                sigma = eig.real
                omega = -eig.imag
                stability = "不稳定" if sigma > 0 else "稳定"
                n_show = f"{n_val:3d}" if k == 1 else "   "
                print(f"  {n_show}  {k:>2d}  {sigma:+10.5f}  "
                      f"{omega:+10.5f}  {stability:>6s}")


def find_stability_boundary(all_results, n_val=1, mode_idx=0):
    """找出哪个模态最先越过稳定性边界 (σ=0)"""
    tracks = track_modes(all_results, n_val, max_modes=2)
    if mode_idx >= len(tracks):
        return None

    track = tracks[mode_idx]
    for i in range(len(track) - 1):
        m1, eig1 = track[i]
        m2, eig2 = track[i + 1]
        if eig1 is None or eig2 is None:
            continue
        if eig1.real <= 0 and eig2.real > 0:
            return (f"Mode {mode_idx+1} (n={n_val}) 在 "
                    f"m ∈ [{m2:.2f}, {m1:.2f}] kg/s 越过稳定边界 "
                    f"(σ: {eig2.real:+.4f} → {eig1.real:+.4f})")
        if eig1.real > 0 and eig2.real <= 0:
            return (f"Mode {mode_idx+1} (n={n_val}) 在 "
                    f"m ∈ [{m1:.2f}, {m2:.2f}] kg/s 回到稳定区 "
                    f"(σ: {eig1.real:+.4f} → {eig2.real:+.4f})")
    return None


# ──────────────────────────────
# 主入口
# ──────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='NCEPU 多工况特征值计算 (基于 CFX 提取数据)'
    )
    parser.add_argument('--csv', type=str, default=DEFAULT_CSV_FILE,
                        help='CFX 气动参数 CSV 文件路径')
    parser.add_argument('--params', type=str, default=DEFAULT_PARAMS_FILE,
                        help='NCEPU 几何参数文件路径')
    parser.add_argument('--rpm', type=float, default=None,
                        help='仅计算指定转速 [rpm]')
    parser.add_argument('--nmax', type=int, default=5,
                        help='最大周向波数 n')
    parser.add_argument('--no-early-stop', action='store_true', default=False,
                        dest='no_early_stop',
                        help='禁用早停: 即使检测到不稳定也扫完所有工况点')
    args = parser.parse_args()

    if not os.path.exists(args.csv):
        print(f"CSV not found: {args.csv}")
        print("请先从 CFX 导出气动参数 CSV, 或使用 --csv 指定路径。")
        print(f"可运行: python -m Geo_and_Base_Flow.NCEPU_CFX_Loader --template "
              f"--csv {args.csv} 生成模板")
        sys.exit(1)

    if not os.path.exists(args.params):
        print(f"Params file not found: {args.params}")
        sys.exit(1)

    print(f"NCEPU 特征值计算: CSV = {args.csv}\n")

    all_results = sweep_eigenvalues_from_cfx(
        args.csv, args.params, rpm_filter=args.rpm,
        n_max=args.nmax, early_stop=not args.no_early_stop)

    if all_results:
        print_eigenvalue_table(all_results)

    print()
    for n_val in range(1, args.nmax + 1):
        for mode_idx in range(2):
            msg = find_stability_boundary(all_results, n_val, mode_idx)
            if msg:
                print(f"  [!] {msg}")

    if all_results:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        rpms_in_result = sorted(set(r['rpm'] for r in all_results))
        rpm_tag = f"{int(rpms_in_result[0])}rpm" if len(rpms_in_result) == 1 else "multi_rpm"
        save_path = os.path.join(script_dir, f"NCEPU_Eigenvalues_{rpm_tag}.csv")
        save_results(all_results, save_path)
