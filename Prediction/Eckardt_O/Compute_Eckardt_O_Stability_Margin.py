"""Compute_Eckardt_O_Eigenvalues.py - Eckardt O 多工况特征值计算

给定流量范围, 逐工况调用均值线程序计算基流参数, 再求解特征值,
将结果保存为 .npz 文件供绘图脚本使用。

用法:
  python Compute_Eckardt_O_Eigenvalues.py --rpm 14000 --m-min 3.0 --m-max 7.0 --points 10
  python Compute_Eckardt_O_Eigenvalues.py --rpm 14000 --m 5.0  (单工况)
"""

import numpy as np
import argparse
import json
import sys
import os

# 将 Stability 根目录加入 sys.path
STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if STABILITY_ROOT not in sys.path:
    sys.path.insert(0, STABILITY_ROOT)

# 均值线包路径
MEANLINE_PKG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))),
    "MeanLIne_Design", "RadCompressor_MeanLine", "radcomp-main"
)
YAML_DATA = os.path.join(MEANLINE_PKG, "data", "known_compressors.yml")
if MEANLINE_PKG not in sys.path:
    sys.path.insert(0, MEANLINE_PKG)

from Geo_and_Base_Flow.Eckardt_O_Geo_Calculator import (
    compute_base_flow, load_yaml, find_record
)
from Model.Model_Eckardt_O import build_cc3_characteristic_matrix
from Matrix.rotor import load_params
from Matrix.impeller import compute_impeller_params
from Eigen_Value_Solver.Eigen_Hunter import hybrid_hunt

# 参数文件路径
PARAMS_FILE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "Geo_and_Base_Flow", "Eckardt_O_params.txt"
)


# ──────────────────────────────
# 参数构建: Geo_Calculator → 稳定性模型
# ──────────────────────────────
def update_params_from_geo(params, geo_data):
    """用 Geo_Calculator 的计算结果更新稳定性模型参数字典

    几何参数保持不变, 只更新流动相关参数并重算派生量。
    """
    U2 = geo_data['U2']

    params['R1_R2'] = params['STA2_radial'] / params['STA3_radial']
    params['r_impeller_exit'] = params['STA3_radial']
    params['r_diffuser_exit'] = params['STA4_radial']

    params['Vx_bar1_imp'] = geo_data['sta2_cm'] / U2
    params['Vtheta_bar1_imp'] = geo_data['sta2_ct'] / U2
    params['Vr_bar2_imp'] = geo_data['sta3_cm'] / U2
    params['Vtheta_bar2_imp'] = geo_data['sta3_ct'] / U2
    params['beta1_imp_deg'] = geo_data['sta2_beta']
    params['alpha2_imp_deg'] = geo_data['sta3_alpha']
    params['lambda_imp'] = geo_data['lambda_imp']
    params['tau_imp'] = geo_data['tau_imp']
    params['dL_dtanbeta1_imp'] = geo_data['dL_dtanbeta1']
    params['AR_imp'] = geo_data['AR_imp_with_density']
    params['Vx_up'] = geo_data['sta1_cm'] / U2
    params['Q'] = geo_data['Q']
    params['G'] = geo_data['G']

    params['tan_beta1_imp'] = np.tan(np.deg2rad(params['beta1_imp_deg']))
    params['tan_beta2_imp'] = np.tan(np.deg2rad(params['beta2_imp_deg']))
    params['tan_alpha1_imp'] = np.tan(np.deg2rad(params['alpha1_imp_deg']))

    params['Vx_plenum'] = params['Q'] / params['r_diffuser_exit']
    params['Vtheta_plenum'] = params['G'] / params['r_diffuser_exit']

    params = compute_impeller_params(params)
    return params

# ──────────────────────────────
# 多工况特征值求解
# ──────────────────────────────
def sweep_eigenvalues(record, rpm, m_values, n_max=5,
                      sigma_range=(-6.0, 2.0), omega_range=(-6.0, 6.0),
                      max_modes=1, true_root_tol=1e-10,
                      early_stop=True):
    """逐工况计算特征值

    参数
    ----
    early_stop : bool
        若为 True, 一旦检测到不稳定特征值 (σ > 0) 即停止后续流量点计算。

    返回
    ----
    all_results : list of dict
        每个元素: {'m': m_flow, 'eigenvalues': {n: [eig1, eig2, ...]},
                   'PR': ..., 'eff': ...}
    """
    base_params = load_params(PARAMS_FILE)

    all_results = []
    n_range = range(1, n_max + 1)

    for i, m_flow in enumerate(m_values):
        print(f"  [{i+1}/{len(m_values)}] m = {m_flow:.3f} kg/s ...",
              end=" ", flush=True)
        try:
            geo_data = compute_base_flow(record, rpm, m_flow)
        except RuntimeError as e:
            print(f"SKIP ({e})")
            continue

        if not geo_data.get('flow_ok', True):
            print(f"SKIP (flow_ok=False) m = {m_flow:.3f} kg/s")
            continue

        dL_dtanbeta1 = geo_data.get('dL_dtanbeta1', None)
        if dL_dtanbeta1 is None or np.isnan(dL_dtanbeta1):
            print(f"SKIP (invalid dL_dtanbeta1) m = {m_flow:.3f} kg/s")
            continue

        params = dict(base_params)
        params = update_params_from_geo(params, geo_data)

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
            eigs_sorted = sorted(found, key=lambda x: x.real,
                                 reverse=True)[:max_modes]
            eigs_n[n_val] = eigs_sorted
            if any(eig.real > 0 for eig in eigs_sorted):
                unstable_detected = True

        status = "UNSTABLE" if unstable_detected else "stable"
        print(f"OK ({status}, PR={geo_data['PR']:.3f}, η={geo_data['eff']:.3f})")
        all_results.append({
            'm': m_flow,
            'eigenvalues': eigs_n,
            'PR': geo_data['PR'],
            'eff': geo_data['eff'],
        })

        if early_stop and unstable_detected:
            remaining = len(m_values) - (i + 1)
            print(f"  [early stop] 检测到不稳定模态 (σ > 0), "
                  f"跳过剩余 {remaining} 个流量点")
            break

    return all_results

# ──────────────────────────────
# 保存 / 加载结果
# ──────────────────────────────
def save_results(all_results, rpm, filepath):
    """将计算结果保存为 CSV 文件"""
    with open(filepath, 'w', encoding='utf-8-sig') as f:
        f.write(f"rpm,{rpm}\n")
        f.write("m,PR,eff,n,sigma,omega\n")
        for res in all_results:
            m = res['m']
            PR = res['PR']
            eff = res['eff']
            for n_val in sorted(res['eigenvalues'].keys()):
                for eig in res['eigenvalues'][n_val]:
                    sigma = eig.real
                    omega = eig.imag
                    f.write(f"{m:.6f},{PR:.6f},{eff:.6f},{n_val},{sigma:.10f},{omega:.10f}\n")
    print(f"[saved] {filepath}")


def load_results(filepath):
    """从 CSV 文件加载计算结果, 返回 all_results 列表"""
    import csv
    all_results = []
    rpm = None
    current_m = None
    current_res = None
    
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == 'rpm':
                rpm = float(row[1])
            elif row[0] == 'm':
                pass
            else:
                m = float(row[0])
                PR = float(row[1])
                eff = float(row[2])
                n_val = int(row[3])
                sigma = float(row[4])
                omega = float(row[5])
                
                if m != current_m:
                    if current_res is not None:
                        all_results.append(current_res)
                    current_m = m
                    current_res = {
                        'm': m,
                        'eigenvalues': {},
                        'PR': PR,
                        'eff': eff,
                    }
                
                if n_val not in current_res['eigenvalues']:
                    current_res['eigenvalues'][n_val] = []
                current_res['eigenvalues'][n_val].append(complex(sigma, omega))
    
    if current_res is not None:
        all_results.append(current_res)
    
    return all_results, rpm


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
        print(f"\n  m = {m:.3f} kg/s  PR = {res['PR']:.4f}  "
              f"\u03b7 = {res['eff']:.4f}")
        print(f"  {'n':>3s}  {'#':>2s}  {'\u03c3':>10s}  "
              f"{'\u03c9':>10s}  {'\u7a33?':>6s}")
        print("  " + "-" * 40)

        for n_val in sorted(res['eigenvalues'].keys()):
            eigs = res['eigenvalues'][n_val]
            for k, eig in enumerate(eigs, start=1):
                sigma = eig.real
                omega = -eig.imag
                stability = "\u4e0d\u7a33\u5b9a" if sigma > 0 else "\u7a33\u5b9a"
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
                    f"m \u2208 [{m2:.2f}, {m1:.2f}] kg/s 越过稳定边界 "
                    f"(\u03c3: {eig2.real:+.4f} \u2192 {eig1.real:+.4f})")
        if eig1.real > 0 and eig2.real <= 0:
            return (f"Mode {mode_idx+1} (n={n_val}) 在 "
                    f"m \u2208 [{m1:.2f}, {m2:.2f}] kg/s 回到稳定区 "
                    f"(\u03c3: {eig1.real:+.4f} \u2192 {eig2.real:+.4f})")
    return None


# ──────────────────────────────
# 主入口
# ──────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Eckardt O 多工况特征值计算 (仅计算)'
    )
    parser.add_argument('--rpm', type=float, default=18000.0,
                        help='转速 [rpm]')
    parser.add_argument('--m', type=float, default=None,
                        help='单工况流量 [kg/s]')
    parser.add_argument('--m-min', type=float, default=3.0,
                        help='最小流量 [kg/s]')
    parser.add_argument('--m-max', type=float, default=6.0,
                        help='最大流量 [kg/s]')
    parser.add_argument('--points', type=int, default=10,
                        help='流量数据点数')
    parser.add_argument('--nmax', type=int, default=5,
                        help='最大周向波数 n')
    parser.add_argument('--no-early-stop', action='store_true', default=False,
                        dest='no_early_stop',
                        help='禁用早停: 即使检测到不稳定也扫完所有流量点')
    args = parser.parse_args()

    # 加载 YAML 数据
    records = load_yaml(YAML_DATA)
    record = find_record(records, "Eckardt O")

    # 确定流量范围
    if args.m is not None:
        m_values = [args.m]
    else:
        m_values = list(np.linspace(args.m_min, args.m_max, args.points))

    print(f"Eckardt O 特征值计算: {args.rpm:.0f} rpm, "
          f"m = {min(m_values):.2f} ~ {max(m_values):.2f} kg/s, "
          f"{len(m_values)} 个工况点\n")

    # 逐工况计算
    all_results = sweep_eigenvalues(
        record, args.rpm, m_values, n_max=args.nmax,
        early_stop=not args.no_early_stop)

    # 打印表格
    if all_results:
        print_eigenvalue_table(all_results)

    # 稳定性边界检测
    print()
    for n_val in range(1, args.nmax + 1):
        for mode_idx in range(2):
            msg = find_stability_boundary(all_results, n_val, mode_idx)
            if msg:
                print(f"  [!] {msg}")

    # 保存结果
    if all_results:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        save_path = os.path.join(
            script_dir,
            f"Eckardt_O_Eigenvalues_{int(args.rpm)}rpm.csv")
        save_results(all_results, args.rpm, save_path)
