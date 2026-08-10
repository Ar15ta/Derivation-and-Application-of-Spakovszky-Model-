# Model_Isolated_Stator.py - 静子孤立模型：上下游无限管道 + 孤立静叶

import os
import sys

# 确保 Stability 根目录在 sys.path 中, 使 Matrix/ 等包可导入
_STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _STABILITY_ROOT not in sys.path:
    sys.path.insert(0, _STABILITY_ROOT)

import numpy as np

# 导入组件模块（每个模块自己读取参数）
from Matrix.stator import B_sta_n_num, compute_stator_params, load_params
from Matrix.axial import Tn_ax_coeff_num, inv_Tn_ax_coeff_num
from Matrix.boundary import IC_upstream_infinite, EC_downstream_infinite


def build_system_matrix(s_val, n_val=1, params=None):
    """
    构建完整系统传递矩阵（静子）

    参数：
    --------
    s_val : complex
        拉普拉斯变量值
    n_val : int
        周向波数
    params : dict or None
        预计算的参数字典, 为 None 时自动从文件加载

    返回：
    --------
    X_sys : np.ndarray (3×3)
        完整系统传递矩阵（上游系数 c^up → 下游系数 c^dn）
    """
    if params is None:
        params = load_params()
        params = compute_stator_params(params)

    Vx = params['Vx']
    mu_up = params.get('mu_up', 1.0)
    mu_dn = params.get('mu_dn', 1.0)

    # 孤立静叶设定：使用静子计算得到的上下游平均周向速度
    Vtheta_up = params.get('Vtheta_bar1_sta', 0.0)
    Vtheta_dn = params.get('Vtheta_bar2_sta', 0.0)

    # 静子传递矩阵
    B_sta = B_sta_n_num(s=s_val, n=n_val, params=params)

    # 上游管道：静子进口 x=0⁻ 处的状态 = T_up(0;0) · c^up
    T_up = Tn_ax_coeff_num(x=0, s=s_val, n=n_val,
                           Vx_bar=Vx, Vtheta_bar=Vtheta_up, x0=0, mu=mu_up)

    # 下游管道：静子出口 x=0⁺ 处的状态 = T_dn(0;0) · c^dn
    inv_T_dn = inv_Tn_ax_coeff_num(x=0, s=s_val, n=n_val,
                                   Vx_bar=Vx, Vtheta_bar=Vtheta_dn, x0=0, mu=mu_dn)

    # 系统传递矩阵
    X_sys = inv_T_dn @ B_sta @ T_up

    return X_sys


def build_characteristic_matrix(s_val, n_val=1, params=None):
    """构建特征矩阵 Y_sys = [EC·X_sys; IC]

    根据边界条件文档，特征值方程为：
        det(Y_sys) = 0，其中 Y_sys = [EC·X_sys; IC]

    上游无限长管道：IC = [[0,1,0],[0,0,1]]  (B_n=0, C_n=0)
    下游无限长管道：EC = [[1,0,0]]            (A_n=0)
    """
    X_sys = build_system_matrix(s_val, n_val, params)

    IC = IC_upstream_infinite(n_val)  # 2×3
    EC = EC_downstream_infinite(n_val)  # 1×3

    Y_sys = np.vstack([
        EC @ X_sys,
        IC
    ])

    return Y_sys


def characteristic_equation(s_val, n_val=1, params=None):
    """计算特征方程的值：det(Y_sys) = det([EC·X_sys; IC])"""
    Y_sys = build_characteristic_matrix(s_val, n_val, params)
    return np.linalg.det(Y_sys)


def make_system_matrix_function(n_val=1):
    """为 solver 提供矩阵接口：预加载参数,避免每次求值重复读文件"""
    params = compute_stator_params(load_params())
    def func(s):
        return build_characteristic_matrix(s, n_val, params)
    return func


if __name__ == '__main__':
    from Eigen_Value_Solver.Eigen_Hunter import hybrid_hunt

    # 加载参数并计算静子派生参数
    params = compute_stator_params(load_params())

    print()
    print("  (=^･ω･^=)  孤立静叶 特征值 (Eigen_Hunter)")
    print()
    header = f"  {'n':>3s}  {'σ':>10s}  {'ω':>10s}  {'s':>24s}  {'稳?':>6s}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    for n_val in range(1, 7):
        Y_func = make_system_matrix_function(n_val)
        eigs = hybrid_hunt(Y_func,
                           sigma_range=(-3.0, 3.0),
                           omega_range=(-6.0, 6.0),
                           verbose=False, true_root_tol=1e-10)
        eigs = sorted(eigs, key=lambda x: x.real, reverse=True)[:1]

        if eigs:
            eig = eigs[0]
            sigma = eig.real
            omega = -eig.imag
            s_str = f"{sigma:+.5f}-{omega:+.5f}j" if omega >= 0 else f"{sigma:+.5f}+{-omega:+.5f}j"
            stability = "不稳定" if sigma > 0 else "稳定"
            print(f"  {n_val:3d}  {sigma:+10.5f}  {omega:+10.5f}  {s_str:>24s}  {stability:>6s}")
        else:
            print(f"  {n_val:3d}  {'--':>10s}  {'--':>10s}  {'未找到特征值':>24s}  {'?':>6s}")
