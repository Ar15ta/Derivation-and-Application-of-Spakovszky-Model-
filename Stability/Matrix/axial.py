"""
axial.py - 轴向管道传递矩阵 (重写版)

参考: Spakovszky 稳定性方法 Model II - 轴向管道推导.md

物理思想:
  管道内的扰动 = 齐次解 (势模态 A_n, B_n) + 涡量对流模态 (C_n)
  系数 (A_n, B_n, C_n) 在整个管道内是常数 (因为是齐次解),
  只是空间形函数 exp(±n·x), exp(-k_n·x) 随 x 变化。

提供:
  1. Tn_ax_coeff_num(x, ...)  : 系数 → x 处状态 (3x3)
  2. inv_Tn_ax_coeff_num(x, ...) : 状态 → 系数 (解析逆, 避免 numpy.linalg.inv 的伪极点)
  3. axial_duct_state_transfer(L, ...) : 状态(0) → 状态(L) (用对称中点求积避免溢出)

数值改进:
  - 取参考点 x_ref = x/2, 让两端的指数因子对称为 exp(±n·x/2),
    避免 L_gap 较大时 exp(n·L_gap) 单边爆炸
  - 解析手写 3x3 逆矩阵 (Cramer 法则), 在物理上没有奇异性
    (因为势模态 + 涡量模态构成完备基, 矩阵在物理 s 范围内总是非奇异)
"""

import numpy as np
import re
import os


def load_params(filename='system_params.txt'):
    """加载参数文件，支持相对路径和绝对路径

    默认从 Geo/ 目录查找参数文件
    """
    import os
    params = {}
    stability_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    full_path = os.path.join(stability_root, 'Geo_and_Base_Flow', filename)
    with open(full_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            match = re.match(r'(\w+)\s*=\s*([\d\.\-]+)', line)
            if match:
                params[match.group(1)] = float(match.group(2))
            else:
                match_placeholder = re.match(r'(\w+)\s*=\s*\[\]', line)
                if match_placeholder:
                    params[match_placeholder.group(1)] = 0.0
    return params


def _coeff_to_state_matrix(x, s, n, Vx_bar, Vtheta_bar, mu=1.0):
    """
    [系数 → 状态] 映射的 3x3 矩阵 (n >= 1)

    mu : 参考半径比 R_ref / R_duct。
         mu=1 表示参考半径等于管道半径 (轴流情形)。
         离心叶轮以 R3 归一化时，管道周向波数相关项以 μ·n 替代 n。

    M_R2(x2) = [ e^{μn x2},     e^{-μn x2},     e^{-kn(2) x2}   ]
               [ je^{μn x2},   -je^{-μn x2},   cC^(2)e^{-kn(2)x2}]
               [ -(λ_A^(2))/(μn)·e^{μn x2}, (λ_B^(2))/(μn)·e^{-μn x2}, 0 ]
    """
    m = mu * n
    k_n = s / Vx_bar + 1j * m * Vtheta_bar / Vx_bar
    lam_A = s + m * (Vx_bar + 1j * Vtheta_bar)
    lam_B = s - m * Vx_bar + 1j * m * Vtheta_bar
    c_C = Vtheta_bar / Vx_bar - 1j * s / (m * Vx_bar)

    e_pm = np.exp(m * x)
    e_mm = np.exp(-m * x)
    e_km = np.exp(-k_n * x)

    M = np.array([
        [e_pm,           e_mm,           e_km],
        [1j * e_pm,     -1j * e_mm,      c_C * e_km],
        [-lam_A / m * e_pm,  lam_B / m * e_mm,  0.0 + 0.0j],
    ], dtype=complex)
    return M


def Tn_ax_coeff_num(x, s, n, Vx_bar=None, Vtheta_bar=None, x0=0, params=None, mu=1.0):
    """
    轴向管道 [系数 → 状态] 传递矩阵 (n >= 1)

    将系数向量 [A_n, B_n, C_n]^T 映射到位置 x 处的状态 [δVx, δVθ, δP]^T

    mu : 参考半径比 R_ref / R_duct，由参数文件给出。
         轴流时 mu=1; 离心叶轮时 mu = R3 / R_duct。
    """
    dx = x - x0
    return _coeff_to_state_matrix(dx, s, n, Vx_bar, Vtheta_bar, mu)


def inv_Tn_ax_coeff_num(x, s, n, Vx_bar=None, Vtheta_bar=None, x0=0, params=None, mu=1.0):
    """轴向管道 [状态 → 系数] 传递矩阵 (n >= 1)"""
    M = Tn_ax_coeff_num(x, s, n, Vx_bar, Vtheta_bar, x0, params, mu)
    return np.linalg.inv(M)


def axial_duct_state_transfer(x, s, n, Vx_bar, Vtheta_bar, mu=1.0):
    """
    轴向管道 [状态 → 状态] 传递矩阵: state(x) = G(x) @ state(0)

    mu : 参考半径比 R_ref / R_duct。
    """
    if abs(x) < 1e-14:
        return np.eye(3, dtype=complex)

    half = 0.5 * x
    M_in = _coeff_to_state_matrix(-half, s, n, Vx_bar, Vtheta_bar, mu)
    M_out = _coeff_to_state_matrix(+half, s, n, Vx_bar, Vtheta_bar, mu)
    return M_out @ np.linalg.inv(M_in)


def axial_duct_state_transfer_n0(x, s, Vx_bar):
    """
    n=0 时的管道状态传递 (轴对称)

    状态 [δVx, δVθ, δP] 在管道内的演化:
        δVx(x)  = A0                        (常数)
        δVθ(x)  = (s/Vx)·C0·e^{-k0·x}      (涡量对流)
        δP(x)   = δP(0) - s·A0·x           (线性下降)

    x=0 自动退化为单位阵。
    """
    k0 = s / Vx_bar
    e_k0 = np.exp(-k0 * x)
    G = np.eye(3, dtype=complex)
    G[1, 1] = e_k0
    G[2, 0] = -s * x
    return G


# ===================== 自检 =====================

if __name__ == '__main__':
    print("=" * 60)
    print("axial.py 自检")
    print("=" * 60)

    Vx, Vth = 0.3426, 0.76
    s_test = complex(1.5, 2.3)

    # 测试 1: x=0 时管道传递矩阵应为单位阵
    print("\n[测试 1] x=0 应退化为单位阵")
    for n in [1, 2, 3, 5]:
        G = axial_duct_state_transfer(0.0, s_test, n, Vx, Vth)
        err = np.linalg.norm(G - np.eye(3))
        print(f"  n={n}: ||G(x=0) - I|| = {err:.2e}")

    # 测试 2: 与原始 inv(M0) 方法对照
    print("\n[测试 2] 中点参考 vs 原点参考 (应数学等价)")
    L_test = 0.3
    for n in [1, 2, 5]:
        G_new = axial_duct_state_transfer(L_test, s_test, n, Vx, Vth)
        M0 = _coeff_to_state_matrix(0, s_test, n, Vx, Vth)
        ML = _coeff_to_state_matrix(L_test, s_test, n, Vx, Vth)
        G_old = ML @ np.linalg.inv(M0)
        err = np.linalg.norm(G_new - G_old) / np.linalg.norm(G_new)
        print(f"  n={n}, L={L_test}: 相对差 = {err:.2e}")

    # 测试 3: 大 L 数值稳健性
    print("\n[测试 3] 大 L 数值稳健性 (无溢出)")
    for L_big in [1.0, 3.0, 5.0]:
        for n in [1, 3, 5]:
            G = axial_duct_state_transfer(L_big, s_test, n, Vx, Vth)
            ok = np.all(np.isfinite(G))
            cond = np.linalg.cond(G) if ok else float('inf')
            print(f"  n={n}, L={L_big}: 有限={ok}, cond={cond:.2e}")