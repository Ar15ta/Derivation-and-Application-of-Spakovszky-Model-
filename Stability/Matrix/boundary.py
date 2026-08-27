import numpy as np


def IC_upstream_infinite(n_val):
    """
    上游无限长管道边界条件（作用于管道系数 [A_n, B_n, C_n]ᵀ）。
    B_n = 0（无向下游衰减的势波逆流传播）
    C_n = 0（涡量只能向下游对流）
    返回 2×3 矩阵 IC。
    """
    return np.array([
        [0, 1, 0],   # B_n = 0
        [0, 0, 1]    # C_n = 0
    ], dtype=complex)



def EC_downstream_infinite(n_val):
    """
    下游无限长管道边界条件。
    令 A_n = 0（避免在 x→+∞ 时发散模式），C_n 在无限管道中允许随流动传播。
    返回 1×3 矩阵 EC。
    """
    return np.array([[1, 0, 0]], dtype=complex)


def B_identity():
    """占位：单位边界矩阵。"""
    return np.eye(3, dtype=complex)

def EC_downstream_plenum(x_exit, s, n, Vx_bar, Vtheta_bar, mu=1.0):
    """
    下游容腔边界条件 (n >= 1)

    mu : 参考半径比 R_ref / R_duct。
         轴流时 mu=1; 离心叶轮时 mu = R3 / R_duct。
    EC = [ (-s/(μn) - Vx - j Vθ) e^{μn x_exit},
           (s/(μn) - Vx + j Vθ) e^{-μn x_exit},
           0 ]
    """
    m = mu * n
    term1 = (-s / m - Vx_bar - 1j * Vtheta_bar) * np.exp(m * x_exit)
    term2 = (s / m - Vx_bar + 1j * Vtheta_bar) * np.exp(-m * x_exit)
    return np.array([[term1, term2, 0.0j]], dtype=complex)