# radial.py - 径向间隙模块
import numpy as np
import re
import os
import matplotlib.pyplot as plt

def load_params(filename='Centrigugal_Params_Generated.txt'):
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


def rad_fun(r, s, n, r0, Q, G):
    """
    Ref 公式：计算 Rn(r), dRn/dr, d2Rn/dr2（梯形法数值积分）。

    辅助积分定义：
      Fn = ∫_{r0}^{r} exp(-sξ²/(2Q) - jnG/Q·ln(ξ)) · ξ^{-n+1} dξ
      Fp = ∫_{r0}^{r} exp(-sξ²/(2Q) - jnG/Q·ln(ξ)) · ξ^{+n+1} dξ

    Rn = r^n · Fn - r^{-n} · Fp
    dRn = n·r^{n-1}·Fn + r^n·fn(r) + n·r^{-n-1}·Fp - r^{-n}·fp(r)
    d2Rn = (n²-n)·r^{n-2}·Fn + 2n·r^{n-1}·fn(r) + r^n·dfn(r)
          - (n²+n)·r^{-n-2}·Fp + 2n·r^{-n-1}·fp(r) - r^{-n}·dfp(r)
    """
    N = 2000  # 积分点数；提高精度以降低 r≈1.2 处静压计算中的数值消失风险
    xi = np.linspace(r0, r, N)
    dxi = (r - r0) / (N - 1)

    exp_arg = -s/(2*Q) * xi**2 - 1j * n * G / Q * np.log(xi)
    exp_vals = np.exp(exp_arg)

    fn_vals = exp_vals * xi**(-n+1)
    fp_vals = exp_vals * xi**(+n+1)

    # 梯形积分 (numpy >= 2.0: trapezoid, 旧版本: trapz)
    _trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
    Fn_int = _trapz(fn_vals, xi)
    Fp_int = _trapz(fp_vals, xi)

    # 末端值（ξ = r）
    exp_r = np.exp(-s/(2*Q) * r**2 - 1j * n * G / Q * np.log(r))
    fn_r = exp_r * r**(-n+1)
    fp_r = exp_r * r**(+n+1)

    # Rn
    Rn = r**n * Fn_int - r**(-n) * Fp_int

    # 一阶导（含 Leibniz 边界项）
    dRn = (n * r**(n-1) * Fn_int + r**n * fn_r
           + n * r**(-n-1) * Fp_int - r**(-n) * fp_r)

    # 被积函数在 r 处的导数（链式法则）
    dphi_dr = -s * r / Q - 1j * n * G / (Q * r)
    dfn_r = exp_r * (dphi_dr * r**(-n+1) + (-n+1) * r**(-n))
    dfp_r = exp_r * (dphi_dr * r**(+n+1) + (+n+1) * r**(+n))

    # 二阶导
    d2Rn = ((n**2 - n) * r**(n-2) * Fn_int
            + 2*n * r**(n-1) * fn_r
            + r**n * dfn_r
            - (n**2 + n) * r**(-n-2) * Fp_int
            + 2*n * r**(-n-1) * fp_r
            - r**(-n) * dfp_r)

    return Rn, dRn, d2Rn


def _T_n(r, s, n, r0, Q, G):
    """
    Ref 公式 (2.64)：径向空间 n>0 谐波系数→状态矩阵 (3×3)。
    映射 [D_n, E_n, F_n]ᵀ → [δṼ_r, δṼ_θ, δP̃]ᵀ。
    """
    Rn, dRn, d2Rn = rad_fun(r, s, n, r0, Q, G)
    T = np.zeros((3, 3), dtype=complex)
    T[0, 0] = 1j * n * r**(n-1)
    T[0, 1] = 1j * n * r**(-n-1)
    T[0, 2] = 1j * n * Rn / r
    T[1, 0] = -n * r**(n-1)
    T[1, 1] = n * r**(-n-1)
    T[1, 2] = -dRn
    T[2, 0] = -1j * s * r**n - 1j * n * (Q + 1j*G) * r**(n-2)
    T[2, 1] = 1j * s * r**(-n) - 1j * n * (Q - 1j*G) * r**(-n-2)
    T[2, 2] = (Q * d2Rn + (s * r + 1j * n * G / r + Q / r) * dRn) / (1j * n)
    return T


def Tn_rad_num(r, s, n, r0, params=None):
    """
    径向无叶空间 n>0 状态→状态传递矩阵 (3×3)。

    状态向量 [δṼ_r, δṼ_θ, δP̃]ᵀ 从进口(r0)到出口(r)。
    M = T(r) · inv(T(r0))
    
    参数：
    --------
    r : float
        出口半径
    s : complex
        拉普拉斯变量
    n : int
        周向波数
    r0 : float
        进口半径
    params : dict or None
        参数字典，为None时自动从文件加载
    
    返回：
    --------
    Tn_rad : np.ndarray (3×3)
        径向传递矩阵
    """
    if params is None:
        params = load_params()
    
    Q = params['Q']
    G = params['G']
    
    T_in = _T_n(r0, s, n, r0, Q, G)
    T_out = _T_n(r,  s, n, r0, Q, G)
    return T_out @ np.linalg.inv(T_in)


def T0_rad_num(r, s, r0, params=None):
    """
    径向空间 n=0 轴对称系数→状态矩阵 (3×2)。
    映射 [D₀, F₀]ᵀ → [δṼ_r, δṼ_θ, δP̃]ᵀ。
    
    参数：
    --------
    r : float
        出口半径
    s : complex
        拉普拉斯变量
    r0 : float
        进口半径
    params : dict or None
        参数字典，为None时自动从文件加载
    
    返回：
    --------
    T0_rad : np.ndarray (3×2)
        径向传递矩阵（n=0）
    """
    if params is None:
        params = load_params()
    
    Q = params['Q']
    G = params['G']
    # 辅助积分 J0 = ∫ exp(-sξ²/(2Q))·ξ dξ (解析)
    J0 = (Q / s) * (np.exp(-s/(2*Q) * r0**2) - np.exp(-s/(2*Q) * r**2))

    # R0 = ∫ ξ^(-3)·exp(-sξ²/(2Q)) dξ (数值)
    N = 2000
    xi = np.linspace(r0, r, N)
    integrand_R0 = np.exp(-s/(2*Q) * xi**2) * xi**(-3)
    _trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz
    R0 = _trapz(integrand_R0, xi)

    T = np.zeros((3, 2), dtype=complex)
    T[0, 0] = 1.0 / r
    T[0, 1] = 0.0
    T[1, 0] = 0.0
    T[1, 1] = (Q / (r * s)) * np.exp(-s/(2*Q) * r**2)
    T[2, 0] = -Q / r**2 - s * np.log(r)
    T[2, 1] = (2 * G / (s * Q)) * R0
    return T


# 保留旧函数别名以保证向后兼容
rad_fun_numeric = rad_fun