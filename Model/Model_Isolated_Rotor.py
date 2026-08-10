# Model.py - 压缩系统稳定性分析模型
# 系统配置：上游无限长 → 孤立转子 → 下游无限长

import os
import sys

# 确保 Stability 根目录在 sys.path 中, 使 Matrix/ 等包可导入
_STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _STABILITY_ROOT not in sys.path:
    sys.path.insert(0, _STABILITY_ROOT)

import numpy as np

# 导入组件模块（每个模块自己读取参数）
from Matrix.rotor import B_rot_n_num
from Matrix.axial import Tn_ax_coeff_num, inv_Tn_ax_coeff_num
from Matrix.boundary import IC_upstream_infinite, EC_downstream_infinite


def check_velocity_triangle(params, tol=5e-3, raise_on_fail=True):
    """检查转子速度三角形是否闭合

    闭合条件 (U=1, 无预旋 α₁=0)：
        V_x · tanα₂ = 1 + V_x · tanβ₂
        V_x · tanβ₁ = V_x · tanα₁ - 1   ⇒  tanβ₁ = -1/V_x (当 α₁=0)

    输入参数中的 (V_x, β₁, β₂, α₁, α₂) 必须同时满足上述关系,
    否则数值系统与解析公式不在同一工作点上。

    参数
    ----------
    params : dict
        含 Vx, beta1_deg, beta2_deg, alpha1_deg, alpha2_deg 的参数字典
    tol : float
        允许残差 (Vx · tanα₂ - 1 - Vx · tanβ₂ 的绝对值)
    raise_on_fail : bool
        True 时不闭合报错, False 时仅警告

    返回
    -------
    dict
        {'closed': bool, 'residual_out': float, 'residual_in': float,
         'Vx_correct': float}
    """
    Vx = params['Vx']
    tan_b1 = np.tan(np.deg2rad(params['beta1_deg']))
    tan_b2 = np.tan(np.deg2rad(params['beta2_deg']))
    tan_a1 = np.tan(np.deg2rad(params['alpha1_deg']))
    tan_a2 = np.tan(np.deg2rad(params['alpha2_deg']))

    # 出口三角形: V_x tanα₂ = 1 + V_x tanβ₂
    res_out = Vx * tan_a2 - (1.0 + Vx * tan_b2)
    # 进口三角形: V_x tanβ₁ = V_x tanα₁ - 1
    res_in = Vx * tan_b1 - (Vx * tan_a1 - 1.0)

    Vx_correct_out = 1.0 / (tan_a2 - tan_b2) if abs(tan_a2 - tan_b2) > 1e-12 else float('nan')
    Vx_correct_in = 1.0 / (tan_a1 - tan_b1) if abs(tan_a1 - tan_b1) > 1e-12 else float('nan')

    closed = (abs(res_out) < tol) and (abs(res_in) < tol)

    if not closed:
        msg = (
            "\n[!] 速度三角形不闭合\n"
            f"    当前 V_x = {Vx:.4f}\n"
            f"    出口残差  V_x·tanα₂ - (1 + V_x·tanβ₂) = {res_out:+.4e}\n"
            f"    进口残差  V_x·tanβ₁ - (V_x·tanα₁ - 1) = {res_in:+.4e}\n"
            f"    出口条件自洽 V_x = 1/(tanα₂ - tanβ₂) = {Vx_correct_out:.4f}\n"
            f"    进口条件自洽 V_x = 1/(tanα₁ - tanβ₁) = {Vx_correct_in:.4f}\n"
            "    数值传递矩阵仅使用 (V_x, β₂, α₁, β₁) 推出 V_θ₂,\n"
            "    解析公式同时使用 α₂, 不自洽时两者会代表不同的物理工作点!\n"
        )
        if raise_on_fail:
            raise ValueError(msg)
        else:
            print(msg)
    else:
        print(f"[OK] 速度三角形闭合 (V_x = {Vx:.4f}, 出口残差 {res_out:+.2e}, 进口残差 {res_in:+.2e})")

    return {
        'closed': closed,
        'residual_out': res_out,
        'residual_in': res_in,
        'Vx_correct_out': Vx_correct_out,
        'Vx_correct_in': Vx_correct_in,
    }


def build_system_matrix(s_val, n_val=1, params=None):
    """
    构建完整系统传递矩阵
    
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
        from Matrix.rotor import load_params, compute_rotor_params
        params = load_params()
        params = compute_rotor_params(params)

    Vx = params['Vx']
    mu_up = params.get('mu_up', 1.0)
    mu_dn = params.get('mu_dn', 1.0)

    # 孤立转子设定：上游远场为无旋（Vθ=0）
    Vtheta_up = 0.0
    # 下游平均周向速度：使用转子出口已计算的周向速度 Vθ2
    Vtheta_dn = params['Vtheta_bar2_rot']

    # 管道参考点统一以执行盘为原点 x=0
    # 模态幅值本身是管道全局常数，管道长度不进入公式
    # （"无限长"仅体现为边界条件 A_n=0 / B_n=0，与 L 无关）

    # 转子传递矩阵
    B_rot = B_rot_n_num(s=s_val, n=n_val, params=params)

    # 上游管道：转子进口 x=0⁻ 处的状态 = T_up(0;0) · c^up
    T_up = Tn_ax_coeff_num(x=0, s=s_val, n=n_val,
                           Vx_bar=Vx, Vtheta_bar=Vtheta_up, x0=0, mu=mu_up)

    # 下游管道：转子出口 x=0⁺ 处的状态 = T_dn(0;0) · c^dn
    # 反解：c^dn = T_dn(0;0)⁻¹ · state(0⁺)
    inv_T_dn = inv_Tn_ax_coeff_num(x=0, s=s_val, n=n_val,
                                   Vx_bar=Vx, Vtheta_bar=Vtheta_dn, x0=0, mu=mu_dn)

    # 系统传递矩阵
    X_sys = inv_T_dn @ B_rot @ T_up

    return X_sys


def build_characteristic_matrix(s_val, n_val=1, params=None):
    """构建特征矩阵 Y_sys = [EC·X_sys; IC]

    根据边界条件文档，特征值方程为：
        det(Y_sys) = 0，其中 Y_sys = [EC·X_sys; IC]
    
    上游无限长管道：IC = [[0,1,0],[0,0,1]]  (B_n=0, C_n=0)
    下游无限长管道：EC = [[1,0,0]]            (A_n=0)
    """
    X_sys = build_system_matrix(s_val, n_val, params)
    
    # 获取边界条件矩阵
    IC = IC_upstream_infinite(n_val)  # 2×3
    EC = EC_downstream_infinite(n_val)  # 1×3
    
    # 构造特征矩阵 Y_sys = [EC·X_sys; IC]，维度 3×3
    Y_sys = np.vstack([
        EC @ X_sys,  # 第一行：下游边界条件约束
        IC           # 第二、三行：上游边界条件约束
    ])
    
    return Y_sys


def characteristic_equation(s_val, n_val=1, params=None):
    """计算特征方程的值：det(Y_sys) = det([EC·X_sys; IC])"""
    Y_sys = build_characteristic_matrix(s_val, n_val, params)
    return np.linalg.det(Y_sys)


def make_system_matrix_function(n_val=1):
    """为 solver 提供矩阵接口：预加载参数,避免每次求值重复读文件"""
    from Matrix.rotor import load_params, compute_rotor_params
    params = load_params()
    params = compute_rotor_params(params)
    def func(s):
        return build_characteristic_matrix(s, n_val, params)
    return func


if __name__ == '__main__':

    from Matrix.rotor import load_params, compute_rotor_params
    from Eigen_Value_Solver.Eigen_Hunter import hybrid_hunt

    # 加载参数并检查速度三角形自洽性
    params = compute_rotor_params(load_params())
    check_velocity_triangle(params, raise_on_fail=True)

    print()
    print("  (ʅ'ω'ʅ)  孤立转子 特征值 (Eigen_Hunter)")
    print()
    header = f"  {'n':>3s}  {'σ':>10s}  {'ω':>10s}  {'s':>24s}  {'稳?':>6s}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    for n_val in range(1, 7):
        Y_func = make_system_matrix_function(n_val)
        eigs = hybrid_hunt(Y_func,
                           sigma_range=(-6.0, 2.0),
                           omega_range=(-6.0, 6.0),
                           verbose=False, true_root_tol=1e-10)
        eigs = sorted(eigs, key=lambda x: x.real, reverse=True)[:1]

        if eigs:
            eig = eigs[0]
            sigma = eig.real
            omega = -eig.imag  # 论文约定 s = σ - jω
            s_str = f"{sigma:+.5f}-{omega:+.5f}j" if omega >= 0 else f"{sigma:+.5f}+{-omega:+.5f}j"
            stability = "不稳定" if sigma > 0 else "稳定"
            print(f"  {n_val:3d}  {sigma:+10.5f}  {omega:+10.5f}  {s_str:>24s}  {stability:>6s}")
        else:
            print(f"  {n_val:3d}  {'--':>10s}  {'--':>10s}  {'未找到特征值':>24s}  {'?':>6s}")