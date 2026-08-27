# Model_NCEPU.py - NCEPU 离心压气机动态系统模型 (无有叶扩压器)
import numpy as np
from Matrix.axial import Tn_ax_coeff_num
from Matrix.boundary import IC_upstream_infinite, EC_downstream_plenum
from Matrix.radial import Tn_rad_num, T0_rad_num, _T_n as _T_n_radial
from Matrix.impeller import B_imp_n_num, compute_impeller_params
from Matrix.rotor import load_params
from Model.Tool.modal_coeffs import extract_modal_coeffs

# NCEPU 专用配置文件
DEFAULT_PARAMS_FILE = 'NCEPU_params.txt'

def compute_cc3_params(params_file=None):
    """加载并计算 NCEPU 所有组件的平均流参数 (无有叶扩压器)

    参数：
    --------
    params_file : str or None
        参数文件路径，为 None 时使用默认文件

    返回：
    --------
    params : dict
        包含所有计算参数的字典
    """
    params = load_params(params_file or DEFAULT_PARAMS_FILE)

    # ==================== 从 STA 几何数据计算派生参数 ====================
    # 站位约定 (Spakovszky CC3, x=0 在叶轮进口执行盘处):
    #   STA2: 叶轮进口
    #   STA3: 叶轮出口 / 无叶扩压器入口
    #   STA4: 无叶扩压器出口 / 容腔入口
    # 长度以叶轮出口轮缘半径 R_3 归一化, 故 radial_sta3 = 1.0。

    # 叶轮参数
    params['R1_R2'] = params['radial_sta2'] / params['radial_sta3']
    # AR_imp 从参数文件读取 (含密度比 ρ₃A₃/(ρ₂A₂), 由 CFX 加载器计算)
    params['r_impeller_exit'] = params['radial_sta3']

    # 无叶扩压器参数 (从叶轮出口 STA3 到 STA4)
    params['r_diffuser_exit'] = params['radial_sta4']

    # 容腔入口速度 (扩压器出口, 由 Q/G 和 r_diffuser_exit 导出)
    params['Vx_plenum'] = params['Q'] / params['r_diffuser_exit']
    params['Vtheta_plenum'] = params['G'] / params['r_diffuser_exit']

    # 计算派生角度参数
    params['tan_beta1_imp'] = np.tan(np.deg2rad(params['beta1_imp_deg']))
    params['tan_beta2_imp'] = np.tan(np.deg2rad(params['beta2_imp_deg']))
    params['tan_alpha1_imp'] = np.tan(np.deg2rad(params['alpha1_imp_deg']))

    # 调试输出计算结果
    print(f"[NCEPU 参数计算] R1_R2 = {params['R1_R2']:.4f}")
    print(f"[NCEPU 参数计算] AR_imp = {params['AR_imp']:.4f}  (含密度)")
    print(f"[NCEPU 参数计算] lambda_imp = {params['lambda_imp']:.4f}")
    print(f"[NCEPU 参数计算] tau_imp = {params['tau_imp']:.4f}")
    print(f"[NCEPU 参数计算] Vx_plenum = {params['Vx_plenum']:.4f}")
    params = compute_impeller_params(params)
    return params


def build_cc3_system_matrix(s, n, params=None, params_file=None):
    """
    构建 NCEPU 系统传递矩阵 X_sys (上游系数 → 下游系数) (无有叶扩压器)

    串联顺序：
        x ∈ (-∞, 0] : 上游轴向导管 (无涡，无预旋)
        x = 0       : 叶轮执行盘
        r ∈ [r2, r3] : 无叶扩压器 (vaneless diffuser)
        x ∈ [0, L_exit] : 下游轴向导管 (连接容腔)
        plenum boundary condition 在 EC 中处理

    参数：
    --------
    s : complex
        拉普拉斯变量
    n : int
        周向波数
    params : dict or None
        预计算的参数字典, 为 None 时自动从文件加载
    params_file : str or None
        参数文件路径，为 None 时使用默认文件

    返回：
    --------
    X_sys : np.ndarray (3×3)
        NCEPU 系统传递矩阵
    """
    if params is None:
        params = compute_cc3_params(params_file)

    # 1. 上游导管 (参考点取 x=0)
    Vx_up = params['Vx_up']
    Vtheta_up = 0.0
    mu_up = params.get('mu_up', 1.0)
    T_up = Tn_ax_coeff_num(x=0, s=s, n=n, Vx_bar=Vx_up, Vtheta_bar=Vtheta_up, x0=0, mu=mu_up)

    # 2. 叶轮
    B_imp = B_imp_n_num(s, n, params)

    # 3. 无叶扩压器 (径向)
    r_in = params['r_impeller_exit']
    r_out = params['r_diffuser_exit']
    if n == 0:
        B_vls = T0_rad_num(r_out, s, r_in, params)
    else:
        B_vls = Tn_rad_num(r_out, s, n, r_in, params)

    # 4. 容腔入口: 扩压器出口状态 → 系数空间 (L=0, 管道退化)
    #    系数空间用于施加 EC_downstream_plenum 边界条件
    Vx_plenum = params['Vx_plenum']
    Vtheta_plenum = params['Vtheta_plenum']
    mu_plenum = params.get('mu_plenum', 1.0)
    T_plenum = Tn_ax_coeff_num(x=0, s=s, n=n, Vx_bar=Vx_plenum, Vtheta_bar=Vtheta_plenum, x0=0, mu=mu_plenum)
    inv_T_plenum = np.linalg.inv(T_plenum)

    X_sys = inv_T_plenum @ B_vls @ B_imp @ T_up
    return X_sys


def build_cc3_characteristic_matrix(s, n, params=None, params_file=None):
    """构建特征矩阵 Y_sys = [EC·X_sys; IC]

    参数：
    --------
    s : complex
        拉普拉斯变量
    n : int
        周向波数
    params : dict or None
        预计算的参数字典, 为 None 时自动从文件加载
    params_file : str or None
        参数文件路径，为 None 时使用默认文件
    """
    if params is None:
        params = compute_cc3_params(params_file)
    X_sys = build_cc3_system_matrix(s, n, params)

    # 边界条件：
    # 上游：无限长管道，无反射 (B_n=0, C_n=0)
    IC = IC_upstream_infinite(n)   # shape (2,3)
    # 下游：容腔边界条件 (L=0, 直接在扩压器出口施加)
    Vx_plenum = params['Vx_plenum']
    Vtheta_plenum = params['Vtheta_plenum']
    mu_plenum = params.get('mu_plenum', 1.0)
    EC = EC_downstream_plenum(0.0, s, n, Vx_plenum, Vtheta_plenum, mu_plenum)

    Y_sys = np.vstack([EC @ X_sys, IC])
    return Y_sys


def characteristic_equation(s, n, params=None):
    Y = build_cc3_characteristic_matrix(s, n, params)
    return np.linalg.det(Y)


def make_cc3_matrix_function(n=1, params_file=None):
    """为 solver 提供矩阵接口：预加载参数,避免每次求值重复读文件

    参数：
    --------
    n : int
        周向波数
    params_file : str or None
        参数文件路径，为 None 时使用默认文件
    """
    params = compute_cc3_params(params_file)
    def func(s):
        return build_cc3_characteristic_matrix(s, n, params)
    return func


# ===================== 模态空间分布 =====================

def compute_modal_shape(s_star, n, params, n_points=50, x_upstream_max=0.2):
    """计算给定特征值 s* 在物理空间的模态形状分布

    物理空间分为三段 (归一化坐标, R_3=0.2m):
      1. 上游轴向管道: x ∈ [-x_upstream_max, 0]   (x=0 = 叶轮进口 STA2)
      2. 叶轮: 执行盘 (零厚度, x=0, 仅状态跳变)
      3. 无叶扩压器: r ∈ [r_impeller_exit, r_diffuser_exit] = [1.0, 1.69]

    状态向量:
      上游/下游管道: [δVx, δVθ, δP]
      扩压器:        [δVr, δVθ, δP]

    返回
    ----
    result : dict
        'segments' : list of dict, 每段含:
            'name'      : 段名 ('Upstream' / 'Impeller' / 'Diffuser')
            'coord'     : np.ndarray, 物理坐标 (归一化)
            'q'         : np.ndarray (N, 3), 复数状态向量
            'amplitude' : np.ndarray (N,), |δV| (速度扰动幅值)
            'pressure'  : np.ndarray (N,), |δP| (压力扰动幅值)
        'eig' : complex, 特征值
        'n'   : int, 周向波数
    """
    # ── 1. SVD 提取上游系数 [A_n, B_n, C_n] (强制上游无限长管道 BC) ──
    # 使用通用工具: SVD 提取零空间向量, 并强制 B_n=C_n=0
    # (上游无限长管道无反射; 消除指数放大导致的远端失真, 不影响特征值 s*)
    # 详见 Model.Tool.modal_coeffs
    Y = build_cc3_characteristic_matrix(s_star, n, params)
    coeffs, _diag = extract_modal_coeffs(Y, force_upstream_bc=True)

    # ── 公共参数 ──
    Vx_up = params['Vx_up']
    Vtheta_up = 0.0
    mu_up = params.get('mu_up', 1.0)
    r_in = params['r_impeller_exit']       # = radial_sta3 = 1.0
    r_out = params['r_diffuser_exit']      # = radial_sta4 = 1.69
    Q = params['Q']
    G = params['G']

    segments = []

    # ── 2. 上游管道: x ∈ [-x_upstream_max, 0] ──
    x_up = np.linspace(-x_upstream_max, 0.0, n_points)
    q_up = np.zeros((n_points, 3), dtype=complex)
    for i, xi in enumerate(x_up):
        M = Tn_ax_coeff_num(x=xi, s=s_star, n=n, Vx_bar=Vx_up,
                            Vtheta_bar=Vtheta_up, x0=0, mu=mu_up)
        q_up[i] = M @ coeffs
    amp_up = np.sqrt(np.abs(q_up[:, 0])**2 + np.abs(q_up[:, 1])**2)
    segments.append({
        'name': 'Upstream',
        'coord': x_up,
        'q': q_up,
        'amplitude': amp_up,
        'pressure': np.abs(q_up[:, 2]),
    })

    # ── 3. 叶轮: x=0 执行盘, 状态跳变 ──
    # 上游 x=0 处状态 → 叶轮进口状态 → B_imp → 叶轮出口状态
    M_up0 = Tn_ax_coeff_num(x=0.0, s=s_star, n=n, Vx_bar=Vx_up,
                            Vtheta_bar=Vtheta_up, x0=0, mu=mu_up)
    q_in_imp = M_up0 @ coeffs                 # 叶轮进口状态 [δVx, δVθ, δP]
    B_imp = B_imp_n_num(s_star, n, params)
    q_out_imp = B_imp @ q_in_imp              # 叶轮出口状态 [δVr, δVθ, δP]
    # 叶轮段用单点表示 (零厚度)
    amp_imp = np.array([np.sqrt(np.abs(q_out_imp[0])**2 + np.abs(q_out_imp[1])**2)])
    segments.append({
        'name': 'Impeller',
        'coord': np.array([0.0]),
        'q': q_out_imp.reshape(1, 3),
        'amplitude': amp_imp,
        'pressure': np.array([np.abs(q_out_imp[2])]),
    })

    # ── 4. 无叶扩压器: r ∈ [r_in, r_out] ──
    r_diff = np.linspace(r_in, r_out, n_points)
    q_diff = np.zeros((n_points, 3), dtype=complex)
    # 用 _T_n 在每个 r 处直接求状态 (以 r_in 为参考)
    for i, ri in enumerate(r_diff):
        T_ri = _T_n_radial(ri, s_star, n, r_in, Q, G)
        T_in = _T_n_radial(r_in, s_star, n, r_in, Q, G)
        T_transfer = T_ri @ np.linalg.inv(T_in)
        q_diff[i] = T_transfer @ q_out_imp
    amp_diff = np.sqrt(np.abs(q_diff[:, 0])**2 + np.abs(q_diff[:, 1])**2)
    segments.append({
        'name': 'Diffuser',
        'coord': r_diff,
        'q': q_diff,
        'amplitude': amp_diff,
        'pressure': np.abs(q_diff[:, 2]),
    })

    return {
        'segments': segments,
        'eig': s_star,
        'n': n,
    }