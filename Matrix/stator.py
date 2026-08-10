# stator.py - 静子模块
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
    return params

def compute_stator_params(params):
    """计算静子派生参数"""
    Vx = params['Vx']
    params['tan_alpha1_sta'] = np.tan(np.deg2rad(params['alpha3_deg']))
    params['tan_alpha2_sta'] = np.tan(np.deg2rad(params['alpha4_deg']))
    params['tan_beta2_rot'] = np.tan(np.deg2rad(params['beta2_deg']))
    params['Vtheta_bar1_sta'] = 1.0 + Vx * params['tan_beta2_rot']
    params['Vtheta_bar2_sta'] = Vx * params['tan_alpha2_sta']
    return params

def B_sta_n_num(s, n, params=None):
    """
    静子执行盘传递矩阵（状态→状态）
    
    参数：
    --------
    s : complex
        拉普拉斯变量
    n : int
        周向波数
    params : dict or None
        参数字典，为None时自动从文件加载
    
    返回：
    --------
    B_sta : np.ndarray (3×3)
        静子传递矩阵
    """
    if params is None:
        params = load_params()
        params = compute_stator_params(params)
    
    tan_alpha2 = params['tan_alpha2_sta']
    tan_alpha1 = params['tan_alpha1_sta']
    Vx_bar = params['Vx']
    Vtheta_bar1 = params['Vtheta_bar1_sta']
    Vtheta_bar2 = params['Vtheta_bar2_sta']
    dL_dtanalpha1 = params['dL_dtanalpha1_sta']
    lambda_sta = params['lambda_sta']
    tau_S = params['tau_S']

    # 损失项（静叶损失迟滞中无 jn 项）
    denom = Vx_bar * (1.0 + s * tau_S)

    B_sta = np.zeros((3, 3), dtype=complex)

    # 速度匹配条件
    B_sta[0, 0] = 1.0
    B_sta[1, 0] = tan_alpha2

    # 压力方程（Spakovszky 静叶执行盘公式）
    B_sta[2, 0] = (-Vtheta_bar2 * tan_alpha2
                   + dL_dtanalpha1 * tan_alpha1 / denom
                   - s * lambda_sta)
    B_sta[2, 1] = Vtheta_bar1 - dL_dtanalpha1 / denom
    B_sta[2, 2] = 1.0

    return B_sta