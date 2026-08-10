# diffuser.py - 扩压器模块
import numpy as np
import re
import os

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
    return params

def B_dif_n_num(s, n, params=None):
    """
    扩压器执行盘传递矩阵（状态→状态）
    
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
    B_dif : np.ndarray (3×3)
        扩压器传递矩阵
    """
    if params is None:
        params = load_params()
    
    ARp_dif = params['ARp_dif']
    tan_alpha2 = params['tan_alpha2_dif']
    Vr_bar1 = params['Vr_bar1_dif']
    Vtheta_bar1 = params['Vtheta_bar1_dif']
    Vr_bar2 = params['Vr_bar2_dif']
    Vtheta_bar2 = params['Vtheta_bar2_dif']
    dL_dtanalpha1 = params['dL_dtanalpha1_dif']
    lambda_dif = params['lambda_dif']
    tau_dif = params['tau_dif']

    denom = Vr_bar1 * (1 + tau_dif * s)
    tan_alpha1 = Vtheta_bar1 / Vr_bar1
    term1 = (-(lambda_dif * s + Vr_bar2 + Vtheta_bar2 * tan_alpha2) / ARp_dif
             + Vr_bar1 + dL_dtanalpha1 * tan_alpha1 / denom)
    term2 = Vtheta_bar1 - dL_dtanalpha1 / denom
    return np.array([
        [1/ARp_dif, 0, 0],
        [tan_alpha2/ARp_dif, 0, 0],
        [term1, term2, 1]
    ], dtype=complex)
