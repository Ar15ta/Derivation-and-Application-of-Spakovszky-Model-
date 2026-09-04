# impeller.py - 叶轮模块
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
            else:
                # 支持 [] 占位符, 解析为 0.0 (待校准参数)
                match_placeholder = re.match(r'(\w+)\s*=\s*\[\]', line)
                if match_placeholder:
                    params[match_placeholder.group(1)] = 0.0
    return params

def compute_impeller_params(params):
    """计算叶轮派生参数"""
    params['tan_beta1_imp'] = np.tan(np.deg2rad(params['beta1_imp_deg']))
    params['tan_beta2_imp'] = np.tan(np.deg2rad(params['beta2_imp_deg']))
    params['tan_alpha1_imp'] = np.tan(np.deg2rad(params['alpha1_imp_deg']))
    return params

def B_imp_n_num(s, n, params=None):
    """
    叶轮执行盘传递矩阵（状态→状态）
    
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
    B_imp : np.ndarray (3×3)
        叶轮传递矩阵
    """
    if params is None:
        params = load_params()
        params = compute_impeller_params(params)
    
    AR_imp = params['AR_imp']
    R1_R2 = params['R1_R2']
    tan_beta1 = params['tan_beta1_imp']
    tan_beta2 = params['tan_beta2_imp']
    tan_alpha1 = params['tan_alpha1_imp']
    Vx_bar1 = params['Vx_bar1_imp']
    Vr_bar2 = params['Vr_bar2_imp']
    Vtheta_bar1 = params['Vtheta_bar1_imp']
    Vtheta_bar2 = params['Vtheta_bar2_imp']
    lambda_imp = params['lambda_imp']
    tau_imp = params['tau_imp']
    dL_dtanbeta1 = params['dL_dtanbeta1_imp']

    denom = 1.0 + tau_imp * (s + 1j * n)
    K_imp = dL_dtanbeta1 / (Vx_bar1 * denom)

    B31 = (tan_beta2 / AR_imp
           - (Vr_bar2 + Vtheta_bar2 * tan_beta2) / AR_imp
           - lambda_imp * (s + 1j * n) / AR_imp
           - R1_R2 * tan_alpha1
           + Vx_bar1
           + K_imp * tan_beta1)

    B32 = Vtheta_bar1 - K_imp

    return np.array([
        [1.0 / AR_imp, 0.0, 0.0],
        [tan_beta2 / AR_imp, 0.0, 0.0],
        [B31, B32, 1.0]
    ], dtype=complex)
