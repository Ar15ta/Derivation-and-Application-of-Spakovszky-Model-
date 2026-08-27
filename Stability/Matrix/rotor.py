# rotor.py - 转子模块
import numpy as np
import re

def load_params(filename='system_params.txt'):
    """加载参数文件，支持相对路径和绝对路径

    默认从 Geo/ 目录查找参数文件
    """
    import os
    params = {}
    # Geo/ 目录位于 Stability 根目录下
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
                match_inf = re.match(r'(\w+)\s*=\s*(inf|-inf)', line, re.IGNORECASE)
                if match_inf:
                    params[match_inf.group(1)] = float(match_inf.group(2))
                else:
                    match_placeholder = re.match(r'(\w+)\s*=\s*\[\]', line)
                    if match_placeholder:
                        params[match_placeholder.group(1)] = 0.0
    return params

def compute_rotor_params(params):
    """计算转子派生参数

    损失斜率识别:
      - 若 params 中已有 dL_dtanbeta1_rot → 直接使用 (不对 phi 导数做转换)
      - 若只有 dL_dtanphi_rot → 按论文公式转换为对 tan(beta1) 的导数:
            dL/d(tan beta1) = dL/d(v_x1) * 1/(tan(alpha1) - tan(beta1))^2
    """
    Vx = params['Vx']
    params['tan_beta1_rot'] = np.tan(np.deg2rad(params['beta1_deg']))
    params['tan_alpha1_rot'] = np.tan(np.deg2rad(params['alpha1_deg']))
    params['tan_beta2_rot'] = np.tan(np.deg2rad(params['beta2_deg']))
    params['Vtheta_bar1_rot'] = Vx * params['tan_alpha1_rot']
    params['Vtheta_bar2_rot'] = 1.0 + Vx * params['tan_beta2_rot']

    if 'dL_dtanbeta1_rot' not in params:
        tan_diff_sq = (params['tan_alpha1_rot'] - params['tan_beta1_rot'])**2
        params['dL_dtanbeta1_rot'] = params['dL_dtanphi_rot'] / tan_diff_sq
    return params

def B_rot_n_num(s, n, params=None):
    """
    转子执行盘传递矩阵（状态→状态）
    
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
    B_rot : np.ndarray (3×3)
        转子传递矩阵
    """
    if params is None:
        params = load_params()
        params = compute_rotor_params(params)
    
    Vx_bar = params['Vx']
    Vtheta_bar1 = params['Vtheta_bar1_rot']
    Vtheta_bar2 = params['Vtheta_bar2_rot']
    tan_beta1 = params['tan_beta1_rot']
    tan_beta2 = params['tan_beta2_rot']
    tan_alpha1 = params['tan_alpha1_rot']
    lambda_rot = params['lambda_rot']
    tau_R = params['tau_R']
    if 'dL_dtanbeta1_rot' not in params:
        tan_diff_sq = (params['tan_alpha1_rot'] - params['tan_beta1_rot'])**2
        params['dL_dtanbeta1_rot'] = params['dL_dtanphi_rot'] / tan_diff_sq
    dL_dtanbeta1 = params['dL_dtanbeta1_rot']

    # 损失项 K（含损失迟滞）
    denom = 1.0 + tau_R * (s + 1j * n)
    K = dL_dtanbeta1 / (Vx_bar * denom)

    B_rot = np.zeros((3, 3), dtype=complex)

    # 速度匹配条件
    B_rot[0, 0] = 1.0
    B_rot[1, 0] = tan_beta2

    # 压力方程（Spakovszky 转子执行盘公式）
    B_rot[2, 0] = (tan_beta2 - tan_alpha1
                   - lambda_rot * (s + 1j * n)
                   - Vtheta_bar2 * tan_beta2
                   + K * tan_beta1)
    B_rot[2, 1] = Vtheta_bar1 - K
    B_rot[2, 2] = 1.0

    return B_rot