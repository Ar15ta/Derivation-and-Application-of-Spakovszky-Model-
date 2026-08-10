import numpy as np

def isolated_rotor_eigenvalue(params, n=1):
    """
    计算孤立转子特征值的解析解 
    
    参数
    ----------
    params : dict
        必须包含以下键值：
        Vx           : 流量系数（无量纲轴向速度）
        tan_beta1    : 转子相对进口气流角的正切值
        tan_beta2    : 转子相对出口气流角的正切值
        tan_alpha2   : 转子绝对出口气流角的正切值
        dL_dtanbeta1 : 转子稳态损失对 tan(beta1) 的导数
        lambda_rot   : 转子叶片排惯性参数
    n : int
        周向谐波数（n≥1）
    
    返回
    -------
    complex
        特征值 s = σ - jω
        其中 σ 为增长率，ω 为旋转率。
        沿用论文 s = σ - jω 约定，与数值模型（Beyn 求解器）保持一致。
    """
    Vx = params['Vx']
    tanb1 = params['tan_beta1']
    tanb2 = params['tan_beta2']
    tana2 = params['tan_alpha2']
    dL = params['dL_dtanbeta1_rot']
    lam = params['lambda_rot']
    
    denom = lam + 2.0 / n
    sigma = (tanb2 + (dL * tanb1 / Vx) - Vx * (1 + tana2**2) + tana2) / denom
    omega = ((dL / Vx) + n * lam + 1) / denom
    
    # 论文定义 s = σ - jω，数值模型沿用此约定
    return sigma - 1j * omega


# ========== 使用示例（基于论文表5.1数据） ==========
if __name__ == '__main__':
    # 从 system_params.txt 读取参数并构造解析解所需的键
    # 注意: 参数文件中给的是 dL/dphi, 需要转换为 dL/d(tan beta1) = dL/dphi / Vx
    try:
        from rotor import load_params as load_rotor_params, compute_rotor_params
        raw = compute_rotor_params(load_rotor_params())
    except Exception:
        # 回退到 stator.load_params 如果 rotor.load_params 不可用
        from stator import load_params as load_rotor_params
        raw = load_rotor_params()
        raw['dL_dtanbeta1_rot'] = raw['dL_dtanphi_rot'] * 1/((params['tan_alpha1_rot']-params['tan_beta1_rot'])**2)

    Vx = raw['Vx']
    tan_beta1 = np.tan(np.deg2rad(raw['beta1_deg']))
    tan_beta2 = np.tan(np.deg2rad(raw['beta2_deg']))
    tan_alpha2 = np.tan(np.deg2rad(raw['alpha2_deg']))
    dL_dtanbeta1_rot = raw['dL_dtanbeta1_rot']
    lambda_rot = raw['lambda_rot']

    params = {
        'Vx': Vx,
        'tan_beta1': tan_beta1,
        'tan_beta2': tan_beta2,
        'tan_alpha2': tan_alpha2,
        'dL_dtanbeta1_rot': dL_dtanbeta1_rot,
        'lambda_rot': lambda_rot
    }

    print()
    print("  (ฅ'ω'ฅ)  孤立转子 解析特征值")
    print()
    header = f"  {'n':>3s}  {'σ':>10s}  {'ω':>10s}  {'s':>24s}  {'稳?':>6s}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for n in [1, 2, 3, 4, 5, 6, 7, 8]:
        eig = isolated_rotor_eigenvalue(params, n)
        sigma = eig.real
        omega = -eig.imag  # 论文约定 s = σ - jω，与 Model_* 系列一致
        s_str = f"{sigma:+.5f}-{omega:+.5f}j" if omega >= 0 else f"{sigma:+.5f}+{-omega:+.5f}j"
        stability = "不稳定" if sigma > 0 else "稳定"
        print(f"  {n:3d}  {sigma:+10.5f}  {omega:+10.5f}  {s_str:>24s}  {stability:>6s}")