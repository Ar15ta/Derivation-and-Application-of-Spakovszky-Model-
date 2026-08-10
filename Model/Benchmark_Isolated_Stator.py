import numpy as np

def isolated_stator_eigenvalue(params, n=1):
    """
    计算孤立静子叶片排的特征值（解析解）

    参数
    ----------
    params : dict
        必须包含以下键值：
        Vx           : 流量系数（无量纲轴向速度）
        tan_alpha1   : 静子进口气流角的正切值
        dL_dtanalpha1: 静子稳态损失对 tan(alpha1) 的导数
        lambda_sta   : 静子叶片排惯性参数
    n : int
        周向谐波数（n ≥ 1）

    返回
    -------
    complex
        特征值 s = σ - jω
        其中 σ 为增长率，ω 为旋转率。
        沿用论文 s = σ - jω 约定，与数值模型（Beyn 求解器）保持一致。
    """
    Vx = params['Vx']
    tana1 = params['tan_alpha1']
    dL = params['dL_dtanalpha1_sta']
    lam = params['lambda_sta']

    denom = lam + 2.0 / n
    sigma = (dL * tana1 / Vx - Vx) / denom
    omega = (dL / Vx) / denom

    # 论文定义 s = σ - jω，数值模型沿用此约定
    # 返回 s = σ - jω 以与数值模型保持一致
    return sigma - 1j * omega


# ========== 使用示例（基于论文表5.1数据） ==========
if __name__ == '__main__':
    # 从 system_params.txt 读取参数并构造解析解所需的键
    from stator import load_params

    raw = load_params()
    Vx = raw['Vx']
    # 静子入口角存在于 system_params.txt 中为 alpha3_deg
    tan_alpha1 = np.tan(np.deg2rad(raw['alpha3_deg']))
    dL_dtanalpha1_sta = raw['dL_dtanalpha1_sta']
    lambda_sta = raw['lambda_sta']

    params = {
        'Vx': Vx,
        'tan_alpha1': tan_alpha1,
        'dL_dtanalpha1_sta': dL_dtanalpha1_sta,
        'lambda_sta': lambda_sta
    }

    print()
    print("  (=^･ω･^=)  孤立静叶 解析特征值")
    print()
    header = f"  {'n':>3s}  {'σ':>10s}  {'ω':>10s}  {'s':>24s}  {'稳?':>6s}"
    print(header)
    print("  " + "-" * (len(header) - 2))
    for n in [1, 2, 3, 4, 5, 6, 7, 8]:
        eig = isolated_stator_eigenvalue(params, n)
        sigma, omega = eig.real, eig.imag
        s_str = f"{sigma:+.5f}{omega:+.5f}j"
        stability = "不稳定" if sigma > 0 else "稳定"
        print(f"  {n:3d}  {sigma:+10.5f}  {omega:+10.5f}  {s_str:>24s}  {stability:>6s}")