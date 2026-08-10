"""
模态系数提取工具 (传递矩阵法后处理)

提供统一的 SVD → coeffs 提取接口, 解决"上游无限长管道边界条件数值残余"
问题. 适用于所有基于 Spakovszky CC3 框架的稳定性模型.

背景
----
传递矩阵法求解特征值 s* 后, 通过 SVD 提取零空间向量得到系数 [A_n, B_n, C_n].
由于 s* 是数值近似根 (det(Y)≈0 但不严格为零), SVD 零空间会混入 ~1e-7 量级
的 B_n, C_n 残余. 这些残余本身很小, 但在上游远端会被指数函数 e^{μn|x|}
放大 (n=4, μ=1.92, |x|=3 时放大 1e10 倍), 严重污染模态形状.

物理上, 上游无限长管道无反射, B_n=C_n=0 应严格成立. 本工具在 SVD 提取后
显式强制该边界条件, 保证模态形状物理正确, 且不影响特征值 (已在求根阶段确定).
"""

import numpy as np


def extract_modal_coeffs(Y, force_upstream_bc=True, normalize=True):
    """从特征矩阵 Y 的零空间提取模态系数 [A_n, B_n, C_n]

    参数
    ----
    Y : ndarray (M, 3) complex
        特征矩阵 Y_sys = [EC·X_sys; IC] (含上游/下游边界条件)
        M 通常 = 3 (2 行 IC + 1 行 EC)
    force_upstream_bc : bool, 默认 True
        是否强制施加上游无限长管道边界条件 (B_n=0, C_n=0).
        - True: 物理正确, 消除指数放大导致的远端失真 (推荐)
        - False: 保留 SVD 原始结果 (仅用于诊断对比)
    normalize : bool, 默认 True
        是否将 coeffs 归一化为单位范数 (|coeffs|=1)

    返回
    ----
    coeffs : ndarray (3,) complex
        模态系数 [A_n, B_n, C_n], 其中
        - A_n: 向上游衰减的势模态 (e^{μn x}, x<0 衰减)
        - B_n: 向叶轮衰减的势模态 (e^{-μn x}, 物理上=0)
        - C_n: 涡量对流模态 (e^{-k_n x}, 物理上=0)
    diag : dict
        诊断信息 (供调试/绘图标注):
        'singular_values' : ndarray, SVD 奇异值 (降序)
        'sigma_min_ratio' : float, σ_min/σ_max (越小越奇异)
        'B_over_A_raw'    : float, 强制清零前的 |B_n/A_n|
        'C_over_A_raw'    : float, 强制清零前的 |C_n/A_n|
        'forced_bc'       : bool, 是否强制了边界条件

    使用示例
    --------
    >>> from Model.Tool.modal_coeffs import extract_modal_coeffs
    >>> Y = build_cc3_characteristic_matrix(s_star, n, params)
    >>> coeffs, diag = extract_modal_coeffs(Y, force_upstream_bc=True)
    >>> q = M(x) @ coeffs   # 计算位置 x 处的状态向量
    """
    Y = np.asarray(Y, dtype=complex)
    if Y.shape[-1] != 3:
        raise ValueError(f"Y 最后一维应为 3 (对应 [A_n,B_n,C_n]), 实际 {Y.shape}")

    # ── SVD 提取零空间向量 (最小奇异值对应的右奇异向量) ──
    U, S, Vh = np.linalg.svd(Y)
    coeffs = Vh[-1].conj()

    # ── 诊断: 强制清零前的残余量 ──
    A_raw = coeffs[0]
    B_over_A_raw = abs(coeffs[1] / A_raw) if abs(A_raw) > 0 else float('inf')
    C_over_A_raw = abs(coeffs[2] / A_raw) if abs(A_raw) > 0 else float('inf')
    sigma_min_ratio = float(S[-1] / S[0]) if S[0] > 0 else float('inf')

    # ── 归一化 ──
    if normalize:
        nrm = np.linalg.norm(coeffs)
        if nrm > 0:
            coeffs = coeffs / nrm

    # ── 强制上游无限长管道边界条件: B_n = 0, C_n = 0 ──
    if force_upstream_bc:
        coeffs = coeffs.copy()   # 避免修改 SVD 内部缓冲
        coeffs[1] = 0.0          # B_n = 0 (无反射势模态)
        coeffs[2] = 0.0          # C_n = 0 (无涡量对流模态)
        # 清零后重新归一化, 使 |coeffs|=1 (仅 A_n 非零 → |A_n|=1)
        if normalize:
            nrm = np.linalg.norm(coeffs)
            if nrm > 0:
                coeffs = coeffs / nrm

    diag = {
        'singular_values': S,
        'sigma_min_ratio': sigma_min_ratio,
        'B_over_A_raw': B_over_A_raw,
        'C_over_A_raw': C_over_A_raw,
        'forced_bc': bool(force_upstream_bc),
    }
    return coeffs, diag


def upstream_amplification_factor(mu, n, x_max):
    """计算 B_n 残余在上游远端的指数放大倍数

    用于评估是否需要强制边界条件. 经验阈值:
        放大倍数 × |B_n/A_n| > 1e-3  →  建议强制清零

    参数
    ----
    mu : float
        参考半径比 R_ref/R_duct (轴流=1, 离心=R2/R_duct)
    n : int
        周向波数
    x_max : float
        上游管道归一化长度 (|x| 最大值)

    返回
    ----
    factor : float
        e^{μn·|x_max|}, B_n 项在远端相对入口的放大倍数
    """
    return float(np.exp(mu * n * x_max))


def should_force_bc(mu, n, x_max, B_over_A_raw, threshold=1e-3):
    """判断是否需要强制上游边界条件

    参数
    ----
    mu, n, x_max : 同 upstream_amplification_factor
    B_over_A_raw : float, SVD 原始 |B_n/A_n| (从 diag 获取)
    threshold : float, 默认 1e-3
        远端残余相对幅值阈值, 超过则建议强制清零

    返回
    ----
    needed : bool
        True 表示需要强制清零
    residual_at_far_end : float
        远端 B_n 残余相对幅值 = |B_n/A_n| × e^{μn|x_max|}
    """
    factor = upstream_amplification_factor(mu, n, x_max)
    residual = abs(B_over_A_raw) * factor
    return residual > threshold, residual


# ===================== 自检 =====================

if __name__ == '__main__':
    print("=" * 60)
    print("modal_coeffs.py 自检")
    print("=" * 60)

    # 构造一个近似奇异矩阵 (真实场景模拟)
    # 真实零空间方向: [1, 1e-7, 3e-8] (A_n 主导, B_n/C_n 微小残余)
    np.random.seed(42)
    true_coeffs = np.array([1.0, 1e-7, 3e-8], dtype=complex)
    true_coeffs /= np.linalg.norm(true_coeffs)

    # 构造 Y 使其零空间近似为 true_coeffs
    # Y @ true_coeffs ≈ 0
    M = 3
    rand = np.random.randn(M, 3) + 1j * np.random.randn(M, 3)
    Y = rand - np.outer(rand @ true_coeffs, true_coeffs.conj())
    # 加一点噪声模拟数值误差
    Y += 1e-10 * (np.random.randn(M, 3) + 1j * np.random.randn(M, 3))

    print("\n[测试 1] 不强制边界条件 (诊断对比)")
    coeffs_raw, diag = extract_modal_coeffs(Y, force_upstream_bc=False)
    print(f"  σ_min/σ_max = {diag['sigma_min_ratio']:.3e}")
    print(f"  |B_n/A_n| = {diag['B_over_A_raw']:.3e}")
    print(f"  |C_n/A_n| = {diag['C_over_A_raw']:.3e}")
    print(f"  coeffs = {coeffs_raw}")

    print("\n[测试 2] 强制边界条件 (默认)")
    coeffs, diag = extract_modal_coeffs(Y, force_upstream_bc=True)
    print(f"  coeffs = {coeffs}")
    print(f"  |A_n| = {abs(coeffs[0]):.6f}, |B_n| = {abs(coeffs[1]):.6e}, |C_n| = {abs(coeffs[2]):.6e}")

    print("\n[测试 3] 放大倍数评估 (Eckardt O: μ=1.923, x_max=3)")
    for n in [1, 2, 3, 4]:
        factor = upstream_amplification_factor(1.923, n, 3.0)
        needed, residual = should_force_bc(1.923, n, 3.0, 1e-7)
        flag = "需强制" if needed else "可忽略"
        print(f"  n={n}: 放大={factor:.3e}, 远端残余={residual:.3e} → {flag}")
