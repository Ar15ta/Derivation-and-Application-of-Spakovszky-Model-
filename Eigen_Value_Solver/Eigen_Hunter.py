"""
Eigen_Hunter.py - 鲁棒的非线性特征值搜索工具

针对压缩系统稳定性分析中"漏检"问题，提供两种互补的全局搜索策略：

  1. Contour 法 (自适应递归围道分割 + Beyn 算法)
     使用 argument principle 计算围道内特征值个数 N_in，
     若 N_in 较大则将围道四分递归，直到每个子围道恰好包裹少数特征值。
     —— 擅长远场强阻尼、密集模态群的批量提取。

  2. Root Locus 法 (Re/Im det = 0 零水平集求交)
     在 (σ, ω) 网格上同时追踪 det Y(s) 的实部零线与虚部零线，
     两线交点即为候选根，再做 Newton 精化。
     —— 擅长原点附近弱奇点 (如 MG 模态)；几何鲁棒、可视化直观。

hybrid_hunt 同时跑两种方法并合并去重，并用最小奇异值 σ_min/σ_max 过滤伪根。

使用示例：
    from Eigen_Hunter import hybrid_hunt
    from Model_Stage import make_system_matrix_function
    Y = make_system_matrix_function(n_val=1)
    eigs = hybrid_hunt(Y, sigma_range=(-6,2), omega_range=(-6,6))
"""

import numpy as np
import warnings
from scipy.linalg import svd

warnings.filterwarnings('ignore', category=RuntimeWarning)


# ===================== 公共工具 =====================

def _log_det_integral(Y_func, contour_pts):
    """
    沿闭合围道计算  N = (1/2πj) ∮ d/ds[log det Y(s)] ds
    使用 log det 的差分而非直接 d/ds(det)，避免行列式溢出。
    contour_pts : 围道节点序列 (必须闭合，最后一点等于第一点)
    返回围道内特征值个数 (近似整数)。
    """
    log_dets = []
    for s in contour_pts:
        sign, logdet = np.linalg.slogdet(Y_func(s))
        log_dets.append(logdet + np.log(sign + 0j))
    # 解决 log 的 2πj 分支跳变
    diffs = np.diff(np.unwrap(np.array(log_dets).imag)) \
            + 1j * 0  # 仅取虚部累加 → 等价 argument principle
    # 等价做法：累计 arg(det) 变化
    args = np.unwrap(np.angle([np.linalg.det(Y_func(s)) for s in contour_pts]))
    N = (args[-1] - args[0]) / (2 * np.pi)
    return N


def _make_rect_contour(sigma_range, omega_range, n_per_side=64):
    """生成矩形围道节点 (闭合)"""
    s0, s1 = sigma_range
    w0, w1 = omega_range
    bot = [complex(s0 + (s1 - s0) * k / n_per_side, w0) for k in range(n_per_side)]
    rgt = [complex(s1, w0 + (w1 - w0) * k / n_per_side) for k in range(n_per_side)]
    top = [complex(s1 - (s1 - s0) * k / n_per_side, w1) for k in range(n_per_side)]
    lft = [complex(s0, w1 - (w1 - w0) * k / n_per_side) for k in range(n_per_side)]
    pts = bot + rgt + top + lft + [complex(s0, w0)]
    return pts


def _safe_det(Y_func, s):
    """安全计算 det，遇异常或非有限值返回 None"""
    try:
        d = np.linalg.det(Y_func(s))
    except Exception:
        return None
    if not np.isfinite(d):
        return None
    return d


def _min_singular(Y_func, s):
    """计算 Y(s) 的最小奇异值（真伪根鉴别量）"""
    try:
        sig = np.linalg.svd(Y_func(s), compute_uv=False)
    except Exception:
        return None
    if not np.all(np.isfinite(sig)):
        return None
    return sig[-1] / sig[0]  # 归一化条件数，避免量级干扰


NAN_C = complex(float('nan'), float('nan'))


def _newton_refine(Y_func, s0, max_iter=60, tol=1e-13):
    """对 det(Y(s))=0 做 Newton 精化（遇 NaN/Inf 早退）"""
    s = complex(s0)
    for _ in range(max_iter):
        det_val = _safe_det(Y_func, s)
        if det_val is None:
            return NAN_C
        if abs(det_val) < tol:
            return s
        h = max(abs(s) * 1e-7, 1e-12)
        dp = _safe_det(Y_func, s + h)
        dm = _safe_det(Y_func, s - h)
        if dp is None or dm is None:
            return NAN_C
        d_det = (dp - dm) / (2 * h)
        if abs(d_det) < 1e-30:
            return s
        ds = det_val / d_det
        if not np.isfinite(ds):
            return NAN_C
        s_new = s - ds
        if abs(ds) < tol * max(abs(s), 1.0):
            return s_new
        s = s_new
    return s


def _dedup(eigs, tol=1e-3):
    """去重 (按欧氏距离, 过滤 NaN/Inf)"""
    unique = []
    for e in eigs:
        if e is None or not np.isfinite(e):
            continue
        if not any(abs(e - u) < tol for u in unique):
            unique.append(e)
    return unique


def _count_eigs_in_rect(Y_func, sigma_range, omega_range, n_per_side=80):
    """利用 argument principle 数围道内特征值个数，遇病态返回 None"""
    pts = _make_rect_contour(sigma_range, omega_range, n_per_side)
    dets = []
    for s in pts:
        d = _safe_det(Y_func, s)
        if d is None or d == 0:
            return None
        dets.append(d)
    args = np.unwrap(np.angle(np.array(dets)))
    N = (args[-1] - args[0]) / (2 * np.pi)
    if not np.isfinite(N):
        return None
    return int(round(N.real))


# ===================== Beyn (单围道) =====================

def _beyn_rect(Y_func, sigma_range, omega_range,
               num_points=192, l=8, eps_rank=1e-10):
    """单一矩形围道 Beyn 法，返回特征值列表"""
    s_min, s_max = sigma_range
    w_min, w_max = omega_range
    width = s_max - s_min
    height = w_max - w_min
    perim = 2 * (width + height)
    nh = max(int(round(num_points * width / perim)), 8)
    nv = max(int(round(num_points * height / perim)), 8)

    s_mid = complex((s_min + s_max) / 2, (w_min + w_max) / 2)
    Y0 = Y_func(s_mid)
    N = Y0.shape[0]
    rng = np.random.default_rng(42)
    V = rng.standard_normal((N, l)) + 1j * rng.standard_normal((N, l))

    A0 = np.zeros((N, l), dtype=complex)
    A1 = np.zeros((N, l), dtype=complex)

    sides = [
        (complex(s_min, w_min), 1.0, nh, width / nh),
        (complex(s_max, w_min), 1j, nv, height / nv),
        (complex(s_max, w_max), -1.0, nh, width / nh),
        (complex(s_min, w_max), -1j, nv, height / nv),
    ]
    for s_start, direction, n_pts, step in sides:
        ds = direction * step
        w = ds / (2j * np.pi)
        for k in range(n_pts):
            s = s_start + direction * (k + 0.5) * step
            try:
                X = np.linalg.solve(Y_func(s), V)
            except np.linalg.LinAlgError:
                continue
            A0 += w * X
            A1 += w * s * X

    U, sig, Vh = svd(A0, full_matrices=False)
    tol = eps_rank * sig[0] if sig[0] > 0 else eps_rank
    rank = int(np.sum(sig > tol))
    if rank == 0:
        return []
    Uk = U[:, :rank]
    Vk = Vh[:rank, :].conj().T
    M = np.diag(1.0 / sig[:rank]) @ (Uk.conj().T @ A1 @ Vk)
    raw = np.linalg.eigvals(M)
    margin = 0.001
    ms = (s_max - s_min) * margin
    mw = (w_max - w_min) * margin
    keep = [e for e in raw
            if s_min + ms < e.real < s_max - ms
            and w_min + mw < e.imag < w_max - mw]
    return keep


# ===================== 方法 1: Contour 法 (自适应) =====================

def contour_hunt(Y_func, sigma_range=(-15.0, 15.0),
                 omega_range=(-15.0, 15.0),
                 max_eigs_per_tile=2, min_tile_size=0.05,
                 num_points=192, l=8, eps_rank=1e-10,
                 verbose=False, depth=0, max_depth=6):
    """
    自适应递归围道分割搜索：
      1. 用 argument principle 数围道内特征值 N_in
      2. 若 N_in <= max_eigs_per_tile，直接 Beyn 提取
      3. 否则将围道一分为四，递归

    这样保证不漏检 (即便特征值很密集)，且每个子围道数值条件好。
    """
    s0, s1 = sigma_range
    w0, w1 = omega_range
    width = s1 - s0
    height = w1 - w0

    if width < min_tile_size or height < min_tile_size or depth >= max_depth:
        # 已经很小，直接 Beyn
        eigs = _beyn_rect(Y_func, sigma_range, omega_range,
                          num_points=num_points, l=l, eps_rank=eps_rank)
        return [_newton_refine(Y_func, e) for e in eigs]

    try:
        N_in = _count_eigs_in_rect(Y_func, sigma_range, omega_range, n_per_side=96)
    except Exception as ex:
        if verbose:
            print(f"  [contour] depth={depth} 计数异常: {ex}")
        N_in = None

    indent = '  ' * depth
    if verbose:
        print(f"{indent}[contour] σ∈[{s0:.2f},{s1:.2f}] ω∈[{w0:.2f},{w1:.2f}] N_in={N_in}")

    if N_in is None:
        # 数值病态：若已较小则直接 Beyn 兜底，否则强制细分
        if width < 2.0 and height < 2.0:
            eigs = _beyn_rect(Y_func, sigma_range, omega_range,
                              num_points=num_points, l=l, eps_rank=eps_rank)
            return _dedup([_newton_refine(Y_func, e) for e in eigs])
        N_in = max_eigs_per_tile + 1

    if N_in <= 0:
        return []

    if N_in <= max_eigs_per_tile:
        eigs = _beyn_rect(Y_func, sigma_range, omega_range,
                          num_points=num_points, l=max(l, N_in + 2),
                          eps_rank=eps_rank)
        refined = [_newton_refine(Y_func, e) for e in eigs]
        if verbose:
            for e in refined:
                print(f"{indent}  → eig = {e:.5f}")
        return refined

    # 一分为四
    sm = 0.5 * (s0 + s1)
    wm = 0.5 * (w0 + w1)
    sub_rects = [
        ((s0, sm), (w0, wm)),
        ((sm, s1), (w0, wm)),
        ((s0, sm), (wm, w1)),
        ((sm, s1), (wm, w1)),
    ]
    all_eigs = []
    for sr, wr in sub_rects:
        all_eigs.extend(contour_hunt(Y_func, sr, wr,
                                     max_eigs_per_tile=max_eigs_per_tile,
                                     min_tile_size=min_tile_size,
                                     num_points=num_points, l=l,
                                     eps_rank=eps_rank, verbose=verbose,
                                     depth=depth + 1, max_depth=max_depth))
    return _dedup(all_eigs, tol=1e-3)


# ===================== 方法 2: 根轨迹法 (Re/Im=0 零水平集求交) =====================

def root_locus_hunt(Y_func, sigma_range=(-6.0, 2.0),
                    omega_range=(-3.0, 3.0),
                    n_sigma=121, n_omega=121,
                    newton_tol=1e-12, dedup_tol=5e-3,
                    true_root_tol=1e-8, verbose=False):
    """
    根轨迹法（Spakovszky 论文原始做法的数值实现）：

      1. 在 (σ, ω) 网格上计算 det Y(s)，其中 s = σ - j·ω（论文约定）。
      2. 找出 Re det = 0 与 Im det = 0 两条零水平集。
         真特征值必同时位于两条曲线上 → 即两线交点。
      3. 在每个交点候选附近做 Newton 精化，再用最小奇异值过滤伪根。

    相比 Beyn (围道积分) + Shotgun (撒点) 的优势:
      - 不依赖 |det| 的大小（det 可能在远离根处也很小，或在根附近异常平坦），
        只看符号变化，鲁棒性强。
      - 对 MG 模态（位于原点附近、|det| 衰减极缓的弱奇点）尤其有效。
      - 几何直观：与论文图相同，可视化即可读出所有根。
    """
    s0, s1 = sigma_range
    w0, w1 = omega_range

    sigmas = np.linspace(s0, s1, n_sigma)
    omegas = np.linspace(w0, w1, n_omega)

    # det 网格 (行=ω, 列=σ)，注意论文 s = σ - jω
    det_grid = np.empty((n_omega, n_sigma), dtype=complex)
    for i, w in enumerate(omegas):
        for j, sg in enumerate(sigmas):
            s = complex(sg, -w)
            try:
                det_grid[i, j] = np.linalg.det(Y_func(s))
            except Exception:
                det_grid[i, j] = complex('nan')

    Re = det_grid.real
    Im = det_grid.imag

    # 找候选交点：单元 (i, j) 的四个角同时存在 Re 符号变化和 Im 符号变化
    seeds = []
    for i in range(n_omega - 1):
        for j in range(n_sigma - 1):
            re_cell = np.array([Re[i, j], Re[i, j + 1], Re[i + 1, j], Re[i + 1, j + 1]])
            im_cell = np.array([Im[i, j], Im[i, j + 1], Im[i + 1, j], Im[i + 1, j + 1]])
            if not (np.all(np.isfinite(re_cell)) and np.all(np.isfinite(im_cell))):
                continue
            re_has_cross = (np.min(re_cell) <= 0) and (np.max(re_cell) >= 0)
            im_has_cross = (np.min(im_cell) <= 0) and (np.max(im_cell) >= 0)
            if re_has_cross and im_has_cross:
                # 单元中心作为种子
                sg_c = 0.5 * (sigmas[j] + sigmas[j + 1])
                wm_c = 0.5 * (omegas[i] + omegas[i + 1])
                seeds.append(complex(sg_c, -wm_c))

    if verbose:
        print(f"[root_locus] 网格 {n_sigma}×{n_omega}, 候选交点 {len(seeds)} 个")

    converged = []
    for s0_seed in seeds:
        try:
            s_fin = _newton_refine(Y_func, s0_seed, max_iter=80, tol=newton_tol)
        except Exception:
            continue
        if not np.isfinite(s_fin):
            continue
        if not (s0 <= s_fin.real <= s1 and w0 <= -s_fin.imag <= w1):
            continue
        cond_inv = _min_singular(Y_func, s_fin)
        if cond_inv is None or cond_inv > true_root_tol:
            continue
        converged.append(s_fin)

    unique = _dedup(converged, tol=dedup_tol)
    if verbose:
        print(f"[root_locus] 收敛 {len(converged)} → 去重 {len(unique)}")
    return sorted(unique, key=lambda x: x.real, reverse=True)


# ===================== 综合搜索 =====================

def hybrid_hunt(Y_func, sigma_range=(-15.0, 15.0),
                omega_range=(-15.0, 15.0), verbose=False,
                true_root_tol=1e-8,
                locus_grid=(121, 121)):
    """
    两策略互补搜索 + 最小奇异值过滤:
      - Contour (Beyn 自适应分割) : 完备性，擅长远场/密集模态群
      - Root Locus (Re/Im=0 求交)  : 几何鲁棒，擅长原点附近弱奇点 (MG 模态)

    所有候选合并后, 用最小奇异值 σ_min/σ_max < true_root_tol 过滤伪根。
    """
    eigs_c = contour_hunt(Y_func, sigma_range, omega_range, verbose=verbose)
    eigs_r = root_locus_hunt(Y_func, sigma_range, omega_range,
                             n_sigma=locus_grid[0], n_omega=locus_grid[1],
                             true_root_tol=true_root_tol, verbose=verbose)
    merged = _dedup(list(eigs_c) + list(eigs_r), tol=1e-3)
    true_roots = []
    for e in merged:
        cond_inv = _min_singular(Y_func, e)
        if cond_inv is not None and cond_inv < true_root_tol:
            true_roots.append(e)
    return sorted(true_roots, key=lambda x: x.real, reverse=True)


# ===================== 特征向量提取 =====================

def extract_eigenvector(Y_func, s_star, tol=1e-8):
    """对特征值 s* 提取对应的右零空间向量 (特征向量)

    利用 SVD: Y(s*) = U·Σ·Vᴴ, 最小奇异值对应的右奇异向量即为
    近似零空间向量 (满足 Y·v ≈ 0)。

    返回归一化向量 v (模为 1)。若 s* 不是真根 (最小奇异值过大) 返回 None。
    """
    try:
        Y = Y_func(s_star)
    except Exception:
        return None
    try:
        U, S, Vh = np.linalg.svd(Y)
    except Exception:
        return None
    if not np.all(np.isfinite(S)):
        return None
    # 最小奇异值 / 最大奇异值 = 条件数倒数, 用于判断是否真根
    if S[0] <= 0 or S[-1] / S[0] > tol:
        return None
    v = Vh[-1].conj()           # 右奇异向量 (Vh 的最后一行共轭)
    return v / np.linalg.norm(v)