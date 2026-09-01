"""Eckardt O 基流计算器 —— 调用 1D 均值线程序, 输出稳定性矩阵所需参数

用法:
  python Eckardt_O_Geo_Calculator.py --rpm 14000 --m 5.0
  python Eckardt_O_Geo_Calculator.py --rpm 14000 --m 5.0 --update-params

输出包括:
  1. 各站位速度三角形、气流角、压力
  2. 无叶扩压器常数 Q 和 G (取扩压器 STA3~STA4 进出口平均值)
  3. 稳定性矩阵所需的全部基流参数
  4. --update-params 选项可自动更新 Eckardt_O_params.txt

角度约定:
  α = 绝对速度与子午面夹角, tan(α) = Vθ/Vm  (正值 = 正周向分量)
  β = 相对速度与子午面夹角, tan(β) = Wθ/Vm  (正值 = 前掠, 负值 = 后掠)
"""

import argparse
import math
import os
import sys
import re

# ──────────────────────────────
# 路径设置
# ──────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Stability 包路径 (向上两级: Geo_and_Base_Flow → Stability)
STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if STABILITY_ROOT not in sys.path:
    sys.path.insert(0, STABILITY_ROOT)

# 均值线包路径
MEANLINE_PKG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
    "MeanLIne_Design", "radcomp-main"
)
YAML_DATA = os.path.join(MEANLINE_PKG, "data", "known_compressors.yml")
if MEANLINE_PKG not in sys.path:
    sys.path.insert(0, MEANLINE_PKG)

import numpy as np

try:
    from ruamel.yaml import YAML
    _yaml_module = None
except ImportError:
    try:
        import yaml as _yaml_module
        YAML = None
    except ImportError as exc:
        raise ImportError("需要 ruamel.yaml 或 PyYAML") from exc

from radcompressor.compressor import Compressor
from radcompressor.condition import OperatingCondition
from radcompressor.geometry import Geometry
from radcompressor.thermo import CoolPropFluid


# ──────────────────────────────
# YAML 加载
# ──────────────────────────────
def load_yaml(path):
    with open(path, "r", encoding="utf-8") as fp:
        if _yaml_module is not None:
            return _yaml_module.safe_load(fp)
        return YAML(typ="safe").load(fp)


def find_record(records, name):
    for r in records:
        if r.get("name") == name:
            return r
    raise ValueError(f"未找到压缩机: {name}")


# ──────────────────────────────
# 单点快速计算 (提取 tanβ₁, μ, ψ, 用于间接法数值求导)
# ──────────────────────────────
def _run_perturbation_point(geom, in0, fld, m_flow, n_rot):
    """运行1D计算, 提取扰动点所需的参数

    返回 (tanβ₁, μ, ψ) 或 (None, None, None):
      tanβ₁: 入口相对流动角的正切
      μ = ΔH/U²: Euler 做功系数
      ψ = (P03-P02)/(ρ₂U²): 叶轮总压升系数
    """
    op = OperatingCondition(in0=in0, fld=fld, m=m_flow, n_rot=n_rot)
    comp = Compressor(geom, op)
    ok = comp.calculate(delta_check=False)
    if not ok or comp.invalid_flag:
        return None, None, None

    imp = comp.imp
    U2 = geom.r4 * n_rot

    beta1 = imp.in2.beta_m if not math.isnan(imp.in2.beta_m) else 0.0
    tan_beta1 = math.tan(math.radians(beta1))

    # Euler 做功系数
    delta_h = imp.out.total.H - comp.in_.total.H
    mu = delta_h / (U2 ** 2)

    # 叶轮总压升系数
    P02 = imp.in2.total.P
    P03 = imp.out.total.P
    rho2 = imp.in2.static.D
    psi = (P03 - P02) / (rho2 * U2 ** 2)

    return tan_beta1, mu, psi


# ──────────────────────────────
# 主计算
# ──────────────────────────────
def compute_base_flow(record, rpm, m_flow):
    """运行均值线计算, 返回各站位流动状态字典

    返回
    ----
    dict with keys:
      flowrate, rpm, U2, tau_u
      sta1_c, sta1_alpha, sta1_cm, sta1_ct, sta1_P0, sta1_P, sta1_M, sta1_a0
      sta2_c, sta2_alpha, sta2_cm, sta2_ct, sta2_beta, sta2_w, sta2_Wtheta, sta2_P0, sta2_P, sta2_M, sta2_Mrel
      sta3_c, sta3_alpha, sta3_cm, sta3_ct, sta3_beta, sta3_w, sta3_Wtheta, sta3_P0, sta3_P, sta3_M, sta3_Mrel
      sta4_c, sta4_alpha, sta4_cm, sta4_ct, sta4_P0, sta4_P, sta4_M
      Q, G                     (无叶扩压器常数, 取 STA3~STA4 进出口平均值)
      mu, loss_imp, omega_imp, omega_vls, PR, eff
    """
    geom = Geometry.from_dict(record["geom"])
    conditions = record["conditions"]
    fld = CoolPropFluid(conditions["fluid"])
    in0 = fld.thermo_prop("PT", float(conditions["in_P"]), float(conditions["in_T"]))
    n_rot = rpm * np.pi / 30.0
    U2 = geom.r4 * n_rot    # 有量纲叶轮出口轮缘速度 [m/s]

    op = OperatingCondition(in0=in0, fld=fld, m=m_flow, n_rot=n_rot)
    comp = Compressor(geom, op)
    ok = comp.calculate(delta_check=False)

    if not ok or comp.invalid_flag:
        raise RuntimeError(f"均值线计算失败 (rpm={rpm}, m={m_flow})")

    imp = comp.imp
    dif = comp.dif

    # ── STA1: 诱导轮入口 (一维 station 1) ──
    s1 = comp.ind.in1               # InducerState: c, alpha, total, static, m_abs
    sta1_c = s1.c
    sta1_alpha = s1.alpha_m if not math.isnan(s1.alpha_m) else 0.0
    sta1_cm = s1.c * math.cos(math.radians(sta1_alpha))
    sta1_ct = s1.c * math.sin(math.radians(sta1_alpha))
    sta1_P0 = s1.total.P
    sta1_P = s1.static.P
    sta1_M = s1.m_abs
    sta1_a0 = s1.total.A              # 总声速

    # ── STA2: 叶轮进口 (一维 station 2) ──
    s2 = imp.in2                      # ImpellerState: + w, beta, m_rel
    sta2_c = s2.c
    sta2_alpha = s2.alpha_m if not math.isnan(s2.alpha_m) else 0.0
    sta2_cm = s2.c * math.cos(math.radians(sta2_alpha))
    sta2_ct = s2.c * math.sin(math.radians(sta2_alpha))
    sta2_w = s2.w
    sta2_beta = s2.beta_m if not math.isnan(s2.beta_m) else 0.0
    sta2_Wtheta = s2.w * math.sin(math.radians(s2.beta_m))
    sta2_P0 = s2.total.P
    sta2_P = s2.static.P
    sta2_M = s2.m_abs
    sta2_Mrel = s2.m_rel
    sta2_rho = s2.static.D            # 入口静态密度 [kg/m³]

    # ── STA3: 叶轮出口 / 无叶扩压器入口 (一维 station 4) ──
    s3 = imp.out
    sta3_c = s3.c
    sta3_alpha = s3.alpha_m if not math.isnan(s3.alpha_m) else 0.0
    sta3_cm = s3.c * math.cos(math.radians(sta3_alpha))
    sta3_ct = s3.c * math.sin(math.radians(sta3_alpha))
    sta3_w = s3.w
    sta3_beta = s3.beta_m if not math.isnan(s3.beta_m) else 0.0
    sta3_Wtheta = s3.w * math.sin(math.radians(s3.beta_m))
    sta3_P0 = s3.total.P
    sta3_P = s3.static.P
    sta3_M = s3.m_abs
    sta3_Mrel = s3.m_rel
    sta3_rho = s3.static.D            # 出口静态密度 [kg/m³]

    # ── STA4: 无叶扩压器出口 (一维 station 5) ──
    s4 = dif.out                      # VanelessState = InducerState
    sta4_c = s4.c
    sta4_alpha = s4.alpha_m if not math.isnan(s4.alpha_m) else 0.0
    sta4_cm = s4.c * math.cos(math.radians(sta4_alpha))
    sta4_ct = s4.c * math.sin(math.radians(sta4_alpha))
    sta4_P0 = s4.total.P
    sta4_P = s4.static.P
    sta4_M = s4.m_abs

    # ── 无叶扩压器常数: 取 STA3~STA4 进出口平均值 ──
    # Q = Vr * r, G = Vtheta * r. 归一化: / (U2 * r4)
    # 入口: Q3 = c3m*r4, G3 = c3t*r4; 出口: Q4 = c4m*r5, G4 = c4t*r5
    Q3 = sta3_cm * geom.r4
    G3 = sta3_ct * geom.r4
    Q4 = sta4_cm * geom.r5
    G4 = sta4_ct * geom.r5
    Q = 0.5 * (Q3 + Q4) / (U2 * geom.r4)
    G_val = 0.5 * (G3 + G4) / (U2 * geom.r4)

    # ── 参考半径比: 用于轴向管道矩阵的周向波数缩放 (离心机以 R2 归一化) ──
    # mu_up = R2 / R_duct_up, mu_plenum = R2 / R_duct_plenum
    mu_up = geom.r4 / geom.r2rms                     # = 1/STA2_radial (叶轮进口均方根半径)
    mu_plenum = geom.r4 / geom.r5                    # = 1/STA4_radial (扩压器出口半径)

    # ── 做功与损失 ──
    delta_h = imp.out.total.H - comp.in_.total.H
    mu = delta_h / (U2 ** 2)                       # 焓升做功系数

    # 叶轮总焓损失 → 无量纲损失系数 (除以 U2^2/2)
    loss_enth = (imp.losses.incidence + imp.losses.skin_friction +
                 imp.losses.blade_loading + imp.losses.clearance +
                 imp.losses.mixing + imp.losses.disc_friction)
    loss_imp = loss_enth / (0.5 * U2 * U2)
    # 无叶扩压器总压损失系数
    omega_vls = 1.0 - sta4_P0 / sta3_P0 if sta3_P0 > 0 else 0.0

    # 归一化转速
    tau_u = U2 / sta1_a0

    # ── 叶轮惯性参数 (Spakovszky 模型Ⅱ公式) ──
    # AR_imp = ρ₃·A₃ / (ρ₂·A₂)  (含密度比的面积比)
    # A₂ = π*(r2s²-r2h²), A₃ = 2π*r4*b4
    A2 = math.pi * (geom.r2s ** 2 - geom.r2h ** 2)
    A3 = 2 * math.pi * geom.r4 * geom.b4
    AR_imp_with_density = (sta3_rho * A3) / (sta2_rho * A2)

    # s_imp = l_comp / R2  (归一化叶轮流道长度)
    s_imp = geom.l_comp / geom.r4

    # λ_imp = s_imp × AR × ln(AR) / (AR - 1)  (公式15)
    # 当 AR→1 时, λ→s_imp (L'Hôpital)
    if abs(AR_imp_with_density - 1.0) < 1e-10:
        lambda_imp = s_imp
    else:
        lambda_imp = s_imp * AR_imp_with_density * math.log(AR_imp_with_density) / (AR_imp_with_density - 1.0)

    # τ_imp = τ_u × 2·s_imp / (Ŵ₁ + Ŵ₂)  (公式20)
    W_hat1 = sta2_w / U2   # 入口归一化相对速度
    W_hat2 = sta3_w / U2   # 出口归一化相对速度
    tau_imp = tau_u * 2.0 * s_imp / (W_hat1 + W_hat2)

    # AR_imp 几何面积比 (不含密度, 供参考)
    AR_imp_geom = A3 / A2

    # ── 损失特性导数: dL/d(tanβ₁) ──
    # 间接法: dL/dtanβ₁ = dμ/dtanβ₁ - dψ/dtanβ₁
    # ψ = μ - L (总压升 = Euler做功 - 损失)
    # 直接差分焓损失会漏掉密度变化和速度三角形改变的附加效应,
    # 间接法从总压升特性线推算, 物理上更完备
    tanb1_base = math.tan(math.radians(sta2_beta))

    # 基点 μ 和 ψ
    delta_h_base = imp.out.total.H - comp.in_.total.H
    mu_base = delta_h_base / (U2 ** 2)
    psi_base = (sta3_P0 - sta2_P0) / (sta2_rho * U2 ** 2)
    loss_norm_base = mu_base - psi_base  # L = μ - ψ

    delta = 0.02
    m_plus = m_flow * (1.0 + delta)
    m_minus = m_flow * (1.0 - delta)

    tanb1_plus, mu_plus, psi_plus = _run_perturbation_point(geom, in0, fld, m_plus, n_rot)
    tanb1_minus, mu_minus, psi_minus = _run_perturbation_point(geom, in0, fld, m_minus, n_rot)

    flow_ok = True
    dmu_dtanbeta1 = 0.0
    dpsi_dtanbeta1 = 0.0
    dL_dtanbeta1 = 0.0

    if (tanb1_plus is not None and tanb1_minus is not None
            and abs(tanb1_plus - tanb1_minus) > 1e-12):
        dpsi_dtanbeta1 = (psi_plus - psi_minus) / (tanb1_plus - tanb1_minus)
        dmu_dtanbeta1 = (mu_plus - mu_minus) / (tanb1_plus - tanb1_minus)
        dL_dtanbeta1 = dmu_dtanbeta1 - dpsi_dtanbeta1
    elif tanb1_plus is not None and mu_plus is not None:
        dpsi_dtanbeta1 = (psi_plus - psi_base) / (tanb1_plus - tanb1_base)
        dmu_dtanbeta1 = (mu_plus - mu_base) / (tanb1_plus - tanb1_base)
        dL_dtanbeta1 = dmu_dtanbeta1 - dpsi_dtanbeta1
    elif tanb1_minus is not None and mu_minus is not None:
        dpsi_dtanbeta1 = (psi_base - psi_minus) / (tanb1_base - tanb1_minus)
        dmu_dtanbeta1 = (mu_base - mu_minus) / (tanb1_base - tanb1_minus)
        dL_dtanbeta1 = dmu_dtanbeta1 - dpsi_dtanbeta1
    else:
        flow_ok = False

    return {
        "flowrate": m_flow,
        "rpm": rpm,
        "U2": U2,
        "tau_u": tau_u,
        # STA1
        "sta1_c": sta1_c, "sta1_alpha": sta1_alpha,
        "sta1_cm": sta1_cm, "sta1_ct": sta1_ct,
        "sta1_P0": sta1_P0, "sta1_P": sta1_P, "sta1_M": sta1_M, "sta1_a0": sta1_a0,
        # STA2
        "sta2_c": sta2_c, "sta2_alpha": sta2_alpha,
        "sta2_cm": sta2_cm, "sta2_ct": sta2_ct,
        "sta2_beta": sta2_beta, "sta2_w": sta2_w, "sta2_Wtheta": sta2_Wtheta,
        "sta2_P0": sta2_P0, "sta2_P": sta2_P, "sta2_M": sta2_M, "sta2_Mrel": sta2_Mrel,
        "sta2_rho": sta2_rho,
        # STA3
        "sta3_c": sta3_c, "sta3_alpha": sta3_alpha,
        "sta3_cm": sta3_cm, "sta3_ct": sta3_ct,
        "sta3_beta": sta3_beta, "sta3_w": sta3_w, "sta3_Wtheta": sta3_Wtheta,
        "sta3_P0": sta3_P0, "sta3_P": sta3_P, "sta3_M": sta3_M, "sta3_Mrel": sta3_Mrel,
        "sta3_rho": sta3_rho,
        # STA4
        "sta4_c": sta4_c, "sta4_alpha": sta4_alpha,
        "sta4_cm": sta4_cm, "sta4_ct": sta4_ct,
        "sta4_P0": sta4_P0, "sta4_P": sta4_P, "sta4_M": sta4_M,
        # 常数
        "Q": Q,
        "G": G_val,
        "mu": mu,
        # 参考半径比 (离心机轴向管道周向波数缩放)
        "mu_up": mu_up,
        "mu_plenum": mu_plenum,
        "loss_imp": loss_imp,
        "loss_enth": loss_enth,
        "omega_vls": omega_vls,
        "PR": imp.out.total.P / imp.in2.total.P,
        "eff": imp.eff,
        # 叶轮惯性参数
        "AR_imp_with_density": AR_imp_with_density,
        "AR_imp_geom": AR_imp_geom,
        "s_imp": s_imp,
        "lambda_imp": lambda_imp,
        "tau_imp": tau_imp,
        # 损失特性导数 (间接法)
        "dL_dtanbeta1": dL_dtanbeta1,
        "dmu_dtanbeta1": dmu_dtanbeta1,
        "dpsi_dtanbeta1": dpsi_dtanbeta1,
        "m_plus": m_plus, "m_minus": m_minus,
        "tanb1_base": tanb1_base,
        "tanb1_plus": tanb1_plus, "tanb1_minus": tanb1_minus,
        # 间接法诊断量
        "mu_base": mu_base, "mu_plus": mu_plus, "mu_minus": mu_minus,
        "psi_base": psi_base, "psi_plus": psi_plus, "psi_minus": psi_minus,
        "loss_base": loss_norm_base,
        "flow_ok": flow_ok,
    }


# ──────────────────────────────
# 格式化输出
# ──────────────────────────────
def print_diagnostics(data):
    """打印完整的诊断信息"""
    U2 = data["U2"]

    print(f"\n{'='*80}")
    print(f"  Eckardt O 基流计算")
    print(f"  转速: {data['rpm']:.0f} rpm  |  流量: {data['flowrate']:.4f} kg/s  |  U2: {U2:.1f} m/s  |  τ_u: {data['tau_u']:.4f}")
    print(f"  叶轮总压比 P04/P02: {data['PR']:.4f}  |  叶轮等熵效率: {data['eff']:.4f}")
    print(f"{'='*80}")

    # ── 各站位速度三角形 ──
    header = f"  {'站位':>6s}  {'α[°]':>8s}  {'β[°]':>8s}  {'c':>8s}  {'cm':>8s}  {'ct':>8s}  {'w':>8s}  {'M_abs':>8s}  {'M_rel':>8s}  {'P0[Pa]':>12s}  {'P[Pa]':>12s}"
    print(f"\n── 各站位速度三角形 (U2={U2:.1f} m/s) ──")
    print(header)
    print("  " + "-" * (len(header) - 2))

    def print_row(name, d, pre, has_rel=False):
        p0 = d.get(f"{pre}_P0", float("nan"))
        ps = d.get(f"{pre}_P", float("nan"))
        alpha = d.get(f"{pre}_alpha", float("nan"))
        beta = d.get(f"{pre}_beta", float("nan")) if has_rel else float("nan")
        c = d.get(f"{pre}_c", float("nan"))
        cm = d.get(f"{pre}_cm", float("nan"))
        ct = d.get(f"{pre}_ct", float("nan"))
        w = d.get(f"{pre}_w", float("nan")) if has_rel else float("nan")
        m_abs = d.get(f"{pre}_M", float("nan"))
        m_rel = d.get(f"{pre}_Mrel", float("nan")) if has_rel else float("nan")
        print(f"  {name:6s}  {alpha:8.2f}  {beta:8.2f}  {c:8.2f}  {cm:8.2f}  {ct:8.2f}  "
              f"{w:8.2f}  {m_abs:8.4f}  {m_rel:8.4f}  {p0:12.1f}  {ps:12.1f}")

    print_row("STA1", data, "sta1", has_rel=False)
    print_row("STA2", data, "sta2", has_rel=True)
    print_row("STA3", data, "sta3", has_rel=True)
    print_row("STA4", data, "sta4", has_rel=False)

    # ── 稳定性矩阵参数 ──
    U2 = data["U2"]
    print(f"\n── 稳定性矩阵归一化参数 (归一化基准: U2={U2:.1f} m/s, r4) ──")
    print(f"  {'参数':<22s}  {'值':>12s}  {'说明':s}")
    print(f"  {'-'*22}  {'-'*12}  {'-'*30}")

    items_stab = [
        ("Vx_bar1_imp",    data["sta2_cm"] / U2,       "叶轮入口子午速度 / U2  (= Vx_up)"),
        ("Vtheta_bar1_imp", data["sta2_ct"] / U2,       "叶轮入口切向速度 / U2  (≈0)"),
        ("Vr_bar2_imp",    data["sta3_cm"] / U2,       "叶轮出口径向速度 / U2"),
        ("Vtheta_bar2_imp", data["sta3_ct"] / U2,       "叶轮出口切向速度 / U2"),
        ("beta1_imp_deg",  data["sta2_beta"],          "入口相对流动角 [°]"),
        ("beta2_imp_deg",  0.0,                        "出口叶片角 [°] (径向=0)"),
        ("alpha1_imp_deg", data["sta2_alpha"],         "入口绝对流动角 [°]"),
        ("alpha2_imp_deg", data["sta3_alpha"],         "出口绝对流动角 [°]"),
        ("AR_imp (geom)",  data["AR_imp_geom"],         "几何面积比 A₃/A₂"),
        ("AR_imp (ρA)",    data["AR_imp_with_density"],  "含密度面积比 ρ₃A₃/(ρ₂A₂)"),
        ("s_imp",          data["s_imp"],               "归一化流道长度 l_comp/R2"),
        ("lambda_imp",     data["lambda_imp"],          "惯性系数 (公式15)"),
        ("tau_imp",        data["tau_imp"],             "损失时间滞后 (公式20)"),
        ("dL_dtanbeta1",   data["dL_dtanbeta1"],        "损失对 tanβ₁ 的导数 (间接法)"),
        ("dmu_dtanbeta1",  data["dmu_dtanbeta1"],       "Euler做功斜率 dμ/dtanβ₁"),
        ("dpsi_dtanbeta1", data["dpsi_dtanbeta1"],      "总压升斜率 dψ/dtanβ₁"),
        ("mu (enthalpy)",  data["mu"],                  "dH/U2^2 (做功系数)"),
        ("Q (进出口均值)", data["Q"],                  "无叶区常数 Vr*r/(U2*r4)"),
        ("G (进出口均值)", data["G"],                  "无叶区常数 Vθ*r/(U2*r4)"),
        ("mu_up",           data["mu_up"],              "上游参考半径比 R2/R_duct_up"),
        ("mu_plenum",       data["mu_plenum"],          "容腔参考半径比 R2/R_duct_plenum"),
        ("Vx_up",          data["sta1_cm"] / U2,       "上游管道轴向速度 / U2"),
        ("Vtheta_up",      0.0,                        "上游无涡"),
        ("Vr_plenum",      data["Q"] / 1.69,           "容腔入口径向速度 = Q/r5"),
        ("Vtheta_plenum",  data["G"] / 1.69,           "容腔入口切向速度 = G/r5"),
    ]

    for name, val, desc in items_stab:
        if isinstance(val, float) and not math.isnan(val):
            print(f"  {name:<22s}  {val:12.6f}  {desc}")

    print(f"\n  * Q, G 取无叶扩压器 STA3~STA4 进出口平均值")
    print(f"  * AR_imp (ρA) 含密度比, 用于 λ_imp 和 τ_imp 计算")
    print(f"  * λ_imp = s_imp × AR×ln(AR)/(AR-1), τ_imp = τ_u × 2s_imp/(Ŵ₁+Ŵ₂)")
    if data.get("flow_ok"):
        method = "中心差分" if (data.get("tanb1_plus") is not None and data.get("tanb1_minus") is not None) else "单侧差分"
        print(f"  * dL_dtanβ₁ = dμ/dtanβ₁ - dψ/dtanβ₁ (间接法, ±2% {method}):")
        print(f"    m⁻={data['m_minus']:.4f}  tanβ₁⁻={data.get('tanb1_minus', 0):.4f}  μ⁻={data.get('mu_minus', 0):.6f}  ψ⁻={data.get('psi_minus', 0):.6f}")
        print(f"    m⁰={data['flowrate']:.4f}  tanβ₁⁰={data['tanb1_base']:.4f}  μ⁰={data.get('mu_base', 0):.6f}  ψ⁰={data.get('psi_base', 0):.6f}")
        print(f"    m⁺={data['m_plus']:.4f}  tanβ₁⁺={data.get('tanb1_plus', 0):.4f}  μ⁺={data.get('mu_plus', 0):.6f}  ψ⁺={data.get('psi_plus', 0):.6f}")
        print(f"    dμ/dtanβ₁ = {data['dmu_dtanbeta1']:+.6f}")
        print(f"    dψ/dtanβ₁ = {data['dpsi_dtanbeta1']:+.6f}")
        print(f"    dL/dtanβ₁ = {data['dL_dtanbeta1']:+.6f}")
    else:
        print(f"  [!] 流量扰动计算失败, dL_dtanβ₁ 保持为 0")


# ──────────────────────────────
# 更新参数文件
# ──────────────────────────────
def update_params_file(data, params_path=None):
    """将计算结果回写到 Eckardt_O_params.txt"""
    if params_path is None:
        params_path = os.path.join(SCRIPT_DIR, "Eckardt_O_params.txt")

    with open(params_path, "r", encoding="utf-8") as f:
        content = f.read()

    U2 = data["U2"]
    replacements = {
        "Vx_bar1_imp":          f"{data['sta2_cm'] / U2:.6f}",
        "Vr_bar2_imp":          f"{data['sta3_cm'] / U2:.6f}",
        "Vtheta_bar2_imp":      f"{data['sta3_ct'] / U2:.6f}",
        "beta1_imp_deg":        f"{data['sta2_beta']:.4f}",
        "alpha2_imp_deg":       f"{data['sta3_alpha']:.4f}",
        "lambda_imp":           f"{data['lambda_imp']:.6f}",
        "tau_imp":              f"{data['tau_imp']:.6f}",
        "dL_dtanbeta1_imp":     f"{data['dL_dtanbeta1']:.6f}",
        "AR_imp":               f"{data['AR_imp_with_density']:.6f}",
        "Vx_up":                f"{data['sta1_cm'] / U2:.6f}",
        "Q":                    f"{data['Q']:.6f}",
        "G":                    f"{data['G']:.6f}",
        "mu_up":                f"{data['mu_up']:.4f}",
        "mu_plenum":            f"{data['mu_plenum']:.4f}",
    }

    for key, value in replacements.items():
        # 只替换 = [] 的行, 跳过已有数值的行
        pattern = rf"^({key}\s*=\s*)\[\]"
        replacement = rf"\g<1>{value}"
        content, count = re.subn(pattern, replacement, content, flags=re.MULTILINE)
        if count:
            print(f"  [OK] {key} = {value}")
        else:
            # 如果已经有数值, 也更新
            pattern2 = rf"^({key}\s*=\s*)[\d\.\-eE]+"
            replacement2 = rf"\g<1>{value}"
            content, count2 = re.subn(pattern2, replacement2, content, flags=re.MULTILINE)
            if count2:
                print(f"  [~] {key} = {value}  (已更新)")

    # 在文件顶部添加工况注释
    header_line = f"# Operating point: {data['rpm']:.0f} rpm, {data['flowrate']:.4f} kg/s, PR_imp(P04/P02)={data['PR']:.4f}\n"
    if not content.startswith("# Operating point:"):
        content = header_line + content
    else:
        content = re.sub(r"^# Operating point:.*\n", header_line, content, count=1, flags=re.MULTILINE)

    with open(params_path, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"\n参数已更新到: {params_path}")


# ──────────────────────────────
# 主入口
# ──────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Eckardt O 基流计算器 —— 输出稳定性矩阵所需参数"
    )
    parser.add_argument("--rpm", type=float, default=14000.0, help="转速 [rpm]")
    parser.add_argument("--m", type=float, default=5.0, dest="m_flow", help="质量流量 [kg/s]")
    parser.add_argument("--update-params", action="store_true",
                        help="自动更新 Eckardt_O_params.txt")
    args = parser.parse_args()

    records = load_yaml(YAML_DATA)
    record = find_record(records, "Eckardt O")

    print(f"运行条件: {args.rpm:.0f} rpm, {args.m_flow:.4f} kg/s")

    data = compute_base_flow(record, args.rpm, args.m_flow)
    print_diagnostics(data)

    if args.update_params:
        update_params_file(data)


if __name__ == "__main__":
    main()
