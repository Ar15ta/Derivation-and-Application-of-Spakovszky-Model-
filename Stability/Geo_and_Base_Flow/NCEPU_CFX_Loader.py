"""NCEPU_CFX_Loader.py - 从 CFX 提取的 CSV 加载 NCEPU 基流参数

CFX 后处理在各站位截面导出: c_m/c_u/rho 用质量流量平均
(massFlowAveraged), p/T/pt/Tt 用面积平均 (areaAveraged)。
本模块读取 CSV, 结合 NCEPU_params.txt 中的几何参数, 计算稳定性模型
所需的全部无量纲参数。

命名约定
────────────────────────────────────────────────────────────────────
速度 (与 Cumpsty/Greitzer 一致):
  c = 绝对速度, w = 相对速度, u = 轮缘(牵连)速度
  分量下标: m=子午, u=周向(circumferential), r=径向, x=轴向
  站位作后缀: _2=叶轮进口(STA2), _3=叶轮出口(STA3), _4=扩压器出口(STA4)
热力学:
  p = 静压, pt = 总压; a = 静声速, at = 总声速
角度:
  alpha = 绝对气流角, beta = 相对气流角, beta_*_blade = 叶片金属角
注意: 子午速度在绝对/相对系相同 (牵连速度只有周向分量)。
      CSV 中的周向速度一律填【绝对(静止系)】值; 相对周向 w_theta 由代码反算。

────────────────────────────────────────────────────────────────────
CSV 文件格式 (首行为表头, 每行一个工况点)
────────────────────────────────────────────────────────────────────
速度分量 c_m/c_u、密度 rho 用质量流量平均 (massFlowAveraged);
压力/温度 p/T/pt/Tt 用面积平均 (areaAveraged)。角度 α/β 与相对速度
w 均由分量计算, 无需提供。切向速度 c_u 一律填 CFX 原始
带符号值, loader 按 rpm 符号自动翻成"旋转方向为正、u3>0"的模型约定。

每个站位 STA2/3/4 统一提供 7 个量:
  c_m   绝对子午速度 [m/s]   (massFlowAveraged)
  c_u   绝对周向速度 [m/s]   (massFlowAveraged)
  rho   静密度 [kg/m³]       (massFlowAveraged)
  p     静压 [Pa]            (areaAveraged)
  T     静温 [K]             (areaAveraged)
  pt    总压 [Pa]            (areaAveraged)
  Tt    总温 [K]             (areaAveraged)

真实气体物性由 CoolProp 按 params.txt 中的 fluid (默认 CO2) 查询:
  静声速 a  = PropsSI('A', T,  p,  fluid)   (能量分布后处理用)
  总焓   h  = PropsSI('H', Tt, pt, fluid)
  叶轮总焓升 dh_t = h3 - h2                  (损失导数自动差分用)
因此 CSV 不再有 a*/dh_t 列。

基本列:
  rpm          转速 [rpm] (带 CFX 原始符号)
  m            质量流量 [kg/s]

各站位:
  c_m2, c_u2, rho2, p2, T2, pt2, Tt2   (叶轮进口)
  c_m3, c_u3, rho3, p3, T3, pt3, Tt3   (叶轮出口; c_u3 缺失走 Wiesner)
  c_m4, c_u4, rho4, p4, T4, pt4, Tt4   (扩压器出口; rho4 缺省=rho3)

派生量 (无需填写, loader 自动算):
  α2/β2/α3/β3  由 c_m、c_u 和牵连速度 u(r) 反算
  w2/w3        相对速度大小 = hypot(c_m, w_theta)
  π_tt         = pt3/pt2

损失导数 dL/dtanβ1 (二选一):
  1) 直接提供列 dL_dtanbeta1
  2) 留空并提供同转速相邻工况 (邻点需含 pt2/pt3/Tt2/Tt3),
     loader 取高/低流量两侧最近点做中心差分; 端点只有单侧时退化为单侧差分
  两者皆无时置 0 (不推荐)。
────────────────────────────────────────────────────────────────────
用法:
  python -m Geo_and_Base_Flow.NCEPU_CFX_Loader --csv NCEPU_CFX.csv
"""

import csv
import math
import os
import re
import numpy as np

from Matrix.rotor import load_params

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PARAMS_FILE = os.path.join(SCRIPT_DIR, "NCEPU_params.txt")
DEFAULT_CSV_FILE = os.path.join(SCRIPT_DIR, "NCEPU_CFX.csv")


# ──────────────────────────────
# 坐标系 / 切向速度符号约定
# ──────────────────────────────
def read_z_axis(params_path):
    """从参数文件读取 z_axis 字符串 (Matrix 的数值解析器会忽略非数字行,
    故这里单独解析)。默认 'down'。

    z_axis = "down": CFD 圆柱坐标 +Z 由叶轮入口指向叶轮内部
    z_axis = "up":   +Z 由叶轮内部指向入口
    """
    val = _read_text_param(params_path, 'z_axis')
    return val.lower() if val else 'down'


def read_fluid(params_path):
    """从参数文件读取 CoolProp 工质名 (数值解析器忽略非数字行, 故单独解析)。
    默认 'CO2'。"""
    val = _read_text_param(params_path, 'fluid')
    return val if val else 'CO2'


def _read_text_param(params_path, key):
    """读取 `key = value` 形式的文本参数 (Matrix 数值解析器不解析此类行)。"""
    try:
        with open(params_path, 'r', encoding='utf-8') as f:
            for line in f:
                m = re.match(r'\s*' + re.escape(key) + r'\s*=\s*([^\s#]+)', line)
                if m:
                    return m.group(1).strip()
    except OSError:
        pass
    return None


def circumferential_sign(rpm, z_axis='down'):
    """计算 c_u_sign: 把 CFX 切向速度翻成"旋转方向为正、轮缘速度 u3>0"的模型约定。

    CFX 圆柱坐标中 +θ 由右手定则 e_θ = e_z × e_r 确定, 而转速 rpm 的符号本身
    就是相对同一个 +z 轴按右手定则定义的。因此叶片(牵连)速度在 CFX +θ 方向上的
    投影符号恒等于 sign(rpm), 与 z 轴指向无关:
      - 若翻转 z 轴定义, 则同一物理旋转的 rpm 符号和 c_u 符号会同时翻转,
        二者比值不变, 故切向翻译只取决于 rpm 符号。
      - z_axis 仅影响轴向速度 c_x 的正负; 这里 c_m 一律按正的子午幅值填入,
        故不参与切向符号判定 (保留参数仅为记录坐标系设置)。

        c_u_sign = sign(rpm) = +1  若 CFX 的 +θ 与旋转方向一致 (rpm≥0)
                               -1  若相反 (rpm<0)
    返回后: u3 = |ω|·R_3,  c_u(model) = c_u_sign · c_u(CFX);
            直填的气流角同样乘 c_u_sign (翻转 θ 轴即翻转角的正负)。
    """
    return 1.0 if rpm >= 0 else -1.0


# ──────────────────────────────
# 真实气体物性 (CoolProp)
# ──────────────────────────────
# 静声速 a = f(T, p); 总焓 h = f(Tt, pt)。
# 工质疑似 sCO2, 不能用 γ=1.4 理想气体回退; CoolProp 不可用或查询失败时返回 0。
try:
    from CoolProp.CoolProp import PropsSI
    _HAS_COOLPROP = True
except Exception:                                   # pragma: no cover
    PropsSI = None
    _HAS_COOLPROP = False


def sound_speed(fluid, T, p):
    """由静温 T[K]、静压 p[Pa] 查静声速 a[m/s]。"""
    if not _HAS_COOLPROP or T <= 0 or p <= 0:
        return 0.0
    try:
        return float(PropsSI('A', 'T', T, 'P', p, fluid))
    except Exception:
        return 0.0


def stagnation_enthalpy(fluid, Tt, pt):
    """由总温 Tt[K]、总压 pt[Pa] 查总焓 h[J/kg]。"""
    if not _HAS_COOLPROP or Tt <= 0 or pt <= 0:
        return 0.0
    try:
        return float(PropsSI('H', 'T', Tt, 'P', pt, fluid))
    except Exception:
        return 0.0


def row_sound_speed(row, sta, fluid):
    """从一行 CSV 取静温/静压查声速 (a<sta> 列已删, 统一由 T,p 查)。"""
    T = float(row[f'T{sta}']) if row.get(f'T{sta}', '') not in ('', None) else 0.0
    p = float(row[f'p{sta}']) if row.get(f'p{sta}', '') not in ('', None) else 0.0
    return sound_speed(fluid, T, p)


def row_enthalpy(row, sta, fluid):
    """从一行 CSV 取总温/总压查总焓。"""
    Tt = float(row[f'Tt{sta}']) if row.get(f'Tt{sta}', '') not in ('', None) else 0.0
    pt = float(row[f'pt{sta}']) if row.get(f'pt{sta}', '') not in ('', None) else 0.0
    return stagnation_enthalpy(fluid, Tt, pt)


# ──────────────────────────────
# 几何派生量
# ──────────────────────────────
def compute_geometry(params):
    """从原始几何 (R_3, r_2h, r_2s, b_3, r_4, l_comp) 计算派生几何量

    若派生量在 params 文件中已有非零值, 则不覆盖 (允许手动指定)。
    """
    R_3 = params.get('R_3', 0.0)
    if R_3 <= 0:
        raise ValueError("R_3 (叶轮出口半径) 必须在 NCEPU_params.txt 中设置为正值")

    r_2h = params.get('r_2h', 0.0)
    r_2s = params.get('r_2s', 0.0)
    b_3 = params.get('b_3', 0.0)
    r_4 = params.get('r_4', 0.0)
    l_comp = params.get('l_comp', 0.0)
    t_imp = params.get('t_imp', 0.0)
    Z_total = params.get('Z_blades', 0) + params.get('Z_split', 0)

    # 叶轮进口均方根半径
    if r_2s > 0 and r_2h >= 0:
        r_2rms = math.sqrt(0.5 * (r_2s ** 2 + r_2h ** 2))
    else:
        r_2rms = 0.0

    # 派生几何 (仅当文件中为 0 时自动计算)
    if r_2rms > 0 and params.get('radial_sta2', 0.0) == 0.0:
        params['radial_sta2'] = r_2rms / R_3
    if r_2s > 0 and r_2h > 0 and params.get('area_sta2', 0.0) == 0.0:
        A2 = math.pi * (r_2s ** 2 - r_2h ** 2)
        params['area_sta2'] = A2 / (R_3 ** 2)
    if b_3 > 0 and params.get('area_sta3', 0.0) == 0.0:
        # 出口总流通面积, 扣除 Z_total 个尾缘厚度的周向堵塞
        # (文档式 A2 = h2·(2πR2/N − t_imp), 总流通面积即 b_3·(2πR3 − Z_total·t_imp))
        A3 = b_3 * max(0.0, 2.0 * math.pi * R_3 - Z_total * t_imp)
        params['area_sta3'] = A3 / (R_3 ** 2)
    if r_4 > 0 and params.get('radial_sta4', 0.0) == 0.0:
        params['radial_sta4'] = r_4 / R_3

    params['R1_R2'] = params['radial_sta2'] / params['radial_sta3']
    params['r_impeller_exit'] = params['radial_sta3']
    params['r_diffuser_exit'] = params['radial_sta4']

    # 参考半径比 (轴向管道周向波数缩放)
    if r_2rms > 0 and params.get('mu_up', 0.0) == 0.0:
        params['mu_up'] = R_3 / r_2rms
    if r_4 > 0 and params.get('mu_plenum', 0.0) == 0.0:
        params['mu_plenum'] = R_3 / r_4

    # 出口叶片金属角仅作非 CFX 路径兜底; CFX 路径下 beta2_imp_deg
    # (Matrix 叶轮出口相对角) 由 cfx_row_to_params 用实测相对气流角覆盖
    if params.get('beta2_imp_deg', 0.0) == 0.0:
        params['beta2_imp_deg'] = params.get('beta_blade_3_deg', 0.0)

    # 归一化叶轮流道长度
    params['s_imp'] = l_comp / R_3 if l_comp > 0 else params.get('s_imp', 0.0)

    return params


# ──────────────────────────────
# Wiesner 滑移因子 (出口角/周向速度的兜底估算)
# ──────────────────────────────
def wiesner_slip_factor(chi_deg, Z_blades, Z_split=0):
    """Wiesner 滑移因子: sigma_W = 1 - sqrt(cos chi)/(Z+Z_split)^0.7

    参数
    ----
    chi_deg : float
        出口叶片金属角 χ3 [°] (径向=0, 后掠<0; 相对子午/径向方向)
    Z_blades : int
        主叶片数
    Z_split : int
        分裂/短叶片数 (无分流叶片填 0)

    返回
    ----
    sigma_W : float  (无有效叶片数时返回 1.0, 即无滑移)
    """
    Z = (Z_blades or 0) + (Z_split or 0)
    if Z <= 0:
        return 1.0
    chi = math.radians(chi_deg)
    cos_chi = max(math.cos(chi), 0.0)
    return 1.0 - math.sqrt(cos_chi) / (Z ** 0.7)


def estimate_exit_flow_angle(chi_deg, c_m3_bar, Z_blades, Z_split=0):
    """由式 (14) 等号关系用 Wiesner 滑移估算出口相对气流角 β3:
        tan β3 = tan χ3 + (σ_W - 1) / c_m3_bar

    c_m3_bar : float
        出口无量纲子午速度 c_m3/u3
    返回 β3 [°]; 若 c_m3_bar 趋零则返回金属角 χ3。
    """
    sigma = wiesner_slip_factor(chi_deg, Z_blades, Z_split)
    if abs(c_m3_bar) < 1e-10:
        return chi_deg
    tan_beta3 = math.tan(math.radians(chi_deg)) + (sigma - 1.0) / c_m3_bar
    return math.degrees(math.atan(tan_beta3))


# ──────────────────────────────
# 单行 CFX 数据 → 参数字典
# ──────────────────────────────
def cfx_row_to_params(row, base_params):
    """将 CFX CSV 的一行 (dict) 转换为完整的稳定性参数字典

    参数
    ----
    row : dict
        CSV 行, 键为列名 (见模块文档)
    base_params : dict
        从 NCEPU_params.txt 加载的基础参数 (含几何)

    返回
    ----
    params : dict
        完整参数字典, 可直接传给 build_cc3_characteristic_matrix
    extra : dict
        额外诊断量 (u3, Q, G, pi_tt, dh_t 等)
    """
    params = dict(base_params)

    R_3 = params['R_3']
    rpm = float(row['rpm'])

    # 根据转速符号和 CFD 的 Z 轴指向, 把切向速度统一翻成
    # "旋转方向为正、u3>0、后掠 χ<0" 的模型约定
    z_axis = params.get('z_axis', 'down')
    c_u_sign = circumferential_sign(rpm, z_axis)
    omega = abs(rpm) * math.pi / 30.0
    u3 = omega * R_3   # 恒为正

    # ── STA2 叶轮进口 ──
    c_m2 = float(row['c_m2'])
    c_u2 = c_u_sign * float(row.get('c_u2', 0.0) or 0.0)
    rho2 = float(row['rho2'])

    # 进口相对周向速度 w_theta = c_u - u(r), r=r_2rms=R_3*radial_sta2
    u_sta2 = u3 * params.get('radial_sta2', 0.0)
    w_theta2 = c_u2 - u_sta2
    w2 = math.hypot(c_m2, w_theta2)            # 相对速度大小(由分量算)
    beta2 = math.degrees(math.atan2(w_theta2, c_m2)) if c_m2 > 0 else 0.0
    alpha2 = math.degrees(math.atan2(c_u2, c_m2)) if c_m2 > 0 else 0.0

    # ── STA3 叶轮出口 ──
    c_m3 = float(row['c_m3'])
    rho3 = float(row['rho3'])

    # c_u3 由 CSV 给出(随切向约定翻转); 缺失时用 Wiesner 滑移因子反推
    # (c_u3/u3 = σ_W + c_m3_bar·tan χ3)
    chi3_deg = params.get('beta_blade_3_deg', 0.0)
    Z_blades = params.get('Z_blades', 0)
    Z_split = params.get('Z_split', 0)
    c_m3_bar = c_m3 / u3
    if row.get('c_u3', '') not in ('', None):
        c_u3 = c_u_sign * float(row['c_u3'])
    else:
        sigma = wiesner_slip_factor(chi3_deg, Z_blades, Z_split)
        c_u3 = u3 * (sigma + c_m3_bar * math.tan(math.radians(chi3_deg)))

    alpha3 = math.degrees(math.atan2(c_u3, c_m3)) if c_m3 > 0 else 0.0

    # 叶轮出口(STA3)相对气流角 β3, 即 Matrix 叶轮出口相对角 beta2_imp。
    # 出口 r=R_3, 牵连速度=u3, 故 w_theta3 = c_u3 - u3。
    w_theta3 = c_u3 - u3
    w3 = math.hypot(c_m3, w_theta3)
    beta3 = math.degrees(math.atan2(w_theta3, c_m3)) if c_m3 > 0 else 0.0

    # ── STA4 扩压器出口 ──
    c_m4 = float(row['c_m4'])
    c_u4 = c_u_sign * float(row['c_u4'])
    rho4 = float(row['rho4']) if row.get('rho4', '') not in ('', None) else rho3

    # ── 热力学量 (静压/声速), 用于能量分布后处理 ──
    # 静声速由静温 T、静压 p 经 CoolProp 查真实气体值 (sCO2 不可用 γ=1.4)。
    fluid = params.get('fluid', 'CO2')
    p2 = float(row['p2']) if row.get('p2', '') not in ('', None) else 0.0
    p3 = float(row['p3']) if row.get('p3', '') not in ('', None) else 0.0
    p4 = float(row['p4']) if row.get('p4', '') not in ('', None) else p3

    a2 = row_sound_speed(row, 2, fluid)
    a3 = row_sound_speed(row, 3, fluid)
    a4 = row_sound_speed(row, 4, fluid) if row.get('T4', '') not in ('', None) else a3

    # 叶轮总焓升 dh_t = h(Tt3,pt3) - h(Tt2,pt2), 供损失导数自动差分之用
    h_t2 = row_enthalpy(row, 2, fluid)
    h_t3 = row_enthalpy(row, 3, fluid)
    dh_t = (h_t3 - h_t2) if (h_t2 > 0 and h_t3 > 0) else 0.0

    # ── 无量纲速度 → Matrix 叶轮接口键
    #    (Matrix 约定: 1=叶轮进口, 2=叶轮出口; 与站位号不同, 保持不动)
    params['Vx_bar1_imp'] = c_m2 / u3
    params['Vtheta_bar1_imp'] = c_u2 / u3
    params['Vr_bar2_imp'] = c_m3 / u3
    params['Vtheta_bar2_imp'] = c_u3 / u3

    params['beta1_imp_deg'] = beta2
    params['beta2_imp_deg'] = beta3
    params['alpha1_imp_deg'] = alpha2
    params['alpha2_imp_deg'] = alpha3

    # ── 上游基流 ──
    if params.get('Vx_up', 0.0) == 0.0:
        params['Vx_up'] = c_m2 / u3

    # ── 无叶扩压器常数 Q, G (进出口平均) ──
    r3_norm = 1.0                      # STA3 = R_3/R_3 = 1
    r4_norm = params['radial_sta4']    # STA4
    Q3 = (c_m3 / u3) * r3_norm
    G3 = (c_u3 / u3) * r3_norm
    Q4 = (c_m4 / u3) * r4_norm
    G4 = (c_u4 / u3) * r4_norm
    Q = 0.5 * (Q3 + Q4)
    G = 0.5 * (G3 + G4)
    params['Q'] = Q
    params['G'] = G
    params['Vx_plenum'] = Q / params['r_diffuser_exit']
    params['Vtheta_plenum'] = G / params['r_diffuser_exit']

    # ── 含密度面积比 AR_imp = ρ3·A3/(ρ2·A2) ──
    if params.get('AR_imp', 0.0) == 0.0 and rho2 > 0 and rho3 > 0:
        A2 = params['area_sta2'] * (R_3 ** 2)
        A3 = params['area_sta3'] * (R_3 ** 2)
        params['AR_imp'] = (rho3 * A3) / (rho2 * A2)

    # ── 惯性系数 lambda_imp ──
    # λ = s·AR·ln(AR)/(AR−1),  s = l_comp/R_3
    if params.get('lambda_imp', 0.0) == 0.0:
        AR = params['AR_imp']
        s_imp = params.get('s_imp', 0.0)
        if s_imp > 0 and AR > 0:
            if abs(AR - 1.0) < 1e-10:
                params['lambda_imp'] = s_imp
            else:
                params['lambda_imp'] = (s_imp * AR * math.log(AR)
                                        / (AR - 1.0))

    # ── 损失时间滞后 tau_imp ──
    # τu 为经验系数 (在 params.txt 中人工给定, 通常 1.0~1.25);
    # τ_imp = τu·2s_imp/(Ŵ1+Ŵ2), Ŵ 为进/出口无量纲相对速度
    if params.get('tau_imp', 0.0) == 0.0 and w2 > 0:
        tau_u = params.get('tau_u', 1.1)
        s_imp = params.get('s_imp', 0.0)
        w1_bar = w2 / u3
        w2_bar = w3 / u3 if w3 > 0 else w1_bar
        if s_imp > 0 and (w1_bar + w2_bar) > 0:
            params['tau_imp'] = tau_u * 2.0 * s_imp / (w1_bar + w2_bar)

    # ── 损失导数 dL/dtanβ1 ──
    if row.get('dL_dtanbeta1', '') not in ('', None):
        params['dL_dtanbeta1_imp'] = float(row['dL_dtanbeta1'])

    # ── 计算 tan 角度 ──
    params['tan_beta1_imp'] = math.tan(math.radians(params['beta1_imp_deg']))
    params['tan_beta2_imp'] = math.tan(math.radians(params['beta2_imp_deg']))
    params['tan_alpha1_imp'] = math.tan(math.radians(params['alpha1_imp_deg']))

    # ── 各站位热力学量 (供能量分布后处理使用) ──
    params['acoustics'] = {
        'u3': u3,
        'rho': {'sta2': rho2, 'sta3': rho3, 'sta4': rho4},
        'a':   {'sta2': a2,   'sta3': a3,   'sta4': a4},
        'p':   {'sta2': p2,   'sta3': p3,   'sta4': p4},
    }

    # π_tt 由进出口总压 pt3/pt2 计算 (仅用于显示, 不进矩阵)
    pt2_val = float(row['pt2']) if row.get('pt2', '') not in ('', None) else 0.0
    pt3_val = float(row['pt3']) if row.get('pt3', '') not in ('', None) else 0.0
    pi_tt = (pt3_val / pt2_val) if pt2_val > 0 and pt3_val > 0 else None

    extra = {
        'u3': u3,
        'c_u_sign': c_u_sign,
        'Q': Q, 'G': G,
        'pi_tt': pi_tt,
        'dh_t': dh_t,
    }
    return params, extra


# ──────────────────────────────
# 损失导数自动差分
# ──────────────────────────────
def compute_dL_dtanbeta1_from_csv(rows, idx, u3, r2_norm=1.0,
                                  c_u_sign=1.0, fluid='CO2'):
    """从 CSV 中查找同转速相邻工况, 自动差分 dL/dtanβ1

    邻点选取: 取同转速下高流量侧和低流量侧各自最近的点 (不设流量容差)。
    两侧都有时做中心差分; 端点工况只有一侧可用时退化为单侧差分
    (最小流量点用高流量侧前向迎风, 最大流量点用低流量侧后向差分)。

    采用 Spakovszky 间接法 (与 Eckardt_O_Geo_Calculator 一致):
        μ = Δh_t / u3²                        (功系数)
        ψ = (pt3 − pt2) / (ρ2 · u3²)          (压力系数)
        L = μ − ψ                             (归一化总压损失)
        dL/dtanβ1 = (Δμ − Δψ) / Δtanβ1

    Δh_t 由总温 Tt、总压 pt 经 CoolProp 查总焓后相减得到。

    参数
    ----
    rows : list[dict]
        所有 CSV 行
    idx : int
        目标行索引
    u3 : float
        叶轮出口轮缘速度 ω·R_3 [m/s] (同转速下恒定)
    r2_norm : float
        STA2 归一化半径 r_2rms/R_3, 用于把周向速度换算成相对气流角
    c_u_sign : float
        切向速度符号翻译 (由 circumferential_sign 给出), u3 为幅值时
        须把 CSV 的 c_u2 同步翻转到模型约定。
    fluid : str
        CoolProp 工质名 (默认 'CO2')。

    需要 CSV 包含以下列用于差分:
      pt2, pt3, Tt2, Tt3  进出口总压 [Pa] / 总温 [K] (查焓)
      rho2                进口静密度 [kg/m³]
      c_m2, c_u2          进口子午/周向速度 (算相对角)

    若数据不足, 返回 None。
    """
    target = rows[idx]
    m0 = float(target['m'])
    rpm0 = float(target['rpm'])

    needed = ['pt2', 'pt3', 'rho2']
    if not all(k in target and target[k] != '' for k in needed):
        return None

    def find_neighbor(side):
        """在同转速中找目标流量侧 (side=+1 大流量/-1 小流量) 距目标最近的点
        (不设流量容差, 取该侧绝对流量差最小者)。"""
        best = None
        best_dm = None
        for r in rows:
            if float(r['rpm']) != rpm0:
                continue
            dm = float(r['m']) - m0
            if side * dm <= 0:          # 必须严格在目标点的指定一侧
                continue
            if best_dm is None or abs(dm) < abs(best_dm):
                best, best_dm = r, dm
        return best

    row_plus = find_neighbor(+1)
    row_minus = find_neighbor(-1)
    if row_plus is None and row_minus is None:
        return None

    u_sta2 = u3 * r2_norm

    def point_vals(r):
        c_m2 = float(r['c_m2'])
        # u3 已取幅值(正), 原始 c_u2 须先翻成模型切向约定;
        # w_theta = c_u2 - u_sta2, β2 = atan2(w_theta, c_m2)
        c_u2 = c_u_sign * float(r.get('c_u2', 0.0) or 0.0)
        w_theta2 = c_u2 - u_sta2
        beta2 = math.atan2(w_theta2, c_m2) if c_m2 > 0 else 0.0
        tanb1 = math.tan(beta2)

        dpt = float(r['pt3']) - float(r['pt2'])
        rho2 = float(r['rho2'])
        psi = dpt / (rho2 * u3 ** 2)

        # 总焓升由 CoolProp 查总焓相减
        delta_h = row_enthalpy(r, 3, fluid) - row_enthalpy(r, 2, fluid)
        mu = delta_h / (u3 ** 2) if delta_h != 0.0 else 0.0
        return tanb1, mu, psi

    t0, mu0, psi0 = point_vals(target)

    if row_plus is not None and row_minus is not None:
        # 中心差分: (高流量点 - 低流量点) / (t+ - t-)
        tp, mup, psip = point_vals(row_plus)
        tm, mum, psim = point_vals(row_minus)
        if abs(tp - tm) < 1e-12:
            return None
        dmu = (mup - mum) / (tp - tm)
        dpsi = (psip - psim) / (tp - tm)
    else:
        # 单侧差分 (端点工况): 只有一个邻点, 相对该邻点求斜率;
        # 最小流量点优先用高流量侧迎风 (前向), 最大流量点用低流量侧 (后向)
        r_nb = row_plus if row_plus is not None else row_minus
        tn, mun, psin = point_vals(r_nb)
        if abs(tn - t0) < 1e-12:
            return None
        dmu = (mun - mu0) / (tn - t0)
        dpsi = (psin - psi0) / (tn - t0)

    return dmu - dpsi


# ──────────────────────────────
# CSV 读取
# ──────────────────────────────
def load_cfx_csv(csv_path):
    """读取 CFX 导出的 CSV, 返回 list[dict]"""
    rows = []
    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for r in reader:
            # 跳过空行
            if not r.get('rpm') or not r.get('m'):
                continue
            rows.append(r)
    return rows


def load_all_operating_points(csv_path=DEFAULT_CSV_FILE,
                              params_file=DEFAULT_PARAMS_FILE):
    """加载 CSV 中全部工况点, 返回 list of (params, extra)

    返回
    ----
    results : list of dict
        每个元素: {'rpm', 'm', 'params', 'extra'}
    """
    base_params = load_params(params_file)
    base_params['z_axis'] = read_z_axis(params_file)
    base_params['fluid'] = read_fluid(params_file)
    base_params = compute_geometry(base_params)

    rows = load_cfx_csv(csv_path)
    results = []
    for row in rows:
        params, extra = cfx_row_to_params(row, base_params)

        # 尝试自动差分损失导数 (同转速 ±2% 流量邻点)
        if params.get('dL_dtanbeta1_imp', 0.0) == 0.0:
            idx = rows.index(row)
            u3 = extra['u3']
            r2_norm = params.get('radial_sta2', 1.0)
            dL = compute_dL_dtanbeta1_from_csv(
                rows, idx, u3, r2_norm,
                c_u_sign=extra['c_u_sign'],
                fluid=base_params.get('fluid', 'CO2'))
            if dL is not None:
                params['dL_dtanbeta1_imp'] = dL

        results.append({
            'rpm': float(row['rpm']),
            'm': float(row['m']),
            'params': params,
            'extra': extra,
        })
    return results


# ──────────────────────────────
# CSV 模板生成
# ──────────────────────────────
CSV_TEMPLATE_HEADER = [
    'rpm', 'm',
    # ── STA2 叶轮进口 ── (c_m/c_u/rho massFlowAveraged; p/T/pt/Tt areaAveraged)
    'c_m2', 'c_u2', 'rho2', 'p2', 'T2', 'pt2', 'Tt2',
    # ── STA3 叶轮出口 ──
    'c_m3', 'c_u3', 'rho3', 'p3', 'T3', 'pt3', 'Tt3',
    # ── STA4 扩压器出口 ──
    'c_m4', 'c_u4', 'rho4', 'p4', 'T4', 'pt4', 'Tt4',
    # 损失导数: 直填 dL_dtanbeta1, 或留空给同转速 ±2% 邻点 (pt/Tt 查焓) 自动差分
    'dL_dtanbeta1',
]


def write_csv_template(path=DEFAULT_CSV_FILE):
    """写出干净的 CSV 模板 (仅表头, 无示例数据; 数据由用户从 CFX 提取后填入)"""
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(CSV_TEMPLATE_HEADER)
    print(f"[template] {path}")


# ──────────────────────────────
# 命令行入口
# ──────────────────────────────
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='NCEPU CFX CSV 加载器')
    parser.add_argument('--csv', default=DEFAULT_CSV_FILE,
                        help='CFX 气动参数 CSV 路径')
    parser.add_argument('--params', default=DEFAULT_PARAMS_FILE,
                        help='几何参数文件路径')
    parser.add_argument('--template', action='store_true',
                        help='生成 CSV 模板后退出')
    args = parser.parse_args()

    if args.template:
        write_csv_template(args.csv)
        raise SystemExit(0)

    if not os.path.exists(args.csv):
        print(f"CSV not found: {args.csv}")
        print("Use --template to generate a template.")
        raise SystemExit(1)

    results = load_all_operating_points(args.csv, args.params)
    print(f"Loaded {len(results)} operating points from {args.csv}\n")
    for r in results:
        p = r['params']
        print(f"  rpm={r['rpm']:.0f}  m={r['m']:.3f} kg/s  "
              f"u3={r['extra']['u3']:.1f} m/s")
        print(f"    c_m2/u3={p['Vx_bar1_imp']:.4f}  "
              f"c_u3/u3={p['Vtheta_bar2_imp']:.4f}  "
              f"β1={p['beta1_imp_deg']:.2f}°  β3(exit)={p['beta2_imp_deg']:.2f}°  "
              f"α3={p['alpha2_imp_deg']:.2f}°")
        print(f"    Q={p['Q']:.4f}  G={p['G']:.4f}  "
              f"AR={p['AR_imp']:.4f}  λ={p['lambda_imp']:.4f}  "
              f"τ={p.get('tau_imp', 0.0):.4f}")
        print(f"    dL/dtanβ1={p.get('dL_dtanbeta1_imp', 0.0):.6f}")
        if r['extra']['pi_tt']:
            print(f"    π_tt={r['extra']['pi_tt']:.4f}")
        print()
