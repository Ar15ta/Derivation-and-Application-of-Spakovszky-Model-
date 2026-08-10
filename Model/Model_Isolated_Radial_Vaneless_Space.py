# Model_Isolated_Radial_Vaneless_Space_using_radial.py
# 使用 radial.py 中的 _T_n 复现论文图 2-17 和 2-18
import os
import sys

# 确保 Stability 根目录在 sys.path 中, 使 Matrix/ 等包可导入
_STABILITY_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _STABILITY_ROOT not in sys.path:
    sys.path.insert(0, _STABILITY_ROOT)

import numpy as np
import matplotlib.pyplot as plt
from Matrix.radial import _T_n   # 直接使用 radial 模块中的原始映射矩阵

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14

# ============================================================
# 参数设置 (完全同版本A)
# ============================================================
Q = 0.215
Gamma = 0.7032
n = 20
R2 = 1.0
R3 = 1.5
s = -1j * n                     # 相对坐标系（原始设置）

Wr_in = 1.0 + 0j
Wtheta_in = (Gamma - 1) / Q + 0j          # = -1.3809
Wtheta_in_abs = np.abs(Wtheta_in)         # 1.3809

# ============================================================
# 求解系数 [D, E, F]
# 注意：这里直接调用 _T_n 代替版本A的 Tn_rad
# 二者接口完全一致: _T_n(r, s, n, r0, Q, G)
# ============================================================
T_in = _T_n(R2, s, n, R2, Q, Gamma)
T_out = _T_n(R3, s, n, R2, Q, Gamma)
A = np.vstack([T_in[0, :], T_in[1, :], T_out[2, :]])
b = np.array([Wr_in, Wtheta_in, 0.0], dtype=complex)
DEF = np.linalg.solve(A, b)

# ============================================================
# 计算各半径处的状态 (仍使用 _T_n)
# ============================================================
r_vals = [1.00, 1.05, 1.10, 1.15, 1.20]
states = []
for r in r_vals:
    T = _T_n(r, s, n, R2, Q, Gamma)
    Wr, Wtheta, P = T @ DEF
    Vr_mean = Q / r
    Vtheta_mean = Gamma / r
    Pt = P + Vr_mean * Wr + Vtheta_mean * Wtheta
    states.append((r, Wr, Wtheta, P, Pt))
    # 调试：打印所有r处的速度值和相位
    data_Wr = Wr / 1.0
    data_Wtheta = Wtheta / Wtheta_in_abs
    print(f"r={r:.2f}: Wr={data_Wr:.6f} (相位={np.degrees(np.angle(data_Wr)):.1f}°), Wtheta={data_Wtheta:.6f} (相位={np.degrees(np.angle(data_Wtheta)):.1f}°)")

# ============================================================
# 绘图 (图2-17) - 归一化
# ============================================================
p_ref = (Q**2 + Gamma**2)            # 入口相对动压 
theta_deg = np.linspace(0, 18, 1000)  # 横坐标范围为0~18°

# 论文中说明：静压和总压扰动以入口处相对动压扰动幅值归一化。
#  inlet relative dynamic head perturbation: δ(½V²) = Vr δWr + Vθ δWθ.
# p_ref = abs(Q * Wr_in + Gamma * Wtheta_in)

# 调试：检查归一化参数
print(f"\n=== 归一化参数 ===")
print(f"Q = {Q}")
print(f"Gamma = {Gamma}")
print(f"Wtheta_in = {Wtheta_in}")
print(f"Wtheta_in_abs = {Wtheta_in_abs}")
print(f"p_ref (inlet dynamic-head perturbation) = {p_ref}")

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
titles = [r'$\delta W_r$', r'$\delta P_t$ (Relative Total Pressure Perturbation)',
          r'$\delta W_\theta$', r'$\delta P$ (Static Pressure Perturbation)']

for ax, title in zip(axes.flatten(), titles):
    for idx, (r, Wr, Wtheta, P, Pt) in enumerate(states):
        if 'W_r' in title:
            data = Wr / 1.0
        elif r'W_\theta' in title:  # 使用原始字符串避免转义
            data = Wtheta / Wtheta_in_abs
        elif 'Total Pressure' in title:
            data = Pt / p_ref
        else:
            data = P / p_ref
        # 扰动波形：幅值 * cos(nθ + φ + offset)
        # 添加90°相位偏移以匹配论文中的波峰位置（波峰从18°调整到14°）
        phase_offset = np.deg2rad(90)  # 90°相位偏移
        pattern = np.abs(data) * np.cos(n * np.radians(theta_deg) + np.angle(data) + phase_offset)
        line_kwargs = {'linewidth': 1.5}
        if abs(r - R2) < 1e-8:
            line_kwargs.update({'linewidth': 2.5, 'linestyle': '-'})
        line = ax.plot(theta_deg, pattern, **line_kwargs)[0]
        if idx < 2:
            label_x = 1.2 + idx * 2.5
            ha = 'left'
        else:
            label_x = 6.0 + (idx - 2) * 3.5
            ha = 'left'

        label_y = np.interp(label_x, theta_deg, pattern)
        ax.text(label_x, label_y, f'r={r:.2f}', color=line.get_color(), fontsize=14,
                verticalalignment='center', horizontalalignment=ha,
                bbox=dict(facecolor='white', edgecolor='white', pad=0.2, alpha=0.75))
    ax.axhline(0.0, color='black', linestyle='--', linewidth=2)
    ax.set_xlim(0, 18)
    ax.set_xlabel('θ (Deg)')
    ax.set_ylabel('Normalized Perturbation')
    ax.set_title(title)
    if 'Static Pressure' in title:
        ax.set_ylim(-0.15, 0.15)
plt.tight_layout()
plt.savefig('fig2_17_from_radial.png', dpi=150)
plt.show()

# ============================================================
# 图2-18: 幅值和相位差 (同样使用 _T_n)
# ============================================================
r_fine = np.linspace(R2, 1.25, 100)
mag_Wr_norm = []
mag_Wtheta_norm = []
phase_diff_deg = []

for r in r_fine:
    T = _T_n(r, s, n, R2, Q, Gamma)
    Wr, Wtheta, _ = T @ DEF
    mag_Wr_norm.append(np.abs(Wr) / 1.0)
    mag_Wtheta_norm.append(np.abs(Wtheta) / Wtheta_in_abs)
    diff = np.angle(Wtheta) - np.angle(Wr)
    diff_deg = np.degrees(diff) % 360.0
    phase_diff_deg.append(diff_deg)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))
ax1.plot(r_fine, mag_Wr_norm, 'o-', label=r'$|\delta W_r|$')
ax1.plot(r_fine, mag_Wtheta_norm, 's-', label=r'$|\delta W_\theta|$')
ax1.set_xlabel('r / R2')
ax1.set_ylabel('Magnitude')
ax1.set_title('Magnitude of Perturbations')
ax1.legend()
ax1.grid(True)
ax1.set_xlim(1.0, 1.25)

ax2.plot(r_fine, phase_diff_deg, 'd-', color='red')
ax2.set_xlabel('r / R2')
ax2.set_ylabel('Phase Difference (Deg)')
ax2.set_title('Phase Difference Between δWθ and δWr')
ax2.grid(True)
ax2.set_xlim(1.0, 1.25)
plt.tight_layout()
plt.savefig('fig2_18_from_radial.png', dpi=300)
plt.show()