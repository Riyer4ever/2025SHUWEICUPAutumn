#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2025 ShuWei-Cup Problem A-3  （1000 点版）
0 → 9.99 s，dt = 0.01 s，共 1000 个采样点
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ---------------- 系统参数 -----------------
m = 2000.0
c = 3600.0
k = 7.2252e6
me = 100.0
r = 0.2
wmax = 100.0          # rad/s
amax = 5000.0         # rad/s^2

# ****** 改为 1000 点 ******
dt = 0.01
T_end = 9.99          # 最大时刻 9.99 s
N = 1000              # 1000 个点
t = np.linspace(0, T_end, N)

# ---------------- 读入最原始扰动（1000 点） -----------------
file_path = r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\#2timeAndForce.xlsx"
F_dist_raw = pd.read_excel(file_path, usecols='B').squeeze().to_numpy()
# 强制 1000 点
if len(F_dist_raw) != N:
    raise ValueError(f"原始数据长度={len(F_dist_raw)}，需要{N}（1000 点）")
F_dist = F_dist_raw

# ---------------- 作动器出力 -----------------
def act_force(theta):
    # theta: (2, N)
    return 2 * me * r * wmax**2 * (np.sin(theta[0]) + np.sin(theta[1]))

# ---------------- RK4 积分 -----------------
def rk4(y0, vy0, theta):
    y, vy, ay = np.zeros(N), np.zeros(N), np.zeros(N)
    y[0], vy[0] = y0, vy0
    F_act = act_force(theta)
    for i in range(N - 1):
        f = F_dist[i] + F_act[i]
        k1y = vy[i]
        k1v = (f - c * vy[i] - k * y[i]) / m
        k2y = vy[i] + 0.5 * dt * k1v
        k2v = (f - c * (vy[i] + 0.5 * dt * k1v) - k * (y[i] + 0.5 * dt * k1y)) / m
        k3y = vy[i] + 0.5 * dt * k2v
        k3v = (f - c * (vy[i] + 0.5 * dt * k2v) - k * (y[i] + 0.5 * dt * k2y)) / m
        k4y = vy[i] + dt * k3v
        k4v = (f - c * (vy[i] + dt * k3v) - k * (y[i] + dt * k3y)) / m
        y[i+1]  = y[i]  + dt * (k1y + 2*k2y + 2*k3y + k4y) / 6
        vy[i+1] = vy[i] + dt * (k1v + 2*k2v + 2*k3v + k4v) / 6
        ay[i] = k1v
    ay[-1] = (F_dist[-1] + F_act[-1] - c * vy[-1] - k * y[-1]) / m
    return y, vy, ay

# ---------------- 速率 & 加速度限幅 -----------------
def rate_limit(theta_raw):
    """保持长度 N 的速率+加速度限幅"""
    theta = np.copy(theta_raw)
    for i in range(2):
        dth = np.diff(theta[i]) / dt
        dth = np.clip(dth, -wmax, wmax)
        ddth = np.diff(dth) / dt
        ddth = np.clip(ddth, -amax, amax)
        # 还原积分（补回长度）
        dth_new = np.concatenate(([0], np.cumsum(ddth * dt))) + dth[0]
        theta[i] = np.concatenate(([theta[i][0]], np.cumsum(dth_new * dt))) + theta[i][0]
    return theta

# ---------------- 在线控制算法 -----------------
Ka_opt = 0.00087   # 离线标定最优增益
alpha  = 0.8       # 低通滤波器

def online_control(ay_meas):
    """返回 theta(2, N) 仅基于 ay_meas"""
    theta = np.zeros((2, N))
    afilt = 0.0
    for k in range(N):
        ay_k = ay_meas[k]
        afilt = alpha * afilt + (1 - alpha) * ay_k
        S     = -np.clip(Ka_opt * afilt, -2.0, 2.0)
        sin_th = np.clip(S / 2.0, -1.0, 1.0)
        theta[0, k] = np.arcsin(sin_th)
        theta[1, k] = np.arcsin(sin_th)
    # 速率 & 加速度限幅
    theta = rate_limit(theta)
    return theta

# ---------------- 主程序 -----------------
# 1. 无控仿真
y_uc, vy_uc, ay_uc = rk4(0.01, 0.0, np.zeros((2, N)))

# 2. 在线控制仿真（仅用 ay_uc 作为“实测”信号）
theta_cl = online_control(ay_uc)
y_cl, vy_cl, ay_cl = rk4(0.01, 0.0, theta_cl)

# 3. 指标
Ih_uc = np.mean(ay_uc**2)
Ih_cl = np.mean(ay_cl**2)
print(f'无控 Ih = {Ih_uc:.3f}')
print(f'在线控制 Ih = {Ih_cl:.3f}  降幅 {100*(Ih_uc-Ih_cl)/Ih_uc:.1f}%')

# ---------------- 绘图 -----------------
plt.figure(figsize=(12, 9))
plt.subplot(3, 1, 1)
plt.plot(t, y_uc, label='Uncontrolled displacement')
plt.plot(t, y_cl, label='Online displacement control')
plt.ylabel('y (m)'); plt.legend(); plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t, ay_uc, label='Uncontrolled acceleration')
plt.plot(t, ay_cl, label='Online control of acceleration')
plt.ylabel('a_y (m/s²)'); plt.legend(); plt.grid()



plt.tight_layout()
plt.savefig('scenario3_result.png', dpi=300)
plt.show()


plt.figure(figsize=(12, 9))

plt.subplot(3, 1, 1)
plt.plot(t, theta_cl[0], label='θ₁(t)')
plt.xlabel('t (s)'); plt.ylabel('Angular velocity (rad)'); plt.legend(); plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t, theta_cl[1], label='θ₂(t)')
plt.xlabel('t (s)'); plt.ylabel('Angular velocity (rad)'); plt.legend(); plt.grid()

plt.tight_layout()
plt.savefig('scenario3_result.png', dpi=300)
plt.show()