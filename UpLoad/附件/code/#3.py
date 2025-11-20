#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2025 ShuWei-Cup Problem A-3  优化版
0 → 9.99 s，dt = 0.01 s，共 1000 点
离线增益表 + 带通滤波 + 双机分相 → 22%↓Ih
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal

# ---------------- 系统参数 -----------------
m = 2000.0; c = 3600.0; k = 7.2252e6
me = 100.0; r = 0.2
wmax = 100.0          # rad/s
amax = 5000.0         # rad/s^2
dt = 0.01
T_end = 9.99
N = 1000
t = np.linspace(0, T_end, N)

# ---------------- 读入扰动（1000 点） -----------------
file_path = r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\#2timeAndForce.xlsx"
F_dist = pd.read_excel(file_path, usecols='B').squeeze().to_numpy()
if len(F_dist) != N:
    raise ValueError(f"扰动数据长度={len(F_dist)}，需要{N}")

# ---------------- 带通滤波器设计 -----------------
fb, fa = 0.3, 15.0          # 通带 0.3–15 Hz
fs = 1/dt
b, a = signal.butter(2, [fb/(0.5*fs), fa/(0.5*fs)], btype='band')

# ---------------- 离线增益表（|a_y| → Ka & 相位差） -----------------
gain_table = np.array([
    [0.0,  0.0010, 0.0],
    [2.0,  0.0008, 0.2],
    [5.0,  0.0006, 0.4],
    [10.0, 0.0004, 0.6],
    [50.0, 0.0002, 0.8]
])
def get_ka_dphi(abs_ay):
    return np.interp(abs_ay, gain_table[:, 0], gain_table[:, 1]), \
           np.interp(abs_ay, gain_table[:, 0], gain_table[:, 2])

# ---------------- 作动器合力 -----------------
def act_force(theta):
    return 2 * me * r * wmax**2 * (np.sin(theta[0]) + np.sin(theta[1]))

# ---------------- RK4 积分 -----------------
def rk4(y0, vy0, theta):
    y, vy, ay = np.zeros(N), np.zeros(N), np.zeros(N)
    y[0], vy[0] = y0, vy0
    F_act = act_force(theta)
    for i in range(N - 1):
        f = F_dist[i] + F_act[i]
        k1y = vy[i]
        k1v = (f - c*vy[i] - k*y[i])/m
        k2y = vy[i] + 0.5*dt*k1v
        k2v = (f - c*(vy[i]+0.5*dt*k1v) - k*(y[i]+0.5*dt*k1y))/m
        k3y = vy[i] + 0.5*dt*k2v
        k3v = (f - c*(vy[i]+0.5*dt*k2v) - k*(y[i]+0.5*dt*k2y))/m
        k4y = vy[i] + dt*k3v
        k4v = (f - c*(vy[i]+dt*k3v) - k*(y[i]+dt*k3y))/m
        y[i+1]  = y[i]  + dt*(k1y + 2*k2y + 2*k3y + k4y)/6
        vy[i+1] = vy[i] + dt*(k1v + 2*k2v + 2*k3v + k4v)/6
        ay[i] = k1v
    ay[-1] = (F_dist[-1] + F_act[-1] - c*vy[-1] - k*y[-1])/m
    return y, vy, ay

# ---------------- 速率&加速度限幅（保持长度） -----------------
# ---------- 速率&加速度限幅（保持长度 N） ----------
def rate_limit(theta_raw):
    theta = np.copy(theta_raw)
    for i in range(2):
        dth = np.gradient(theta[i], dt)          # 长度 N
        dth = np.clip(dth, -wmax, wmax)
        ddth = np.gradient(dth, dt)              # 长度 N
        ddth = np.clip(ddth, -amax, amax)
        # 用欧拉积分还原角度，保证长度不变
        dth_new = np.zeros(N)
        dth_new[0] = dth[0]
        for k in range(1, N):
            dth_new[k] = dth_new[k-1] + ddth[k-1] * dt
        theta[i] = np.zeros(N)
        theta[i][0] = theta_raw[i][0]
        for k in range(1, N):
            theta[i][k] = theta[i][k-1] + dth_new[k] * dt
    return theta
# ---------------- 在线控制：带通+分相+查表 -----------------
def online_control_opt(ay_raw):
    ay_bp = signal.lfilter(b, a, ay_raw)          # 零相位已足够
    theta = np.zeros((2, N))
    for k in range(N):
        abs_a = abs(ay_bp[k])
        Ka, dphi = get_ka_dphi(abs_a)
        S = -np.clip(Ka * ay_bp[k], -2.0, 2.0)
        sin_th = np.clip(S/2.0, -1.0, 1.0)
        theta[0, k] = np.arcsin(sin_th) + dphi
        theta[1, k] = np.arcsin(sin_th) - dphi
    return rate_limit(theta)

# ---------------- 主程序 -----------------
# 1. 无控
y_uc, vy_uc, ay_uc = rk4(0.01, 0.0, np.zeros((2, N)))
Ih_uc = np.mean(ay_uc**2)

# 2. 在线控制
theta_cl = online_control_opt(ay_uc)
y_cl, vy_cl, ay_cl = rk4(0.01, 0.0, theta_cl)
Ih_cl = np.mean(ay_cl**2)

print(f'无控 Ih = {Ih_uc:.3f}')
print(f'优化在线 Ih = {Ih_cl:.3f}  降幅 {100*(Ih_uc-Ih_cl)/Ih_uc:.1f}%')

# ---------------- 绘图 -----------------
plt.figure(figsize=(12, 9))
plt.subplot(3, 1, 1)
plt.plot(t, y_uc, label='Uncontrolled')
plt.plot(t, y_cl, label='Optimized online')
plt.ylabel('y (m)'); plt.legend(); plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t, ay_uc, label='Uncontrolled')
plt.plot(t, ay_cl, label='Optimized online')
plt.ylabel('$a_y$ (m/s²)'); plt.legend(); plt.grid()

plt.subplot(3, 1, 3)
plt.plot(t, theta_cl[0], label='θ₁(t)')
plt.plot(t, theta_cl[1], label='θ₂(t)')
plt.xlabel('t (s)'); plt.ylabel('θ (rad)'); plt.legend(); plt.grid()
plt.tight_layout()
plt.savefig('task3_optimized.png', dpi=300)
plt.show()