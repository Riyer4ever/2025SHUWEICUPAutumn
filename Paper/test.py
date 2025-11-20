import numpy as np
import pandas as pd

# 系统参数
m_body = 2000          # 车体质量（不含作动器） kg
m_ecc = 100            # 单个偏心块质量 kg
r_ecc = 0.2            # 偏心块旋转半径 m
k = 7225200            # 等效刚度系数 N/m
c = 3600               # 等效阻尼系数 N·s/m
dt = 0.01              # 采样步长 s
T = 10                 # 总时长 s
t = np.arange(0, T+dt, dt)

# 读取场景1数据：横向扰动力
# 假设 Excel 文件名为 data.xlsx，Sheet 名为“场景1”
df1 = pd.read_excel(r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\data.xlsx", sheet_name='场景1')
F_disturb = df1['横向扰动力(单位：N)'].values

# 初始化数组
y = np.zeros_like(t)      # 横向位移
v = np.zeros_like(t)      # 横向速度
a = np.zeros_like(t)      # 横向加速度（相对轮子的）

# 初始条件
y[0] = 0.01
v[0] = 0

# 无作动器情况
for i in range(1, len(t)):
    F_net = F_disturb[i] - k*y[i-1] - c*v[i-1]
    a[i-1] = F_net / m_body
    v[i] = v[i-1] + a[i-1]*dt
    y[i] = y[i-1] + v[i]*dt
# 补上最后一个加速度
F_net = F_disturb[-1] - k*y[-1] - c*v[-1]
a[-1] = F_net / m_body
ih_no_act = np.mean(a**2)

# 读取作动器角度
theta11 = df1['离心式作动器1的第1组偏心块角度(单位：rad)'].values
theta12 = df1['离心式作动器1的第2组偏心块角度(单位：rad)'].values
theta21 = df1['离心式作动器2的第1组偏心块角度(单位：rad)'].values
theta22 = df1['离心式作动器2的第2组偏心块角度(单位：rad)'].values

# 重置状态
y[:] = 0; v[:] = 0; a[:] = 0
y[0] = 0.01

# 有作动器情况
for i in range(1, len(t)):
    # 计算作动器合力（两组对转，y方向叠加）
    F_act = m_ecc * r_ecc * (
        (np.cos(theta11[i]) + np.cos(theta12[i])) +
        (np.cos(theta21[i]) + np.cos(theta22[i]))
    ) * (100**2)   # 假设最大角速度 100 rad/s
    F_net = F_disturb[i] + F_act - k*y[i-1] - c*v[i-1]
    a[i-1] = F_net / m_body
    v[i] = v[i-1] + a[i-1]*dt
    y[i] = y[i-1] + v[i]*dt
F_net = F_disturb[-1] + F_act - k*y[-1] - c*v[-1]
a[-1] = F_net / m_body
ih_with_act = np.mean(a**2)

print(ih_no_act, ih_with_act)