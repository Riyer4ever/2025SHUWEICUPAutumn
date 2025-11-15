import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 定义变量
m   = 2000          # 车体等效质量 kg
k   = 7225200       # 等效横向刚度 N/m
c   = 3600          # 等效横向阻尼 N·s/m
dt  = 0.01          # 步长 s
t   = np.arange(0, 10, dt)          # 0–10 s

m_block = 100       # 单块偏心质量 kg
r       = 0.2       # 旋转半径 m

# 1. 扰动力：单列 B 1000 点
F_dist = pd.read_excel(r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\timeAndForce.xlsx", usecols='B').squeeze().to_numpy()

# 2. 四块独立角度：五列 [time, theta1, theta2, theta3, theta4]（rad）
ang_df = pd.read_excel(r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\theta.xlsx")
theta = ang_df[['theta1', 'theta2', 'theta3', 'theta4']].to_numpy().T  # 4×N
omega = np.gradient(theta, dt, axis=1)                                # 4×N

# 总横向离心力 = 四块分量求和
F_act = np.sum(m_block * r * omega**2 * np.sin(theta), axis=0)  # (N,)
F_total = F_dist + F_act

# RK4方法做积分
y = np.zeros_like(t)
v = np.zeros_like(t)
y[0], v[0] = 0.01, 0          # 初始条件

def acc(y, v, F_val):
    return (F_val - c*v - k*y) / m

for i in range(len(t)-1):
    Fi = F_total[i]
    k1y = v[i]
    k1v = acc(y[i], v[i], Fi)

    k2y = v[i] + 0.5*dt*k1v
    k2v = acc(y[i] + 0.5*dt*k1y, v[i] + 0.5*dt*k1v, Fi)

    k3y = v[i] + 0.5*dt*k2v
    k3v = acc(y[i] + 0.5*dt*k2y, v[i] + 0.5*dt*k2v, Fi)

    k4y = v[i] + dt*k3v
    k4v = acc(y[i] + dt*k3y, v[i] + dt*k3v, Fi)

    y[i+1] = y[i] + (dt/6)*(k1y + 2*k2y + 2*k3y + k4y)
    v[i+1] = v[i] + (dt/6)*(k1v + 2*k2v + 2*k3v + k4v)

# 加速度（公式）
a = (F_total - c*v - k*y) / m

# 横向振动指标
I_h = np.mean(a**2)

# 画图
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(t, y)
plt.xlabel('time / s'); plt.ylabel('y / m')
plt.title('Lateral displacement (4-independent blocks)')

plt.subplot(1, 2, 2)
plt.plot(t, a)
plt.xlabel('time / s'); plt.ylabel('ÿ / m·s⁻²')
plt.title('Lateral acceleration (4-independent blocks)')
plt.tight_layout()
plt.show()

# 输出横向振动指数
print(f'四块独立角度 — 有控横向振动指标  I_h = {I_h:.6f} (m²/s⁴)')