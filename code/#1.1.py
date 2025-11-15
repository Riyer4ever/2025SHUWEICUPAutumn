import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 定义变量
m = 2000          # kg
k = 7225200         # N/m
c = 3600          # N·s/m
dt = 0.01        # 步长
t = np.arange(0, 10, dt)
F = pd.read_excel(r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\#1timeAndForce.xlsx", usecols='B').squeeze().to_numpy()  # 扰动力

# 初值定义
y = np.zeros_like(t)
v = np.zeros_like(t)
y[0], v[0] = 0.01, 0          # y(0)=0.01 m, v(0)=0

# 加速度
def acc(y, v, F_val):
    return (F_val - c*v - k*y) / m

# RK4 主循环
for i in range(len(t)-1):
    Fi = F[i]
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

# 计算加速度
a = (F - c*v - k*y) / m

# 振动指标
I_h = np.mean(a**2)

# 制图
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(t, y)
plt.xlabel('time / s'); plt.ylabel('y / m'); plt.title('Lateral displacement (RK4)')

plt.subplot(1, 2, 2)
plt.plot(t, a)
plt.xlabel('time / s'); plt.ylabel('ÿ / m·s⁻²'); plt.title('Lateral acceleration (RK4)')
plt.tight_layout()
plt.show()

# 输出横向振动指标
print(f'RK4 无控振动指标  I_h = {I_h:.4f} (m²/s⁴)')