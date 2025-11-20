import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 原始参数值
m = 2000  # 质量 (kg)
k = 7225200  # 刚度 (N/m)
c = 3600  # 阻尼 (N·s/m)

# 蒙特卡洛模拟次数
num_simulations = 100

# 存储每次模拟后的振动指标 I_h
I_h_values = []

# 定义系统模型
def system_model(m, c, k, F_dist, dt, T):
    t = np.arange(0, 10, dt)
    y = np.zeros_like(t)
    v = np.zeros_like(t)
    y[0], v[0] = 0.01, 0  # 初始条件
    for i in range(len(t) - 1):
        a = (F_dist[i] - c * v[i] - k * y[i]) / m
        v[i + 1] = v[i] + a * dt
        y[i + 1] = y[i] + v[i + 1] * dt
    # 确保最后一个点也被更新
    a_last = (F_dist[-1] - c * v[-1] - k * y[-1]) / m
    v[-1] = v[-2] + a_last * dt
    y[-1] = y[-2] + v[-1] * dt
    a = (F_dist - c * v - k * y) / m
    I_h = np.mean(a**2)
    return I_h

# 读取扰动力 F_dist
file_path = r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\#2timeAndForce.xlsx"
F_dist_df = pd.read_excel(file_path, usecols='B')  # 假设扰动力在第 B 列
F_dist = F_dist_df.squeeze().to_numpy()

# 检查数据长度
T = 9.99  # 总时间 (s)
dt = 0.01  # 时间步长 (s)
expected_length = int(T / dt) + 1  # 时间数组长度应该是 1000
if len(F_dist) != expected_length:
    raise ValueError(f"扰动力数据长度 {len(F_dist)} 不匹配，应为 {expected_length}")

# 进行蒙特卡洛模拟
for _ in range(num_simulations):
    # 随机扰动参数
    m_perturbed = m * (1 + np.random.uniform(-0.1, 0.1))
    k_perturbed = k * (1 + np.random.uniform(-0.15, 0.15))
    c_perturbed = c * (1 + np.random.uniform(-0.2, 0.2))
    
    # 计算振动指标 I_h
    I_h = system_model(m_perturbed, c_perturbed, k_perturbed, F_dist, dt, T)
    I_h_values.append(I_h)

# 绘制箱线图
plt.figure(figsize=(10, 6))
plt.boxplot(I_h_values, labels=[r'Vibration Index  $I_h$ '])
plt.title(r'Boxplot of Vibration Index $I_h$ under Parameter Perturbations')
plt.ylabel(r'Vibration Index $I_h$ (m²/s⁴)')
plt.grid(True)
plt.show()