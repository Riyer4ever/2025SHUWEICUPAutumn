import numpy as np
import matplotlib.pyplot as plt
import cma
import pandas as pd

# ---------------- 系统参数 -----------------
m = 2000.0
c = 3600.0
k = 7.2252e6
me = 100.0
r = 0.2
wmax = 100.0          # 最大角速度 rad/s
amax = 5000.0         # 最大角加速度 rad/s^2
T = 10.0
dt = 0.01
N = int(T/dt)+1
t = np.linspace(0, T, N)

# ---------------- 读入扰动 -----------------
F_dist = pd.read_excel(r"D:\workspace\SHUWEICUP2544531_2025.11.14\ProblemA\#2timeAndForce.xlsx", usecols="B").squeeze().to_numpy()   # shape (1001,)

# ---------------- 作动器出力 -----------------
def act_force(theta):
    """
    theta: (2,N) 两个作动器角位移
    返回 (N,) 总 y 向力
    """
    # 角速度/加速度限幅（后处理）
    return 2*me*r*wmax**2 * (np.sin(theta[0]) + np.sin(theta[1]))

# ---------------- RK4 积分 -----------------
def rk4(y0, vy0, theta):
    y = np.zeros(N)
    vy = np.zeros(N)
    ay = np.zeros(N)
    y[0], vy[0] = y0, vy0
    F_act = act_force(theta)
    for i in range(N-1):
        fdist = F_dist[i]
        f1v = (fdist + F_act[i] - c*vy[i] - k*y[i])/m
        f1y = vy[i]
        f2v = (fdist + F_act[i] - c*(vy[i]+0.5*dt*f1v) - k*(y[i]+0.5*dt*f1y))/m
        f2y = vy[i] + 0.5*dt*f1v
        f3v = (fdist + F_act[i] - c*(vy[i]+0.5*dt*f2v) - k*(y[i]+0.5*dt*f2y))/m
        f3y = vy[i] + 0.5*dt*f2v
        f4v = (fdist + F_act[i] - c*(vy[i]+dt*f3v) - k*(y[i]+dt*f3y))/m
        f4y = vy[i] + dt*f3v
        vy[i+1] = vy[i] + dt*(f1v + 2*f2v + 2*f3v + f4v)/6
        y[i+1]  = y[i]  + dt*(f1y + 2*f2y + 2*f3y + f4y)/6
        ay[i] = f1v
    ay[-1] = (F_dist[-1] + F_act[-1] - c*vy[-1] - k*y[-1])/m
    return y, vy, ay

# ---------------- 参数化 theta -----------------
nF = 6          # Fourier 阶数
def theta_from_coef(coef):
    """
    coef: (2, 2*nF+1)  每个作动器 [a0,a1,b1,...,aF,bF]
    返回 (2,N)
    """
    th = np.zeros((2, N))
    for i in range(2):
        c = coef[i]
        a0 = c[0]
        ab = c[1:].reshape(-1,2)   # [[a1,b1],...,[aF,bF]]
        tmp = a0
        for j in range(nF):
            tmp += ab[j,0]*np.cos((j+1)*2*np.pi*t/T) + ab[j,1]*np.sin((j+1)*2*np.pi*t/T)
        th[i] = tmp
    return th

# ---------------- 目标函数 -----------------
def objective(coef_flat):
    coef = coef_flat.reshape(2, -1)
    th = theta_from_coef(coef)
    # 简单限幅（可再加速度约束）
    y, vy, ay = rk4(0.01, 0.0, th)
    Ih = np.mean(ay**2)
    return Ih

# ---------------- 优化 -----------------
dim = 2*(2*nF+1)
x0 = np.zeros(dim)
es = cma.CMAEvolutionStrategy(x0, 0.5)
es.optimize(objective, iterations=300)
best_coef = es.result.xbest.reshape(2, -1)

# ---------------- 最优结果重算 -----------------
theta_opt = theta_from_coef(best_coef)
y_uc, vy_uc, ay_uc = rk4(0.01, 0.0, np.zeros((2,N)))   # 无控
y_c,  vy_c,  ay_c  = rk4(0.01, 0.0, theta_opt)         # 有控

Ih_uc = np.mean(ay_uc**2)
Ih_c  = np.mean(ay_c**2)

print('无控 Ih =', Ih_uc)
print('有控 Ih =', Ih_c, '降幅 {:.1f}%'.format(100*(Ih_uc-Ih_c)/Ih_uc))

# ---------------- 绘图 -----------------
plt.figure(figsize=(12,8))
plt.subplot(3,1,1)
plt.plot(t, y_uc, label='无控位移')
plt.plot(t, y_c, label='有控位移')
plt.ylabel('y (m)'); plt.legend(); plt.grid()
plt.subplot(3,1,2)
plt.plot(t, ay_uc, label='无控加速度')
plt.plot(t, ay_c, label='有控加速度')
plt.ylabel('a_y (m/s²)'); plt.legend(); plt.grid()
plt.subplot(3,1,3)
plt.plot(t, theta_opt[0], label='θ₁(t)')
plt.plot(t, theta_opt[1], label='θ₂(t)')
plt.xlabel('t (s)'); plt.ylabel('角位移 (rad)'); plt.legend(); plt.grid()
plt.tight_layout()
plt.show()