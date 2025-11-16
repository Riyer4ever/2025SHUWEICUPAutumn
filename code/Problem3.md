## 1 问题目标（现实场景）
- **扰动力未知**（无法测量）  
- **仅可实时采集**车体横向加速度信号　$a_y(t)$  
- 设计**在线控制策略**：仅用加速度反馈生成离心作动器角位移指令　$\theta(t)$  
- 目标：最小化横向振动指标　$I_h=\dfrac{1}{T}\int_0^T a_y^2(t)\,dt$

---

## 2 系统与测量模型

| 符号 | 含义 | 值 |
|---|---|---|
| $m$ | 车体质量 | 2000 kg |
| $c$ | 等效阻尼 | 3600 N·s/m |
| $k$ | 等效刚度 | 7.2252×10⁶ N/m |
| $m_e$ | 单偏心块质量 | 100 kg |
| $r$ | 旋转半径 | 0.2 m |
| $\omega_{\max}$ | 最大角速度 | 100 rad/s |
| $\alpha_{\max}$ | 最大角加速度 | 5000 rad/s² |

**运动方程**  
$$
m\ddot{y}+c\dot{y}+ky=F_{\text{dist}}(t)+2m_e r\omega^2\bigl[\sin\theta_1(t)+\sin\theta_2(t)\bigr]
$$
**可测信号**：仅　$a_y(t)=\ddot{y}(t)$　（采样频率　$f_s=100\ \text{Hz}$）

---

## 3 设计思路：加速度反馈 + 在线估计

### 3.1 两步策略
1. **在线估计**  
   用加速度计数据实时估计　$\dot{y}$　与　$y$　（积分器 + 高通滤波消除漂移）。
2. **反馈控制**  
   采用**加速度负反馈 + 阻尼注入**思想，令作动器产生与加速度反向的力：
   $$
   \sin\theta_1(t)+\sin\theta_2(t)=-K_a\cdot a_y(t)
   $$
   其中　$K_a$　为反馈增益，通过离线扫频或在线优化确定，并满足饱和限幅：
   $$
   |\sin\theta_i|\le 1\quad\Rightarrow\quad|K_a\cdot a_y|\le 2
   $$

### 3.2 抗饱和与平滑
- 对　$a_y$　做一阶低通滤波（截止　$f_c\approx 15\ \text{Hz}$）抑制高频噪声
- 滤波后信号再经积分器链（HP→LP）得到　$\hat{y},\ \hat{\dot{y}}$　用于监控与增益调度
- 角速度/加速度约束：在　$\theta$　指令后加**速率限制器**（≤100 rad/s，≤5000 rad/s²）

### 3.3 在线实现流程（每 0.01 s）
1. 读取　$a_y(k)$
2. 滤波：　$a_f(k)=\alpha\, a_f(k-1)+(1-\alpha)\,a_y(k)$
3. 计算总正弦指令：　$S(k)=-\text{sat}\bigl(K_a\cdot a_f(k),\ -2,\ 2\bigr)$
4. 分配：　$\sin\theta_1(k)=\sin\theta_2(k)=S(k)/2$
5. 反三角：　$\theta_i(k)=\arcsin\bigl(\text{sat}(S(k)/2,\ -1,\ 1)\bigr)$
6. 速率/加速度限幅输出

---

## 4 离线调参与验证
- 用 Scenario 2 数据离线扫频　$K_a\in[0,\ 0.002]$　寻找最小　$I_h$
- 最优　$K_a^*$　代入在线算法，运行 RK4 仿真
- 输出：无控 vs 仅加速度反馈的　$y(t),\ a_y(t),\ \theta(t),\ I_h$

---

## 5 输出结果
- 车身横向位移曲线　$y(t)$　（无控 vs 加速度反馈）
- 车身横向加速度曲线　$a_y(t)$
- 在线生成的角位移指令　$\theta_1(t),\ \theta_2(t)$
- 横向振动指标　$I_h$　对比与降幅百分比

---

## 6 代码框架（Python 伪代码）
```python
# 离线最优 Ka 已标定
Ka_opt = 0.00087

# 初始化
theta = np.zeros((2, N))
afilt = 0
for k in range(1, N):
    ay_k = ay_meas[k]          # 来自传感器或仿真
    afilt = 0.8*afilt + 0.2*ay_k
    S = -np.clip(Ka_opt*afilt, -2.0, 2.0)
    sin_th = np.clip(S/2.0, -1.0, 1.0)
    theta[:, k] = np.arcsin(sin_th)

# 限速率/加速度（简单后处理）
theta = rate_limit(theta, dt, wmax, amax)

# RK4 闭环仿真
y_cl, vy_cl, ay_cl = rk4(0.01, 0.0, theta)
Ih_cl = np.mean(ay_cl**2)