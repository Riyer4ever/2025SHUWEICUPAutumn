## 1 问题目标
在**附录 1-Scenario 2**给出的横向扰动力作用下，为两台对角安装的离心作动器设计**角位移曲线**  
$\theta_1(t),\;\theta_2(t)$，使车体横向振动指标

$$
I_h=\frac{1}{T}\int_0^T a_y^2(t)\,dt
$$

最小，并给出无控/有控对比曲线及指标数值。

---

## 2 系统模型

| 符号 | 含义 | 数值 |
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
m\ddot{y}+c\dot{y}+ky=F_{\text{dist}}(t)+F_{\text{act}}(t)
$$

**作动器合力**（两台，各含 4 块对称反向旋转，$\omega\equiv\omega_{\max}$

$$
F_{\text{act}}(t)=2m_e r\omega^2\bigl[\sin\theta_1(t)+\sin\theta_2(t)\bigr]
$$

---

## 3 求解步骤

1. **读入扰动力**  
   时间列 $t\in[0,9.99]$ s，步长$\Delta t=0.01$ s，共 1001 点。

2. **参数化角位移**  
   采用 6 阶 Fourier 级数（可调整）：
   $$
   \theta_i(t)=a_0^{(i)}+\sum_{k=1}^{6}\Bigl[a_k^{(i)}\cos\!\tfrac{2\pi k t}{T}+b_k^{(i)}\sin\!\tfrac{2\pi k t}{T}\Bigr],\quad i=1,2
   $$
   优化变量：$2\times(2\times 6+1)=26$ 维实向量。

3. **数值积分**  
   - 固定步长 RK4 解二阶微分方程  
   - 初始条件：$y(0)=0.01\,\text{m},\;\dot{y}(0)=0\,\text{m/s}$  
   - 加速度直接由方程给出：
  $$
   a_y(t)=\ddot{y}(t)=\frac{F_{\text{dist}}(t)+F_{\text{act}}(t)-c\dot{y}(t)-ky(t)}{m}
  $$

1. **目标函数**
   $$
   J=I_h=\frac{1}{N}\sum_{i=0}^{N-1}a_{y,i}^2
   $$
   后处理检查角速度$/$加速度超限则加 penalty。

1. **优化算法**  
   使用 CMA-ES 离线求解，收敛后获得最优 Fourier 系数 $→$ $\theta_1(t),\;\theta_2(t)$。

2. **对比输出**  
   - 无控：$F_{\text{act}}\equiv 0$
   - 有控：采用最优 $\theta(t)$ 
   - 生成位移、加速度曲线及 $I_h$ 数值与降幅百分比。

---

## 4 输出结果
- 车身横向位移曲线 $y(t)$（无控 vs 有控）
- 车身横向加速度曲线 $a_y(t)$（无控 vs 有控）
- 最优角位移指令 $\theta_1(t),\;\theta_2(t)$
- 横向振动指标 $I_h$ 对比（数值 + 降幅）