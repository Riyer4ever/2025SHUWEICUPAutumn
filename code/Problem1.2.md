# 第一题第二问求解思路（有控工况）

## 1. 问题目标
- **开启离心作动器**
- 按附录 1-Scenario 1 给定的四块偏心角度 $\theta(t)$生成控制力
- 输出：有控横向位移曲线、横向加速度曲线、角位移曲线及横向振动指标 $I_h$，并与无控对比

## 2. 物理模型

| 项目 | 公式 |
| --- | --- |
| 单块离心力大小 | $F_0 = m_{\text{block}} r \omega^2$ |
| 同组两块横向分量 | $2 m_{\text{block}} r \omega^2 \sin\theta$ |
| 两台对角共四块 | $F_{\text{act}}(t) = 4 m_{\text{block}} r \omega^2(t) \sin\theta(t)$ |

**总运动方程**：
$$
m\ddot{y} + c\dot{y} + ky = F_{\text{dist}}(t) + F_{\text{act}}(t)
$$

## 3. 求解步骤

1. **读入角度数据**
   - 时间列：$t \in [0, 10]$ s，步长 $\Delta t = 0.01$ s
   - 角度：从 `Scenario1_theta.txt` 读取 $\theta(t)$

2. **计算角速度与作动力**
   - 角速度：\( \omega(t) = \dfrac{d\theta}{dt} \)（用 `np.gradient`）
   - 作动力：\( F_{\text{act}}(t) = 4 \cdot m_{\text{block}} \cdot r \cdot \omega^2(t) \cdot \sin\theta(t) \)

3. **数值积分**
   - 用固定步长 RK4 解总运动方程
   - 初始条件同第一问：
     \[
     y(0) = 0.01\ \text{m}, \quad \dot{y}(0) = 0\ \text{m/s}
     \]

4. **计算加速度**
   - 由方程直接得到：
     \[
     \ddot{y}(t) = \frac{F_{\text{dist}}(t) + F_{\text{act}}(t) - c\dot{y}(t) - ky(t)}{m}
     \]

5. **振动指标**
   \[
   I_h = \frac{1}{T}\int_0^T \ddot{y}^2(t)\,dt \approx \frac{1}{N}\sum_{i=0}^{N-1} \ddot{y}_i^2
   \]

## 4. 输出
- 有控横向位移曲线 \( y(t) \)
- 有控横向加速度曲线 \( \ddot{y}(t) \)
- 角位移曲线 \( \theta(t) \)（与附录1一致）
- 有控横向振动指标 \( I_h \)（数值），与无控并排对比