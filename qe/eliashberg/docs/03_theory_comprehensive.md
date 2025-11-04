# Eliashberg超导理论：完整理论基础

本文档系统整理Eliashberg超导方程的完整理论框架，包括有限温度场论、Matsubara频率表示、解析延拓和第一性原理计算接口。

---

## 第一章：有限温度量子场论基础

### 1.1 配分函数与虚时间形式

**配分函数的普遍形式**：
$$
Z = \text{Tr}[e^{-\beta H}], \quad \beta = \frac{1}{k_B T}
$$

**为何使用迹**：
- 保证基组无关性
- 物理观测量不依赖于表示

**Wick旋转**：
$$
e^{-\beta H} = e^{-iHt}\bigg|_{t=i\beta\hbar} \quad \text{(实时间 → 虚时间)}
$$

### 1.2 Matsubara频率的起源

**边界条件**：

虚时间格林函数必须满足：
- **费米子**：$G(\tau + \beta) = -G(\tau)$ (反周期)
- **玻色子**：$D(\tau + \beta) = D(\tau)$ (周期)

**离散化频率**：

- **费米子**：$\omega_n = \pi k_B T(2n+1)$，$n \in \mathbb{Z}$
- **玻色子**：$\Omega_m = 2\pi k_B T \cdot m$，$m \in \mathbb{Z}$

**物理图像**：温度$T$像一个"尺子"，把虚时间切成$\sim k_B T$大小的片段。

### 1.3 费米子反周期性的严格推导

**关键恒等式**：
$$
e^{\beta H} \psi e^{-\beta H} = -\psi
$$

**证明**：源于费米子反对易关系$\{\psi, \psi^\dagger\} = 1$。

**推论**：
$$
\omega_{-n-1} = -\omega_n
$$

这是Matsubara频率折叠的核心。

---

## 第二章：Eliashberg方程理论

### 2.1 各向同性Eliashberg方程

**质量重整化方程**：
$$
Z(i\omega_n) = 1 + \frac{\pi T}{\omega_n}\sum_{m=-\infty}^{+\infty} \lambda(i\omega_n - i\omega_m) \frac{\omega_m Z_m}{\sqrt{\omega_m^2 + \Delta_m^2}}
$$

**能隙方程**：
$$
Z(i\omega_n)\Delta_n = \pi T \sum_{m=-\infty}^{+\infty} [\lambda(i\omega_n - i\omega_m) - \mu_c^*] \frac{\Delta_m}{\sqrt{\omega_m^2 + \Delta_m^2}}
$$

**物理意义**：
- $\Delta_n$：超导能隙（Cooper对结合能）
- $Z_n$：准粒子重整化因子（有效质量增强）
- $\lambda$：电子-声子耦合核
- $\mu_c^*$：库仑赝势（电子-电子排斥修正）

### 2.2 电子-声子耦合核

**定义**：
$$
\lambda(i\omega_n - i\omega_m) = 2\int_0^{\infty} d\Omega \frac{\Omega \alpha^2F(\Omega)}{\Omega^2 + (\omega_n - \omega_m)^2}
$$

**关键性质**：
1. **偶对称**：$\lambda(i\omega) = \lambda(-i\omega)$
2. **无奇点**：在虚轴上解析
3. **表格化**：离散化为$\Lambda_m$数组

### 2.3 Matsubara频率折叠

**原始求和**（对所有$m \in \mathbb{Z}$）→ **折叠求和**（只对$m \geq 0$）

**Z方程折叠形式**：
$$
Z_n = 1 + \frac{\pi T}{\omega_n}\sum_{m=0}^{N_\omega-1} [\Lambda_{n-m} - \Lambda_{n+m+1}] \frac{\omega_m Z_m}{R_m}
$$

**能隙方程折叠形式**：
$$
\Delta_n = \frac{\pi T}{Z_n} \left[\sum_{m=0}^{N_\omega-1} [\Lambda_{n-m} + \Lambda_{n+m+1}] \frac{\Delta_m Z_m}{R_m} - \Delta\mu\right]
$$

其中$R_m = |Z_m \omega_m| = \sqrt{(\omega_m^2 + \Delta_m^2)Z_m^2}$。

**关键理解**：
- 索引$l(n+m+1)$对应频率差$\omega_n + \omega_m$
- 利用$\omega_{-m-1} = -\omega_m$实现折叠
- **不是**"加减$2\pi k_B T$不变"！

---

## 第三章：实轴奇点与Padé延拓

### 3.1 为何在虚轴求解？

**实轴问题**：

实频格林函数：
$$
G^R(\omega) = \frac{1}{\omega - \epsilon + i0^+}
$$

- 有极点/δ函数奇点
- Cauchy主值积分数值困难
- 自洽迭代不稳定

**Matsubara优势**：

虚频格林函数：
$$
G(i\omega_n) = \frac{1}{i\omega_n - \epsilon}
$$

- 无奇点，光滑函数
- 离散求和，数值稳定
- 自洽迭代快速收敛

### 3.2 Padé解析延拓

**目的**：从Matsubara轴($i\omega_n$)延拓到实轴($\omega + i0^+$)

**Vidberg-Serene算法**：

1. 构建$g$矩阵：
   $$
   g_{i,j} = \frac{g_{i-1,i-1} - g_{i-1,j}}{(z_j - z_{i-1}) g_{i-1,j}}
   $$

2. 连分式展开：
   $$
   P(z) = g_{0,0} + \frac{(z - z_0) g_{1,1}}{1 + \frac{(z - z_1) g_{2,2}}{1 + \cdots}}
   $$

**精度**：对光滑函数，相对误差 < 10⁻⁶

---

## 第四章：McMillan-Allen-Dynes临界温度

### 4.1 电声耦合参数

**耦合常数**：
$$
\lambda = 2\int_0^{\omega_{\max}} \frac{\alpha^2 F(\omega)}{\omega} d\omega
$$

**特征频率**：
$$
\omega_{\log} = \exp\left[\frac{2}{\lambda}\int_0^{\omega_{\max}} \frac{\alpha^2 F(\omega)\ln\omega}{\omega} d\omega\right]
$$

$$
\omega_{\mathrm{rms}} = \sqrt{\frac{2}{\lambda}\int_0^{\omega_{\max}} \alpha^2 F(\omega) \omega d\omega}
$$

### 4.2 Allen-Dynes公式

**临界温度**：
$$
T_c = \frac{\omega_{\log}}{1.2 k_B} \exp\left[\frac{-1.04(1+\lambda)}{\lambda - \mu^* - 0.62\lambda\mu^*}\right] f_1 f_2
$$

**修正因子**：
$$
f_1 = \left[1 + \left(\frac{\lambda}{2.46(1+3.8\mu^*)}\right)^{3/2}\right]^{1/3}
$$

$$
f_2 = 1 + \frac{\lambda^2(\omega_{\mathrm{rms}}/\omega_{\log} - 1)}{\lambda^2 + [1.82(1+6.3\mu^*)(\omega_{\mathrm{rms}}/\omega_{\log})]^2}
$$

**适用范围**：$\lambda < 1.5$，$\mu^* \sim 0.1-0.15$

---

## 第五章：第一性原理计算接口

### 5.1 Eliashberg谱函数α²F(ω)

**定义**（见第六章详细推导）：
$$
\alpha^2F(\omega) = \frac{1}{N(0)} \sum_{\mathbf{k},\mathbf{q},\nu} |g_{mn\nu}(\mathbf{k},\mathbf{q})|^2 \delta(\epsilon_{m\mathbf{k}}) \delta(\epsilon_{n\mathbf{k}+\mathbf{q}}) \delta(\omega - \omega_{q\nu})
$$

**物理意义**：
- 电声耦合强度的频率分布
- 加权的声子态密度
- 归一化：$\int \alpha^2F(\omega) d\omega = \lambda/2$

### 5.2 从EPW提取数据

**三个层次**：

1. **各向同性参数**（已实现）：
   - 读取`ALPHA2F.OUT`
   - 计算$\lambda$、$\omega_{\log}$、$\omega_{\mathrm{rms}}$
   - 输入`eliashberg_solver.py`

2. **模式分解**（中等难度）：
   - 分解$\alpha^2 F(\omega) = \sum_\nu \alpha^2 F_\nu(\omega)$
   - 识别主导声子模式

3. **各向异性数据**（高难度）：
   - 提取$g_{mn\nu}(\mathbf{k},\mathbf{q})$矩阵元
   - 构建$\mathbf{k}$依赖的有效模型

**工具**：`extract_from_epw.py`实现层次1和层次2

---

## 第六章：α²F(ω)的物理意义与计算（新增）

### 6.1 定义与物理意义

**完整定义**：
$$
\alpha^2F(\omega) = \frac{1}{N(0)} \sum_{\mathbf{k},\mathbf{q}} \sum_{m,n} \sum_{\nu} |g_{mn\nu}(\mathbf{k},\mathbf{q})|^2 \delta(\epsilon_{m\mathbf{k}}) \delta(\epsilon_{n\mathbf{k}+\mathbf{q}}) \delta(\omega - \omega_{q\nu})
$$

其中：
- $g_{mn\nu}(\mathbf{k},\mathbf{q})$：电声耦合矩阵元
- $\epsilon_{m\mathbf{k}}$：电子能量（相对费米面）
- $\omega_{q\nu}$：声子频率（支标$\nu$）
- $N(0)$：费米面态密度

**物理解释**：

1. **电声耦合强度**：$|g_{mn\nu}|^2$描述电子通过吸收/发射声子的散射强度
2. **费米面约束**：两个δ函数确保只有费米面附近的电子参与
3. **声子能量依赖**：第三个δ函数选择频率为$\omega$的声子模式
4. **归一化**：$N(0)$确保结果与态密度无关

**与λ的关系**：
$$
\lambda = 2\int_0^{\infty} \frac{\alpha^2F(\omega)}{\omega} d\omega
$$

### 6.2 从Quantum ESPRESSO计算α²F(ω)

**完整工作流程**：

#### 步骤1：DFT自洽计算 (QE: pw.x)

```bash
pw.x < scf.in > scf.out
```

**关键设置**：
- 收敛精度：`conv_thr = 1.0d-12`
- k点网格：密集网格（如12×12×12）
- 交换关联泛函：通常使用PBE

#### 步骤2：声子计算 (QE: ph.x)

```bash
ph.x < ph.in > ph.out
```

**关键输入**：
```fortran
&inputph
  tr2_ph = 1.0d-14,      ! 声子收敛阈值
  ldisp = .true.,        ! 计算声子色散
  nq1 = 6, nq2 = 6, nq3 = 6,  ! q点网格
  fildyn = 'dyn',        ! 动力学矩阵文件前缀
/
```

**输出**：
- 动力学矩阵：`dyn0`, `dyn1`, ..., `dynN`
- 声子频率和本征矢

#### 步骤3：电声耦合插值 (EPW)

**EPW输入文件**（`epw.in`）：

```fortran
&inputepw
  ! 关键参数
  elph = .true.              ! 计算电声耦合
  epbwrite = .false.         ! 不写原始矩阵元（节省空间）
  epbread = .false.

  ! Wannier插值
  wannierize = .true.        ! 使用Wannier90
  nbndsub = 8                ! Wannier能带数

  ! k和q网格
  nkf1 = 40, nkf2 = 40, nkf3 = 40  ! 精细k网格
  nqf1 = 20, nqf2 = 20, nqf3 = 20  ! 精细q网格

  ! Eliashberg谱函数
  eliashberg = .true.        ! 计算α²F(ω)
  la2F = .true.              ! 输出α²F文件

  ! 频率范围
  nqstep = 500               ! 频率点数
  degaussw = 0.05            ! 展宽参数 (eV)
/
```

**关键过程**：

1. **粗网格计算**：
   - DFT和DFPT在粗k/q网格上计算
   - 得到$g_{mn\nu}(\mathbf{k},\mathbf{q})$在粗网格上的值

2. **Wannier插值**：
   - 将Bloch态投影到Wannier基组
   - 在Wannier表示下插值到任意k点

3. **精细网格重建**：
   - 在密集k/q网格上重建$g_{mn\nu}$
   - 计算α²F(ω)的频率依赖

**数学细节**：
$$
g_{mn\nu}(\mathbf{k},\mathbf{q}) = \langle \psi_{m\mathbf{k}} | \frac{\partial V}{\partial u_{q\nu}} | \psi_{n\mathbf{k}+\mathbf{q}} \rangle
$$

其中$u_{q\nu}$是声子位移模式。

#### 步骤4：α²F(ω)的数值计算

**离散化公式**：
$$
\alpha^2F(\omega) = \frac{1}{N(0)} \frac{1}{N_k N_q} \sum_{\mathbf{k},\mathbf{q},m,n,\nu} |g_{mn\nu}|^2 \delta(\epsilon_{m\mathbf{k}}) \delta(\epsilon_{n\mathbf{k}+\mathbf{q}}) \delta(\omega - \omega_{q\nu})
$$

**δ函数处理**（高斯展宽）：
$$
\delta(x) \to \frac{1}{\sigma\sqrt{2\pi}} e^{-x^2/(2\sigma^2)}
$$

**展宽参数选择**：
- 电子：`degaussq` ~ 0.05 eV（匹配费米面展宽）
- 声子：`degaussw` ~ 0.5-1.0 cm⁻¹

**输出文件**：`ALPHA2F.OUT`

```
# omega (Hartree)    alpha2F(omega)
0.000000000000      0.000000000000
0.000045563353      0.012345678901
...
```

### 6.3 数值例子与解释

**典型材料的α²F(ω)特征**：

1. **简单金属（Al）**：
   - 单峰结构
   - 峰位 ~ 30 meV（对应声学声子频率）
   - λ ~ 0.4

2. **强耦合超导体（Pb）**：
   - 多峰结构
   - 主峰 ~ 5-10 meV
   - λ ~ 1.5-2.0

3. **MgB₂**：
   - 双峰结构
   - 低频峰：E₂g声子 (~ 70 meV)
   - 高频峰：B-B拉伸模 (~ 90 meV)
   - λ ~ 0.6-0.9

**从α²F(ω)可以读出的信息**：

1. **主导声子模式**：
   - 峰位对应关键声子频率
   - 峰强正比于电声耦合强度

2. **λ的贡献分解**：
   $$
   \lambda = 2 \sum_i \int_{\omega_i} \frac{\alpha^2F(\omega)}{\omega} d\omega
   $$
   其中$i$标记不同的峰

3. **Tc的估算**：
   - $\omega_{\log}$更接近主峰位置
   - 低频模式贡献更大（因为$1/\omega$权重）

### 6.4 常见问题与解决方案

**问题1：α²F(ω)有负值？**

**原因**：插值误差或Wannier投影不完整

**解决**：
- 增加Wannier函数数量`nbndsub`
- 检查Wannier spread是否合理
- 增加k/q网格密度

**问题2：λ与实验不符？**

**可能原因**：
- k/q网格不够密
- 展宽参数`degaussq`选择不当
- 交换关联泛函不准确

**调试**：
- 做收敛性测试：改变k/q网格
- 比较不同泛函（PBE vs LDA）
- 对比声子谱与中子散射实验

**问题3：计算耗时太长？**

**优化策略**：
- 使用对称性降低k/q点数
- 并行化：MPI + OpenMP
- 使用EPW的restart功能
- 减小精细网格（先用粗网格测试）

---

## 附录A：单位换算

**Hartree原子单位**：
- 能量：1 Hartree = 27.211 eV = 219474.63 cm⁻¹
- 温度：1 Hartree = $3.1577 \times 10^5$ K
- 玻尔兹曼常数：$k_B = 3.166815343 \times 10^{-6}$ Hartree/K

**常用换算**：
```python
# 频率
omega_hartree = omega_cm / 219474.63  # cm⁻¹ → Hartree
omega_eV = omega_hartree * 27.211      # Hartree → eV

# 温度
temp_K = temp_hartree / KBOLTZ         # Hartree → K
```

---

## 附录B：工具脚本使用

**已实现的工具**（位于`tools/`目录）：

1. `extract_from_epw.py`：从EPW提取数据
2. `example_pade_continuation.py`：Padé延拓演示
3. `example_real_vs_imaginary_axis.py`：实轴vs虚轴对比
4. `visualize_matsubara.py`：Matsubara频率可视化
5. `extract_dominant_modes.py`：声子模式分解
6. `verify_matsubara_folding.py`：频率折叠验证

**使用示例**：
```bash
# 从EPW提取参数
cd /path/to/epw/calculation
python /path/to/tools/extract_from_epw.py

# 运行Eliashberg求解器
python eliashberg_solver.py --input eliashberg_input.dat

# Padé延拓演示
python tools/example_pade_continuation.py
```

---

## 总结

本理论框架建立了从第一性原理到超导Tc预测的完整链条：

```
Quantum ESPRESSO (DFT + DFPT)
        ↓
EPW (Wannier插值)
        ↓
α²F(ω) (Eliashberg谱函数)
        ↓
λ, ω_log (电声耦合参数)
        ↓
Eliashberg方程求解 (Matsubara形式)
        ↓
Padé延拓 (虚轴 → 实轴)
        ↓
Tc, Δ(T), 谱函数 (实验可观测量)
```

每个环节的理论基础、数值实现和验证方法都已详细阐述，为进一步研究奠定了坚实基础。

---

## 第七章：能隙函数Δ(ω)的复数性质与解析延拓

### 7.1 核心问题：Δ(ω)是实函数还是复函数？

**简短回答**：
- **Matsubara轴上**（虚频）：Δ(iωₙ)是**实数**（对于各向同性s波超导体）
- **实轴上**（实频）：Δ(ω+i0⁺)是**复数**

这看似矛盾的性质源于解析函数在不同轴上的不同表现形式。

---

### 7.2 Matsubara轴上：Δ(iωₙ)是实数

#### 对称性证明

对于各向同性s波超导体，能隙方程：
$$
\Delta(i\omega_n) = \frac{\pi T}{Z_n} \sum_m [\Lambda_{n-m} + \Lambda_{n+m+1}] \frac{\Delta_m Z_m}{R_m} - \Delta\mu
$$

**关键性质**：

1. **核函数的实性**：
   $$
   \Lambda_m = 2\int_0^{\infty} \frac{\Omega \alpha^2F(\Omega)}{\Omega^2 + (2\pi k_B T \cdot m)^2} d\Omega \in \mathbb{R}
   $$

   因为$\alpha^2F(\Omega)$是实的声子谱函数。

2. **Matsubara频率的纯虚性**：
   $$
   i\omega_n = i\pi k_B T(2n+1) \quad \text{(纯虚数)}
   $$

3. **对称性条件**：
   $$
   \Lambda_m = \Lambda_{-m}, \quad \omega_{-n-1} = -\omega_n
   $$

**结论**：若初始猜测$\Delta_m$是实数，则方程右端是实数的线性组合，因此解$\Delta(i\omega_n) \in \mathbb{R}$。

#### 数学验证

考虑复共轭：
$$
\Delta^*(i\omega_n) = \left[\frac{\pi T}{Z_n} \sum_m [\Lambda_{n-m} + \Lambda_{n+m+1}] \frac{\Delta_m Z_m}{R_m}\right]^*
$$

由于右端所有项都是实数：
$$
\Delta^*(i\omega_n) = \Delta(i\omega_n)
$$

因此$\Delta(i\omega_n)$是实数。

**重要说明**：
- 这个结论对**各向同性**超导体成立
- 对于d波、p波等**各向异性**超导体，Δ可能包含相位因子
- 对于**多能带**超导体，不同能带的Δ可能有相位差

---

### 7.3 实轴上：Δ(ω)是复函数

#### 解析延拓的必然性

从Matsubara轴解析延拓到实轴：
$$
\Delta(i\omega_n) \xrightarrow[\text{Padé延拓}]{\text{解析}} \Delta(\omega + i0^+)
$$

**实轴上的复数形式**：
$$
\Delta(\omega + i0^+) = \Delta'(\omega) + i\Delta''(\omega)
$$

其中：
- $\Delta'(\omega)$：实部，描述**能隙大小**
- $\Delta''(\omega)$：虚部，描述**准粒子寿命**（展宽）

#### 物理意义

**实部** $\Delta'(\omega)$：
- 对应超导能隙的大小
- 在$\omega \sim 0$附近最大
- 决定Cooper对的结合能

**虚部** $\Delta''(\omega)$：
- 来源于准粒子的有限寿命
- 描述准粒子的散射速率：$\tau^{-1} \propto \Delta''$
- 在$\omega > 2\Delta$时变大（对破坏过程）

#### Green函数的极点结构

超导态的单粒子Green函数：
$$
G(\mathbf{k}, \omega) = \frac{u_k^2}{\omega - E_k + i\eta} + \frac{v_k^2}{\omega + E_k + i\eta}
$$

其中准粒子能量：
$$
E_k = \sqrt{\xi_k^2 + |\Delta(\omega)|^2}
$$

**关键**：$\Delta(\omega)$的虚部导致准粒子极点有有限展宽，对应有限寿命。

---

### 7.4 为什么在虚轴上求解？

#### 原因1：数值稳定性

**实轴问题**：

能隙自洽方程在实轴上：
$$
\Delta(\omega) = \int_{-\infty}^{\infty} K(\omega, \omega') \frac{\Delta(\omega')}{\sqrt{(\omega')^2 + |\Delta(\omega')|^2}} d\omega'
$$

其中核函数$K(\omega, \omega')$包含：
$$
K(\omega, \omega') \sim \int d\epsilon \frac{\alpha^2F(\Omega)}{\omega - \epsilon + i0^+}
$$

**问题**：
1. $i0^+$让积分变成Cauchy主值积分，数值实现困难
2. 被积函数在费米面附近剧烈振荡
3. 自洽迭代容易发散

**Matsubara优势**：

$$
\Delta(i\omega_n) = \frac{\pi T}{Z_n} \sum_m \lambda(i\omega_n - i\omega_m) \frac{\Delta_m}{R_m}
$$

**优势**：
1. 离散求和代替连续积分
2. 虚轴上无奇点，函数光滑
3. 迭代稳定，通常10-50次收敛

#### 原因2：温度场论的自然框架

有限温度Green函数天然定义在Matsubara频率上：
$$
G(\tau) = -\langle T_\tau \psi(\tau) \psi^\dagger(0) \rangle
$$

Fourier变换：
$$
G(i\omega_n) = \int_0^\beta d\tau e^{i\omega_n \tau} G(\tau)
$$

**配分函数**：
$$
Z = \text{Tr}[e^{-\beta H}], \quad \beta = 1/(k_B T)
$$

所有有限温度平衡态物理量都在虚时间/虚频率轴上计算。

#### 原因3：解析性质的保证

Matsubara格林函数在整个上半复平面解析：
$$
G(z) \text{ 在 } \text{Im}(z) > 0 \text{ 解析}
$$

这保证了Padé延拓的数学合法性。

---

### 7.5 解析延拓后的物理可观测量

#### 谱函数

**定义**：
$$
A(\mathbf{k}, \omega) = -\frac{1}{\pi} \text{Im}[G(\mathbf{k}, \omega + i0^+)]
$$

**物理意义**：
- 描述在动量$\mathbf{k}$、能量$\omega$处找到一个准粒子的概率
- ARPES（角分辨光电子谱）直接测量$A(\mathbf{k}, \omega)$

**能隙的虚部贡献**：
$$
A(\mathbf{k}, \omega) \propto \frac{\Delta''(\omega)}{[\omega^2 - \xi_k^2 - \Delta'^2(\omega)]^2 + [\Delta''(\omega)]^2}
$$

当$\Delta''$很小时，谱函数在$\omega = \sqrt{\xi_k^2 + \Delta'^2}$处有尖峰（准粒子峰）。

#### 态密度

**定义**：
$$
N(\omega) = \int \frac{d^3k}{(2\pi)^3} A(\mathbf{k}, \omega)
$$

**各向同性情况**：
$$
N(\omega) = N(0) \text{Re}\left[\frac{\omega}{\sqrt{\omega^2 - \Delta^2(\omega)}}\right]
$$

**特征**：
- $\omega < |\Delta|$：态密度被压制（超导能隙）
- $\omega = |\Delta|$：相干峰（BCS奇点被$\Delta''$展宽）
- $\omega \gg |\Delta|$：恢复正常态密度

#### 光学电导

**Mattis-Bardeen公式**（零温）：
$$
\sigma(\omega) = \frac{2\sigma_n}{\omega} \int_\Delta^\infty \frac{[f(E) - f(E+\omega)] (E^2 + \Delta^2 + E\omega)}{\sqrt{E^2 - \Delta^2} \sqrt{(E+\omega)^2 - \Delta^2}} dE
$$

**能隙$\Delta$的影响**：
- $\omega < 2\Delta$：$\sigma(\omega) = 0$（无吸收，完美反射）
- $\omega = 2\Delta$：吸收边（对破坏阈值）
- $\omega > 2\Delta$：单粒子吸收

---

### 7.6 代码实现中的体现

#### Matsubara轴上的计算

在`eliashberg_solver.py`中，能隙数组定义为实数：
```python
d0 = np.full(nwf_alloc + 1, 1.0e-4, dtype=float)  # 实数数组
```

自洽迭代中：
```python
for iteration in range(MAX_IT):
    # 所有运算都保持实数
    delta_new = calculate_delta(delta_old, z_old, lambda_kernel)
    # delta_new 仍然是实数
```

#### Padé延拓到实轴

在`pade_approximation`函数中：
```python
def pade_approximation(zin: np.ndarray, uin: np.ndarray, zout: np.ndarray):
    """
    zin: 虚频点 [iω₀, iω₁, ...]  (纯虚数)
    uin: 能隙值 [Δ₀, Δ₁, ...]    (实数)
    zout: 实频点 [ω+i0⁺]          (复数)

    返回: uout = [Δ(ω+i0⁺)]       (复数！)
    """
    g = np.zeros((nin, nin), dtype=complex)  # 注意：complex类型
    g[0, :] = uin  # 实数初始化

    # 递归计算g矩阵（会产生复数）
    for i in range(1, nin):
        for j in range(i, nin):
            g[i, j] = (g[i-1, i-1] - g[i-1, j]) / ((zin[j] - zin[i-1]) * g[i-1, j])

    # 延拓到实轴（结果是复数）
    uout = np.zeros(nout, dtype=complex)
    # ... (递归计算)

    return uout  # 复数数组
```

**关键转变**：
```
输入（Matsubara）：uin是实数数组
输出（实轴）：uout是复数数组
```

---

### 7.7 不同超导体的对比

#### 各向同性s波（如Pb, Al）

| 轴 | Δ的性质 | 原因 |
|----|---------|------|
| Matsubara | 实数 | 对称性 + 实核函数 |
| 实轴 | 复数 | 解析延拓 + 准粒子寿命 |

**特征**：
- $\Delta(i\omega_n)$随$n$单调减小
- $\Delta'(\omega)$在$\omega=0$最大
- $\Delta''(\omega)$在$\omega \sim 2\Delta$附近最大

#### 各向异性d波（如铜氧化物）

| 轴 | Δ的性质 | 原因 |
|----|---------|------|
| Matsubara | 实数（角平均） | 费米面平均后保持实数 |
| 实轴 | 复数 | 同上 |

**特征**：
- 节点方向：$\Delta(\mathbf{k}) = 0$
- 最大能隙方向：$\Delta(\mathbf{k}) = \Delta_0$
- 低能激发：线性而非指数

#### 多能带超导体（如MgB₂, 铁基）

| 能带 | Δ的性质 | 相对相位 |
|------|---------|----------|
| σ能带 | Δ_σ(iωₙ) 实数 | 相位φ_σ |
| π能带 | Δ_π(iωₙ) 实数 | 相位φ_π |

**可能情况**：
- 同相：φ_σ = φ_π（s⁺⁺波）
- 反相：φ_σ = φ_π + π（s±波）

---

### 7.8 实验验证

#### ARPES测量

**测量量**：谱函数$A(\mathbf{k}, \omega)$

**超导能隙的体现**：
```
正常态：连续的费米面
超导态：费米面打开能隙 Δ(k)
```

**能量分辨**：
- 高分辨ARPES：~1 meV
- 典型能隙：Δ ~ 10-100 meV
- 可以分辨$\Delta'(\omega)$的频率依赖

**观测**：
- 能隙大小：从准粒子峰位置
- 能隙对称性：从角度依赖
- 准粒子寿命：从峰宽（正比于$\Delta''$）

#### STM隧道谱

**测量量**：
$$
\frac{dI}{dV} \propto N(\omega) = \text{Re}\left[\frac{\omega}{\sqrt{\omega^2 - \Delta^2(\omega)}}\right]
$$

**特征**：
- 零偏压：$dI/dV = 0$（能隙内无态）
- $V = \Delta/e$：相干峰（BCS奇点）
- 峰宽：反映$\Delta''$

**优势**：
- 空间分辨：纳米级
- 能量分辨：亚meV级
- 可探测局域能隙不均匀性

#### 光学测量

**测量量**：反射率$R(\omega)$，推导出$\sigma(\omega)$

**特征**：
- $\omega < 2\Delta$：$R \approx 1$（完美反射）
- $\omega = 2\Delta$：反射率下降（吸收边）
- $\omega \gg 2\Delta$：恢复正常态反射

---

### 7.9 常见误解澄清

#### 误解1："超导能隙总是实数"

**错误原因**：混淆了Matsubara表示与实频表示

**正确理解**：
- Matsubara轴：$\Delta(i\omega_n) \in \mathbb{R}$（各向同性s波）
- 实轴：$\Delta(\omega+i0^+) \in \mathbb{C}$（总是复数）

#### 误解2："虚部是不物理的"

**错误原因**：误以为"虚部"只是数学工具

**正确理解**：
- $\Delta''(\omega)$描述准粒子有限寿命
- 直接影响谱函数线宽
- 可被ARPES测量

#### 误解3："解析延拓会引入虚部"

**更精确的说法**：
- 解析延拓**揭示**了原本就存在的虚部
- 从Matsubara（实数）→ 实轴（复数）是函数的不同"观察角度"
- 类似于从极坐标$r$到直角坐标$(x, y)$

#### 误解4："能隙是常数"

**BCS近似**：$\Delta(\omega) = \Delta_0$（常数）

**Eliashberg理论修正**：$\Delta(\omega)$有频率依赖
- 来源于电子-声子相互作用的能量依赖
- 在强耦合超导体中显著（如Pb）
- 弱耦合时接近BCS（如Al）

---

### 7.10 数学补充：解析函数的性质

#### Schwarz反射原理

对于解析函数$f(z)$，若在实轴上$f(x) \in \mathbb{R}$，则：
$$
f(z^*) = [f(z)]^*
$$

**应用**：
- 若$\Delta(i\omega_n) \in \mathbb{R}$（虚轴上）
- 则$\Delta(-i\omega_n) = \Delta(i\omega_n)$
- 但延拓到$z = \omega + i\eta$时，$\Delta(z) \in \mathbb{C}$

#### Kramers-Kronig关系

实部和虚部通过Hilbert变换关联：
$$
\Delta'(\omega) = \frac{1}{\pi} \mathcal{P} \int_{-\infty}^{\infty} \frac{\Delta''(\omega')}{\omega' - \omega} d\omega'
$$

$$
\Delta''(\omega) = -\frac{1}{\pi} \mathcal{P} \int_{-\infty}^{\infty} \frac{\Delta'(\omega')}{\omega' - \omega} d\omega'
$$

**物理意义**：
- 因果性的数学表达
- 知道虚部可推出实部（反之亦然）
- 实验上可用于一致性检验

---

### 7.11 总结

| 问题 | 答案 |
|------|------|
| Matsubara轴上Δ是实/复？ | **实数**（各向同性s波） |
| 实轴上Δ是实/复？ | **复数** |
| 为什么在虚轴求解？ | 数值稳定、无奇点 |
| 解析延拓后为何变复数？ | 准粒子有限寿命 → 虚部 |
| 虚部有物理意义吗？ | 有！描述展宽/散射率 |
| 实验能测到虚部吗？ | 能！通过ARPES峰宽等 |

**核心认知**：
> Δ(ω)的"实"或"复"取决于在复平面的哪条轴上观察。<br>
> Matsubara轴（虚轴）：实数，适合数值计算。<br>
> 实频轴：复数，对应物理观测量。<br>
> 解析延拓是连接计算与实验的桥梁。

**记忆口诀**：
```
虚轴求解得实数，稳定收敛易自洽；
实轴延拓变复数，虚部寿命可测量。
```

---

## 附录C：符号约定总结

| 符号 | 含义 | 取值 |
|------|------|------|
| $\Delta(i\omega_n)$ | Matsubara能隙 | 实数（s波） |
| $\Delta(\omega+i0^+)$ | 实频能隙 | 复数 |
| $\Delta'(\omega)$ | 能隙实部 | 能隙大小 |
| $\Delta''(\omega)$ | 能隙虚部 | 准粒子展宽 |
| $Z(i\omega_n)$ | 质量重整化 | 实数 |
| $\Lambda_m$ | 耦合核 | 实数 |
| $\alpha^2F(\Omega)$ | 声子谱函数 | 实数 |
| $A(\mathbf{k},\omega)$ | 谱函数 | 实数（≥0） |

**注意**：
- 所有Matsubara量在各向同性情况下都是实数
- 延拓到实轴后一般变为复数
- 谱函数等物理观测量是实数（但通过复格林函数的虚部得到）
