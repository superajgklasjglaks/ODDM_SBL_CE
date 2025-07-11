# Massive MIMO 方向角估计算法

本文件夹包含了基于ODDM信道估计算法改编的Massive MIMO方向角估计算法实现。

## 算法概述

本实现模仿了ODDM系统中时延和多普勒频移的特性，将其应用于Massive MIMO系统的发射方向角(AoD)和接收方向角(AoA)估计。主要包含以下算法：

1. **1D SBL** - 一维稀疏贝叶斯学习算法
2. **2D SBL** - 二维稀疏贝叶斯学习算法
3. **分层1D SBL** - 分层一维稀疏贝叶斯学习算法
4. **分层2D SBL** - 分层二维稀疏贝叶斯学习算法
5. **OMP** - 正交匹配追踪算法
6. **Traditional Beamforming** - 传统波束成形算法

## 系统参数

- **发射天线数**: 32
- **接收天线数**: 32
- **天线间距**: 0.5λ (半波长)
- **载波频率**: 28 GHz
- **角度范围**: ±60°

## 文件说明

### 主要文件
- `massive_mimo_angle_estimation.m` - 主仿真脚本，比较所有算法性能
- `test_mimo_angle_estimation.m` - 测试脚本，验证算法正确性

### 核心算法函数
- `MIMO_CE_1D_SBL.m` - 1D SBL角度估计算法
- `MIMO_CE_2D_SBL.m` - 2D SBL角度估计算法
- `MIMO_hierarchical_SBL_refined_1D.m` - 分层1D SBL细化算法
- `MIMO_hierarchical_SBL_refined_2D.m` - 分层2D SBL细化算法
- `MIMO_OMP.m` - OMP角度估计算法
- `MIMO_traditional_beamforming.m` - 传统波束成形算法

### 支持函数
- `MIMO_Array_Response_Tx.m` - 发射阵列响应函数
- `MIMO_Array_Response_Rx.m` - 接收阵列响应函数
- `MIMO_First_Order_Linear_Approximation.m` - 一阶线性近似(1D)
- `MIMO_First_Order_Linear_Approximation_2D.m` - 一阶线性近似(2D)
- `MIMO_Gen_measurement_matrix.m` - 测量矩阵生成函数
- `MIMO_wt_derivation_aod.m` - 发射角导数函数
- `MIMO_wt_derivation_aoa.m` - 接收角导数函数

## 使用方法

### 快速测试
```matlab
% 运行测试脚本验证算法
run('test_mimo_angle_estimation.m')
```

### 完整仿真
```matlab
% 运行完整的性能比较仿真
run('massive_mimo_angle_estimation.m')
```

### 自定义参数
可以在主脚本中修改以下参数：
- 角度分辨率 (`r_aod`, `r_aoa`)
- 天线数量 (`Nt`, `Nr`)
- 路径数 (`P`)
- 仿真帧数 (`N_fram`)
- SNR范围 (`SNR_dB`)

## 算法特点

### 与ODDM的对应关系
- **时延 ↔ 发射方向角(AoD)**
- **多普勒频移 ↔ 接收方向角(AoA)**
- **采样函数 ↔ 阵列响应函数**
- **测量矩阵 ↔ 导向矩阵**

### 性能指标
- **NMSE** - 归一化均方误差
- **CRLB** - 克拉美-罗下界
- **算法复杂度** - 网格点数比较

## 输出结果

仿真完成后将显示：
1. 各算法在不同SNR下的NMSE性能曲线
2. 详细的数值结果表格
3. 算法复杂度分析

## 注意事项

1. 确保MATLAB路径包含此文件夹
2. 建议先运行测试脚本验证环境
3. 完整仿真可能需要较长时间，可适当减少仿真参数
4. 如遇到内存不足，可减少天线数或网格分辨率

## 参考

本实现基于ODDM信道估计算法，将其核心思想应用于Massive MIMO方向角估计问题。