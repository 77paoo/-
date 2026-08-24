# 电磁波墙体探测快速入门 Demo

## 目标
用 Python 仿真探地雷达（GPR）检测墙内导线 —— 从 A-scan 到 B-scan 再到 SAR 聚焦成像，完整走通。

## 前置条件
```bash
pip install numpy matplotlib scipy
```

## 运行方法
```bash
cd 04-em-wall-detection
python demo_gpr.py
```

## 预期输出
1. 终端打印各导线检测位置和误差
2. 生成 `gpr_bscan.png`（原始 B-scan，可见双曲线）
3. 生成 `sar_image.png`（SAR 聚焦后，双曲线→点目标）

## 核心流程

```
Ricker脉冲 → A-scan单道 → 多道B-scan → 后向投影SAR → 导线检测
```

## 你可以尝试

1. 修改导线深度和位置，观察 B-scan 双曲线形态变化
2. 比较不同阵列孔径对成像分辨率的影响
3. 添加噪声后观察对检测精度的影响
