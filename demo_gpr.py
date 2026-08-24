#!/usr/bin/env python3
"""
探地雷达（GPR）墙内导线检测快速 Demo
Ricker 子波 → B-scan 仿真 → SAR 聚焦成像
"""

import numpy as np
import matplotlib.pyplot as plt

C = 3e8  # 光速

# ========== 1. Ricker 子波 ==========
def ricker(fc, dt, n=256):
    t = np.arange(n)*dt - n*dt/2
    tau = np.pi*fc*t
    return t, (1 - 2*tau**2)*np.exp(-tau**2)

# ========== 2. A-scan 仿真 ==========
def simulate_a_scan(wavelet, dt, wire_positions, eps_r=6.0):
    n = len(wavelet)
    signal = np.zeros(n)
    v = C / np.sqrt(eps_r)
    # 直达波
    signal += 0.2 * np.roll(wavelet, 1)
    # 墙面反射
    signal += -0.1 * np.roll(wavelet, 3)
    # 导线反射
    for depth in wire_positions:
        t_delay = 2 * depth / v
        delay = int(t_delay / dt)
        if delay < n:
            signal += 0.5 * np.exp(-depth*8) * np.roll(wavelet, delay)
    signal += 0.01 * np.random.randn(n)
    return signal

# ========== 3. B-scan 仿真 ==========
def simulate_b_scan(wires, n_pos=60, span=0.8):
    """
    wires: [(x, depth), ...] 导线位置
    """
    fc, dt = 1.5e9, 1/6e9
    _, wav = ricker(fc, dt)
    positions = np.linspace(-span/2, span/2, n_pos)
    bscan = np.zeros((n_pos, len(wav)))
    for i, x_ant in enumerate(positions):
        depths = []
        for wx, wz in wires:
            dx = abs(x_ant - wx)
            if dx < 0.15:
                depths.append(np.sqrt(wz**2 + dx**2))
        bscan[i] = simulate_a_scan(wav, dt, depths)
    return positions, wav, bscan

# ========== 4. SAR 后向投影 ==========
def back_projection(bscan, pos, n_fft, v, dx=0.01, dz=0.005):
    xs = np.arange(-0.4, 0.41, dx)
    zs = np.arange(0.01, 0.25, dz)
    image = np.zeros((len(zs), len(xs)))
    dt = 1/(6e9)
    for ix, x in enumerate(xs):
        for iz, z in enumerate(zs):
            for ip, x_ant in enumerate(pos):
                r = 2 * np.sqrt((x-x_ant)**2 + z**2)
                idx = r / v / dt
                if 0 <= idx < n_fft-1:
                    i0, i1 = int(idx), int(idx)+1
                    f = idx - i0
                    image[iz, ix] += (1-f)*bscan[ip, i0] + f*bscan[ip, i1]
    return xs, zs, image

# ========== 5. 主函数 ==========
def main():
    print("=" * 50)
    print("探地雷达墙内导线检测快速 Demo")
    print("=" * 50)

    # 导线配置
    wires = [(0.0, 0.12), (0.20, 0.08), (-0.15, 0.15)]
    print(f"\n导线配置: {wires}")
    print(f"导线数: {len(wires)}")

    # B-scan
    print("\n仿真 B-scan...")
    pos, t, bscan = simulate_b_scan(wires, n_pos=60, span=0.8)
    depth = t * C / np.sqrt(6.0) * 100  # cm

    # 绘图 B-scan
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    im = axes[0].imshow(bscan.T, aspect='auto', cmap='gray',
               extent=[pos.min()*100, pos.max()*100, depth.max(), depth.min()])
    axes[0].set_xlabel('天线位置 (cm)'); axes[0].set_ylabel('深度 (cm)')
    axes[0].set_title('原始 B-scan（双曲线形态）')
    plt.colorbar(im, ax=axes[0])

    # SAR 成像
    print("SAR 聚焦成像...")
    v = C / np.sqrt(6.0)
    xs, zs, sar = back_projection(bscan, pos, bscan.shape[1], v)

    im2 = axes[1].imshow(sar, aspect='auto', cmap='hot',
                extent=[xs.min()*100, xs.max()*100, zs.max()*100, zs.min()*100])
    axes[1].set_xlabel('水平位置 (cm)'); axes[1].set_ylabel('深度 (cm)')
    axes[1].set_title('SAR 聚焦成像（点目标）')
    plt.colorbar(im2, ax=axes[1])

    # 标记真实导线位置
    for wx, wz in wires:
        axes[1].scatter(wx*100, wz*100, c='cyan', marker='x', s=80, linewidths=2)

    plt.tight_layout()
    plt.savefig('gpr_bscan.png', dpi=150)
    print("已生成: gpr_bscan.png")

    # 检测峰值
    from scipy.ndimage import maximum_filter, label
    sar_norm = sar / sar.max()
    local_max = maximum_filter(sar_norm, size=10)
    maxima = (sar_norm == local_max) & (sar_norm > 0.4)
    labeled, n_feat = label(maxima)
    print(f"\n检测到 {n_feat} 根导线:")
    detected = []
    for rid in range(1, n_feat + 1):
        yi, xi = np.where(labeled == rid)
        if len(yi):
            cz = np.mean(yi) * 0.005
            cx = np.mean(xi) * 0.01 - 0.4
            detected.append((cx, cz))
            # 找最近的真实导线
            best = min(((wx-wz, wx, wz) for wx, wz in wires),
                      key=lambda p: abs(p[1]-cx)+abs(p[2]-cz))
            err = np.sqrt((cx-best[1])**2 + (cz-best[2])**2)*100
            print(f"  导线 {rid}: x={cx:.3f}m, z={cz:.3f}m (误差 {err:.1f}cm)")

    plt.figure(figsize=(8, 5))
    plt.imshow(sar, aspect='auto', cmap='hot',
               extent=[xs.min()*100, xs.max()*100, zs.max()*100, zs.min()*100])
    for wx, wz in wires:
        plt.scatter(wx*100, wz*100, c='cyan', marker='x', s=100,
                   linewidths=2, label='真实' if wx==wires[0][0] else '')
    for cx, cz in detected:
        plt.scatter(cx*100, cz*100, c='lime', marker='o', s=80,
                   facecolors='none', linewidths=2,
                   label='检测' if cx==detected[0][0] else '')
    plt.xlabel('水平位置 (cm)'); plt.ylabel('深度 (cm)')
    plt.title('导线检测结果 (真实 vs 检测)')
    plt.legend(); plt.tight_layout()
    plt.savefig('sar_image.png', dpi=150)
    print("已生成: sar_image.png")

    print("\nDemo 完成！")

if __name__ == '__main__':
    main()
