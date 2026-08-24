%% 探地雷达墙内导线检测快速 Demo - MATLAB 版
% Ricker 子波 → B-scan → SAR 聚焦成像 → 导线检测
% 纯 MATLAB 实现

clear; clc; close all;
rng(42);
fprintf('========================================\n');
fprintf('探地雷达墙内导线检测 Demo (MATLAB)\n');
fprintf('========================================\n');

%% 参数
c = 3e8;            % 光速
eps_r = 6.0;        % 墙体介电常数
v = c / sqrt(eps_r); % 波速
fc = 1.5e9;         % 中心频率
dt = 1 / 6e9;       % 时间步长
n_fft = 256;        % FFT点数

%% 1. Ricker 子波
fprintf('\n[1/3] Ricker 子波生成...\n');
t_ricker = ((0:n_fft-1) * dt) - (n_fft * dt / 2);
tau = pi * fc * t_ricker;
wavelet = (1 - 2 * tau.^2) .* exp(-tau.^2);

%% 2. B-scan 仿真
fprintf('\n[2/3] 仿真 B-scan...\n');
wires = [0.0, 0.12; 0.20, 0.08; -0.15, 0.15];  % [x, depth]
fprintf('  导线数: %d\n', size(wires, 1));
fprintf('  导线位置: [%.2f, %.2f], [%.2f, %.2f], [%.2f, %.2f]\n', wires');

n_pos = 60;
span = 0.8;
positions = linspace(-span/2, span/2, n_pos);
bscan = zeros(n_pos, n_fft);

for i = 1:n_pos
    x_ant = positions(i);
    depths = [];
    for w = 1:size(wires, 1)
        dx = abs(x_ant - wires(w, 1));
        if dx < 0.15
            depths(end+1) = sqrt(wires(w, 2)^2 + dx^2);
        end
    end
    
    % A-scan
    signal = zeros(1, n_fft);
    signal = signal + 0.2 * circshift(wavelet, 1);
    signal = signal - 0.1 * circshift(wavelet, 3);
    
    for d = depths
        t_delay = 2 * d / v;
        delay = round(t_delay / dt);
        if delay < n_fft
            signal = signal + 0.5 * exp(-d * 8) * circshift(wavelet, delay);
        end
    end
    bscan(i, :) = signal + 0.01 * randn(1, n_fft);
end

depth_axis = t_ricker * v / 2 * 100;  % cm

%% 3. SAR 后向投影
fprintf('\n[3/3] SAR 聚焦成像...\n');
dx_sar = 0.01;
dz_sar = 0.005;
xs_sar = -0.4:dx_sar:0.4;
zs_sar = 0.01:dz_sar:0.25;
sar_img = zeros(length(zs_sar), length(xs_sar));

for ix = 1:length(xs_sar)
    for iz = 1:length(zs_sar)
        for ip = 1:length(positions)
            r = 2 * sqrt((xs_sar(ix)-positions(ip))^2 + zs_sar(iz)^2);
            idx = r / v / dt;
            if idx >= 1 && idx < n_fft-1
                i0 = floor(idx);
                i1 = i0 + 1;
                f = idx - i0;
                sar_img(iz, ix) = sar_img(iz, ix) + (1-f)*bscan(ip, i0) + f*bscan(ip, i1);
            end
        end
    end
end
fprintf('  SAR 图像分辨率: %d x %d\n', size(sar_img));

%% 4. 绘图
figure('Position', [100 100 1400 500]);

subplot(1, 2, 1);
imagesc(positions*100, depth_axis, bscan');
axis xy; colormap(gray); xlabel('天线位置 (cm)'); ylabel('深度 (cm)');
title('原始 B-scan（双曲线形态）'); colorbar;

subplot(1, 2, 2);
imagesc(xs_sar*100, zs_sar*100, sar_img);
axis xy; colormap(hot); hold on;
scatter(wires(:,1)*100, wires(:,2)*100, 100, 'c', 'x', 'LineWidth', 2);
xlabel('水平位置 (cm)'); ylabel('深度 (cm)');
title('SAR 聚焦成像（点目标）'); colorbar;

saveas(gcf, 'gpr_bscan.png');
fprintf('\n已生成: gpr_bscan.png\n');

% 峰值检测
sar_norm = sar_img / max(sar_img(:));
threshold = 0.4;
peaks = imregionalmax(sar_norm);
[peak_rows, peak_cols] = find(peaks & sar_norm > threshold);

fprintf('\n检测到 %d 根导线:\n', length(peak_rows));
for k = 1:length(peak_rows)
    cx = xs_sar(peak_cols(k));
    cy = zs_sar(peak_rows(k));
    fprintf('  导线 %d: x=%.3f m, z=%.3f m\n', k, cx, cy);
end

figure('Position', [100 100 600 450]);
imagesc(xs_sar*100, zs_sar*100, sar_img);
axis xy; colormap(hot); hold on;
scatter(wires(:,1)*100, wires(:,2)*100, 120, 'c', 'x', 'LineWidth', 3, 'DisplayName', '真实');
for k = 1:min(3, length(peak_rows))
    scatter(xs_sar(peak_cols(k))*100, zs_sar(peak_rows(k))*100, 80, 'g', 'o', ...
        'LineWidth', 2, 'DisplayName', sprintf('检测 #%d', k));
end
xlabel('水平位置 (cm)'); ylabel('深度 (cm)'); title('导线检测结果');
legend('Location', 'best');

saveas(gcf, 'sar_image.png');
fprintf('已生成: sar_image.png\n');

fprintf('\n========================================\n');
fprintf('Demo 完成！\n');
fprintf('========================================\n');
