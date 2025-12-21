% 读取图像并转换为灰度
I = imread("9.jpg");
I_coloured = I;

if size(I, 3) == 3
    I = rgb2gray(I);
else
    I_coloured = cat(3, I, I, I);
end

I = im2double(I);

% 计算局部平均值
window_width = 71;
sigma = (window_width - 1) / 6;
Is = imgaussfilt(I, sigma, 'FilterSize', window_width); fname = '高斯滤波器';

% 多尺度滤波融合
w = [5 71 351];
s = (w - 1) / 6;
Imf = zeros([size(I) 3]);
Igf = zeros([size(I) 3]);

for i = 1:3
    Igf(:, :, i) = imgaussfilt(I, s(i), 'FilterSize', w(i));
    tmp_m = fspecial("average", w(i));
    Imf(:, :, i) = imfilter(I, tmp_m, 'replicate');
end

figure; imshow(Imf); title('均值滤波器融合结果（5、71、351分别对应R、G、B三通道）');
figure; imshow(Igf); title('高斯滤波器融合结果（5、71、351分别对应R、G、B三通道）');

% 计算Weber对比度
Cw = (I - Is) ./ (Is + eps);

stats = summary(Cw(:));
disp('Weber 对比度统计特征：');
disp(stats);

% 显示直方图
figure('Name', sprintf('Weber 对比度直方图（%s，%d）', fname, window_width));
histogram(Cw, 'FaceColor', [.4 .8 .6], 'EdgeColor', 'none'); xlabel('Weber 对比度'); ylabel('频率');
title(sprintf('Weber 对比度分布情况（%s，%d）', fname, window_width));

% 显示热力图
figure;
imagesc(reshape(Cw, size(I)));
colormap pink; colorbar; title(sprintf('Weber 对比度热力图（%s，%d）', fname, window_width)); axis image;

% 对比度可视化与变换
gamma = .4;
Cw_disp = mat2gray(Cw);

I_hsv = rgb2hsv(I_coloured);
I_hsv(:, :, 3) = Cw_disp;
Iw = hsv2rgb(I_hsv);

I_hsv(:, :, 3) = Cw_disp .^ gamma;
Iwg = hsv2rgb(I_hsv);

I_hsv(:, :, 3) = histeq(Cw_disp);
Iwh = hsv2rgb(I_hsv);

figure; imshow(Iw); title(sprintf('Weber 对比度重建结果（%s，%d）', fname, window_width));
figure; imshow(Iwg); title(sprintf('伽马变换后的 Weber 对比度重建结果（%s，%d，gamma=%0.1f）', fname, window_width, gamma));
figure; imshow(Iwh); title(sprintf('直方图均衡化后的 Weber 对比度重建结果（%s，%d）', fname, window_width));

% 基于Weber对比度进行图像增强
alpha = .8;
beta = .85;
gamma = .4;

Cw_new = abs(Cw) .^ alpha;
Cw_new = beta .* sign(Cw) .* Cw_new;
Is_new = Is .^ gamma;

I_out = Is_new .* (1 + Cw_new);
I_out(I_out > 1) = 1.;
I_out(I_out < 0) = 0.;

figure;
f = restoreColour(I_coloured, I_out); imshow(f); title(sprintf('Weber 对比度图像增强结果（%s，%d）', fname, window_width));