close all; clc; clear;

% 读取图像并转换为灰度
I = imread("9.jpg");
I_coloured = I;

if size(I, 3) == 3
    I = rgb2gray(I);
else
    I_coloured = cat(3, I, I, I);
end

I = im2double(I);

% 计算RMS对比度
wsize = [7 7];
h = ones(wsize) / prod(wsize);
Iave = imfilter(I, h, 'replicate');

normDiff = (I - Iave) ./ (Iave + eps);
normDiff = normDiff .^ 2;
Cr = sqrt(imfilter(normDiff, h, 'replicate'));

stats = summary(Cr(:));
disp('RMS 对比度统计特征：');
disp(stats);

% 显示直方图
figure;
histogram(Cr, 'FaceColor', [.4 .8 .6], 'EdgeColor', 'none'); xlabel('RMS 对比度'); ylabel('频率');
title(sprintf('RMS 对比度分布情况（wsize=[%d %d]）', wsize(1), wsize(2)));

% 显示热力图
figure;
imagesc(reshape(Cr, size(I)));
colormap pink; colorbar; title(sprintf('RMS 对比度热力图（wsize=[%d %d]）', wsize(1), wsize(2))); axis image;

% 显示对比度可视化
figure; imshow(Cr); title('RMS 对比度重建图像结果');

% 基于RMS对比度进行图像增强
Cr = rescale(Cr, 0, 1);

alpha = 1.2;
beta = 1.6;
smooth_sigma = 1.0;

G = 1 + alpha * (1 - Cr) .^ beta;
G = imgaussfilt(G, smooth_sigma);

Ie = Iave + G .* (I - Iave);
Ie = min(max(Ie, 0), 1);

f = restoreColour(I_coloured, Ie .^ .6);
figure; imshow(f); title('RMS 对比度图像增强结果');