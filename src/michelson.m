close all; clc; clear;

%% 读取图像
I = imread("9.jpg");
I_coloured = I;

if size(I, 3) == 3
    I = rgb2gray(I);
else
    I_coloured = cat(3, I, I, I);
end

I = im2double(I);

%% 计算 Michelson 对比度
wsize = [7 7];  % 领域大小
Imax = ordfilt2(I, wsize(1)*wsize(2), true(wsize));
Imin = ordfilt2(I, 1, true(wsize));

Cm = (Imax - Imin) ./ (Imax + Imin + eps);

% 获取统计特征
stats = summary(Cm(:));
disp('Michelson 对比度统计特征：');
disp(stats);

% 1. 使用直方图显示分布情况
figure('Name', 'Michelson 对比度直方图');
histogram(Cm, 'FaceColor', [.4 .8 .6], 'EdgeColor', 'none'); xlabel('Michelson 对比度'); ylabel('频率');
title(sprintf('Michelson 对比度分布情况（wsize=[%d %d]）', wsize(1), wsize(2)));

% 2. 显示热力图
figure;
imagesc(reshape(Cm, size(I)));
colormap pink; colorbar; title(sprintf('Michelson 对比度热力图（wsize=[%d %d]）', wsize(1), wsize(2))); axis image;

% 3. 直接显示 Michelson 对比度
figure; imshow(Cm); title(sprintf('Michelson 对比度可视化（wsize=[%d %d]）', wsize(1), wsize(2)));

%% 使用 Michelson 对比度进行图像增强（对比度大小正比于 Imax/Imin）
target_cm = .8;
gamma = .6;

% Cm = 1 - 2 / (k + 1)
target_ratio = 2 / (1 - target_cm) - 1;

[rows, cols] = size(I);
m = wsize(1);
n = wsize(2);
pad_m = floor(m / 2);
pad_n = floor(n / 2);

I_padded = padarray(I, [pad_m, pad_n], 'replicate', 'both'); 
B_matrix = im2col(I_padded, [m n], 'sliding');

I_min_vec = min(B_matrix, [], 1); 
I_max_vec = max(B_matrix, [], 1); 
Cm_vec = Cm(:)'; 

% 确定需要增强的窗口索引
valid_indices = find(I_max_vec ~= I_min_vec & Cm_vec < 0.6);
B_enhanced_matrix = B_matrix; 

% 预计算增强块
for k = valid_indices
    I_min_k = I_min_vec(k);
    I_max_k = I_max_vec(k);
    
    uplimit = min(1, I_min_k * target_ratio);
    if I_max_k < .3, uplimit = max(uplimit, I_max_k * 2.4); end
    
    B_k = B_matrix(:, k);
    
    % 线性放大
    B_new_k = rescale(B_k, I_min_k, uplimit);
    B_enhanced_matrix(:, k) = B_new_k;
end


I_enhanced = zeros(rows, cols, 'double');
% I_enhanced_sum = zeros(rows, cols, 'double');
% I_count = zeros(rows, cols, 'double');

k_index = 0;
for j = 1:cols
    for i = 1:rows
        k_index = k_index + 1;
        [r_start, r_end, c_start, c_end] = get_valid_range(i, j, m, n, rows, cols);
        r_offset = r_start - i + 1;
        r_offset_end = r_end - i + 1;
        c_offset = c_start - j + 1;
        c_offset_end = c_end - j + 1;
        
        B_new = reshape(B_enhanced_matrix(:, k_index), [m, n]);
        
        % 不采用累加平均
        I_enhanced(r_start:r_end, c_start:c_end) = B_new(r_offset:r_offset_end, c_offset:c_offset_end);

        % 简单累加（权重=1）
        % I_enhanced_sum(r_start:r_end, c_start:c_end) = I_enhanced_sum(r_start:r_end, c_start:c_end) + B_new(r_offset:r_offset_end, c_offset:c_offset_end);
        % I_count(r_start:r_end, c_start:c_end) = I_count(r_start:r_end, c_start:c_end) + 1;
    end
end

% I_enhanced = I_enhanced_sum ./ I_count;
I_eg = I_enhanced .^ gamma;
%%
% figure; 
% subplot(211); imshow(I); title('原图');
% subplot(212); imshow(I_enhanced); title('增强后');

gamma = .6;

I_ec = restoreColour(I_coloured, I_enhanced .^ gamma);
figure;
imshow(I_ec); title('Michelson 对比度图像增强结果');
% subplot(211); imshow(I_coloured);
% subplot(212); imshow(I_ec);

%% 判断是否越界
function [r_start, r_end, c_start, c_end] = get_valid_range(i, j, m, n, rows, cols)
    r_start = i;
    r_end = i + m - 1;
    c_start = j;
    c_end = j + n - 1;
    
    r_start = max(r_start, 1);
    r_end = min(r_end, rows);
    c_start = max(c_start, 1);
    c_end = min(c_end, cols);
end
