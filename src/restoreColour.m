function recoverd_img = restoreColour(img, weber)
    recoverd_img = rgb2hsv(img);
    recoverd_img(:, :, 3) = weber;
    recoverd_img(:, :, 2) = recoverd_img(:, :, 2) * 1.2;
    recoverd_img = hsv2rgb(recoverd_img);
end