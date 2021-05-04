function J = bil_filt(input_img, hsize, sigma_d, sigma_r)
% evaluate bilateral filter on input image for Durand tonemapping
% hsize - size of filter
% sigmda_d - stdev of Gaussian filter for spatial domain
% sigmda_r - stdev of Gaussian filter for intensity domain

% spatial domain weighting
f = Gaussian_filter(hsize, sigma_d);
offset = (hsize-1)/2;

rows = size(input_img, 1);
cols = size(input_img, 2);
J = zeros(rows, cols);
for y = 1:rows
   for x = 1:cols
         y_min = max(y-offset, 1);
         y_max = min(y+offset, rows);
         x_min = max(x-offset, 1);
         x_max = min(x+offset, cols);
         I = input_img(y_min:y_max, x_min:x_max); % current working window
      
         % intensity domain weighting
         g = exp(-(I - input_img(y,x)).^2 / (2 * sigma_r^2));

         fg = f((y_min:y_max) - y + offset + 1, (x_min:x_max) - x + offset + 1) .* g;
         norm_term = sum(fg(:));
         J(y,x) = 1 / norm_term * sum(fg(:) .* I(:));    
   end
end
end
