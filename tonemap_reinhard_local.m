function [rgb_im, bad_points_map] = tonemap_reinhard_local(rad_map, a, sat)
% http://www.cmap.polytechnique.fr/~peyre/cours/x2005signal/hdr_photographic.pdf
% tonemap radiance map to viewable hdr image
% rad_map - radiance map
% a - key value
% sat - saturation parameter, usually between 0.4 and 0.6

num_scales = 8;
kernel_scale = 1.6;
alpha1 = .35;
alpha2 = 1.6 * alpha1;
phi = 15; % sharpening
eps = .08;
N = size(rad_map,1) * size(rad_map,2); % number of pixels
d = 1e-4; % small offset so log(0) does not occur


lum_map = compute_luminance_map(rad_map);
L_avg = exp(1 / N * sum(log(lum_map(:) + d)));
L = a / L_avg * lum_map;
for s_idx = 1:num_scales
    
    s_R1 = kernel_scale ^ (s_idx-1);
    kernel_size_R1 = ceil(s_R1);
    if mod(kernel_size_R1,2) == 0
        kernel_size_R1 = kernel_size_R1 + 1;
    end
    
    s_R2 = kernel_scale ^ s_idx;
    kernel_size_R2 = ceil(s_R2);
    if mod(kernel_size_R2,2) == 0
        kernel_size_R2 = kernel_size_R2 + 1;
    end
    
    R1 = Gauss_filt_reinhard(kernel_size_R1, alpha1, s_R1);
    R2 = Gauss_filt_reinhard(kernel_size_R2, alpha2, s_R2);

    V1_temp = conv2(L, R1, 'same');
    V2_temp = conv2(L, R2, 'same');
    V(:,:,s_idx) = (V1_temp - V2_temp) ./ (2^phi * a ./ s_R1.^2 + V1_temp);
    V1(:,:,s_idx) = V1_temp;
end

rows = size(lum_map,1);
cols = size(lum_map,2);
V1_all_scales = zeros(rows, cols);
% keep track of points where threshold is surpassed
bad_points_map = zeros(rows,cols);
for y = 1:rows
    for x = 1:cols
        % get first scale where V entry < eps
        entry_vals = V(y,x,:);
        entry_vals = entry_vals(:);
        s_idx = find(entry_vals < eps, 1, 'first');
        if isempty(s_idx)
            % if no V's are small enough, select smallest
            [~, s_idx] = min(entry_vals);
            bad_points_map(y,x) = 1;
        end
        % use that scale
        V1_all_scales(y,x) = V1(y,x,s_idx);
    end
end
L_d = L ./ (1 + V1_all_scales);

% Preserve color while compressing dynamic range
[r_out, g_out, b_out] = perform_range_reduction(rad_map, lum_map, sat, L_d);

rgb_im(:,:,1) = r_out;
rgb_im(:,:,2) = g_out;
rgb_im(:,:,3) = b_out;

end
