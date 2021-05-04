function rgb_im = tonemap_reinhard_global(rad_map, a, sat)
% http://www.cmap.polytechnique.fr/~peyre/cours/x2005signal/hdr_photographic.pdf
% tonemap radiance map to viewable hdr image
% rad_map - radiance map
% a - key value
% sat - saturation parameter, usually between 0.4 and 0.6

d = 1e-4; % small offset so log(0) does not occur
N = size(rad_map,1) * size(rad_map,2); % number of pixels

lum_map = compute_luminance_map(rad_map);
L_avg = exp(1 / N * sum(log(lum_map(:) + d)));
L = a / L_avg * lum_map;
L_d = L .* (1 + L ./ max(L(:))^2) ./ (1 + L);

% Preserve color while compressing dynamic range
[r_out, g_out, b_out] = perform_range_reduction(rad_map, lum_map, sat, L_d);

rgb_im(:,:,1) = r_out;
rgb_im(:,:,2) = g_out;
rgb_im(:,:,3) = b_out;

end
