function rgb_im = tonemap_durand_local(rad_map, target_contrast, hsize, sigma_d, sigma_r, sat)
% target_contrast - for computing compression factor, 3 - 20 seems appropriate
% hsize - size of filter (must be odd)
% sigmda_d - stdev of Gaussian filter for spatial domain
% sigmda_r - stdev of Gaussian filter for intensity domain

intensity = compute_luminance_map(rad_map);

log_intensity = log(intensity);
log_base = bil_filt(log_intensity, hsize, sigma_d, sigma_r);
log_detail = log_intensity - log_base;
target_contrast = log(target_contrast);
compression_factor = target_contrast / (max(log_base(:)) - min(log_base(:)));
log_absolute_scale = max(log_base(:)) * compression_factor;
log_output_intensity = log_base * compression_factor + log_detail - log_absolute_scale;

L_d = exp(log_output_intensity);

% Preserve color while compressing dynamic range
[r_out, g_out, b_out] = perform_range_reduction(rad_map, intensity, sat, L_d);

rgb_im(:,:,1) = r_out;
rgb_im(:,:,2) = g_out;
rgb_im(:,:,3) = b_out;

end
