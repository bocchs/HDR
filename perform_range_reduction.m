function [r_out, g_out, b_out] = perform_range_reduction(rad_map, lum_map, s, L_d)
% Scale HDR image to LDR image while preserving color
%
% rad_map - radiance map
% lum_map - luminance map
% s - saturation parameter, usually between 0.4 and 0.6
% L_d - luminance map after HDR compression
%
% https://www.cs.huji.ac.il/~danix/hdr/hdrc.pdf
% C_out = (C_in / L_in)^s * L_out
% s is usually between 0.4 and 0.6
% For preserving color while reducing dynamic range
% https://www.cl.cam.ac.uk/~rkm38/pdfs/mantiuk09cctm.pdf

r = rad_map(:,:,1);
b = rad_map(:,:,2);
g = rad_map(:,:,3);

r_out = (r ./ lum_map) .^s .* L_d;
g_out = (b ./ lum_map) .^s .* L_d;
b_out = (g ./ lum_map) .^s .* L_d;

% clamp channels to 1
% view histogram of a channel to visualize about how many pixels are clipped
% histogram(r_out(:));
r_out(r_out > 1) = 1;
g_out(g_out > 1) = 1;
b_out(b_out > 1) = 1;

end