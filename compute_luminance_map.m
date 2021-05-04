function L = compute_luminance_map(rad_map)
% https://en.wikipedia.org/wiki/Relative_luminance
% rad_map - 3-channel radiance map
R = rad_map(:,:,1);
G = rad_map(:,:,2);
B = rad_map(:,:,3);
L = .2126 * R + .7152 * G + .0722 * B;
end
