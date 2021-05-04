clear

% base_folder = 'Garage_02/';
% base_folder = 'Garage_03/';
% base_folder = 'Garage_04_foil/';
% base_folder = 'Garage_04_less/';
% base_folder = 'Garage_04/';
base_folder = 'Garage_05b/';
% base_folder = 'Porch_02/';
% base_folder = 'Study_01/';
% base_folder = 'GolfSun_01a/';
% base_folder = 'Living_01/';
% base_folder = 'BookPainting/';
% base_folder = 'CameraBook/';
base_folder = strcat('images/',base_folder);

num_points = 500;
lambda = 1;
% get radiance map and response functions for each channel
[rad_map, g1, g2, g3, exposures, etimes] = compute_radiance_map(base_folder, num_points, lambda);

num_exposures = length(etimes);
[etimes,inds] = sort(etimes); % order exposures from fastest to slowest

for k = 1:num_exposures
    subplot(ceil(sqrt(num_exposures)),ceil(sqrt(num_exposures)),k)
    imshow(exposures{inds(k)})
    title(num2str(etimes(k)))
end
sgtitle('Exposure Times (sec)') 

% plot response functions
pixel_range = 0:255;
figure(2)
plot(pixel_range,g1,'r')
hold on
plot(pixel_range,g2,'g')
plot(pixel_range,g3,'b')
hold off
xlabel('Pixel Value')
ylabel('log Exposure')
title('Response Function')
xlim([0 255])

figure(3)
plot(g1,pixel_range,'r')
hold on
plot(g2,pixel_range,'g')
plot(g3,pixel_range,'b')
hold off
xlabel('log Exposure')
ylabel('Pixel Value')
title('Response Function')
ylim([0 255])

% rad maps used in HDR/tonemapping papers
% rad_map = hdrread('memorial.hdr');
% rad_map = hdrread('nave.hdr');


a = 1;
sat = .5;
hdr_im_reinhard_global = tonemap_reinhard_global(rad_map, a, sat);
% subplot(2,2,1)
f = figure();
ax = axes(f);
imshow(hdr_im_reinhard_global, [], 'Parent', ax)
title(ax, 'Reinhard Global')

a = 1;
sat = .5;
[hdr_im_reinhard_local, bad_points_map] = tonemap_reinhard_local(rad_map, a, sat);
% subplot(2,2,2)
f = figure();
ax = axes(f);
imshow(hdr_im_reinhard_local, [], 'Parent', ax)
title(ax, 'Reinhard Local')


sat = .6;
target_contrast = 5;
hsize = 5;
sigma_d = 2; % 5 % spatial kernel sigma 
sigma_r = 2; % 1 % intensity kernel sigma
hdr_im_durand = tonemap_durand_local(rad_map, target_contrast, hsize, sigma_d, sigma_r, sat);
f = figure();
ax = axes(f);
imshow(hdr_im_durand, [], 'Parent', ax)
title(ax, 'Durand Local')

hdr_im_matlab_global = tonemap(rad_map);
% convert from uint8 to double [0-1]
hdr_im_matlab_global = im2double(hdr_im_matlab_global);
% subplot(2,2,3)
f = figure();
ax = axes(f);
imshow(hdr_im_matlab_global, [], 'Parent', ax)
title(ax, 'Matlab Global')

hdr_im_matlab_local = localtonemap(single(rad_map));
hdr_im_matlab_local = double(hdr_im_matlab_local);
% subplot(2,2,4)
f = figure();
ax = axes(f);
imshow(hdr_im_matlab_local, [], 'Parent', ax)
title(ax, 'Matlab Local')

figure()
img_array(:,:,:,1) = hdr_im_reinhard_global;
img_array(:,:,:,2) = hdr_im_reinhard_local;
img_array(:,:,:,3) = hdr_im_durand;
img_array(:,:,:,4) = hdr_im_matlab_global;
img_array(:,:,:,5) = hdr_im_matlab_local;

montage(img_array, 'size', [1 NaN])
