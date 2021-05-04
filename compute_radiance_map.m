function [E, g1, g2, g3, exposures, etimes] = compute_radiance_map(base_folder, num_points, l)
% https://www.pauldebevec.com/Research/HDR/debevec-siggraph97.pdf
% returns radiance map E given a base folder containing multiple exposures
% of an image
% inputs:
% base_folder - contains multiple exposures for one image
% num_points - number of points for computing response curve
% l - lambda: smoothing parameter in optimization problem

% outputs:
% E - radiance map
% g1 - response function for channel 1 (R)
% g2 - response function for channel 2 (G)
% g3 - response function for channel 3 (B)

img_folder = dir(strcat(base_folder,'*.jpg'));
num_images = size(img_folder,1);

fname = img_folder(1).name;
path = strcat(base_folder,fname);
img = imread(path);
% randomly select points to use for computing response curve
y_points = randsample(size(img,1), num_points);
x_points = randsample(size(img,2), num_points);
points = [y_points x_points];
Z1 = zeros(num_points, num_images);
Z2 = zeros(num_points, num_images);
Z3 = zeros(num_points, num_images);

etimes = zeros(num_images,1);
for j = 1:num_images
    fname = img_folder(j).name;
    path = strcat(base_folder,fname);
    img = imread(path);
    exposures{j} = img; % save each exposure for viewing later
    c1 = img(:,:,1); % R channel
    c2 = img(:,:,2); % G channel
    c3 = img(:,:,3); % B channel
    for idx = 1:num_points
        Z1(idx,j) = c1(points(idx,1), points(idx,2));
        Z2(idx,j) = c2(points(idx,1), points(idx,2));
        Z3(idx,j) = c3(points(idx,1), points(idx,2));
    end
    info = imfinfo(path);
    etime = info.DigitalCamera.ExposureTime;
    etimes(j) = etime;
end
Z1_min = min(Z1(:)); % should be (close to) 0
Z1_max = max(Z1(:)); % should be (close to) 255

Z2_min = min(Z2(:));
Z2_max = max(Z2(:));

Z3_min = min(Z3(:));
Z3_max = max(Z3(:));

% obtain weighting hat function
w1 = weighting_func(Z1_min, Z1_max);
w2 = weighting_func(Z2_min, Z2_max);
w3 = weighting_func(Z3_min, Z3_max);

% compute response function for each channel
B = log(etimes);
[g1,lE1] = gsolve(Z1,B,l,w1);
[g2,lE2] = gsolve(Z2,B,l,w2);
[g3,lE3] = gsolve(Z3,B,l,w3);


% Eqn 6 of paper
% compute radiance map
rows = size(img,1);
cols = size(img,2);
log_E1 = zeros(rows, cols); % log radiance map
w1_zij_sum = zeros(rows, cols);
log_E2 = zeros(rows, cols);
w2_zij_sum = zeros(rows, cols);
log_E3 = zeros(rows, cols);
w3_zij_sum = zeros(rows, cols);
for j=1:num_images
    fname = img_folder(j).name;
    path = strcat(base_folder,fname);
    info = imfinfo(path);
    etime = info.DigitalCamera.ExposureTime;
    img = imread(path); % select jth exposure image
    c1 = img(:,:,1);
    c2 = img(:,:,2);
    c3 = img(:,:,3);
    for y = 1:size(img,1)
        for x = 1:size(img,2)
            z1ij = c1(y,x); % pixel intensity
            log_E1(y,x) = log_E1(y,x) + w1(z1ij+1) * (g1(z1ij+1) - log(etime));
            w1_zij_sum(y,x) = w1_zij_sum(y,x) + w1(z1ij+1);
            
            z2ij = c2(y,x);
            log_E2(y,x) = log_E2(y,x) + w2(z2ij+1) * (g2(z2ij+1) - log(etime));
            w2_zij_sum(y,x) = w2_zij_sum(y,x) + w2(z2ij+1);
            
            z3ij = c3(y,x);
            log_E3(y,x) = log_E3(y,x) + w3(z3ij+1) * (g3(z3ij+1) - log(etime));
            w3_zij_sum(y,x) = w3_zij_sum(y,x) + w3(z3ij+1);
        end
    end
end
% after summing over all images, a weight might still be 0
% for most image sequences, a weight's final value will be > 0
% replace any 0's with 1 (next smallest number)
% leaving any 0's will cause NaNs
w1_zij_sum(w1_zij_sum == 0) = 1;
w2_zij_sum(w2_zij_sum == 0) = 1;
w3_zij_sum(w3_zij_sum == 0) = 1;
log_E1 = log_E1 ./ w1_zij_sum;
log_E2 = log_E2 ./ w2_zij_sum;
log_E3 = log_E3 ./ w3_zij_sum;

E1 = exp(log_E1);
E2 = exp(log_E2);
E3 = exp(log_E3);
E = zeros(rows,cols,3);
E(:,:,1) = E1;
E(:,:,2) = E2;
E(:,:,3) = E3;

end
