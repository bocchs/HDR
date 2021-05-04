function w = weighting_func(Z_min, Z_max)
% get hat weighting function for use in computing response function
w = zeros(Z_max,1);
for z = 0:Z_max
    if z <= .5*(Z_min + Z_max)
        w(z+1) = z - Z_min;
    else
        w(z+1) = Z_max - z;
    end
end
end
