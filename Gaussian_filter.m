% https://www.mathworks.com/matlabcentral/fileexchange/46189-generating-gaussian-filter-2d-matrix
% create gaussian filter
function g = Gaussian_filter(size, sigma)
%filter size, odd number
g=zeros(size,size); 
x0=(size+1)/2; %center
y0=(size+1)/2; %center
for i=-(size-1)/2:(size-1)/2
    for j=-(size-1)/2:(size-1)/2
        x=i+x0; %row
        y=j+y0; %col
        g(y,x)= 1/(2*pi*sigma^2) * exp(-((x-x0)^2+(y-y0)^2)/(2*sigma^2));
    end
end
%normalize gaussian filter
% sum1=sum(g);
% sum2=sum(sum1);
% g=g/sum2;
end
