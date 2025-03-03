function XConv = convolve3DGaussian(XMat,hsize,sigma)
[X, Y, Z] = ndgrid(-floor(hsize/2):floor(hsize/2));

% Calculate the Gaussian kernel
gaussian_kernel = exp(-(X.^2 + Y.^2 + Z.^2) / (2 * sigma^2));
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));  % Normalize the kernel

A = ones(size(XMat));
A(isnan(XMat)) = 0;

M = XMat;
index=find(A==0);
M(index) = 0;

counts = convn(A, gaussian_kernel, 'same');
sums = convn(M, gaussian_kernel, 'same');
XConv= sums ./counts .* A;
