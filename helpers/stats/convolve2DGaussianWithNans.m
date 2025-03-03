function XConv = convolve2DGaussianWithNans(X, hsize, sigma, gaussianKernel)
    if nargin < 4
        gaussianKernel = fspecial('gaussian', hsize, sigma);
    end
    
    % Normalize the Gaussian kernel
    gaussianKernel = gaussianKernel / sum(gaussianKernel(:));
    
    % Get the size of the data and the kernel
    [rows, cols] = size(X);
    [kRows, kCols] = size(gaussianKernel);
    
    % Initialize the output matrix
    XConv = nan(size(X));
    
    % Loop through each element of the input data
    for i = 1:rows
        for j = 1:cols
            % Initialize variables to calculate the weighted sum and the weight sum
            weightedSum = 0;
            weightSum = 0;
            
            % Loop through the elements of the kernel
            for m = 1:kRows
                for n = 1:kCols
                    % Calculate the coordinates in the input data
                    x = i - floor(kRows / 2) + m;
                    y = j - floor(kCols / 2) + n;
                    
                    % Check if the coordinates are within the data bounds
                    if x >= 1 && x <= rows && y >= 1 && y <= cols
                        % Check if the element is not NaN or -Inf
                        if ~isnan(X(x, y)) && X(x, y) > -Inf && X(x, y) < Inf
                            % Update the weighted sum and weight sum
                            weightedSum = weightedSum + X(x, y) * gaussianKernel(m, n);
                            weightSum = weightSum + gaussianKernel(m, n);
                        end
                    end
                end
            end
            
            % Update the output matrix with the weighted average
            if weightSum > 0
                XConv(i, j) = weightedSum / weightSum;
            end
        end
    end
end
