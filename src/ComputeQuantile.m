%*******************************************************************************
% ComputeQuantile
%
% Compute quantiles from a vector list
% Note that MATLAB's built-in does interpolation while this function does not
%*******************************************************************************
function Q = ComputeQuantile(Y, T)
    Y = sort(Y);
    N = length(Y);
    Q = zeros(length(T), 1);
    for t = 1:1:length(T)
        idx = ceil(T(t)*N);
        Q(t) = Y(idx);
    end
end
