function y = safe_mean(x, dim)
%SAFE_MEAN Mean with NaN omission using base MATLAB syntax.

if nargin < 2
    y = mean(x, 'omitnan');
else
    y = mean(x, dim, 'omitnan');
end
end
