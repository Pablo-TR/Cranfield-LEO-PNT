function y = safe_std(x, dim)
%SAFE_STD Standard deviation with NaN omission using base MATLAB syntax.

if nargin < 2
    y = std(x, 'omitnan');
else
    y = std(x, 0, dim, 'omitnan');
end
end
