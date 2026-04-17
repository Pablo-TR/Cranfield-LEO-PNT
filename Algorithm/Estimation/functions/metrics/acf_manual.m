function [acf, lags] = acf_manual(x, maxLag)
%ACF_MANUAL Toolbox-free autocorrelation normalised by the zero-lag value.

x = x(:);
x = x(~isnan(x));
if isempty(x)
    lags = (-maxLag:maxLag).';
    acf = nan(size(lags));
    return;
end

x = x - mean(x);
den = sum(x.^2);
lags = (-maxLag:maxLag).';
acf = nan(size(lags));

if den <= eps
    return;
end

n = numel(x);
for ii = 1:numel(lags)
    lag = lags(ii);
    if lag >= 0
        a = x(1+lag:n);
        b = x(1:n-lag);
    else
        L = -lag;
        a = x(1:n-L);
        b = x(1+L:n);
    end
    acf(ii) = sum(a .* b) / den;
end
end
