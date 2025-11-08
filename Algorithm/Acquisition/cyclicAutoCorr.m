function [R, alpha] = cyclicAutoCorr(Z, k, Fsr, alpha_hat)
    Z = circshift(Z,-k);
    alpha = alpha_hat./Fsr;
    N = length(Z);
    n = 0:1:N-1; 
    R = zeros(1,length(alpha));
    for j= 1:1:length(alpha)
        sum_arr = Z.*conj(Z).*exp(-2*pi*alpha(j).*n.*1i);
        R(j) = (1/N)*sum(sum_arr);
    end
end

