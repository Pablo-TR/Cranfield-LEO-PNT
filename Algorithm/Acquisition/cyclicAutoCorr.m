%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            cyclicAutoCorr.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: -
% Author: Pablo, Fredo, Marti & Dan
% Date: 9th November 2025
% Last Update: 21st November 2025
% Version: 1.0.0
%
% Description: Computes the cyclic autocorrelation function of a complex-valued 
% signal Z for a given set of cyclic frequencies. The input signal is first
% circularly shifted by k samples and then processed to estimate the cyclic
% autocorrelation at each specified cyclic frequency.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

