%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            plotFFT.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Plots
% Author: Pablo, Marti, Dan and Fredo
% Date: 9th November 2025
% Last Update: 9th November 2025
% Version: 1.0.0
%
% Description: Plots the Fast Fourier Transform of the signal in dB.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotFFT(fShifted, signalFFT)
    figure(1);
    plot(fShifted, 10*log10(abs(fftshift(signalFFT))));
    title('"Noisy" FFT');
    xlabel('Frequency [MHz]');
    ylabel('Amplitude [dB]');
    grid on;
    axis tight;
end
