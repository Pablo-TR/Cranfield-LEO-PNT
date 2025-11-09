%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            plotPSD.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Plots
% Author: Pablo, Marti, Dan and Fredo
% Date: 9th November 2025
% Last Update: 9th November 2025
% Version: 1.0.0
%
% Description: Plots the Power Spectral Density of the signal in dB.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotPSD(fpsd, pxx)
    plot(fpsd, 10*log10(pxx)); 
    title('"Smooth" PSD (pwelch with averaging)');
    xlabel('Frequency [MHz]');
    ylabel('Power [dB]');
    grid on;
    axis tight;
end
