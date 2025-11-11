%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            plotAutocorrelation.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Plots
% Author: Pablo, Marti, Dan and Fredo
% Date: 9th November 2025
% Last Update: 11th November 2025
% Version: 1.0.0
%
% Description: Plots the autocorrelation of the signal.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotAutocorrelation(alphaHat, Rnormalized)
    figure();
    plot(alphaHat .* 1e-6, Rnormalized, 'LineWidth', 1.5); % X-axis in MHz
    title('Normalised Cyclic Spectrum');
    xlabel('Cyclic Frequency (alpha) [MHz]');
    ylabel('Normalised Amplitude');
    grid on;
    axis tight;
    ylim([0, 1.1]); % Force Y-axis from 0 to ~1
end
