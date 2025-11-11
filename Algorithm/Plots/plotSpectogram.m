%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            plotSpectogram.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Plots
% Author: Pablo, Marti, Dan and Fredo
% Date: 9th November 2025
% Last Update: 9th November 2025
% Version: 1.0.0
%
% Description: Plots the 3D Spectogram of the signal.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSpectogram(t,f,s)
    figure()
    surf(t, f, s);
    shading interp;
    box on;
    axis tight;
    view(120, 30);
    title('3D Spectrogram (Surf Plot)');
    xlabel('Time [ms]');
    ylabel('Frequency [MHz]');
    zlabel('Power [dB]');
    fprintf('Displaying plots.\n');
end
