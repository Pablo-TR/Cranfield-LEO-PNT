%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            peaksCorrelation.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: [Module Name, e.g. Tracking]
% Author: Alfredo
% Date: 7th March 2026
% Last Update: 7th March 2026
% Version: 1.0.0
%
% Description: Outputs the product of the correlation between the signal
% and the local replica against time. Peaks should be seen every [x]ms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

%%
tStart = 0;
tEnd = 10;

signalID = 2;
signalData = getSignalDataByID(signalID);
Fs = signalData.SamplingFreq;
Fsym = 2.304e8;

[signalDetrended, ~, timeArray, binData] = readSignal(signalID, tStart, tEnd);

[~, local] = localReplica(Fs, Fsym);


% Correlation between incoming signal and synchronisation sequence
[rcorr, lags] = xcorr(signalDetrended, local);

time = tStart + (1000*lags/Fs);
rcorrNorm = (abs(rcorr) - min(abs(rcorr)))/(max(abs(rcorr))-min(abs(rcorr)));

figure
plot(time, rcorrNorm)
xlabel('Time (ms)','interpreter','latex','fontsize',14)
ylabel('$|$R($\tau$)$|$','interpreter','latex','fontsize',14);
xlim([tStart time(end)])
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;


