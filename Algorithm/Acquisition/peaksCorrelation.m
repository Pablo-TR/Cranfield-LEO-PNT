%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            peaksCorrelation.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Acquisition
% Author: Alfredo
% Date: 7th March 2026
% Last Update: 7th March 2026
% Version: 1.0.0
%
% Description: Outputs the product of the correlation between the signal
% and the local replica against time. Peaks should be seen every 1 ms.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

%%
tStart = 0;
tEnd = 5;

signalID = 0;
signalData = getSignalDataByID(signalID);
Fs = signalData.SamplingFreq;
Fsym = 2.304e8;

[r, ~, ~, ~] = readSignal(signalID, tStart, tEnd);
% r = signalGenerator(Fs, 5, 0, 1);
[~, local] = localReplica(Fs, Fsym);
% local = localExtendedReplica(Fs, tEnd-tStart);


% Correlation between incoming signal and synchronisation sequence
[rcorr, lags] = xcorr(r, local);

mu = mean(rcorr);
sigma = std(rcorr);

time = tStart + (1000*lags/Fs);
rcorrNorm = (abs(rcorr) - min(abs(rcorr)))/(max(abs(rcorr))-min(abs(rcorr)));

rcorrVar = (rcorr - mu)/sigma;
%
figure
plot(time, abs(rcorr))
xlabel('Time (ms)','interpreter','latex','fontsize',18)
ylabel('Likelihood','interpreter','latex','fontsize',18);
title('Correlation')
xlim([tStart time(end)])
%ylim([0 16])
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;


%%
% Peak value
peakVal = max(abs(rcorr));

% Exclude region around the peak to estimate noise
[~, peakIdx] = max(abs(rcorr));
guard = round(0.001 * Fs); % e.g. 1 ms guard (adjust if needed)

noiseIdx = true(size(rcorr));
noiseIdx(max(1,peakIdx-guard):min(length(rcorr),peakIdx+guard)) = false;

% Noise RMS (recommended): Interference floor
noiseRMS = sqrt(mean(abs(rcorr(noiseIdx)).^2));

% SNR in dB (amplitude ratio)
SNR_dB = 20 * log10(peakVal / noiseRMS)
%%
figure
plot(real(r(peakIdx-500000:peakIdx-500000+600)), imag(r(peakIdx-500000:peakIdx-500000+600)), '.')

replicaFound = r(peakIdx-500000 : 1 : peakIdx-500000+600);
txFilter = rcosdesign(0.1, 40, 30, 'sqrt');
txFilter = resample(txFilter, 625, 288);
txSig = upfirdn(replicaFound, txFilter, 30);           % Pulse shaping

rxFiltered = conv(txSig, txFilter, 'same');

rxSamples = rxFiltered(40*30+1 : 30 : end-40*30);

figure
plot(real(rxSamples), imag(rxSamples), '.')




%%
tStart = 0;
tEnd = 1;

signalID = 0;
signalData = getSignalDataByID(signalID);
Fs = signalData.SamplingFreq;
Fsym = 2.304e8;

[signalDetrended, ~, ~, ~] = readSignal(signalID, tStart, tEnd);

[~, local] = localReplica(Fs, Fsym);
signalDetrended = signalDetrended(:);
local = local(:);

% Correlation between incoming signal and synchronisation sequence
% Correlation via convolution
rcorr = conv(signalDetrended, conj(flipud(local)), 'same');

% For 'valid', output index k corresponds to local aligned with
lags = 0:length(rcorr)-1;

time = tStart + (1000*lags/Fs);
% Find index of maximum peak
[peakValue, peakIdx] = max(abs(rcorr));

% Corresponding lag and time
peakLag = lags(peakIdx);
peakTime = time(peakIdx);

rcorrNorm = (abs(rcorr) - min(abs(rcorr)))/(max(abs(rcorr))-min(abs(rcorr)));

figure
plot(time, abs(rcorr))
xlabel('Time (ms)','interpreter','latex','fontsize',14)
ylabel('$|$R($\tau$)$|$','interpreter','latex','fontsize',14);
%ylim([0 18])
xlim([tStart time(end)])
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;

%%
preIdx  = 4*1+1 : 1 : length(signalDetrended)-4*1;
preSym  = signalDetrended(preIdx);
figure
plot(real(preSym(203579:203779)), imag(preSym(203579:203779)), 'x')
grid on

rxSamples = rcorr(4*1+1 : 1 : end-4*1);
figure
plot(real(rxSamples(203579:203779)), imag(rxSamples(203579:203779)), 'x')
grid on
