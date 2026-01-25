%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            mainSignalValidation.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: -
% Author: Pablo, Fredo, Marti & Dan
% Date: 9th November 2025
% Last Update: 9th November 2025
% Version: 1.0.0
%
% Description: Compares and validates the given OneWeb signal to the
% results in (Komodromos, Z; Humphreys, T.E)
%
% Komodromos, Z. M., & Humphreys, T. E. (2025). Signal parameter estimation 
% and demodulation of the OneWeb Ku-band downlink [Preprint]. 
% The University of Texas at Austin. 
% https://radionavlab.ae.utexas.edu/wp-content/uploads/komodromos_oneweb_parameter_estimation.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear; close all;  % Reset MATLAB environment

%% Read Signal
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
Fs = 3.20e8; % Sampling frequency of the signal [MHz].
tStart = 400; % Start time of signal to read [ms].
tEnd = 420; % Stop time [ms]
[signalDetrended, N, timeArray] = readSignal(filename, Fs, tStart, tEnd); % Returns DETRENDED signal

%% FFT Calculation
signalFFT = fft(signalDetrended);
fShifted = (-N/2 : N/2-1) * (Fs / N) * 1e-6; % in MHz
fprintf('FFT calculation complete.\n');
plotFFT(fShifted, abs(signalFFT))

%% PSD calculation
fprintf('Calculating PSD...\n');

% Small window segment size.
windowLength = floor(0.025*N);
window = hamming(windowLength);
nOverlap = windowLength / 2; % 50% overlap between windows
nfft = windowLength; % Same number of FFT points as the window

% This function will now average (N / win_len) * 2 segments.
[pxx, fpsd] = pwelch(signalDetrended, window, nOverlap, nfft, Fs, 'centered');
w = abs(fpsd(1)-fpsd(2));
fpsd = fpsd * 1e-6; % Convert frequency to MHz
pxx = pxx.*w; % Remove normalisation with frequency.
fprintf('PSD (pwelch) calculation complete.\n');
plotPSD(fpsd, pxx);

%% Cyclic Autocorrelation
fprintf('Calculating cyclic autocorrelation...\n');

% Use only a percentage of the processed signal.
tChunk = 0.1*t; % 10% of the processed signal 
lChunk = floor(tChunk * Fs / 1000);
signalChunk = signalDetrended(1:lChunk); 

% The frequency range to study (based on paper)
alphaHat = (230.39: 0.0001: 230.41) .* 1e6;

% Autocorrelation
[R, alpha] = cyclicAutoCorr(signalChunk, 0, Fs, alphaHat);

% Normalise
Rnormalized = abs(R) / max(abs(R));

plotAutocorrelation(alphaHat, Rnormalized)


%% Spectogram
fprintf('Calculating data for 3D spectrogram...\n');
winSpec = hamming(4096);
overlapSpec = 2048;
nfftSpec = 4096;
[sSpec, fSpec, tSpec] = spectrogram(signalDetrended, winSpec, overlapSpec, nfftSpec, Fs, 'centered');
tSpec = tSpec * 1000; % tspec to ms
fSpec = fSpec * 1e-6; % fspec to MHz
sSpec = 10*log10(abs(sSpec)); % Output of spectogram in db
plotSpectogram(tSpec,fSpec,sSpec);