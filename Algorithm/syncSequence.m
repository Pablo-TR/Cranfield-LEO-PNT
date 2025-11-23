%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            syncSequence.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: -
% Author: Pablo, Fredo, Marti & Dan
% Date: 21st November 2025
% Last Update: 21st November 2025
% Version: 1.0.0
%
% Description: Recreating the synchronisation sequence from the OneWeb signal
% seen in (Komodromos, Z; Humphreys, T.E).
% Komodromos, Z. M., & Humphreys, T. E. (2025). Signal parameter estimation 
% and demodulation of the OneWeb Ku-band downlink [Preprint]. 
% The University of Texas at Austin. 
% https://radionavlab.ae.utexas.edu/wp-content/uploads/komodromos_oneweb_parameter_estimation.pdf
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% System parameters
filename = 'OW-sample500-250-5s-ram-29.bin';
t   = 0.01; % Interval of time of signal to read [ms]. Reads from t = 0ms to t 

M = 4;              % QPSK modulation
k = log2(M);        % Bits per symbol
Fsym = 2.304e8;     % Symbol rate (Hz)
Fs = 3.2e8;         % Sampling frequency (Hz)
L_up = 25;          % Samples per symbol
M_dn = 18;
span = 20;          % Filter span (symbols), number of symbols covered
rollOff = 0.1;      % Obtained from paper

% Doppler search
fD_min   = -200e3;     % Min Doppler [Hz]
fD_max   =  200e3;     % Max Doppler [Hz]
fD_step  =  2e3;       % Doppler step [Hz]
fD_vec   = fD_min:fD_step:fD_max;
%% 1. Hexadecimal sequence of symbols (SS) conversion to Binary
qSSHex = ['B5D0 CDB5 66F9 5A93 F90B 0060 834E 073C 9EC3 EAAA D425 C677' ...
    '93B0 EE1F 993C 5CF5 2FFE 5839 CC7E 5170 FE09 31EF 33CD 3E13 16F4 ' ...
    '3E9E 2A17 5D4B 2D9B E629 2E62 6386 B994 6849 7811 5074 5930 417E ' ...
    '3338 E497 3A3A 5B05 CFBD 5A8F 669D 9D31 EEB8 B48C B7E2 2DBA'];

qSSHex = erase(qSSHex, " "); % Remove spaces
qSSBin = "";
for hex = 1:length(qSSHex)
    binByte = dec2bin(hex2dec(qSSHex(hex)), 4);
    qSSBin = qSSBin + binByte;
end

qSSBin = char(qSSBin) - '0';

%% 2. QPSK Modulation (bi2de is binary to decimal conv)
dataSym = bi2de(reshape(qSSBin, k, []).', 'left-msb').';   % a_m symbol phase value 0, 1, 2 or 3 (left-msb treats left bit as highest)
txSym = pskmod(dataSym, M, 0);  % π/4 QPSK, converts integer symbol values into complex points in the in-phase quadrature diagram

%% 3. SRRC Pulse Shaping Filter

txFilter = rcosdesign(rollOff, span, L_up, 'sqrt'); % Tx SRRC filter, span*L_up+1 samples of the filter
dataSymReSampled = upsample(dataSym, L_up);
txSig = upfirdn(txSym, txFilter, L_up);  % Pulse shaping

figure
plot(real(txSig(1:2000)), 'b+-', 'LineWidth', 1.2); hold on;
plot(imag(txSig(1:2000)), 'r+-', 'LineWidth', 1.2);
title('Transmitted Baseband (SRRC Pulse Shaped)');
xlabel('Sample Index'); ylabel('Amplitude'); grid on;
legend('I','Q');


groupDelay = (span * L_up) / 2;   % group delay in samples at 25x rate

% Trim the transient (centre-align the pulse train)
txSigAligned = txSig(groupDelay+1 : end-groupDelay);
txSigSam = txSigAligned(1:M_dn:end);    % take every 18th sample

txSigSam = txSigSam / max(txSigSam);

[signalDetrended, N, timeArray, binData] = readSignal(filename, Fs, t); % Returns DETRENDED signal






