%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            syncSequence.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: -
% Author: Pablo, Fredo, Marti & Dan
% Date: 21st November 2025
% Last Update: 26th November 2025
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
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
tStart   = 0; % Start time of signal to read [ms].
tEnd = 20; % End time of signal to read [ms].

M = 4;              % QPSK modulation
k = log2(M);        % Bits per symbol
Fsym = 2.304e8;     % Symbol rate (Hz)
Fs = 3.2e8;         % Sampling frequency (Hz)
span = 15;          % Filter span (symbols), number of symbols covered
sps = 18;           % Samples per symbol
rollOff = 0.1;      % Obtained from paper

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
dataSym = fliplr(dataSym); % Flips entire sequence
txSym = pskmod(dataSym, M, 0);  % QPSK, converts integer symbol values into complex points in the in-phase quadrature diagram

Nsym = length(txSym); % Number of symbols

t_in  = (0:Nsym-1)/Fsym; % Time array containing the time at which each symbol has been sent, duration of SS = 1.73e-06s
t_out = 0:1/Fs:t_in(end); % Time axis of samples at 320 MHz

% Upsampling, so that symbols frequency matches the income signal
% frequency, from Fsym to Fs
txSymResampled = resample(txSym, Fs, Fsym);

[signalDetrended, ~, timeArray, binData] = readSignal(filename, Fs, tStart, tEnd);

% Correlation between incoming signal and synchronisation sequence
[rcorr, lags] = xcorr(signalDetrended, txSymResampled);

time = tStart + (1000*lags/Fs);

%% 3. SRRC Pulse Shaping Filter, maximise correlation by varying span and sps?

txFilter = rcosdesign(rollOff, span, sps, 'sqrt'); % Tx SRRC filter, span*L_up+1 samples of the filter
txSig = upfirdn(txSym, txFilter, 25, 18);  % Pulse shaping

[rcorr2, lags2] = xcorr(signalDetrended, txSig);
r2 = abs(rcorr2) / (norm(signalDetrended)*norm(txSig));
time2 = tStart + (1000*lags2/Fs);

%% Shape of SRRC Pulse
figure
plot(1:length(txFilter), txFilter)
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;

%%
figure
plot(real(txSig), 'b+-', 'LineWidth', 1.2); hold on;
plot(imag(txSig), 'r+-', 'LineWidth', 1.2);
title('Transmitted Baseband (SRRC Pulse Shaped)');
xlabel('Sample Index'); ylabel('Amplitude'); grid on;
legend('I','Q');

%%
figure
plot(time, abs(rcorr))
xlabel('Time (ms)','interpreter','latex','fontsize',14)
ylabel('$|$R($\tau$)$|$','interpreter','latex','fontsize',14);
xlim([tStart time(end)])
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;

%%
figure
plot(real(txSymResampled), imag(txSymResampled), '.', 'MarkerSize', 8);
axis equal; grid on;
title('Constellation Before Matched Filter');
xlabel('In-phase'); ylabel('Quadrature');

%%
fSig = conv(signalDetrended, txFilter, 'same');
figure
plot(real(fSig), imag(fSig), '.', 'MarkerSize', 8);
axis equal; grid on;
title('Constellation Before Matched Filter');
xlabel('In-phase'); ylabel('Quadrature');

%% Circular correlation (In progress)
N = length(signalDetrended);
X = fft(signalDetrended);

r_circ = ifft( X .* conj(X) );

figure;
plot(timeArray*1000, abs(r_circ));
xlabel('Time (ms)', 'interpreter','latex','fontsize',14);
ylabel('|R_{xx}(k)|', 'interpreter','latex','fontsize',14);
title('Circular Autocorrelation Magnitude');
grid on;


