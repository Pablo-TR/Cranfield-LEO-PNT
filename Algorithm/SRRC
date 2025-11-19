clc; clear; close all;

%% System Parameters
M = 4;              % QPSK modulation
k = log2(M);        % Bits per symbol
numSymbols = 2000;  % Number of transmitted symbols
sps = 8;            % Samples per symbol
alpha = 0.1;        % Roll-off factor
span = 6;           % Filter span (symbols)
EbNo_dB = 0;        % SNR (dB)

%% 1. Generate Random Bitstream
dataIn = randi([0 1], numSymbols * k, 1);

%% 2. QPSK Modulation
dataSym = bi2de(reshape(dataIn, k, []).').';
txSym = pskmod(dataSym, M, pi/4);  % π/4 QPSK

%% 3. SRRC Pulse Shaping Filter
txFilter = rcosdesign(alpha, span, sps, 'sqrt'); % Tx SRRC filter
txSig = upfirdn(txSym, txFilter, sps);           % Pulse shaping

%% 4. AWGN Channel
rxSig = awgn(txSig, EbNo_dB, 'measured');

%% 5. Receiver Matched Filtering
rxFilter = txFilter;                             % Rx SRRC filter
rxFiltered = conv(rxSig, rxFilter, 'same');      % Matched filtering

%% 6. Symbol Sampling
rxSamples = rxFiltered(span*sps+1 : sps : end-span*sps);

%% =========================Frequency and ISI Analysis=========================
Nfft = 2048;
[H_srrc, f] = freqz(txFilter, 1, Nfft, sps);
[H_rc, f2]  = freqz(conv(txFilter, txFilter), 1, Nfft, sps);


%% ======================Plot Results==============================


figure('Color','w','Position',[100 100 1100 800]);

% Random Bitstream
subplot(3,3,1);
stairs(dataIn(1:400), 'LineWidth', 1.5);
title('Random Binary Bitstream (dataIn)');
xlabel('Bit Index');
ylabel('Bit Value (0 or 1)');
ylim([-0.2 1.2]);
grid on;

% Transmitted SRRC waveform
subplot(3,3,2);
plot(real(txSig(1:400)), 'b+-', 'LineWidth', 1.2); hold on;
plot(imag(txSig(1:400)), 'r+-', 'LineWidth', 1.2);
title('Transmitted Baseband (SRRC Pulse Shaped)');
xlabel('Sample Index'); ylabel('Amplitude'); grid on;
legend('I','Q');

% Transmitted SRRC waveform
subplot(3,3,2);
plot(real(txSig(1:400)), 'b+-', 'LineWidth', 1.2); hold on;
plot(imag(txSig(1:400)), 'r+-', 'LineWidth', 1.2);
title('Transmitted Baseband (SRRC Pulse Shaped)');
xlabel('Sample Index'); ylabel('Amplitude'); grid on;
legend('I','Q');

% Received signal with noise
subplot(3,3,3);
plot(real(rxSig(1:400)), 'b+-', 'LineWidth', 1.2); hold on;
plot(imag(rxSig(1:400)), 'r+-', 'LineWidth', 1.2);
title('Received Signal (with AWGN)');
xlabel('Sample Index'); ylabel('Amplitude'); grid on;

% After SRRC matched filtering
subplot(3,3,4);
plot(real(rxFiltered(1:400)), 'b+-', 'LineWidth', 1.2); hold on;
plot(imag(rxFiltered(1:400)), 'r+-', 'LineWidth', 1.2);
title('Output After Matched Filtering');
xlabel('Sample Index'); ylabel('Amplitude'); grid on;

% SRRC filter impulse response
subplot(3,3,5);
plot(txFilter, 'k', 'LineWidth', 1.5);
title('SRRC Filter Impulse Response');
xlabel('Samples'); ylabel('Amplitude'); grid on;

% SRRC frequency response
subplot(3,3,6);
plot(f, 20*log10(abs(H_srrc)/max(abs(H_srrc))), 'b', 'LineWidth', 1.5);
title('Frequency Response of SRRC Filter');
xlabel('Normalized Frequency (×π rad/sample)');
ylabel('Magnitude (dB)'); grid on; ylim([-80 5]);

% RC (SRRC*SRRC) frequency response
subplot(3,3,7);
plot(f2, 20*log10(abs(H_rc)/max(abs(H_rc))), 'r', 'LineWidth', 1.5);
title('Overall Response (Raised Cosine)');
xlabel('Normalized Frequency (×π rad/sample)');
ylabel('Magnitude (dB)'); grid on; ylim([-80 5]);


%  Constellation before matched filter
subplot(3,3,8);
% use the SAME timing offset as rxSamples to make a fair comparison
preIdx  = span*sps+1 : sps : length(rxSig)-span*sps;
preSym  = rxSig(preIdx);
plot(real(preSym), imag(preSym), '.', 'MarkerSize', 8);
axis equal; grid on;
title('Constellation Before Matched Filter');
xlabel('In-phase'); ylabel('Quadrature');


% Constellation after matched filter
subplot(3,3,9);
plot(real(rxSamples), imag(rxSamples), '.', 'MarkerSize', 8);
axis equal; grid on;
title('Constellation After Matched Filter');
xlabel('In-phase'); ylabel('Quadrature');




