%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            pulseShapingOptimisation.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Algorithm
% Author: Alfredo
% Date: 8th March 2026
% Last Update: 8th March 2026
% Version: 1.0.0
%
% Description: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all

%% Pulse-shaping optimisation

tStart = 0;
tEnd = 1;

signalID = 0;
signalData = getSignalDataByID(signalID);
Fs = signalData.SamplingFreq;
Fsym = signalData.SymbolFreq;
[signalDetrended, ~, timeArray, binData] = readSignal(signalID, tStart, tEnd);

M = 4;              % QPSK modulation
k = log2(M);        % Bits per symbol
span = 1:1:40;          % Filter span (symbols), number of symbols covered
sps = 1:1:40;           % Samples per symbol
maxCorr = zeros(length(span), length(sps));
rollOff = 0.05:0.05:0.1;

% 1. Hexadecimal sequence of symbols (SS) conversion to Binary
qSSHex = ['B5D0 CDB5 66F9 5A93 F90B 0060 834E 073C 9EC3 EAAA D425 C677' ...
        '93B0 EE1F 993C 5CF5 2FFE 5839 CC7E 5170 FE09 31EF 33CD 3E13 16F4 ' ...
        '3E9E 2A17 5D4B 2D9B E629 2E62 6386 B994 6849 7811 5074 5930 417E ' ...
        '3338 E497 3A3A 5B05 CFBD 5A8F 669D 9D31 EEB8 B48C B7E2 2DBA'];

% qSSHex = ['D03A B3D0 44F1 051E F1AD AA4A 9E87 ACEB 17BE 7555 3860 B4CC' ...
%       '1EDA 772F 11EB 0BF0 6FF7 09E1 BBC7 02CA F7A1 E27F EEB3 E72E 24F8 ' ...
%       'E717 652C 038D 631D 7461 6746 4E94 D118 4981 C922 0AC8 01EA 82C7 ' ...
%       'EEE9 781C E5E5 0DA0 BFD3 059F 4413 13E2 77D9 D89B DC76 63D5'];

qSSHex = erase(qSSHex, " "); % Remove spaces
qSSBin = "";
for hex = 1:length(qSSHex)
    binByte = dec2bin(hex2dec(qSSHex(hex)), 4);
    qSSBin = qSSBin + binByte;
end
    
qSSBin = char(qSSBin) - '0';

% 2. QPSK Modulation (bi2de is binary to decimal conv)
dataSym = fliplr(bi2de(reshape(qSSBin, k, []).', 'left-msb').');   % a_m symbol phase value 0, 1, 2 or 3 (left-msb treats left bit as highest) and flips entire sequence
rawLocalReplica = exp(1i.*dataSym*pi/2);  % QPSK, converts integer symbol values into complex points in the in-phase quadrature diagram

% 3. SRRC Pulse Shaping Filter, maximise correlation by varying span and sps?

for r = 1:length(rollOff)
    disp(['rollOff index = ', num2str(r), ', value = ', num2str(rollOff(r))])

    for i = 1:length(span)
        disp(['span index = ', num2str(i), ', value = ', num2str(span(i))])

        for j = 1:length(sps)
            disp(['sps index = ', num2str(j), ', value = ', num2str(sps(j))])

            % rcosdesign requires span*sps to be even
            if mod(span(i) * sps(j), 2) ~= 0
                continue
            end

            txFilter = rcosdesign(rollOff(r), span(i), sps(j), 'sqrt');
            localReplica = upfirdn(rawLocalReplica, txFilter, sps(j), 1);

            [p_num, q_den] = rat(Fs/Fsym);
            localReplica = resample(localReplica, p_num, q_den);
            gd = span(i)*sps(j)/2;

            localReplica = localReplica(gd+1 : sps(j) : end-gd);

            % Correlation between incoming signal and synchronisation sequence
            [rcorr, lags] = xcorr(signalDetrended, localReplica);

            maxCorr(i, j, r) = max(abs(rcorr));
        end
    end
end

% Best value for each span across sps and rollOff
maxmax = zeros(length(span), 1);
for i = 1:length(span)
    maxmax(i) = max(max(maxCorr(i, :, :)));
end

% Global optimum
[bestCorr, linearIdx] = max(maxCorr(:));
[iBest, jBest, rBest] = ind2sub(size(maxCorr), linearIdx);

bestSpan = span(iBest);
bestSps = sps(jBest);
bestRollOff = rollOff(rBest);

disp('Best parameters found:')
disp(['Max correlation = ', num2str(bestCorr)])
disp(['Best span       = ', num2str(bestSpan)])
disp(['Best sps        = ', num2str(bestSps)])
disp(['Best rollOff    = ', num2str(bestRollOff)])
