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
span = 2:1:100;          % Filter span (symbols), number of symbols covered
sps = 2:1:100;           % Samples per symbol
maxCorr = zeros(length(span), length(sps));
rollOff = 0.05:0.05:0.5;

% 1. Hexadecimal sequence of symbols (SS) conversion to Binary
qSSHex = ['D03A B3D0 44F1 051E F1AD AA4A 9E87 ACEB 17BE 7555 3860 B4CC' ...
        '1EDA 772F 11EB 0BF0 6FF7 09E1 BBC7 02CA F7A1 E27F EEB3 E72E 24F8 ' ...
        'E717 652C 038D 631D 7461 6746 4E94 D118 4981 C922 0AC8 01EA 82C7 ' ...
        'EEE9 781C E5E5 0DA0 BFD3 059F 4413 13E2 77D9 D89B DC76 63D5'];

qSSHex = erase(qSSHex, " "); % Remove spaces
qSSBin = "";
for hex = 1:length(qSSHex)
    binByte = dec2bin(hex2dec(qSSHex(hex)), 4);
    qSSBin = qSSBin + binByte;
end
    
qSSBin = char(qSSBin) - '0';

% 2. QPSK Modulation (bi2de is binary to decimal conv)
dataSym = fliplr(bi2de(reshape(qSSBin, k, []).', 'left-msb').');   % a_m symbol phase value 0, 1, 2 or 3 (left-msb treats left bit as highest) and flips entire sequence
rawLocalReplica = pskmod(dataSym, M, 0);  % QPSK, converts integer symbol values into complex points in the in-phase quadrature diagram
    
[p_num, q_den] = rat(Fs/Fsym);
rawLocalReplica = resample(rawLocalReplica, p_num, q_den);

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
            localReplica = upfirdn(rawLocalReplica, txFilter, 1, 1);

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
