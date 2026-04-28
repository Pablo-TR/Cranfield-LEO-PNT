%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            localReplica.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Acquisition
% Author: Alfredo
% Date: 7th March 2026
% Last Update: 7th March 2026
% Version: 1.1.0
%
% Description: Creation of local replica in desired sampling frequency (Fs)
%
% Komodromos, Z. M., & Humphreys, T. E. (2025). Signal parameter estimation 
% and demodulation of the OneWeb Ku-band downlink [Preprint]. 
% The University of Texas at Austin. 
% https://radionavlab.ae.utexas.edu/wp-content/uploads/komodromos_oneweb_parameter_estimation.pdf
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rawLocalReplica, localReplica] = localReplica(Fs, Fsym)

    M = 4;              % QPSK modulation
    k = log2(M);        % Bits per symbol
    span = 10;          % Filter span (symbols), number of symbols covered
    sps = 5;           % Samples per symbol
    rollOff = 0.1;      % Obtained from paper

    [p_num, q_den] = rat(Fs/Fsym);

    % 1. Hexadecimal sequence of symbols (SS) conversion to Binary
    qSSHex = ['B5D0 CDB5 66F9 5A93 F90B 0060 834E 073C 9EC3 EAAA D425 C677' ...
        '93B0 EE1F 993C 5CF5 2FFE 5839 CC7E 5170 FE09 31EF 33CD 3E13 16F4 ' ...
        '3E9E 2A17 5D4B 2D9B E629 2E62 6386 B994 6849 7811 5074 5930 417E ' ...
        '3338 E497 3A3A 5B05 CFBD 5A8F 669D 9D31 EEB8 B48C B7E2 2DBA'];

    % qSSHex = ['D03A B3D0 44F1 051E F1AD AA4A 9E87 ACEB 17BE 7555 3860 B4CC' ...
    %     '1EDA 772F 11EB 0BF0 6FF7 09E1 BBC7 02CA F7A1 E27F EEB3 E72E 24F8 ' ...
    %     'E717 652C 038D 631D 7461 6746 4E94 D118 4981 C922 0AC8 01EA 82C7 ' ...
    %     'EEE9 781C E5E5 0DA0 BFD3 059F 4413 13E2 77D9 D89B DC76 63D5'];
 

    qSSHex = erase(qSSHex, " "); % Remove spaces
    qSSBin = "";
    for hex = 1:length(qSSHex)
        binByte = dec2bin(hex2dec(qSSHex(hex)), 4);
        qSSBin = qSSBin + binByte;
    end
    
    qSSBin = char(qSSBin) - '0';

    % 2. QPSK Modulation (bi2de is binary to decimal conv)
    dataSym = fliplr(bi2de(reshape(qSSBin, k, []).', 'left-msb').');   % a_m symbol phase value 0, 1, 2 or 3 (left-msb treats left bit as highest) and flips entire sequence
    rawLocalReplica = exp(1i.*dataSym*pi/2); % QPSK, converts integer symbol values into complex points in the in-phase quadrature diagram
    
    
    
    % localReplica = resample(rawLocalReplica, p_num, q_den);
    % figure
    % plot(real(rawLocalReplica))
    % xlabel('Sample Index')
    % ylabel('Amplitude')


    % 3. SRRC Pulse Shaping Filter, maximise correlation by varying span and sps?
    txFilter = rcosdesign(rollOff, span, sps, 'sqrt'); % Tx SRRC filter, span*L_up+1 samples of the filter
    localReplica = upfirdn(rawLocalReplica, txFilter, sps, 1);  % Pulse shaping

    figure
    plot(txFilter, 'k', 'LineWidth', 1.5);
    title('SRRC Filter Impulse Response');
    xlabel('Samples'); ylabel('Amplitude'); grid on;
    
    
    % figure
    % plot(real(rawLocalReplica))

    localReplica = resample(localReplica, p_num, q_den);

    localReplica = localReplica(span*sps+1 : sps : end-span*sps);

    rawLocalReplica = resample(rawLocalReplica, p_num, q_den);
    % figure
    % plot(real(localReplica))

    % localReplica = resample(localReplica(floor(0.5*(span*sps))+1:end-ceil(0.5*(span*sps))),p_num, q_den);
    
    
end
