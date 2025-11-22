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
M = 4;              % QPSK modulation
k = log2(M);        % Bits per symbol

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



