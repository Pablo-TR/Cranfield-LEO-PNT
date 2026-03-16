%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            dopplerCorrelation.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Algorithm
% Author: Alfredo
% Date: 8th March 2026
% Last Update: 8th March 2026
% Version: 1.0.0
%
% Description: Time-domain Doppler estimation via Doppler wipe-off and correlation.
% A grid of candidate Doppler frequencies is applied to the received
% signal to wipe off the carrier offset. Each wiped signal is
% cross-correlated with a local replica, and the maximum correlation
% peak determines the Doppler estimate.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fd, correlation] = dopplerCorrelation(r, signalData, tCoh)

    Fs = signalData.SamplingFreq;
    Fsym = signalData.SymbolFreq; % Symbol frequency [Hz]
    
    tArray = (0:length(r)-1).'/Fs; % Time vector (Nx1)
    
    if tCoh <= 1
        [~, r_local] = localReplica(Fs, Fsym);
    else
        r_local = localExtendedReplica(Fs, tCoh);
    end
    

    r = r(:);
    r_local = r_local(:);
    
    fd = (-1000:100:1000)'; % Doppler frequency
    
    
    correlation = zeros(length(fd),1);
    argmax_idx = zeros(length(fd),1);
    
    for i=1:length(fd)
        disp(fd(i))
        r_wiped = r .* exp(-1j*2*pi*fd(i)*tArray);
        % Correlate with local replica
        rcorr = xcorr(r_wiped, r_local);
        %rcorr = conv(r_wiped, conj(flipud(r_local)), "valid");
        
        % Detection metric: peak magnitude
        [correlation(i), argmax_idx(i)] = max(abs(rcorr));
    end
