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

function [fd, correlation] = dopplerCorrelation(signalID, tStart, tCoh)

    signalData = getSignalByID(signalID);
    Fs = signalData.SamplingFreq;
    Fsym = signalData.SymbolFreq; % Symbol frequency [Hz]
    
    nSamples = round(Fs*tCoh*1e-3); % Number of samples
    tArray    = (0:nSamples-1).' / Fs; % Time vector (Nx1)
    
    [~, r_local] = localReplica(Fs, Fsym);
    %r_local = localExtendedReplica(Fs, tCoh);
    
    [r, ~, ~, ~] = readSignal(signalID, tStart, tStart+tCoh);
    r = r(:);
    
    fd = -15000:1000:15000; % Doppler frequency
    
    
    correlation = zeros(length(fd),1);
    
    for i=1:length(fd)
        disp(i)
        r_wiped = r .* exp(-1j*2*pi*fd(i)*tArray);
        % Correlate with local replica
        rcorr = xcorr(r_wiped, r_local);
        
        % Detection metric: peak magnitude
        correlation(i) = max(abs(rcorr));
    end
