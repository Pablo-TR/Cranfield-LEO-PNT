%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            dopplerSearch.m
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

function [doppler, timeArr] = dopplerSearch(signalID, tInterval, tCoh)

signalData = getSignalByID(signalID);
Fs = signalData.SamplingFreq;
Fsym = signalData.SymbolFreq; % Symbol frequency [Hz]

nSamples = round(Fs*tCoh*1e-3); % Number of samples
tArray    = (0:nSamples-1).' / Fs; % Time vector (Nx1)

[~, r_local] = localReplica(Fs, Fsym);

numTimeBlocks = (tInterval(2)-tInterval(1))/tCoh;

doppler = zeros(numTimeBlocks,1);
timeArr = zeros(numTimeBlocks,1);
fd = 0;

for n = 1:numTimeBlocks
    disp(n);
    disp(numTimeBlocks);    

    tStart = (n-1)*tCoh;
    tEnd     = tStart + tCoh; % End time [ms]
    [r, ~, ~, ~] = readSignal(signalID, tStart, tEnd);
    
    r = r(:);
    fd_limit = 10000;

    % Wipe-off candidate Doppler
    r_wiped = r .* exp(-1j*2*pi*fd*tArray);
    
    while fd_limit > 100
        
        % Correlate with local replica
        rcorr = xcorr(r_wiped, r_local);
    
        % Detection metric: peak magnitude
        likelihood = max(abs(rcorr));
        
        fd_upper = fd + fd_limit;
        r_wiped_upper = r .* exp(-1j*2*pi*fd_upper*tArray);
        rcorr_upper= xcorr(r_wiped_upper, r_local);
        likelihood_upper= max(abs(rcorr_upper));

        fd_lower = fd - fd_limit;
        r_wiped_lower= r .* exp(-1j*2*pi*fd_lower*tArray);
        rcorr_lower= xcorr(r_wiped_lower, r_local);
        likelihood_lower= max(abs(rcorr_lower));

        if likelihood_lower < likelihood && likelihood_upper < likelihood
            fd_limit = fd_limit/2;
            continue
        elseif likelihood_lower > likelihood
            fd = fd_lower;
            r_wiped = r .* exp(-1j*2*pi*fd*tArray);
            continue
        elseif likelihood_upper > likelihood
            fd = fd_upper;
            r_wiped = r .* exp(-1j*2*pi*fd*tArray);
        end
    end

    doppler(n) = fd;
    timeArr(n) = tStart;

end


return