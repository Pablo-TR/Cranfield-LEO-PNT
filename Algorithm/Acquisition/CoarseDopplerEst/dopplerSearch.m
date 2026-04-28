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

function [doppler, timeArr] = dopplerSearch(r, signalData, tInterval, tCoh, fd, fdRange)

r = r(:);

Fs = signalData.SamplingFreq;
Fsym = signalData.SymbolFreq; % Symbol frequency [Hz]

nSamples = round(Fs*tCoh*1e-3); % Number of samples
tArray    = (0:nSamples-1).' / Fs; % Time vector (Nx1)

tStart = tInterval(1);

if tCoh <= 1
        [~, r_local] = localReplica(Fs, Fsym);
    else
        r_local = localExtendedReplica(Fs, tCoh);
end

numTimeBlocks = (tInterval(2)-tInterval(1))/tCoh;

doppler = zeros(numTimeBlocks,1);
timeArr = zeros(numTimeBlocks,1);

for n = 1:numTimeBlocks
    disp(n);
    disp(numTimeBlocks);    

    Nx = numel(r);
    Nr = numel(r_local);
    Nfft = Nr + Nx - 1;
    selectInd = floor(Nr/2) : Nfft - ceil(Nr/2);

    fftLocal = conj(fft(r_local, Nfft));

    % Wipe-off candidate Doppler
    r_wiped = r .* exp(-1j*2*pi*fd*tArray);
    
    while fdRange > 0.1
        disp(fdRange)
        
        % Correlate with local replica
        fftWiped = fft(r_wiped, Nfft);
        inverse = ifft(fftWiped.*fftLocal, Nfft);
        likelihood = max(abs(inverse(selectInd)));
        
        fd_upper = fd + fdRange;
        r_wiped_upper = r .* exp(-1j*2*pi*fd_upper*tArray);
        fftWiped_upper = fft(r_wiped_upper, Nfft);    
        inverse_upper = ifft(fftWiped_upper.*fftLocal, Nfft);
        likelihood_upper = max(abs(inverse_upper(selectInd)));

        fd_lower = fd - fdRange;
        r_wiped_lower= r .* exp(-1j*2*pi*fd_lower*tArray);
        fftWiped_lower = fft(r_wiped_lower, Nfft);    
        inverse_lower = ifft(fftWiped_lower.*fftLocal, Nfft);
        likelihood_lower = max(abs(inverse_lower(selectInd)));

        if likelihood_lower < likelihood && likelihood_upper < likelihood
            fdRange = fdRange/2;
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