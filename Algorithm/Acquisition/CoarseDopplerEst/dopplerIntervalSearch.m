%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            dopplerIntervalSearch.m
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

function doppler = dopplerIntervalSearch(r, signalData, tCoh, fd)

r = r(:);

Fs = signalData.SamplingFreq;
Fsym = signalData.SymbolFreq; % Symbol frequency [Hz]

nSamples = round(Fs*tCoh*1e-3); % Number of samples
tArray    = (0:nSamples-1).' / Fs; % Time vector (Nx1)
isNotFound = 1;


if tCoh <= 1
        [~, r_local] = localReplica(Fs, Fsym);
    else
        r_local = localExtendedReplica(Fs, tCoh);
end

 

Nx = numel(r);
Nr = numel(r_local);
Nfft = Nr + Nx - 1;
selectInd = Nr : Nfft;

fftLocal = conj(fft(r_local, Nfft));

% Wipe-off candidate Doppler
r_wiped = r .* exp(-1j*2*pi*fd*tArray);
    
while isNotFound
        
    % Correlate with local replica
    fftWiped = fft(r_wiped, Nfft);
    inverse = ifft(fftWiped.*fftLocal, Nfft);
    likelihood = max(abs(inverse(selectInd)));
        
    fd_upper = fd + 1000;
    r_wiped_upper = r .* exp(-1j*2*pi*fd_upper*tArray);
    fftWiped_upper = fft(r_wiped_upper, Nfft);    
    inverse_upper = ifft(fftWiped_upper.*fftLocal, Nfft);
    likelihood_upper = max(abs(inverse_upper(selectInd)));

    fd_lower = fd - 1000;
    r_wiped_lower= r .* exp(-1j*2*pi*fd_lower*tArray);
    fftWiped_lower = fft(r_wiped_lower, Nfft);    
    inverse_lower = ifft(fftWiped_lower.*fftLocal, Nfft);
    likelihood_lower = max(abs(inverse_lower(selectInd)));

    if likelihood_lower < likelihood && likelihood_upper < likelihood
        isNotFound = 0;
    elseif likelihood_lower > likelihood
        fd = fd_lower;
        r_wiped = r .* exp(-1j*2*pi*fd*tArray);
        continue
    elseif likelihood_upper > likelihood
        fd = fd_upper;
        r_wiped = r .* exp(-1j*2*pi*fd*tArray);
    end
end

isNotNotFound = 1;

while isNotNotFound
        
    % Correlate with local replica
    fftWiped = fft(r_wiped, Nfft);
    inverse = ifft(fftWiped.*fftLocal, Nfft);
    likelihood = max(abs(inverse(selectInd)));
        
    fd_upper = fd + 2;
    r_wiped_upper = r .* exp(-1j*2*pi*fd_upper*tArray);
    fftWiped_upper = fft(r_wiped_upper, Nfft);    
    inverse_upper = ifft(fftWiped_upper.*fftLocal, Nfft);
    likelihood_upper = max(abs(inverse_upper(selectInd)));

    fd_lower = fd - 2;
    r_wiped_lower= r .* exp(-1j*2*pi*fd_lower*tArray);
    fftWiped_lower = fft(r_wiped_lower, Nfft);    
    inverse_lower = ifft(fftWiped_lower.*fftLocal, Nfft);
    likelihood_lower = max(abs(inverse_lower(selectInd)));

    if likelihood_lower < likelihood && likelihood_upper < likelihood
        isNotNotFound = 0;
    elseif likelihood_lower > likelihood
        fd = fd_lower;
        r_wiped = r .* exp(-1j*2*pi*fd*tArray);
        continue
    elseif likelihood_upper > likelihood
        fd = fd_upper;
        r_wiped = r .* exp(-1j*2*pi*fd*tArray);
    end
end


doppler = fd;



return