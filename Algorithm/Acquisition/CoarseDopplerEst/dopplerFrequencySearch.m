%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            dopplerFrequencySearch.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Acquisition
% Author: Alfredo
% Date: 12th March 2026
% Last Update: 12th March 2026
% Version: 1.0.0
%
% Description: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [result, tArray] = dopplerFrequencySearch(r, signalData, tCoh, fDoppler)
    
    Fs = signalData.SamplingFreq;
    Fsym = signalData.SymbolFreq;

    r = r(:);

    if tCoh <= 1
        [~, r_local] = localReplica(Fs, Fsym);
    else
        r_local = localExtendedReplica(Fs, tCoh);
    end
    
    r_local = r_local(:);
    
    tArray = (0:length(r)-1)'/Fs;
    
    Nx = numel(r);
    Nr = numel(r_local);
    Nfft = Nr + Nx - 1;
    selectInd = Nr : Nfft;
    
    result = zeros(length(fDoppler), length(selectInd));
    
    fftLocal = conj(fft(r_local, Nfft));
    
    for i = 1:length(fDoppler)
        disp(i)
        disp(length(fDoppler))
    
        % Wipe-off candidate Doppler
        r_wiped = r .* exp(-1j*2*pi*fDoppler(i)*tArray);
    
        % Compute the FFT of the wiped signal
        fftWiped = fft(r_wiped, Nfft);

        inverse = ifft(fftWiped.*fftLocal, Nfft);

        % rcorr = xcorr(r_wiped, r_local);
    
        % Store the result of the inverse FFT
        result(i, :) = abs(inverse(selectInd));

    end

