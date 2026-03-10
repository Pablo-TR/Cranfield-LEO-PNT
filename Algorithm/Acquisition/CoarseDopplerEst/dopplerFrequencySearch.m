function [result, tArray, fDoppler] = dopplerFrequencySearch(signalID, tCoh)

Fsym = 2.304e8; % Symbol frequency [Hz]

signalData = getSignalByID(signalID);
Fs = signalData.SamplingFreq;

[r_local, ~] = localReplica(Fs, Fsym);
r_local = r_local(:);

fDoppler = -15000:100:15000; % Doppler frequency

tStart = 0;
tEnd     = tStart + tCoh; % End time [ms]
[r, ~, ~, ~] = readSignal(signalID, tStart, tEnd);
r = r(:);

tArray = (0:length(r)-1)'/Fs;

Nx = numel(r);
Nr = numel(r_local);
Nfft = Nr + Nx - 1;
selectInd = floor(Nr/2) : Nfft - ceil(Nr/2);

result = zeros(length(fDoppler),length(selectInd));

fftLocal = conj(fft(r_local, Nfft));

for i = 1:length(fDoppler)    

    % Wipe-off candidate Doppler
    r_wiped = r .* exp(-1j*2*pi*fDoppler(i)*tArray);

    % Compute the FFT of the wiped signal
    fftWiped = fft(r_wiped, Nfft);

    inverse = ifft(fftWiped.*fftLocal, Nfft);

    % Store the result of the inverse FFT
    result(i, :) = abs(inverse(selectInd));

end

