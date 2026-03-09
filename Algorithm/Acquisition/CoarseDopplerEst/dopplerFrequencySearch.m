function [result, tArray, fDoppler] = dopplerFrequencySearch(signalID, tCoh)

Fsym = 2.304e8; % Symbol frequency [Hz]

signalData = getSignalByID(signalID);
Fs = signalData.SamplingFreq;

nSamples = round(Fs*tCoh*1e-3); % Number of samples
tArray   = (0:nSamples-1).' / Fs; % Time vector (Nx1)


[r_local, ~] = localReplica(Fs, Fsym);

fDoppler = -50000:1000:50000; % Doppler frequency

tStart = 0;
tEnd     = tStart + tCoh; % End time [ms]
[r, ~, ~, ~] = readSignal(signalID, tStart, tEnd);
r = r(:);

Nx = numel(r);
Nr = numel(r_local);
Nfft = Nr + Nx - 1;
selectInd = floor(Nr/2) : Nfft - ceil(Nr/2);

result = zeros(length(fDoppler),length(selectInd));

fftLocal = conj(fft(local, Nfft)');

for i = 1:length(fDoppler)    

    % Wipe-off candidate Doppler
    r_wiped = r .* exp(-1j*2*pi*fDoppler(i)*tArray);

    % Compute the FFT of the wiped signal
    fftWiped = fft(r_wiped, Nfft);

    inverse = ifft(fftWiped.*fftLocal, Nfft);

    % Store the result of the inverse FFT
    result(i, :) = abs(inverse(selectInd));

end


