function result = dopplerFrequencySearch(signalID, Tcoh)

signalData = getSignalByID(signalID);
Fs = signalData.SamplingFreq;

nSamples = round(Fs*Tcoh*1e-3); % Number of samples
tArray    = (0:nSamples-1).' / Fs; % Time vector (Nx1)

Fsym = 2.304e8;
[r_local, ~] = localReplica(Fs, Fsym);

fd = -50000:100:50000;

tStart = 0;
tEnd     = tStart + 1; % End time [ms]
[r, ~, ~, ~] = readSignal(signalID, tStart, tEnd);
r = r(:);

Nx = numel(r);
Nr = numel(r_local);
Nfft = Nr + Nx - 1;
selectInd = floor(Nr/2) : Nfft - ceil(Nr/2);

result = zeros(length(fd),length(selectInd));

fftLocal = fft(conj(r_local), Nfft)';

for i = 1:length(fd)    

    % Wipe-off candidate Doppler
    r_wiped = r .* exp(-1j*2*pi*fd(i)*tArray);

    % Compute the FFT of the wiped signal
    fftWiped = fft(r_wiped, Nfft);

    inverse = ifft(fftWiped.*fftLocal, Nfft);

    % Store the result of the inverse FFT
    result(i, :) = abs(inverse(selectInd));

end


