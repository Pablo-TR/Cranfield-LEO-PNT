function result = dopplerFrequencySearch(filename, Fs, timeSignal, t)

fc = 2.304e8;
[r_local, ~] = localReplica(Fs, fc);

%fd = -50000:100:50000;

fd=12600;
tStart = 0;
tEnd     = tStart + 1; % End time [ms]
[r, ~, ~, ~] = readSignal(filename, Fs, tStart, tEnd);
r = r(:);

Nx = numel(r);
Nr = numel(r_local);
Nfft = Nr + Nx - 1;
selectInd = floor(Nr/2) : Nfft - ceil(Nr/2);

result = zeros(length(fd),length(selectInd));

fftLocal = fft(conj(r_local), Nfft)';

for i = 1:length(fd)    

    % Wipe-off candidate Doppler
    r_wiped = r .* exp(-1j*2*pi*fd(i)*t);

    % Compute the FFT of the wiped signal
    fftWiped = fft(r_wiped, Nfft);

    inverse = ifft(fftWiped.*fftLocal, Nfft);

    % Store the result of the inverse FFT
    result(i, :) = abs(inverse(selectInd));

end


