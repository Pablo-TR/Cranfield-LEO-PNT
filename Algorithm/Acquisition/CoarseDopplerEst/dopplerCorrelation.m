function [doppler, correlation] = dopplerCorrelation(filename, Fs, Tcoh, tStart, tEnd)

N    = round(Fs*Tcoh*1e-3); % Number of samples
t    = (0:N-1).' / Fs; % Time vector (Nx1)

fc = 2.304e8;        % Carrier frequency [Hz] (passed to replica)
[r_local, ~] = localReplica(Fs, fc);

[r, ~, ~, ~] = readSignal(filename, Fs, tStart, tEnd);
r = r(:);

doppler = -100000:500:100000;


correlation = zeros(length(doppler),1);

for i=1:length(doppler)
    disp(i)
    r_wiped = r .* exp(-1j*2*pi*doppler(i)*t);
    % Correlate with local replica
    rcorr = xcorr(r_wiped, r_local);
    
    % Detection metric: peak magnitude
    correlation(i) = max(abs(rcorr));
end
