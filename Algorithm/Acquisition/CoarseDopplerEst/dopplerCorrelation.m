function [doppler, correlation] = dopplerCorrelation(signalID, tInterval, Tcoh)

signalData = getSignalByID(signalID);
Fs = signalData.SamplingFreq;

nSamples = round(Fs*Tcoh*1e-3); % Number of samples
tArray    = (0:nSamples-1).' / Fs; % Time vector (Nx1)

Fsym = 2.304e8;        % Carrier frequency [Hz] (passed to replica)
[r_local, ~] = localReplica(Fs, Fsym);

[r, ~, ~, ~] = readSignal(signalID, tInterval(1), tInterval(2));
r = r(:);

doppler = -100000:500:100000;


correlation = zeros(length(doppler),1);

for i=1:length(doppler)
    disp(i)
    r_wiped = r .* exp(-1j*2*pi*doppler(i)*tArray);
    % Correlate with local replica
    rcorr = xcorr(r_wiped, r_local);
    
    % Detection metric: peak magnitude
    correlation(i) = max(abs(rcorr));
end
