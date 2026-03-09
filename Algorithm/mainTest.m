clc; clear all; close all;

Fsym = 2.304e8; Fs = 5e8; Fd = 4200; SNR = 100; %dBW

timeSignal = 1;
[synSig, timePeaks] = signalGenerator(Fs, timeSignal, Fd, SNR);

local = localReplica(Fs, Fsym);
localExtended = localExtendedReplica(Fs, timeSignal);

% Correlation between incoming signal and synchronisation sequence
[rcorr, lags] = xcorr(synSig, local);

time = lags'/Fs;
% Keep only time >= 0
idx = time >= 0;
time = time(idx);
rcorrNorm = (abs(rcorr) - min(abs(rcorr)))/(max(abs(rcorr))-min(abs(rcorr)));
rcorrNorm = rcorrNorm(idx);
[maxcorr, idxMax] = max(rcorrNorm);

tol = 1e-8;
assert(abs(timePeaks - (time(idxMax)*1e3-floor(time(idxMax)*1e3)))<tol, 'Correlation peak does not match ')

%%
figure
plot(time, rcorrNorm)
xlabel('Time (ms)','interpreter','latex','fontsize',14)
ylabel('$|$R($\tau$)$|$','interpreter','latex','fontsize',14);
xlim([0 time(end)])
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;


%% dopplerFrequencySearch

fDoppler = -5000:100:5000; % Doppler frequency

localExtended = localExtended(:);

Nx = numel(synSig);
Nr = numel(localExtended);
Nfft = Nr + Nx - 1;
selectInd = floor(Nr/2) : Nfft - ceil(Nr/2);

result = zeros(length(fDoppler),length(selectInd));

fftLocal = conj(fft(localExtended, Nfft));

t = (0:length(synSig)-1)'/Fs;

for i = 1:length(fDoppler)
    % Wipe-off candidate Doppler
    r_wiped = synSig .* exp(-1j*2*pi*fDoppler(i)*t);

    % Compute the FFT of the wiped signal
    fftWiped = fft(r_wiped, Nfft);

    inverse = ifft(fftWiped.*fftLocal, Nfft);

    % Store the result of the inverse FFT
    result(i, :) = abs(inverse(selectInd));

end

%%
[T,F] = meshgrid(t,fDoppler);

% Find maximum correlation value
[maxVal, idx] = max(result(:));

% Get row and column indices
[row, col] = ind2sub(size(result), idx);

% Corresponding time and frequency
tMax = T(row, col);
fMax = F(row, col)

%%
surf(T, F, result, 'EdgeColor','none');
xlabel('Time'); ylabel('Frequency'); zlabel('Correlation');
view(3); colorbar; shading interp;


%%
figure
plot(fDoppler, result(:,col))
xlabel('Freq')
ylabel('Corr')

%% dopplerSearch

numTimeBlocks = 1;
doppler = zeros(numTimeBlocks,1);
timeArr = zeros(numTimeBlocks,1);
fd = 0;

t = (0:length(synSig)-1)'/Fs;

for n = 1:numTimeBlocks
    disp(n);
    disp(numTimeBlocks);    

    tStart = 0*1;
    tEnd     = tStart + 1; % End time [ms]
    
    synSig = synSig(:);
    fd_limit = 10000;

    % Wipe-off candidate Doppler
    r_wiped = synSig .* exp(-1i*2*pi*fd*t);
    
    while fd_limit > 10
        
        % Correlate with local replica
        rcorr = xcorr(r_wiped, local);
    
        % Detection metric: peak magnitude
        likelihood = max(abs(rcorr));
        
        fd_upper = fd + fd_limit;
        r_wiped_upper= synSig .* exp(-1i*2*pi*fd_upper*t);
        rcorr_upper= xcorr(r_wiped_upper, local);
        likelihood_upper= max(abs(rcorr_upper));

        fd_lower = fd - fd_limit;
        r_wiped_lower = synSig .* exp(-1i*2*pi*fd_lower*t);
        rcorr_lower= xcorr(r_wiped_lower, local);
        likelihood_lower= max(abs(rcorr_lower));

        if likelihood_lower < likelihood && likelihood_upper < likelihood
            fd_limit = fd_limit/2;
            continue
        elseif likelihood_lower > likelihood
            fd = fd_lower;
            r_wiped = synSig .* exp(-1i*2*pi*fd*t);
            continue
        elseif likelihood_upper > likelihood
            fd = fd_upper;
            r_wiped = synSig .* exp(-1i*2*pi*fd*t);
        end
    end

    doppler(n) = fd
    timeArr(n) = tStart;

end


%% dopplerCorrelation

fd = -10000:50:10000; % Doppler frequency

t = (0:length(synSig)-1)'/Fs;

correlation = zeros(length(fd),1);

for i=1:length(fd)
    disp(i)
    r_wiped = synSig .* exp(-1i*2*pi*fd(i)*t);
    % Correlate with local replica
    rcorr = xcorr(r_wiped, local);
    
    % Detection metric: peak magnitude
    correlation(i) = max(abs(rcorr));
end


%%
[a,b] = max(correlation);
fd(b)
figure
plot(fd, correlation)
grid on;
xlabel('Frequency')
ylabel('Correlation')