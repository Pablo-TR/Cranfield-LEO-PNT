clc; clear; close all;

%%
filename = 'OW-sample500-250-5s-ram-29.bin';

tDur = 5; % [ms]
tSignal = 2500; % [ms]
Fs = 3.2e8;
Fsym = 2.304e8;

nPeaks = 5;
nBlocks = tSignal/tDur;
tPeaksHighest = nan(1, nBlocks);

%%
[rawLocalReplica, localReplica] = localReplica(Fs, Fsym);
for i = 1:nBlocks
    disp([i, nBlocks])
    tStart = 0+(i-1)*tDur;
    tEnd = tStart + tDur;
    [signalDetrended, ~, timeArray, binData] = readSignal(filename, Fs, tStart, tEnd);
    [rcorr, lags] = xcorr(signalDetrended, rawLocalReplica);

    time = (1000*lags/Fs); % Starts at 0 secs as lags does not have prev knowledge
    rcorrNorm = (abs(rcorr) - min(abs(rcorr)))/(max(abs(rcorr))-min(abs(rcorr)));

    [allPeaks, allLocs] = sort(rcorrNorm, 'descend');
    selectedPeaks = [];
    selectedLocs  = [];
    
    % Select the 5 highest peaks and make sure they are at a distance greater
    % than the length of the local replica.
    j = 1;
    while length(selectedLocs) < nPeaks && j <= length(allLocs)
        if isempty(selectedLocs) || all(abs(allLocs(j) - selectedLocs) >= length(localReplica))
            selectedPeaks(end+1) = allPeaks(j);
            selectedLocs(end+1)  = allLocs(j);
        end
        j = j + 1;
    end
    
    peakTimes = time(selectedLocs);
    tPeaksHighestArray = peakTimes(1:5);
    tPeaksHighest(i) = min(tPeaksHighestArray);

end

%%

time = tDur/2:tDur:tDur*nBlocks; % Assume time delay within the block is constant

% Plotting time delay relative to the first measurement
figure
plot(time, tPeaksHighest-tPeaksHighest(1))
xlabel('Time (ms)','interpreter','latex','fontsize',14)
ylabel('$\tau_D$ (ms)','interpreter','latex','fontsize',14);
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1.2);
grid on;
