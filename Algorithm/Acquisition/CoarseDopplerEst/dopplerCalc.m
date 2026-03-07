clear all; clc; close all

%% Params
FsRec   = 3.2e8;          % Sampling frequency [Hz]
FsOsc = 5e8; 
Tcoh = 1;           % Coherent integration time [ms]
NRec    = round(FsRec*Tcoh*1e-3); % Number of samples
tRec    = (0:NRec-1).' / FsRec; % Time vector (Nx1)

NOsc    = round(FsOsc*Tcoh*1e-3); % Number of samples
tOsc    = (0:NOsc-1).' / FsOsc; % Time vector (Nx1)



% Read received signal
filenameRec = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
filenameOsc = 'OW-sample500-250-5s-ram-29.bin';

timeSignal = 1; % Time of the signal read [ms]
numTimeBlocks = timeSignal/Tcoh;
%%
[doppler, time] = dopplerSearch(filenameOsc, FsOsc, timeSignal, Tcoh, tOsc);


%%
[doppler, correlation] = dopplerCorrelation(filenameOsc, FsOsc, Tcoh, 0, 1);

%%
result = dopplerFrequencySearch(filenameOsc, FsOsc, timeSignal, tOsc);
%%
figure
plot(tOsc,result)
grid on;
%%
fd = -50000:100:50000;
[maxVal, linearIndex] = max(result(:));
[row, col] = ind2sub(size(result), linearIndex);

disp(maxVal)
disp(fd(row))   % frequency index
disp(tOsc(col))   % time index


%%
figure
plot(tOsc, result(500,:))
grid on;
xlabel('Time')
ylabel('Correlation')

%%

[T,F] = meshgrid(tOsc,fd);

surf(T, F, result, 'EdgeColor','none');
xlabel('Time'); ylabel('Frequency'); zlabel('Correlation');
view(3); colorbar; shading interp;
%%
figure
plot(doppler/1e3, correlation)
grid on;
xlabel('Doppler frequency [kHz]')
ylabel('Correlation')

%% Plot time vs Doppler
figure
hold on
plot(time, doppler/1e3, 'LineWidth', 1.5)
grid on
%ylim([0 80])
xlabel('Time [ms]')
ylabel('Doppler frequency [kHz] (that gave most likelihood)')

