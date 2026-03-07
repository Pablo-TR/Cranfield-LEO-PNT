clear all; clc; close all

%% Params

signalID = 0;
Tcoh = 1;           % Coherent integration time [ms]

tStart = 0;
tEnd = 10;
tInterval = [tStart tEnd]; % Time interval of the signal to read [ms]

%%
[doppler, time] = dopplerSearch(signalID, tInterval, Tcoh);

%%
[fDoppler, correlation] = dopplerCorrelation(signalID, tInterval, Tcoh);

%%
result = dopplerFrequencySearch(signalID, Tcoh);

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
ylim([0 20])
xlabel('Time [ms]')
ylabel('Doppler frequency [kHz] (that gave most likelihood)')

