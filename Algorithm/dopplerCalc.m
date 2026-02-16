clear; clc; close all

%% Params
Fs   = 3.2e8;          % Sampling frequency [Hz]
Tcoh = 10;           % Coherent integration time [ms]
N    = round(Fs*Tcoh*1e-3); % Number of samples
t    = (0:N-1).' / Fs; % Time vector (Nx1)
fc   = 2.304e8;        % Carrier frequency [Hz] (passed to replica)

% Load replica
[r_local, ~] = localReplica(Fs, fc);

% Ensure replica is Nx1 for correlation
r_local = r_local(:);

% Read received signal
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
timeSignal = 10; % Time of the signal read [ms]
numTimeBlocks = timeSignal/Tcoh;

% Doppler search (wipe-off then correlate)
fd_proposed = -80000:1000:80000;         % Search grid [Hz]

doppler = zeros(length(numTimeBlocks),2);
timeArr = zeros(length(numTimeBlocks),1);

for n = 1:numTimeBlocks

    disp(n);
    disp(numTimeBlocks);
    tStart = 10 + (n-1)*Tcoh;
    tEnd     = 10 + n*Tcoh; % End time [ms]
    [r, ~, timeArray, ~] = readSignal(filename, Fs, tStart, tEnd);
    
    r = r(:);
    likelihood  = zeros(numel(fd_proposed),1);
    
    for i = 1:numel(fd_proposed)
        disp(fd_proposed(i))
        % Wipe-off candidate Doppler
        r_wiped = r .* exp(-1j*2*pi*fd_proposed(i)*t);
    
        % Correlate with local replica
        rcorr = xcorr(r_wiped, r_local);
    
        % Detection metric: peak magnitude
        likelihood(i) = max(abs(rcorr));
    end

    [a, b] = max(likelihood);

    doppler(n,:) = [fd_proposed(b), a];
    timeArr(n) = tStart;

end

%%

figure; plot(fd_proposed/1e3, likelihood, 'b-')
ylabel('Likelihood')
xlabel('Doppler frequency [kHz]')
grid on

%% Plot likelihood vs Doppler
figure
hold on
plot(timeArr, doppler(:,1)/1e3, 'LineWidth', 1.5)
grid on
ylim([0 80])
xlabel('Time [ms]')
ylabel('Doppler frequency [kHz] (that gave most likelihood)')



%%

% Adding Doppler
fd_added = 1000; % [Hz] Doppler you artificially add (for testing)
r_doppler = r .* exp(1j*2*pi*fd_added*t);

r_wipedDoppler = r_doppler .* exp(-1j*2*pi*fd_added*t);

[rcorr, lags] = xcorr(r, r_local);
time = tStart + (1000*lags/Fs);
%
figure
plot(time, abs(rcorr))
xlabel('Time (ms)','interpreter','latex','fontsize',14)
ylabel('$|$R($\tau$)$|$','interpreter','latex','fontsize',14);
xlim([tStart time(end)])
title('No Doppler')
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;

[rcorr2, lags2] = xcorr(r_doppler, r_local);
time2 = tStart + (1000*lags2/Fs);
%
figure
plot(time2, abs(rcorr2))
xlabel('Time (ms)','interpreter','latex','fontsize',14)
ylabel('$|$R($\tau$)$|$','interpreter','latex','fontsize',14);
xlim([tStart time2(end)])
title('Doppler')
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;

[rcorr3, lags3] = xcorr(r_wipedDoppler, r_local);
time3 = tStart + (1000*lags/Fs);
%
figure
plot(time3, abs(rcorr3))
xlabel('Time (ms)','interpreter','latex','fontsize',14)
ylabel('$|$R($\tau$)$|$','interpreter','latex','fontsize',14);
xlim([tStart time3(end)])
title('Dopplered No doppler')
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 1.2);
grid on;
