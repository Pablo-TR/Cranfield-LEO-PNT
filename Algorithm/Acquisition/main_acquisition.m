%% Params
clear all
Fs = 3.2e8; %Sampling frequency [Hz]
Tcoh = 5e-3;% Coherent sampling period [s]
N = round(Fs*Tcoh); % Number of samples per coherent subaccumulation
F = Fs/N*(-N/2:N/2-1);
t = (0:N-1)' / Fs; % Time vector for N samples
fc = 2.304e8; % Carrier frequency [Hz] not used by now 
fd0 = -10000;% Doppler added to local oscillator signal [Hz]
j = 0;
n = 0:1:N-1;

%% Load replica
[r_local, ~] = localReplica(Fs, fc);

%% Read local oscilator signal
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
tStart   = 0; % Start time of signal to read [ms].
tEnd = (tStart + Tcoh)*1e3; % End time of signal to read [ms].
[r, ~, timeArray, ~] = readSignal(filename, Fs, tStart, tEnd);
%% Visualise fft of local oscilator signal
X = abs(fftshift(fft(r)));
figure
plot(F,X)
%% Add Doppler to local oscilator signal
fd = dopplerCurve(j, fd0, N); 
theta_doppler = 2*pi*cumsum(fd).*(1/Fs);
r_doppler = r.*exp(1j*theta_doppler);
X_doppler = abs(fftshift(fft(r_doppler)));
%% Visualise Doppler in local oscilator
figure
plot(F,X_doppler)


%% Doppler search
% Right now unactive for debugging, we change fd_propsed manually
fd_proposed = [-100 -10000]; %Example of Dopppler search. Second one is the right one

for i=1:1:length(fd_proposed)
    %Wipe-off
    r_wiped = r_doppler.*exp(-1j*2*pi*(fd_proposed(i).*n).*(1/Fs));
    X_wiped = fftshift(fft(r_wiped));
  
    %Compute detection metric
    [rcorr, lags] = xcorr(r_wiped, r_local);
    %% Plot result from the 'search'
    time = tStart + (1000*lags/Fs);
    figure
    plot(time, abs(rcorr));
    title(['Correlated signal for Doppler ', num2str(fd_proposed(i)), 'Hz'])
    xlabel('Time [ms]')
    ylabel('Magnitude')
    xlim([tStart, time(end)])
end


%%
function [f_d, fd0] = dopplerCurve(j, fd0, N)
%{
v_s = 7500;% m/s (LEO satellite)
c = 3e8;
f_d_max = (f_c/c)*v_s;
t0 = 50;                % seconds (controls curvature)

% Doppler
f_d = -f_d_max * (t-t(end)/2)./sqrt((t-t(end)/2).^2+t0^2);
%}
f_d = fd0 +j*5; 
f_d = f_d.*ones(1,N);
end