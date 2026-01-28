Fs = 3.2e8;
Fsym = 2.304e8;
t = 20; %ms
[x_l, ~]  = localReplica(Fs, Fsym);
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
x_r = readSignal(filename,Fs, t);
[tau, fd] = performGridSearch(x_l, x_r, Fs);
%% Create sinusoid without Doppler and plot it in frequency domain
clear all
Fs = 10000;            % Sampling frequency                    
Ts = 1/Fs;            % Sampling time
L = 1500000;             % Length of signal
t = (0:L-1)*Ts;        % Time vector
fc = 1e8;
r = exp(1j*2*pi*2500*t); % Baseband sinusoid
Y = abs(fftshift(fft(r)));
F = Fs/L*(-L/2:L/2-1);
figure
plot(F,Y,"LineWidth",3)
figure
plot(t, real(r))
%% Create doppler and add it to sinusoid
k = 1:1:length(r);
fd = dopplerCurve(t, fc); %100.*t;%;
theta_doppler = 2*pi*cumsum(fd).*Ts;%2.*pi*(50.*t.^2); % [rad]
r = r.*exp(1j*theta_doppler);
Y_doppler = abs(fftshift(fft(r)));
figure
plot(t, real(r))
figure
plot(t,fd)
figure
plot(F,Y_doppler)

%% Wipeoff doppler (numerically) and plot
theta_doppler_num = 2*pi*cumsum(fd).*Ts;
r_wiped = r.*exp(-1j*theta_doppler_num);
Y_doppler_wiped = abs(fftshift(fft(r_wiped)));
figure
plot(F,Y_doppler_wiped)

%% What if we don't know doppler values? -> Doppler search
% This approach considers doppler constant at small intervals of the signal
Tcoh = 0.01; % Coherent integration time [s]
K = round(Tcoh * Fs);       % samples per block
N = floor(length(r) / K);  % number of blocks
delta_fd = 1/Tcoh;
fd_search = -2000:100:2000;
for i = 0:1:N-1
    fprintf('######Subaccumulation %d/%d######\n',i,N)
    %Create subaccumulation
    x = r(i*K+1:K*(i+1));
    n = 0:1:length(x)-1;
    for j = 1:1:length(fd_search)
        fprintf('Doppler %d/%d\n',j,length(fd_search))
        %Search for each doppler bin
        fd_proposed = fd_search(j);
        %Wipe-off
        x_wiped = x.*exp(-1j*2*pi*(fd_proposed.*n).*Ts);
        %Compute detection metric
        M(j) = abs(sum(x_wiped));
    end
    idx = find(M == max(M));
    fd_detected(i+1) = fd_search(idx);
end

%%
t_mid = Tcoh/2:Tcoh:t(end);
scatter(t_mid, fd_detected)
hold on
plot(t, fd, 'LineWidth',3)
title('True Doppler vs acquired Doppler')
legend('Acquired doppler','True doppler')
xlabel('Time [s]')
ylabel('Frequency [Hz]')



%%
function f_d = dopplerCurve(t, f_c)
v_s = 7500;% m/s (LEO satellite)
c = 3e8;
f_d_max = (f_c/c)*v_s;
t0 = 50;                % seconds (controls curvature)

% Doppler
f_d = -f_d_max * (t-t(end)/2)./sqrt((t-t(end)/2).^2+t0^2);
end
%%
function [tau_found, fd_found] = performGridSearch(x_l, x_r, Fs)
Ts = 1/Fs;
fd_max = 614*10^3;
tau_max = 4e-3;
fd = -fd_max:3070:fd_max;
tau = -tau_max:1e-4:tau_max;
% Repeat local replica to be the same length as the received signal
K = length(x_r);
x_l = repmat(x_l, ceil(K/length(x_l)), 1);
x_l = x_l(1:K);k = 0:K-1;
for i = 1:1:length(fd)
    fprintf('Searching for fd %d/%d\n', i, length(fd));
    for j = 1:1:length(tau)
    fprintf('Searching for tau %d/%d\n', j, length(tau));
% Convert tau delay into samples
        N = round(tau(j)*Fs);
% Circshift (account for time delay) of local replica
    x_local_shifted = circshift(x_l,N);
%FFT of shifted local replica
    x_local_fft_shifted = fft(x_local_shifted, K);
% Perform wipe-off to received signal
    x_wiped = x_r .* exp(-1j*2*pi*fd(i).*k.*Ts);
% Compute FFT of wiped received signal
    x_wiped_fft = fft(x_wiped, K);
% Cross-correlate with time-shifted local replica
    Y = x_wiped_fft .* conj(x_local_fft_shifted);
% Compute ambiguity function
    R = (1/K^2) * ifft(Y, K);
% Obtain maximum and indices
    S(i,j)  = max(R.^2);
    end
end
    [z_max, idx] = max(S(:));
    [m, n] = ind2sub(size(S), idx);
    fd_found = fd(m);
    tau_found = tau(n);
% Plot ambiguity function
figure
surf(fd, tau,S)
end

%% Steps
% 1. Wipe-off Doppler (Doppler search)
% 2. Correlate with local replica
% 3. Obtain Doppler (max correlation)
% 4. Obtain delay. 