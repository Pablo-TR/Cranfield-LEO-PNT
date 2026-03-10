clc;
clear;
close all;

%% ========================================================================
%  READ INPUT SIGNAL
%  ========================================================================
filename  = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
tStart_ms = 0;
tEnd_ms   = 100;
[r, ~, ~, ~] = readSignal(filename, Fs, tStart_ms, tEnd_ms);
r = r(:);
K = floor(length(r)/N);

%% ========================================================================
%  ACQUISITION STAGE
%  ========================================================================



%% ========================================================================
%  NOISE / MEASUREMENT MODEL
%  ========================================================================

SNRdB = 40;
SNRin = 10^(SNRdB/10);

sigma_phi2 = 1 / (2 * SNRin);
R_fd       = sigma_phi2 / (2*pi*T)^2;

%% ========================================================================
%  TRUE DOPPLER GENERATION
%  ========================================================================

fd = generateDoppler(fd0_gen, fddot, t, 0);

figure;
plot(t, fd, 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Doppler frequency [Hz]');
title('True Doppler evolution');
grid on;

%% ========================================================================
%  BLOCK-BASED MEASUREMENTS
%  ========================================================================

N_coh = round(Fs * T);          % Samples per block
K     = floor(length(r) / N_coh);

fd_1ms      = zeros(K, 1);
fd_1ms_true = zeros(K, 1);

for k = 1:K
    idx = 1 + (k-1)*Nrep;
    fd_1ms_true(k) = fd(idx);
    fd_1ms(k)      = fd(idx) + sqrt(R_fd) * randn;
end

t_blk = (0:N_coh-1).' / Fs;  % Time axis within one block

% Time axis at block centers
t_k = ((0:K-1) * N_coh + (N_coh-1)/2).' / Fs;

%% ========================================================================
%  KALMAN FILTER SETUP (2-STATE MODEL)
%  State: x = [fd; fd_dot]
%  ========================================================================

Fmat = [1, T;
        0, 1];

Hmat = eye(2);

% Initial state
fd0     = fd(1);
fd_dot0 = (fd_1ms(2) - fd_1ms(1)) / T;
x       = [fd0; fd_dot0];

% Initial covariance
P = diag([R_fd, 10000^2]);

% Process noise covariance
q = 1e2;
Q = q * [T^3/3, T^2/2;
         T^2/2, T];

% Measurement noise covariance
sigma_fd     = R_fd;  
sigma_fd_dot = 10000;     % [Hz/s]

R = diag([R_fd, sigma_fd_dot^2]);

%% ========================================================================
%  STORAGE VARIABLES
%  ========================================================================

x_out = zeros(2, K);
z_out = zeros(2, K-1);

x_out(:,1) = x;

%% ========================================================================
%  KALMAN FILTER LOOP
%  ========================================================================

for k = 1:K-1

    % -----------------------------
    % Prediction
    % -----------------------------
    x = Fmat * x;
    P = Fmat * P * Fmat' + Q;

    % -----------------------------
    % Measurements
    % -----------------------------
    fd_meas     = fd_1ms(k+1);
    fd_dot_meas = (fd_1ms(k+1) - fd_1ms(k)) / T;

    y = [fd_meas;
         fd_dot_meas];

    % -----------------------------
    % Update
    % -----------------------------
    S  = Hmat * P * Hmat' + R;
    Kk = P * Hmat' / S;
    z  = y - Hmat * x;

    x = x + Kk * z;
    P = (eye(2) - Kk * Hmat) * P;

    % Save outputs
    x_out(:,k+1) = x;
    z_out(:,k)   = z;
end

%% ========================================================================
%  PERFORMANCE METRICS
%  ========================================================================

RMSE_est  = rmse(fd_1ms_true, x_out(1,:).');
RMSE_meas = rmse(fd_1ms_true, fd_1ms);

%% ========================================================================
%  PLOTS
%  ========================================================================

figure;
plot(t_k, x_out(1,:), 'r--', 'LineWidth', 2); hold on;
plot(t_k, fd_1ms, 'k', 'LineWidth', 1);
plot(t_k, fd_1ms_true, 'g', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Doppler [Hz]');
title('Kalman filter Doppler estimation');
legend('Estimated Doppler', 'Measured Doppler', 'True Doppler', 'Location', 'best');
grid on;

figure;
plot(t_k, x_out(2,:), 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Doppler rate [Hz/s]');
title('Estimated Doppler derivative');
grid on;

figure;
plot(t_k(2:end), z_out(1,:), 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Residual on f_d [Hz]');
title('Innovation / residual for Doppler');
grid on;

figure;
plot(t_k(2:end), z_out(2,:), 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Residual on \it{df_d}/dt [Hz/s]');
title('Innovation / residual for Doppler derivative');
grid on;

figure;
plot(t_k, fd_1ms, 'LineWidth', 1.1); hold on;
plot(t_k, fd_1ms_true, 'LineWidth', 1.1);
xlabel('Time [s]');
ylabel('Doppler [Hz]');
title('Measured Doppler vs true Doppler');
legend('Doppler + noise', 'True Doppler', 'Location', 'best');
grid on;

figure;
plot(t_k, RMSE_meas * ones(size(t_k)), 'LineWidth', 1.5); hold on;
plot(t_k, RMSE_est  * ones(size(t_k)), 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('RMSE [Hz]');
title('RMSE comparison');
legend('Measurement RMSE', 'KF RMSE', 'Location', 'best');
grid on;

%% ========================================================================
%  EXTRA: SIGNAL GENERATION WITH LOCAL REPLICA + NOISE + DOPPLER
%  ========================================================================

[r_loc, loc] = localReplica(Fs, Fsym); 
r_loc = r_loc.';
Npss  = length(r_loc);

for k = 0:floor((N_gen - 1) / Nrep)

    i0 = k * Nrep + 1;
    i1 = min(i0 + Npss - 1, N_gen);

    r(i0:i1) = r(i0:i1) + r_loc(1:(i1 - i0 + 1));

    Ps     = mean(abs(r(i0:i1)).^2);
    sigma2 = Ps / SNRin;

    w = sqrt(sigma2/2) * ...
        (randn(size(r(i0:i1))) + 1j * randn(size(r(i0:i1))));

    r(i0:i1) = r(i0:i1) + w;
end

figure;
plot(t, real(r));
xlabel('Time [s]');
ylabel('Real part of signal');
title('Signal before Doppler modulation');
grid on;

phi = 2*pi*cumsum(fd) / Fs;
r   = r .* exp(1j * phi);

figure;
plot(t, real(r));
xlabel('Time [s]');
ylabel('Real part of signal');
title('Signal after Doppler modulation');
grid on;

%% ========================================================================
%  EXTRA TEST
%  ========================================================================

time_test = 0:600;
fd_test   = generateDoppler(0, 0, time_test, 0);

figure;
plot(time_test, fd_test, 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Doppler frequency [Hz]');
title('Generated Doppler test');
grid on;

%% ========================================================================
%  FUNCTIONS
%  ========================================================================

function fd = generateDoppler(fd0, fddot, t, doptype)

    if doptype == 1
        idx = ones(1, length(t));
        fd  = fd0 .* idx;   % + fddot .* t;
        fd  = fd(:);

    else
        tpass = 600;
        fc    = 2.9e9;
        c     = 3e8;
        v     = 7500;
        d     = 400e3;

        t  = t(:);
        tc = t - tpass/2;

        rho  = sqrt(d.^2 + (v .* tc).^2);
        rdot = (v.^2 .* tc) ./ rho;   % [m/s]
        fd   = -(fc/c) .* rdot;       % [Hz]
    end
end