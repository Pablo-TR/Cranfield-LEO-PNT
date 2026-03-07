clc; clear;

%% -- Read signal --
%{
filename = 'rx_stream-320.00-1975.00-2-nvme2.56-10db.bin';
tStart_ms = 0;
tEnd_ms   = 100;
[r, ~, ~, ~] = readSignal(filename, Fs, tStart_ms, tEnd_ms);
r = r(:);
K = floor(length(r)/N);
%}

              


%% -- or Generate dummy signal --
Fs = 3.2e8;
Fsym = 2.304e8;
Tsim = 0.05;
T = 10^-3;
Tpss = 400/Fsym;
Trep = 1e-3;
Nrep = round(Fs*Trep);
fd0_gen = 100000;
fddot = 2000;
N_gen = round(Fs*Tsim);
r = zeros(N_gen,1);
t = (0:N_gen-1).'/Fs;

SNRdB = 30;
SNRin = 10^(SNRdB/10);

sigma_phi2 = 1/(2*SNRin);             
R_fd  = sigma_phi2 / (2*pi*T)^2;

fd = generateDoppler(fd0_gen, fddot,t,0);

figure()
plot(t, fd)
xlabel('Time [s]')
ylabel('Doppler frequency [Hz]')




%% --- KF parameters (block-based) ---
                 
N_coh = round(Fs*T);       % samples per block
K    = floor(length(r)/N_coh);

for k = 1:1:K
    fd_1ms(k) = fd(1+(k-1)*Nrep) + sqrt(R_fd)*randn;
    fd_1ms_true(k) = fd(1+(k-1)*Nrep);
end


tblk = (0:N_coh-1).'/Fs;      % time inside one block [s]

% State transition / measurement
Fmat = [1 T; 0 1];
Hmat = eye(2);

% --- KF initialization ---
fd0     = fd(1);
fd_dot0 = (fd_1ms(2)-fd_1ms(1))/(T);
x  = [fd0; fd_dot0];
P  = diag([R_fd, 10^2]);

q = 100000;
Q = q * [T^3/3  T^2/2;
         T^2/2  T];

sigma_fd    = R_fd;     % Hz
sigma_fddot = 10;       % Hz/s
R = diag([R_fd, sigma_fddot^2]);


% --- Block time axis (center of each block) ---
t_k = ((0:K-1)*N_coh + (N_coh-1)/2).'/Fs;  % seconds

x_out = [];
z_out = [];

x_out = [x_out x];

%% --- KF loop ---
for k = 1:K-1
    % --- Predict ---
    x = Fmat*x;
    P = Fmat*P*Fmat' + Q;

    % --- Measurements  ---
    fd_meas = fd_1ms(k+1);
    fddot_meas = (fd_1ms(k+1) - fd_1ms(k))/T;

    y = [fd_meas; fddot_meas];
    

    % --- Update ---
    S  = Hmat*P*Hmat' + R;
    Kk = P*Hmat'/S;
    z = y - Hmat*x;
    x  = x + Kk*z;
    z_out = [z_out z];
    P  = (eye(2) - Kk*Hmat)*P;

    x_out = [x_out x];
end

RMSE_est  = rmse(fd_1ms_true,x_out(1,:));
RMSE_true = rmse(fd_1ms_true, fd_1ms);

%% --- Plots ---
figure;
plot(t_k, x_out(1,:), 'r--', 'LineWidth',2); hold on;
plot(t_k, fd_1ms, 'k-');
plot(t_k, fd_1ms_true, 'g-');
legend('Estimated Doppler','Measured Doppler (True + noise)', 'True Doppler');
xlabel('Time (s)'); ylabel('Doppler (Hz)'); grid on;

figure;
plot(t_k, x_out(2,:), 'LineWidth',2);
xlabel('Time (s)'); ylabel('Measured Doppler derivative (Hz/s)'); grid on;

figure;
plot(t_k(2:end), z_out(1,:));
xlabel('Time (s)'); ylabel('Residual f_d (Hz)'); grid on;

figure;
plot(t_k(2:end), z_out(2,:));
xlabel('Time (s)'); ylabel('Residual f_ddot (Hz/s)'); grid on;

figure;
plot(t_k, fd_1ms);
hold on
plot(t_k, fd_1ms_true);
xlabel('Time (s)'); ylabel('Doppler frequency (Hz)'); grid on;
legend('Doppler + noise', 'True Doppler')

figure;
plot(t_k, RMSE_true)
hold on
plot(t_k, RMSE_est)
xlabel('Time (s)'); ylabel('RMSE (Hz)'); grid on;
legend('RMSE measurements', 'RMSE KF')





%% Extra
[r_loc, loc] = localReplica(Fs, Fsym);
r_loc = r_loc';

Npss = length(r_loc);

for k = 0:floor((N_gen-1)/Nrep)
    i0 = k*Nrep + 1;
    i1 = min(i0 + Npss - 1, N_gen);
    r(i0:i1) = r(i0:i1) + r_loc(1:(i1-i0+1));
    Ps = mean(abs(r(i0:i1)).^2);          
    sigma2 = Ps / (SNRin);
    w = sqrt(sigma2/2) * (randn(size(r(i0:i1))) + 1j*randn(size(r(i0:i1))));
    r(i0:i1) = r(i0:i1) + w;
end

figure()
plot(t,real(r))
xlabel('Time [s]')
ylabel('Real part of signal')

phi = 2*pi*cumsum(fd)/Fs;
r = r.*exp(1j*phi);

figure()
plot(t,real(r))
xlabel('Time [s]')
ylabel('Real part of signal')

time = 0:600;
fd = generateDoppler(0, 0, time, 0);
figure()
plot(time, fd)
xlabel('Time [s]')
ylabel('Doppler Frequency [Hz]')

%% Functions

function fd = generateDoppler(fd0, fddot, t, doptype)
    if doptype == 1
        idx = ones(1,length(t));
        fd = fd0.*idx; %+ fddot.*t;   % Hz
        fd = fd';
    else
        tpass = 600;
        fc = 2.9e9;
        c = 3*10^8;
        t = t(:);
        v = 7500;
        d = 400e3;
        tc=t - tpass/2;
        rho  = sqrt(d.^2 + (v.*tc).^2);
        rdot = (v.^2 .* tc) ./ rho;      % [m/s],
        fd   = -(fc/c) .* rdot;         % [Hz]
    end
end