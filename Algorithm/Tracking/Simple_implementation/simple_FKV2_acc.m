clc; clear;

%% -- Generate dummy signal --
Fs = 3.2e8;
Fsym = 2.304e8;
Tsim = 0.05;
T = 1e-3;                 % KF update period
Tpss = 400/Fsym;
Trep = 1e-3;
Nrep = round(Fs*Trep);
fd0_gen = 100000;
fddot_gen = 2000;         % this is fd_dot slope if doptype==1
N_gen = round(Fs*Tsim);
r = zeros(N_gen,1);
t = (0:N_gen-1).'/Fs;

SNRdB = 5;
SNRin = 10^(SNRdB/10);

sigma_phi2 = 1/(2*SNRin);
R_fd  = sigma_phi2 / (2*pi*T)^2;

fd = generateDoppler(fd0_gen, fddot_gen, t, 0);

figure()
plot(t, fd)
xlabel('Time [s]')
ylabel('Doppler frequency [Hz]')
grid on

%% --- KF parameters ---
N_coh = round(Fs*T);
K = floor(length(r)/N_coh);

fd_1ms = zeros(K,1);
fd_1ms_true = zeros(K,1);

for k = 1:K
    fd_1ms(k) = fd(1+(k-1)*Nrep) + sqrt(R_fd)*randn;
    fd_1ms_true(k) = fd(1+(k-1)*Nrep);
end

%% ------------------------------------------------------------------------
%  3-state KF: [fd ; fd_dot ; fd_ddot]
%  Constant-acceleration model

Fmat = [1, T, 0.5*T^2;
        0, 1, T;
        0, 0, 1];

Hmat = [1, 0, 0];     % measure only fd

% --- Initialisation ---
fd0      = fd_1ms(1);
fd_dot0  = (fd_1ms(2) - fd_1ms(1))/T;
fd_ddot0 = 0;         

x = [fd0; fd_dot0; fd_ddot0];

P = diag([R_fd, (1e4)^2, (0e4)^2]);  

% Process noise:
qj = 1e-5;   

Q = qj * [T^5/20, T^4/8,  T^3/6;
          T^4/8,  T^3/3,  T^2/2;
          T^3/6,  T^2/2,  T];

R = R_fd;    

% --- Time axis ---
t_k = ((0:K-1)*N_coh + (N_coh-1)/2).'/Fs;

x_out = zeros(3,K);
z_out = zeros(1,K-1);
x_out(:,1) = x;

%% --- KF loop ---
for k = 1:K-1

    % Predict
    x = Fmat*x;
    P = Fmat*P*Fmat' + Q;

    % Measurement
    y = fd_1ms(k+1);

    % Update
    S  = Hmat*P*Hmat' + R;     % scalar
    Kk = P*Hmat'/S;            % 3x1
    z  = y - Hmat*x;           % innovation scalar

    x = x + Kk*z;
    P = (eye(3) - Kk*Hmat)*P;

    x_out(:,k+1) = x;
    z_out(k) = z;
end

RMSE_est  = rmse(fd_1ms_true, x_out(1,:)');
RMSE_true = rmse(fd_1ms_true, fd_1ms);

%% --- Plots ---
figure;
plot(t_k, x_out(1,:), 'r--', 'LineWidth', 2); hold on;
plot(t_k, fd_1ms, 'k-');
plot(t_k, fd_1ms_true, 'g-', 'LineWidth', 1.2);
legend('Estimated Doppler','Measured Doppler', 'True Doppler');
xlabel('Time (s)');
ylabel('Doppler (Hz)');
grid on;

figure;
plot(t_k, x_out(2,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Estimated Doppler rate (Hz/s)');
grid on;

figure;
plot(t_k, x_out(3,:), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Estimated Doppler acceleration (Hz/s^2)');
grid on;

figure;
plot(t_k(2:end), z_out, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Innovation on f_d (Hz)');
grid on;

figure;
plot(t_k, fd_1ms, 'k'); hold on;
plot(t_k, fd_1ms_true, 'g');
xlabel('Time (s)');
ylabel('Doppler frequency (Hz)');
legend('Doppler + noise', 'True Doppler');
grid on;

figure;
plot(t_k, RMSE_true*ones(size(t_k))); hold on;
plot(t_k, RMSE_est*ones(size(t_k)));
xlabel('Time (s)');
ylabel('RMSE (Hz)');
legend('RMSE measurements', 'RMSE KF');
grid on;

%% Functions

function fd = generateDoppler(fd0, fddot, t, doptype)
    if doptype == 1
        t = t(:);
        fd = fd0 + fddot.*t;   % linear Doppler evolution
    else
        tpass = 600;
        fc = 2.9e9;
        c = 3e8;
        t = t(:);
        v = 7500;
        d = 400e3;
        tc = t - tpass/2;
        rho  = sqrt(d.^2 + (v.*tc).^2);
        rdot = (v.^2 .* tc) ./ rho;
        fd   = -(fc/c) .* rdot;
    end
end