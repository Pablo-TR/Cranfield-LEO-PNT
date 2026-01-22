clc; clear;

%% --- Simulation parameters ---
T = 0.0001;          % 1 ms sample period
N = 2000000;          % 20 s
t = (0:N-1)*T;

f_c = 1575.42e6;    % L1 GNSS carrier (Hz)
c = 3e8;

%% --- True Doppler ---
v_s = 7500;
t0 = 50;
f_d_max = f_c/c*v_s;
doppler_Hz = f_d_max * ((t - t(end)/2) ./ sqrt((t - t(end)/2).^2 + t0^2));

f0 = 10e3;
f_true = f0 + doppler_Hz;

phi_true = cumsum(f_true*T);

%% --- IQ measurements ---
sigma = 0.05;
I = cos(2*pi*phi_true) + sigma*randn(1,N);
Q = sin(2*pi*phi_true) + sigma*randn(1,N);

%% --- EKF initialization (4th-order) ---
x = [0; f0; 0; 0];     % [phase; freq; freq_rate; doppler]
P = eye(4);
Qk = diag([1e-8 10 1 50]);
Rk = sigma^2 * eye(2);

phi_est = zeros(1,N);
freq_est = zeros(1,N);
doppler_est = zeros(1,N);

%% --- Subaccumulation parameters ---
sub_len = 10;  % accumulate 10 samples
num_sub = floor(N/sub_len);

%% --- EKF loop with subaccumulation & wipeoff ---
for k = 1:num_sub
    idx = (k-1)*sub_len + (1:sub_len);

    % --- Subaccumulated IQ ---
    I_sub = sum(I(idx));
    Q_sub = sum(Q(idx));
    
    % --- Carrier wipeoff using previous phase estimate ---
    phi_hat = x(1);
    z = (I_sub + Q_sub*1j) * exp(-1j*2*pi*phi_hat);  % complex rotation

    % Separate back to real/imag
    z_real = real(z);
    z_imag = imag(z);

    % --- EKF Prediction ---
    F = [1 T*sub_len 0.5*(T*sub_len)^2 T*sub_len;
         0 1 T*sub_len 0;
         0 0 1 0;
         0 0 0 1];

    x = F*x;
    P = F*P*F' + Qk;

    % --- Measurement model ---
    z_hat = [1; 0];  % after wipeoff, expected residual is (I=1,Q=0)
    H = [0 0 0 0; 0 0 0 0]; % simplified linearization (small residual)
    H(1,1)=1; H(2,1)=0;    % small-angle approx
    K = P*H'/(H*P*H' + Rk);
    x = x + K*([z_real; z_imag]-z_hat);
    P = (eye(4) - K*H)*P;

    % --- Save estimates ---
    phi_est(idx) = x(1);
    freq_est(idx) = x(2);
    doppler_est(idx) = x(4);
end

%% --- Plot results ---
figure;
subplot(3,1,1);
plot(t, phi_true, 'k', t, unwrap(phi_est), 'r--');
legend('True phase','Estimated phase'); xlabel('Time (s)'); ylabel('Phase (cycles)'); grid on;

subplot(3,1,2);
plot(t, f_true, 'k', t, freq_est + doppler_est, 'r--');
legend('True freq','Estimated freq + Doppler'); xlabel('Time (s)'); ylabel('Freq (Hz)'); grid on;

subplot(3,1,3);
plot(t, doppler_Hz, 'k', t, doppler_est, 'r--');
legend('True Doppler','Estimated Doppler'); xlabel('Time (s)'); ylabel('Doppler (Hz)'); grid on;
