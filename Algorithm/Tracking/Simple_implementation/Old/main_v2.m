%% Info
% Implementation of Doppler tracking based on EKF for a simple sinusoid
% from: 
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=32088 [1]
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=42091 [2]

%% Generate sinusoid with Doppler using kinematics defined by paper
clear all;
% Define constants
f_c = 2.9e9;     % Carrier frequency in Hz. If 0 no carrier, baseband
c = 3e8;     % Speed of light in m/s
Fs = 200e3;   % Sampling frequency (0.1 megasamples/s )
T = 100;       % Total sampling time
t = 0:1/Fs:T;% Time vector for the duration of the signal
phi_0 = 0;   % Initial phase [rad]
T_sub = 1e-4;% Subaccumulation interval [s]
A = 1; % Signal amplitude gain.
SNR_dB = 40; % Signal to noise ratio

% Generate Doppler frequency experienced by signal
f_d = compute_doppler(f_c,t,c);
f_d_coarse = f_d(1) + rand*2; % Simulate acquisition
plot_doppler(f_d,t)
%%
% Generate complex sinusoid with doppler
 % Compute doppler phase
 phi_d = 2*pi*cumsum(f_d)/Fs;
 % Add initial phase to doppler phase
 phi = phi_0 + phi_d;
 % Compute signal: signal = doppler effect + baseband signal
 x = A.*exp(1j*(phi));
 % Add noise
 x = awgn(x, SNR_dB, 'measured');
 SNR = 10^(SNR_dB/10);
 % Compute noise variance for later
 sigma_sqr_noise = A^2/SNR;


%% Signal subaccumulations

K = floor(T/T_sub); % Number of subaccumulations
N = floor(length(x)/K); %Number of samples in a subaccumulation
r_complex = zeros(1,K); %Subaccumulated signal

%for k = 1:1:K
%    r_complex(k) = (1/N) * sum(x(((k-1)*N+1):k*N)); 
%end
%r_IQ = [real(r_complex); imag(r_complex)]';


%% Compute discriminators: only if the model from [2] is used
%{
Z_I = zeros(1,K);
Z_Q = zeros(1,K);
for k=1:1:K-1
    Z_I(k) = r_IQ(k+1,1)*r_IQ(k,2)-r_IQ(k+1,2)*r_IQ(k,1);
    Z_Q(k) = r_IQ(k+1,1)*r_IQ(k,1)+r_IQ(k+1,2)*r_IQ(k,2);
end
%}
%% Obtain phase from discriminators

% Compute differential phase from discriminators
%dPhi = atan2(Z_Q, Z_I); uncomment if model from [2] is used

%% Extended Kalman Filter [1]

%Initialise
sigma_phase_sqr = 0.5/SNR; % Variance of phase
F = [1 T_sub (T_sub^2)/2 (T_sub^3)/6; 
    0 1 T_sub (T_sub^2)/2;
    0 0 1 T_sub
    0 0 0 1];
    % Initial model values [theta, theta_dot], ideally from acquisition
x_prev = [0 2*pi*f_d_coarse 0 0]';
   % Initial P values
sigma_theta = pi;                 % rad (phase totally unknown)
sigma_w0    = 2*pi*50;           % rad/s (e.g. ±500 Hz uncertainty)
sigma_w1    = 2*pi*5;            % rad/s^2
sigma_w2    = 2*pi*10;             % rad/s^3

P_prev = diag([sigma_theta^2, ...
               sigma_w0^2, ...
               sigma_w1^2, ...
               sigma_w2^2]);
C = [1 0 0 0];
I4 = eye(4);
I2 = eye(2);
R = (sigma_sqr_noise / N) * I2;%sigma_sqr_noise.* I2;
Q = (2*pi*200)^2.*[T_sub^6/252  T_sub^5/72   T_sub^4/30  T_sub^3/24;
                      T_sub^5/72,  T_sub^4/20,  T_sub^3/8,  T_sub^2/6;
                      T_sub^4/30,  T_sub^3/8,   T_sub^2/3,  T_sub/2;
                      T_sub^3/24,  T_sub^2/6,   T_sub/2,    1]; %Multiply by T_sub [1]

for k = 1:1:K
    fprintf('Iteration %d/%d\n', k, K);
    %Predict
    x_pred = F*x_prev;
    P_pred= F*P_prev*F' + Q;

    %Wipe-off Doppler
    idx = (k-1)*N + (1:N);
    t_local = (0:N-1).' / Fs;

    phi_hat = x_pred(1) + x_pred(2) * t_local + 0.5 * x_pred(3) * t_local.^2 + (1/6) * x_pred(4) * t_local.^3;

    x_wiped = x(idx)' .* exp(-1j * phi_hat);
    r_k = mean(x_wiped);
    
    x_meas = [real(r_k); imag(r_k)];
    
    %Update/Correct
    h = A*[cos(C*x_pred); sin(C*x_pred)]; % Measurement model evaluated with what we predict to have (x_pred)
    H = A*[-sin(C*x_pred) 0 0 0; 
        cos(C*x_pred) 0 0 0];
    G = P_pred*H'*inv((H*P_pred*H' + R));
    x_corr = x_pred + G*(x_meas-h);
    x_corr(1) = wrapToPi(x_corr(1));
    P_corr = (I4-G*H)*P_pred;
    
    
    % Store the corrected state and covariance for the next iteration
    x_prev = x_corr;
    P_prev = P_corr;
    
    x_estimation(:,k) = x_corr;
end


%%
%Real phase vs kalman estimated phase
t_sub = T_sub/2:T_sub:T;
true_data_sub = zeros(1,K);

for k = 1:K
    idx = (k-1)*N + (1:N);
    aux = f_d(idx);%atan2(imag(x(idx)),real(x(idx)));
    true_data_sub(k) = mean(aux);
end

figure
plot(t_sub, true_data_sub)
hold on
scatter(t_sub, f_d_coarse+x_estimation(2,:)/(2*pi))


%% Auxiliary functions
function [f_d] = compute_doppler(f_c,t,c)

% Parameters
v_s = 7500;              % m/s (LEO satellite)
f_d_max = (f_c/c)*v_s;
t0 = 50;                % seconds (controls curvature)

% Doppler
f_d = f_d_max * (t-t(end)/2)./sqrt((t-t(end)/2).^2+t0^2);
end

%% Plots
function plot_vel(t,v)
    figure
    plot(t, v)
    hold on
    title('Velocity')
    xlabel('time [s]')
    ylabel('V [m/s]')
    hold off
end

function plot_doppler(fd, t)
    figure
    plot(t, fd)
    hold on
    title('Doppler frequency')
    xlabel('time [s]')
    ylabel('Fd [Hz]')
    hold off
end

function plot_phases(t, p_real, p_estimated)
plot(t,p_real)
hold on
plot(t,p_estimated);
legend('Real phase', 'EKF phase')
title('Signal phase comparison')
xlabel('Time [s]')
ylabel('Phase [rad]')
hold off
end