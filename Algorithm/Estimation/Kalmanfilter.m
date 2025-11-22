%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            ${}.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: [Estimation]
% Author: [Marti Trilla]
% Date: ${09/11/2025}
% Last Update: []
% Version: 1.0.0
%
% Description: [Kalman filter example]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;  % Reset MATLAB environment

%% 1. Definitiions
% Define the time step
dt = 1;

% State vector (x)
% x = [position; velocity] (our desired variables)

% State Transition Matrix (A)
%    x_k     =     A   *     x_k-1
% [pos_k, vel_k] = [1  dt;  
%                   0   1]*[pos_k-1, velk_1]
%
% pos_k = pos_k-1 + vel_k-1 * dt
% vel_k = vel_k-1
A = [1 dt; 
     0 1];

% Observation Matrix (H)
% We only measure the position
% z = [1 0] * [pos; vel]
H = [1 0];

% Measurement Noise Covariance (R) 
% How noisy is our sensor (real data).
% Ideally we know from datasheet of hardware.
R = 100; %(higher = more noise)

% Process Noise Covariance (Q)
% How reliable is our model.
% We add a little noise to acceleration.
Q = [0.01 0; 
     0 0.01]; %(high=bad model, we assume our object is gonna have big and fast changes)

%% 2. Simulate Example Data (Ground Truth and Noisy Measurement) 
% Simulation duration
n_steps = 100;

% Store the ground truth (for comparison)
ground_truth = zeros(n_steps, 1);
% Real (constant) velocity
v_real = 2; 

for t = 2:n_steps
    ground_truth(t) = ground_truth(t-1) + v_real * dt;
end

% Create noisy measurements (z)
% randn generates Gaussian (normal) noise
noisy_measurements = ground_truth + sqrt(R) * randn(n_steps, 1);

%% 3. Kalman Filter 

% Initial state estimate (we start at 0, 0)
x_estimate = [0; 0]; % [initial_pos; initial_vel]

% Initial error covariance (P) (P: Total uncertainty of hte model (hardware
% and our model)
% We start with a high uncertainty (need to initialize)
P_estimate = [0.1 0; 
              0 0.1];

% Arrays to store the results
pos_estimate_hist = zeros(n_steps, 1);
vel_estimate_hist = zeros(n_steps, 1);
pos_predict_hist = zeros(n_steps, 1);

%% 4. The Kalman Filter Loop

for t = 1:n_steps
    
    % --- PHASE 1: PREDICTION ---
    
    % Predict the next state
    % x_pred = A * x_estimate_previous
    x_pred = A * x_estimate;
    
    % Predict the error covariance
    % P_pred = A * P_estimate_previous * A' + Q
    P_pred = A * P_estimate * A' + Q;
    
    % Save the prediction (just for plotting)
    pos_predict_hist(t) = x_pred(1);

    % --- PHASE 2: UPDATE (CORRECTION) ---
    
    % Calculate the Kalman Gain (K)
    % K = P_pred * H' * inv(H * P_pred * H' + R)
    K = P_pred * H' / (H * P_pred * H' + R);
    
    % Get the current measurement
    z = noisy_measurements(t);
    
    % Correct the state estimate
    % x_estimate = x_pred + K * (z - H * x_pred)
    x_estimate = x_pred + K * (z - H * x_pred);
    
    % Correct the error covariance
    % P_estimate = (I - K * H) * P_pred
    P_estimate = (eye(2) - K * H) * P_pred;
    
    % Save the final estimate for this step
    pos_estimate_hist(t) = x_estimate(1);
    vel_estimate_hist(t) = x_estimate(2);
end

%% 5. Visualize the Results 
figure;
hold on;
grid on;

% Plot the ground truth
plot(1:n_steps, ground_truth, 'k', 'LineWidth', 2, 'DisplayName', 'Ground Truth');

% Plot the noisy measurements
plot(1:n_steps, noisy_measurements, 'rx', 'DisplayName', 'Noisy Measurements (Sensor)');

% Plot the Kalman filter estimate
plot(1:n_steps, pos_estimate_hist, 'b-', 'LineWidth', 2, 'DisplayName', 'Kalman Filter Estimate');

% (Optional) Plot the prediction before correction
% plot(1:n_steps, pos_predict_hist, 'g--', 'DisplayName', 'Prediction (Before Correction)');

legend('Location', 'northwest');
title('Kalman Filter: Tracking an Object in 1D');
xlabel('Time Step');
ylabel('Position');
hold off;

% Plot the velocity estimate
figure;
plot(1:n_steps, vel_estimate_hist, 'm-', 'LineWidth', 2, 'DisplayName', 'Estimated Velocity');
yline(v_real, 'k--', 'DisplayName', 'Real Velocity');
title('Velocity Estimate (Hidden State)');
xlabel('Time Step');
ylabel('Velocity');
legend;
grid on;