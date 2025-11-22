%% 1. Definitions
% Define the time step
dt = 1;
% State vector (x)
% x = [position; velocity] (our desired variables)
% State Transition Matrix (A)
%    x_k     =     A   * x_k-1
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

% --- We define the "base" values for simulation ---
% This R is for the *simulation* (the "true" sensor noise)
R_base_for_simulation = 1000; 
% This Q is the "base" model for comparison
Q_base_for_comparison = [0.01 0; 0 0.01]; 


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
% We create this ONCE using our base R, so all tests use the same data
noisy_measurements = ground_truth + sqrt(R_base_for_simulation) * randn(n_steps, 1);


%% 3. (Original) Kalman Filter - NO LONGER USED, MOVED INTO LOOPS
% This section is now inside the loops below to reset the filter


%% 4. (Original) The Kalman Filter Loop - NO LONGER USED, MOVED INTO LOOPS
% This section is now inside the loops below


%% 5. EXPERIMENT 1: Compare 10 different R values

% --- Define the vector of 10 R values to test (logarithmic scale) ---
% 10 values spaced logaritmically from 10^0 (1) to 10^4 (10000)
R_values_to_test = logspace(0, 4, 10); 
% --- Use a colormap for 10 distinct colors ---
colors = parula(10); 

figure;
hold on;
grid on;
title('Experiment 1: Tuning R (10 values)');
xlabel('Time Step');
ylabel('Position');

% Plot the common data first
plot(1:n_steps, ground_truth, 'k', 'LineWidth', 2, 'DisplayName', 'Ground Truth');
plot(1:n_steps, noisy_measurements, 'rx', 'MarkerSize', 8, 'DisplayName', 'Noisy Measurements');

% Loop for each R value
for i = 1:length(R_values_to_test)
    
    R_current = R_values_to_test(i); % This is the R the *filter* will use
    Q_current = Q_base_for_comparison; % Q is fixed
    
    % --- Copied from original Section 3 ---
    x_estimate = [0; 0];
    P_estimate = [0.1 0; 0 0.1];
    pos_estimate_hist = zeros(n_steps, 1);
    vel_estimate_hist = zeros(n_steps, 1);
    pos_predict_hist = zeros(n_steps, 1);
    % --- End of copied Section 3 ---

    % --- Copied from original Section 4 ---
    for t = 1:n_steps
        % PREDICTION
        x_pred = A * x_estimate;
        P_pred = A * P_estimate * A' + Q_current; % Use Q_current
        
        pos_predict_hist(t) = x_pred(1);
        
        % UPDATE (CORRECTION)
        % --- This is the key change for this experiment ---
        K = P_pred * H' / (H * P_pred * H' + R_current); % Use R_current
        % ---
        
        z = noisy_measurements(t);
        x_estimate = x_pred + K * (z - H * x_pred);
        P_estimate = (eye(2) - K * H) * P_pred;
        
        pos_estimate_hist(t) = x_estimate(1);
        vel_estimate_hist(t) = x_estimate(2);
    end
    % --- End of copied Section 4 ---
    
    % Plot the result for THIS R value
    % We use colors(i, :) because parula(10) returns a 10x3 matrix of RGB values
    plot(1:n_steps, pos_estimate_hist, 'Color', colors(i, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Estimate (R = %.0f)', R_current));
end
legend('Location', 'northwest');
hold off;


%% 6. EXPERIMENT 2: Compare 10 different Q values

% --- Define the vector of 10 Q diagonal values to test (logarithmic scale) ---
% 10 values spaced logaritmically from 10^-5 (0.00001) to 10^0 (1)
q_diag_values_to_test = logspace(-5, 2, 10);
% --- Use a colormap for 10 distinct colors ---
colors = parula(10);

figure;
hold on;
grid on;
title('Experiment 2: Tuning Q (10 values)');
xlabel('Time Step');
ylabel('Position');

% Plot the common data first
plot(1:n_steps, ground_truth, 'k', 'LineWidth', 2, 'DisplayName', 'Ground Truth');
plot(1:n_steps, noisy_measurements, 'rx', 'MarkerSize', 8, 'DisplayName', 'Noisy Measurements');

% Loop for each Q value
for i = 1:length(q_diag_values_to_test)
    
    R_current = R_base_for_simulation; % R is fixed
    
    % --- This is the key change for this experiment ---
    q_val = q_diag_values_to_test(i);
    Q_current = [q_val 0; 0 q_val]; % Construct the Q matrix
    % ---
    
    % --- Copied from original Section 3 ---
    x_estimate = [0; 0];
    P_estimate = [0.1 0; 0 0.1];
    pos_estimate_hist = zeros(n_steps, 1);
    vel_estimate_hist = zeros(n_steps, 1);
    pos_predict_hist = zeros(n_steps, 1);
    % --- End of copied Section 3 ---

    % --- Copied from original Section 4 ---
    for t = 1:n_steps
        % PREDICTION
        % --- This is the key change for this experiment ---
        x_pred = A * x_estimate;
        P_pred = A * P_estimate * A' + Q_current; % Use Q_current
        % ---

        pos_predict_hist(t) = x_pred(1);
        
        % UPDATE (CORRECTION)
        K = P_pred * H' / (H * P_pred * H' + R_current); % Use R_current (which is fixed)
        
        z = noisy_measurements(t);
        x_estimate = x_pred + K * (z - H * x_pred);
        P_estimate = (eye(2) - K * H) * P_pred;
        
        pos_estimate_hist(t) = x_estimate(1);
        vel_estimate_hist(t) = x_estimate(2);
    end
    % --- End of copied Section 4 ---
    
    % Plot the result for THIS Q value
    plot(1:n_steps, pos_estimate_hist, 'Color', colors(i, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('Estimate (Q diag = %.5f)', q_val));
end
legend('Location', 'northwest');
hold off;