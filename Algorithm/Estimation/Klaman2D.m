%% 1. Definitions (2D - Modelo CV)
% Define the time step
dt = 1;

% State vector (x)
% x = [pos_x; vel_x; pos_y; vel_y] (4x1 vector)
% State Transition Matrix (A) - Modelo de Velocidad Constante
A = [1 dt 0  0; 
     0 1  0  0;
     0 0  1  dt;
     0 0  0  1];

% Observation Matrix (H) - Solo medimos posicion
H = [1 0 0 0;
     0 0 1 0];

% --- Parámetros Base para la SIMULACIÓN ---
n_steps = 200; % Más pasos para ver mejor el recorrido

% "VERDADERO" Ruido del Proceso (Para crear el camino aleatorio)
% Qué tan errático es el movimiento *real* del objeto.
real_accel_std = 0.5; % m/s^2

% "VERDADERO" Ruido del Sensor (Nivel de ruido GPS)
% >>> CAMBIO: Aumentado 10x para ser "MAS RUIDOSO" <<<
R_val_sim = 20000; % m^2 (Varianza EXTREMADAMENTE alta)
R_sim_matrix = [R_val_sim 0; 0 R_val_sim];

% --- Parámetros Base para el FILTRO ---
% Un R y Q "base" para los experimentos.
R_base_filter = R_sim_matrix; % Asumimos que CONOCEMOS el ruido del sensor
Q_base_filter_val = 0.1;
Q_base_filter = [0 0 0 0; 0 Q_base_filter_val 0 0; 0 0 0 0; 0 0 0 Q_base_filter_val];


%% 2. Simulate Example Data (Recorrido Aleatorio + Ruido GPS)

% --- Crear el Recorrido "Verdadero" (Ground Truth) ---
ground_truth_state = zeros(n_steps, 4); % [px, vx, py, vy]
% Matriz de entrada de aceleración (para la simulación)
B = [0.5*dt^2 0;
     dt       0;
     0        0.5*dt^2;
     0        dt];

for t = 2:n_steps
    % Calcular aceleración aleatoria "real"
    real_accel = randn(2, 1) * real_accel_std;
    % Propagar el estado verdadero (el objeto se mueve)
    ground_truth_state(t, :) = A * ground_truth_state(t-1, :)' + B * real_accel;
end
% Extraer solo las posiciones [x, y]
ground_truth_pos = ground_truth_state(:, [1, 3]);

% --- Crear Mediciones Ruidosas (Señal GPS) ---
% Generar ruido masivo
noise = randn(n_steps, 2) * sqrt(R_val_sim);
noisy_measurements = ground_truth_pos + noise;


%% 3. & 4. Kalman Filter Loops (Movidos a los experimentos)

% >>> CAMBIO: Vectores de 3 elementos para guardar los errores <<<
rmse_results_R = zeros(1, 3);
rmse_results_Q = zeros(1, 3);

% --- Colormap para 3 valores ---
colors = parula(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. EXPERIMENTO 1: Sintonizando R (Confianza en el Sensor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% >>> CAMBIO: 3 valores de R para probar <<<
R_values_to_test = [100, 20000, 500000]; 
% 100 = Subestimado (confía demasiado en el sensor)
% 20000 = Correcto (cercano a la verdad)
% 500000 = Sobreestimado (ignora demasiado el sensor)

figure('Name', 'Experimento 1: Sintonizando R (Trayectorias)');
hold on; grid on; axis equal;
title('Sintonizando R: Trayectorias Resultantes (Ruido Extremo)');

% Graficar la verdad y las mediciones (una sola vez)
plot(ground_truth_pos(:, 1), ground_truth_pos(:, 2), 'k', 'LineWidth', 3, 'DisplayName', 'Verdad (Oculta)');
% Graficamos solo el 10% de las mediciones, si no, es ilegible
plot(noisy_measurements(1:10:end, 1), noisy_measurements(1:10:end, 2), 'rx', 'MarkerSize', 8, 'DisplayName', 'Mediciones (GPS)');

for i = 1:length(R_values_to_test)
    
    R_val_current = R_values_to_test(i);
    R_current = [R_val_current 0; 0 R_val_current]; % R del filtro
    Q_current = Q_base_filter; % Q del filtro (fijo)
    
    % --- Resetear el filtro ---
    x_estimate = [noisy_measurements(1, 1); 0; noisy_measurements(1, 2); 0];
    P_estimate = eye(4) * 1e4; % Incertidumbre inicial muy alta
    pos_estimate_hist = zeros(n_steps, 2);
    
    % --- Bucle del Filtro ---
    for t = 1:n_steps
        % PREDICT
        x_pred = A * x_estimate;
        P_pred = A * P_estimate * A' + Q_current;
        % UPDATE
        z = noisy_measurements(t, :)';
        K = P_pred * H' / (H * P_pred * H' + R_current);
        x_estimate = x_pred + K * (z - H * x_pred);
        P_estimate = (eye(4) - K * H) * P_pred;
        % Guardar
        pos_estimate_hist(t, :) = x_estimate([1, 3])';
    end
    
    % --- Graficar la trayectoria resultante ---
    plot(pos_estimate_hist(:, 1), pos_estimate_hist(:, 2), 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', sprintf('R = %.0e', R_val_current));

    % --- Calcular y guardar el error (RMSE) ---
    error = pos_estimate_hist - ground_truth_pos;
    rmse_results_R(i) = sqrt(mean(sum(error.^2, 2)));
end
legend('Location', 'northeast');
xlabel('Posición X'); ylabel('Posición Y');
hold off;

% --- Gráfica de Error vs R ---
figure('Name', 'Experimento 1: Error de R');
% Usamos 'semilogx' para espaciar los puntos correctamente
semilogx(R_values_to_test, rmse_results_R, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on; grid on;
% Marcar el valor "verdadero" de R
xline(R_val_sim, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('R Verdadera (%.0e)', R_val_sim));
title('Error (RMSE) vs. Sintonización de R');
xlabel('Valor de R (log)');
ylabel('Error Total (RMSE)');
legend('show');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. EXPERIMENTO 2: Sintonizando Q (Confianza en el Modelo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% >>> CAMBIO: 3 valores de Q para probar <<<
q_diag_values_to_test = [1e-4, 0.1, 10];
% 1e-4 = Demasiado rígido (confía demasiado en el modelo)
% 0.1  = Equilibrado (cercano a la verdad)
% 10   = Demasiado nervioso (desconfía demasiado del modelo)


figure('Name', 'Experimento 2: Sintonizando Q (Trayectorias)');
hold on; grid on; axis equal;
title('Sintonizando Q: Trayectorias Resultantes (Ruido Extremo)');

% Graficar la verdad y las mediciones (una sola vez)
plot(ground_truth_pos(:, 1), ground_truth_pos(:, 2), 'k', 'LineWidth', 3, 'DisplayName', 'Verdad (Oculta)');
plot(noisy_measurements(1:10:end, 1), noisy_measurements(1:10:end, 2), 'rx', 'MarkerSize', 8, 'DisplayName', 'Mediciones (GPS)');

for i = 1:length(q_diag_values_to_test)
    
    R_current = R_base_filter; % R del filtro (fijo y "correcto")
    q_val_current = q_diag_values_to_test(i);
    Q_current = [0 0 0 0; 0 q_val_current 0 0; 0 0 0 0; 0 0 0 q_val_current]; % Q del filtro
    
    % --- Resetear el filtro ---
    x_estimate = [noisy_measurements(1, 1); 0; noisy_measurements(1, 2); 0];
    P_estimate = eye(4) * 1e4;
    pos_estimate_hist = zeros(n_steps, 2);
    
    % --- Bucle del Filtro ---
    for t = 1:n_steps
        % PREDICT
        x_pred = A * x_estimate;
        P_pred = A * P_estimate * A' + Q_current;
        % UPDATE
        z = noisy_measurements(t, :)';
        K = P_pred * H' / (H * P_pred * H' + R_current);
        x_estimate = x_pred + K * (z - H * x_pred);
        P_estimate = (eye(4) - K * H) * P_pred;
        % Guardar
        pos_estimate_hist(t, :) = x_estimate([1, 3])';
    end
    
    % --- Graficar la trayectoria resultante ---
    plot(pos_estimate_hist(:, 1), pos_estimate_hist(:, 2), 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', sprintf('Q = %.0e', q_val_current));

    % --- Calcular y guardar el error (RMSE) ---
    error = pos_estimate_hist - ground_truth_pos;
    rmse_results_Q(i) = sqrt(mean(sum(error.^2, 2)));
end
legend('Location', 'northeast');
xlabel('Posición X'); ylabel('Posición Y');
hold off;

% --- Gráfica de Error vs Q ---
figure('Name', 'Experimento 2: Error de Q');
semilogx(q_diag_values_to_test, rmse_results_Q, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Error (RMSE) vs. Sintonización de Q');
xlabel('Valor de Q (diagonal, log)');
ylabel('Error Total (RMSE)');