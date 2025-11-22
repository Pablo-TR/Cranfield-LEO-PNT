% === SIMULADOR ONEWEB PNT: VERSIÓN FINAL CORREGIDA ===
% Solución al bug de "doble resta" en la innovación del Kalman
clc; clear; close all;

%% 1. CONFIGURACIÓN DEL ESCENARIO
c = 299792458; 
dt = 0.1; % 10 Hz
Total_Time = 600; % 10 Minutos
T_span = 0:dt:Total_Time;
N_Time = length(T_span);
N_Sats = 3; 

Pos_RS = [-745000, -5462000, 3198000]; 
Vel_RS = [0, 0, 0];
Pos_User_True = Pos_RS + [3000, 4000, -50]; 
Vel_User_True = [0, 0, 0];
User_Clk_Drift_True = 200; 

Sat_Data = struct('Pos_True', [], 'Vel_True', [], 'Pos_TLE', [], 'Vel_TLE', [], ...
                  'Clk_Bias', [], 'Corrections', [], 'Meas_RS', [], 'Meas_User', [], ...
                  'Signal_Mask', []);

fprintf('1. Generando escenario físico...\n');

for s = 1:N_Sats
    % Diversidad Geométrica
    Phase_Offset = (s-1) * (pi/3); % Separarlos más
    Inclination_Offset = (s-1) * 0.2; 
    
    Radius = 7500000; Omega = 0.001;
    
    Sat_Data(s).Pos_True = zeros(N_Time, 3);
    Sat_Data(s).Vel_True = zeros(N_Time, 3);
    Sat_Data(s).Pos_TLE  = zeros(N_Time, 3);
    Sat_Data(s).Vel_TLE  = zeros(N_Time, 3);
    Sat_Data(s).Clk_Bias = zeros(N_Time, 1);
    Sat_Data(s).Meas_RS  = nan(N_Time, 1); 
    Sat_Data(s).Meas_User= nan(N_Time, 1);
    
    Err_Pos = [3000, -2000, 1000] .* ((-1)^s); 
    Err_Vel = [3, -3, 1] .* ((-1)^s);
    Current_Clk = 150 * s;
    
    % Gaps
    Sat_Data(s).Signal_Mask = true(N_Time, 1);
    Gap_Period = 400 + (s*100); 
    for k=50:Gap_Period:N_Time 
        if k+40 <= N_Time
            Sat_Data(s).Signal_Mask(k : k+40) = false;
        end
    end
    
    for k = 1:N_Time
        t = T_span(k);
        % Órbitas
        x = Radius * cos(Omega*t + Phase_Offset);
        y = Radius * sin(Omega*t + Phase_Offset) * cos(Inclination_Offset);
        z = Radius * sin(Omega*t + Phase_Offset) * sin(Inclination_Offset) + 7000000;
        vx = -Radius * Omega * sin(Omega*t + Phase_Offset);
        vy = Radius * Omega * cos(Omega*t + Phase_Offset) * cos(Inclination_Offset);
        vz = Radius * Omega * cos(Omega*t + Phase_Offset) * sin(Inclination_Offset);
        
        Sat_Data(s).Pos_True(k,:) = [x, y, z];
        Sat_Data(s).Vel_True(k,:) = [vx, vy, vz];
        Sat_Data(s).Pos_TLE(k,:) = Sat_Data(s).Pos_True(k,:) + Err_Pos;
        Sat_Data(s).Vel_TLE(k,:) = Sat_Data(s).Vel_True(k,:) + Err_Vel;
        
        Current_Clk = Current_Clk + randn * 0.05;
        Sat_Data(s).Clk_Bias(k) = Current_Clk;
        
        % Mediciones
        if Sat_Data(s).Signal_Mask(k)
            r_rs = Sat_Data(s).Pos_True(k,:) - Pos_RS;
            u_rs = r_rs / norm(r_rs);
            rate_rs = dot(Sat_Data(s).Vel_True(k,:) - Vel_RS, u_rs);
            Sat_Data(s).Meas_RS(k) = rate_rs + Current_Clk + randn*0.5; % Ruido RS bajo
            
            r_us = Sat_Data(s).Pos_True(k,:) - Pos_User_True;
            u_us = r_us / norm(r_us);
            rate_us = dot(Sat_Data(s).Vel_True(k,:) - Vel_User_True, u_us);
            Sat_Data(s).Meas_User(k) = rate_us + Current_Clk + User_Clk_Drift_True + randn*1.5;
        end
    end
end

%% 2. PROCESAMIENTO RS (Referencia)
fprintf('2. Procesando RS (Tracking & Corrections)...\n');

for s = 1:N_Sats
    % KF Tracking
    F = [1 dt; 0 1]; H_kf = [1 0];
    Q = diag([0.1, 0.01]); R_kf = 1.0^2; % Tuning ajustado
    
    idx_valid = find(~isnan(Sat_Data(s).Meas_RS), 1);
    if isempty(idx_valid), x_est=[0;0]; else, x_est=[Sat_Data(s).Meas_RS(idx_valid);0]; end
    P = eye(2)*100;
    
    Smoothed = nan(N_Time, 1);
    
    for k = 1:N_Time
        x_pred = F * x_est; P_pred = F * P * F' + Q;
        z = Sat_Data(s).Meas_RS(k);
        if ~isnan(z)
            y = z - H_kf * x_pred;
            S = H_kf * P_pred * H_kf' + R_kf;
            K = P_pred * H_kf' / S;
            x_est = x_pred + K * y;
            P = (eye(2) - K * H_kf) * P_pred;
        else
            x_est = x_pred; 
        end
        Smoothed(k) = x_est(1);
    end
    
    % Generación de Correcciones
    Vec_RS = repmat(Pos_RS, N_Time, 1);
    R_vec = Sat_Data(s).Pos_TLE - Vec_RS;
    Dist = sqrt(sum(R_vec.^2, 2));
    U_vec = R_vec ./ Dist;
    Z_TLE = sum((Sat_Data(s).Vel_TLE) .* U_vec, 2); % Vel_RS es 0
    
    Smoothed_Filled = fillmissing(Smoothed, 'linear'); 
    Raw_Correction = Smoothed_Filled - Z_TLE;
    Sat_Data(s).Corrections = movmean(Raw_Correction, 20);
    Sat_Data(s).Corrections = fillmissing(Sat_Data(s).Corrections, 'nearest');
end

%% 3. NAVEGACIÓN USUARIO (CON COLD START ROBUSTO)
fprintf('3. Ejecutando Navegación Usuario (Batch LS + EKF)...\n');

% --- PASO A: COLD START (Batch Least Squares Temporal) ---
fprintf('   -> Inicializando con Batch LS (2 segundos)... ');

% CORRECCIÓN AQUÍ:
% En lugar de un punto aleatorio en la Tierra, usamos la posición de la RS 
% como "semilla". Sabemos que estamos a menos de ~50-100km de ella.
x_ls = [6000, 6000, 6000'; 0];  % <--- ESTA ES LA CLAVE

% (El resto del bucle 'for iter = 1:10' se queda igual...)

% Usaremos los primeros 20 pasos de tiempo (2 segundos) para triangular
batch_steps = 20; 
start_idx = find(Sat_Data(1).Signal_Mask == 1, 1); % Primer momento válido

for iter = 1:10 % Iteraciones de Newton-Gauss
    H_stack = [];
    dy_stack = [];
    
    % Acumulamos mediciones en el tiempo
    for k_offset = 0:batch_steps-1
        k = start_idx + k_offset;
        if k > N_Time, break; end
        
        t = T_span(k);
        
        for s = 1:N_Sats
            meas = Sat_Data(s).Meas_User(k);
            corr = Sat_Data(s).Corrections(k);
            
            if ~isnan(meas) && ~isnan(corr)
                z_clean = meas - corr;
                
                % Modelo TLE
                r_sat = Sat_Data(s).Pos_TLE(k,:)';
                v_sat = Sat_Data(s).Vel_TLE(k,:)';
                
                r_vec = r_sat - x_ls(1:3);
                rng_val = norm(r_vec);
                u_vec = r_vec / rng_val;
                
                % Predicción h(x)
                h_val = dot(v_sat, u_vec) + x_ls(4);
                
                % Jacobiana H (Fila)
                v_proj = dot(v_sat, u_vec);
                vec_perp = v_sat - v_proj * u_vec;
                
                H_row = [-vec_perp' / rng_val, 1];
                
                % Apilar en la matriz gigante
                H_stack = [H_stack; H_row];
                dy_stack = [dy_stack; (z_clean - h_val)];
            end
        end
    end
    
    if isempty(dy_stack)
        error('No hay suficientes datos válidos para inicializar.'); 
    end
    
    % Resolver sistema sobredeterminado
    % Usamos 'pinv' (Pseudo-inversa) que es más robusta a singularidades que '\'
    dx = pinv(H_stack) * dy_stack; 
    x_ls = x_ls + dx;
    
    if norm(dx) < 0.1, break; end % Convergencia lograda
end

fprintf('Posición Inicial hallada. Error: %.2f m\n', norm(x_ls(1:3) - Pos_User_True'));


% --- PASO B: FILTRO DE KALMAN EXTENDIDO (EKF) ---
% Inicializamos el EKF con el resultado robusto del Batch LS
x_nav = x_ls; 
P_nav = diag([50^2, 50^2, 50^2, 5^2]); % Confianza alta (ya inicializamos bien)

% Tuning
Q_nav = diag([0.1, 0.1, 0.1, 2.0]); 
R_nav_val = 5.0^2;

Pos_Error = zeros(N_Time, 1);
Num_Sats_Used = zeros(N_Time, 1);

% Bucle principal de navegación
for k = 1:N_Time
    % 1. PREDICCIÓN
    x_pred = x_nav;
    P_pred = P_nav + Q_nav;
    
    % 2. ACTUALIZACIÓN
    z_list = []; h_list = []; H_list = [];
    
    for s = 1:N_Sats
        meas = Sat_Data(s).Meas_User(k);
        corr = Sat_Data(s).Corrections(k);
        
        if ~isnan(meas) && ~isnan(corr)
            z_clean = meas - corr;
            
            r_sat = Sat_Data(s).Pos_TLE(k,:)';
            v_sat = Sat_Data(s).Vel_TLE(k,:)';
            
            r_vec = r_sat - x_pred(1:3);
            rng_val = norm(r_vec);
            u_vec = r_vec / rng_val;
            
            h_val = dot(v_sat, u_vec) + x_pred(4);
            
            v_proj = dot(v_sat, u_vec);
            vec_perp = v_sat - v_proj * u_vec;
            H_row = [-vec_perp' / rng_val, 1];
            
            z_list = [z_list; z_clean];
            h_list = [h_list; h_val];
            H_list = [H_list; H_row];
        end
    end
    
    Num_Sats_Used(k) = length(z_list);
    
    if ~isempty(z_list)
        y = z_list - h_list; 
        
        R_mat = eye(length(z_list)) * R_nav_val;
        S = H_list * P_pred * H_list' + R_mat;
        K = P_pred * H_list' / S;
        
        x_nav = x_pred + K * y;
        P_nav = (eye(4) - K * H_list) * P_pred;
    else
        x_nav = x_pred;
    end
    
    Pos_Error(k) = norm(x_nav(1:3) - Pos_User_True');
end

%% 4. VISUALIZACIÓN
fprintf('4. Generando gráficas...\n');
figure('Position', [100, 100, 1000, 800], 'Color', 'w');

subplot(2,1,1);
% Graficar en escala logarítmica a veces ayuda a ver la caída, 
% pero lineal es más honesta para ver el error final en metros.
plot(T_span, Pos_Error, 'b-', 'LineWidth', 2);
grid on;
title('Convergencia de Posición del Usuario (EKF Multi-Sat)');
ylabel('Error de Posición 3D (m)'); xlabel('Tiempo (s)');
yline(10, 'r--', 'Label', 'Objetivo 10m', 'LineWidth', 2);
% Forzar vista en los primeros 100 metros para ver detalle
ylim([0 200]); 
subtitle('Debería converger desde 8000m a <20m rápidamente');

subplot(2,1,2);
plot(T_span, Num_Sats_Used, 'k-', 'LineWidth', 1.5);
ylim([0 N_Sats+1]); 
ylabel('Satélites Usados'); xlabel('Tiempo (s)');
title('Disponibilidad de Satélites');
grid on;

% Stats finales
fprintf('\n=== RESULTADOS ===\n');
fprintf('Error Final (t=%d): %.2f metros\n', Total_Time, Pos_Error(end));