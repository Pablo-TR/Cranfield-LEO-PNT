%% ════════════════════════════════════════════════════════════════════════
% ANALYSIS_Phase6C_Tesis.m
% Comparativa de geometrías LEO con Handover forzado a < 40º de Elevación
% ════════════════════════════════════════════════════════════════════════

clear; clc; close all;

fprintf('══════════════════════════════════════════════════════════\n');
fprintf('  PHASE 6C — Sequential Handovers (Strict 40 deg drop)\n');
fprintf('══════════════════════════════════════════════════════════\n');

BASE_SEED   = 1000;
P6_SIGMA_D  = 50.0; 
test_DR     = [0.1, 1.0, 5.0]; 
policies    = {'stickyBest', 'stickyDiverse'};
P6_N_MC     = 10;   % Semillas de Monte Carlo

% Tiempo asegurado para ver el salto de 3 satélites
SIM_DURATION = 1500; 
dt = 0.1;
K_len = SIM_DURATION / dt + 1; 

res_mu_East  = nan(length(policies), length(test_DR), K_len);
res_std_East = nan(length(policies), length(test_DR), K_len);

bo_starts_cell = cell(length(policies), 1);
bo_ends_cell   = cell(length(policies), 1);

max_t_plot = 0; % Para recortar la gráfica al final del 3er pase

for pi_idx = 1:length(policies)
    pol = policies{pi_idx};
    fprintf('\n  --- Evaluando Política: %s ---\n', pol);
    
    for dri = 1:length(test_DR)
        sigDR = test_DR(dri);
        fprintf('    Ejecutando σ_DR = %.1f Hz/s ... ', sigDR);
        
        c = build_base_cfg_local(); 
        c.duration_sec             = SIM_DURATION;
        c.dt                       = dt;
        c.sigma_doppler_Hz         = P6_SIGMA_D;
        c.sigma_dopplerRate_Hzps   = sigDR;
        c.meas.useDoppler          = true;
        c.meas.useDopplerRate      = true;
        c.tropo.enable             = false;
        c.capture.mode             = pol;
        c.capture.switchBlackout_s = 10;
        c.capture.minDwell_s       = 30;
        c.clk.truth.initBias_mps   = 0.10;
        c.P0_clkBias               = (0.15)^2;

        % -- 1. Run reference case to capture tracking --
        c.seedMeas = BASE_SEED + 600;
        out_ref = real_strtg_sequential_EKF_geom(c);
        t_axis = out_ref.t(:)';
        
        if dri == 1
            % Store blackout intervals.
            trk_mask = out_ref.trk.blackout(:)';
            b_start = []; b_end = []; in_bo = false;
            for k = 1:length(trk_mask)
                if trk_mask(k) && ~in_bo
                    b_start(end+1) = t_axis(k); in_bo = true;
                elseif ~trk_mask(k) && in_bo
                    b_end(end+1) = t_axis(k-1); in_bo = false;
                end
            end
            if in_bo, b_end(end+1) = t_axis(end); end
            
            bo_starts_cell{pi_idx} = b_start;
            bo_ends_cell{pi_idx}   = b_end;
            
            % Encontrar instante de muerte del 3er satélite
            sid_v = out_ref.trk.sid(:)';
            switch_epochs = find(diff(sid_v) ~= 0 & sid_v(1:end-1) ~= 0);
            
            if length(switch_epochs) >= 3
                t_cut = t_axis(switch_epochs(3)) + 15; 
            else
                t_cut = t_axis(end);
            end
            max_t_plot = max(max_t_plot, t_cut);
            
            if pi_idx == 1
                sky_data_best = out_ref;
            else
                sky_data_div = out_ref;
            end
        end
        
        % -- 2. Monte Carlo --
        % -- 2. Monte Carlo --
        mc_errE = nan(P6_N_MC, length(t_axis));
        
        c_base = c; 
        
        parfor n = 1:P6_N_MC
            c_local = c_base; 
            c_local.seedMeas = BASE_SEED + 650 + pi_idx*100 + dri*20 + n;
            
            r_full = real_strtg_sequential_EKF_geom(c_local);
            
            mc_errE(n,:) = r_full.errENU(1,:); 
        end
        
        res_mu_East(pi_idx, dri, 1:length(t_axis))  = mean(mc_errE, 1, 'omitnan');
        res_std_East(pi_idx, dri, 1:length(t_axis)) = std(mc_errE, 0, 1, 'omitnan');
        fprintf('Completado.\n');
    end
end

%% ── 3. PLOT ERROR RESULTS ─────────────────────────────────────────────
fh6c_horiz = figure('Name','P6C RAAN Geometry', 'Color','w', 'Position',[50 150 1500 450]);

color_line = {[0.85 0.15 0.15], [0.15 0.65 0.15]}; 
color_band = {[0.95 0.30 0.30], [0.30 0.85 0.30]}; 
color_bo   = {[1.0 0.85 0.85], [0.85 1.0 0.85]};

for dri = 1:length(test_DR)
    subplot(1, 3, dri); hold on; grid on; box on;
    
    for pi_idx = 1:length(policies)
        b_start = bo_starts_cell{pi_idx};
        b_end   = bo_ends_cell{pi_idx};
        for bi = 1:length(b_start)
            if b_start(bi) <= max_t_plot
                patch([b_start(bi), b_end(bi), b_end(bi), b_start(bi)], ...
                      [-1000, -1000, 1000, 1000], ...
                      color_bo{pi_idx}, 'EdgeColor', 'none', ...
                      'FaceAlpha', 0.6, 'HandleVisibility', 'off');
            end
        end
    end
    
    for pi_idx = 1:length(policies)
        mu_curve  = squeeze(res_mu_East(pi_idx, dri, 1:length(t_axis)));
        std_curve = squeeze(res_std_East(pi_idx, dri, 1:length(t_axis)));
        mu_curve = mu_curve(:)'; std_curve = std_curve(:)';
        
        fill([t_axis, fliplr(t_axis)], [mu_curve + std_curve, fliplr(mu_curve - std_curve)], ...
             color_band{pi_idx}, 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off'); 
        plot(t_axis, mu_curve, '-', 'Color', color_line{pi_idx}, 'LineWidth', 2.0, ...
             'DisplayName', sprintf('Policy: %s', policies{pi_idx}));
    end
    
    yline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    ylim([-50, 500]); xlim([0, max_t_plot]);
    xlabel('Time (s)', 'FontWeight', 'bold');
    
    if dri == 1
        ylabel('East Error (m)', 'FontWeight', 'bold');
        legend('Location', 'northeast');
    end
    title(sprintf('Noise \\sigma_{DR} = %.1f Hz/s', test_DR(dri)), 'FontSize', 12);
end
sgtitle({'Phase 6C: Cross-Track Error Recovery (3 Satellites Pass)', 'Switch strictly triggered when Elev < 40^{\circ}'}, 'FontSize', 15, 'FontWeight', 'bold');

%% ── 4. PLOT THREE-SATELLITE SKYPLOT ──────────────────────────────────
fh_sky = figure('Name', 'Skyview Comparison', 'Color', 'w', 'Position', [100 100 1100 550]);
draw_sky_local(sky_data_best, 'Política: stickyBest (Same Plane)', max_t_plot);
subplot(1, 2, 2, polaraxes); 
draw_sky_local(sky_data_div, 'Política: stickyDiverse (Crossing RAAN)', max_t_plot);
sgtitle('Geometric Diversity: 3-Satellite Forced Handover', 'FontSize', 15, 'FontWeight', 'bold');


%% ═══════════════════════ FUNCIONES DEL SCRIPT ════════════════════════════
function draw_sky_local(out_data, title_str, t_limit)
    pax = gca; hold on;
    pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';    
    pax.RLim = [0 90]; pax.RTick = [0 30 50 60 90];
    pax.RTickLabel = {'90° (Zenith)','60°','40° (Drop)','30°','0° (Horizon)'}; 
    pax.GridAlpha = 0.3;
    
    th = linspace(0, 2*pi, 100);
    polarplot(th, ones(size(th))*50, 'r--', 'LineWidth', 1.5, 'DisplayName', '40° Handover Boundary');
    
    t_ax = out_data.t(:)';
    valid_time_idx = find(t_ax <= t_limit);
    trk_sid = out_data.trk.sid(valid_time_idx)';
    tracked_sats = unique(trk_sid(trk_sid > 0));
    colors = lines(length(tracked_sats));
    
    for i = 1:length(tracked_sats)
        s = tracked_sats(i);
        az_rad = deg2rad(out_data.sat.az(valid_time_idx, s));
        el_deg = out_data.sat.elev(valid_time_idx, s); 
        r_plot = 90 - el_deg;
        valid = el_deg > 0;
        polarplot(az_rad(valid), r_plot(valid), '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1.0, 'HandleVisibility', 'off');
        
        active_idx = find(trk_sid == s);
        if ~isempty(active_idx)
            polarplot(az_rad(active_idx), r_plot(active_idx), '.', 'Color', colors(i,:), 'MarkerSize', 10, 'HandleVisibility', 'off');
            start_idx = active_idx(1);
            polarplot(az_rad(start_idx), r_plot(start_idx), 'p', 'MarkerSize', 14, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('Sat %d Start', s));
        end
    end
    title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'bestoutside');
end

function cfg = build_base_cfg_local()
    cfg.tleFile = fullfile('InputData', 'oneweb.tle');
    cfg.duration_sec  = 100; cfg.dt = 0.1;
    cfg.startUTC = []; cfg.stopUTC = []; cfg.startFromTLE = 'first'; cfg.maxTLEAgeDays = 7;
    cfg.rx.lat_deg = 52.06869; cfg.rx.lon_deg = -0.60884; cfg.rx.alt_m = 100;
    cfg.rx.vel_ecef = [0;0;0]; cfg.rx.acc_ecef = [0;0;0];
    cfg.fc_Hz  = 11.325e9; cfg.c = 299792458; cfg.lambda = cfg.c / cfg.fc_Hz;
    cfg.minElevDeg = 30; cfg.coarseStep = 10; cfg.minMaxElev = 30; cfg.maxSats = 12; cfg.constellationSelectionWindow_s = 120;
    cfg.useSGP4 = false; cfg.allowJ2Fallback = true;
    cfg.init.errENU_m = [400; 400; 0]; cfg.init.sigENU_m = [500; 500; 10];
    cfg.P0_clkBias  = (0.1*2)^2; cfg.P0_clkDrift = (2.5e-7*2)^2;
    cfg.sigma_pos_rw = 1e-4;
    cfg.clk.truth.initBias_mps = 0.10; cfg.clk.truth.initDrift_mps2 = 5e-5;
    cfg.clk.truth.sigmaBias_rw_mps_sqrt = 5e-3 * 2; cfg.clk.truth.sigmaDrift_rw_mps2_sqrt = 1.0e-4 * 2;
    cfg.clk.filt.sigmaBias_rw_mps_sqrt = 1.0e-2; cfg.clk.filt.sigmaDrift_rw_mps2_sqrt = 5.0e-4;
    cfg.clk.filter_init_bias_mps = 0.0; cfg.clk.filter_init_drift_mps2 = 0.0;
    cfg.tropo.enable = false; cfg.tropo.el_floor_deg = 5; cfg.tropo.map_epsilon_deg = 0.5;
    cfg.sigma_doppler_Hz = 2.0; cfg.sigma_dopplerRate_Hzps = 0.5;
    cfg.meas.enable_lighttime = false; cfg.meas.colored_tau_s = 0.5; cfg.meas.sigma_col0_mps = 0.03;
    cfg.meas.Rinflate = 2.0; cfg.meas.useDoppler = true; cfg.meas.useDopplerRate = false; cfg.meas.fdJacStep_m = 0.10;
    cfg.filter.useNISGate = true; cfg.filter.nisGate = 16.0; cfg.filter.maxIEKFIter = 10; cfg.filter.iekfTol = 1e-5;
    cfg.filter.maxStep_m = 25; cfg.filter.startupInflateEpochs = 50; cfg.filter.startupRmult = 25;
    cfg.height.enable = true; cfg.height.sigma_m = 10.0; cfg.height.z0_m = 0.0;
    cfg.capture.mode = 'stickyBest'; cfg.capture.minDwell_s = 30; cfg.capture.switchBlackout_s = 3; cfg.capture.elev_margin_deg = 20;
    cfg.seedMeas = 42;
end
