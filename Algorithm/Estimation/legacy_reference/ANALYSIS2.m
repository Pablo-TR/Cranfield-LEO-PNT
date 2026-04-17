%% ═══════════════════════════════════════════════════════════════════════════
%  RUN_FULL_ANALYSIS.m
%
%  PURPOSE:
%    Automated Monte Carlo simulation campaign for the LEO Doppler Sequential
%    IEKF (real_strtg_sequential_EKF_fastrealistic.m).
%
%    Executes five analysis phases, each averaging over N_MC random seeds to
%    isolate geometry effects from noise variance.
%
%  PHASES:
%    Phase 1 — Baseline Doppler-only & coherent clock oscillator tuning
%    Phase 2 — Physical realism: tropospheric mismatch bleed-through
%    Phase 3 — Observability degradation: Doppler noise sweep
%    Phase 4 — Doppler-rate measurements
%    Phase 5 — Deep IEKF parameter optimisation (Pareto front)
%
%  KEY DESIGN RULES (enforced throughout):
%    • cfg.duration_sec = 100 s  &  cfg.dt = 0.01 s  — fixed timeline
%    • Vertical / Up axis is EXCLUDED from all metrics (height constrained)
%    • All KPIs separate East, North, and joint 2D Horizontal error
%    • Steady-state RMS = last SS_WINDOW_S seconds
%    • Convergence = first time 2D Horiz < CONV_THR_M  AND stays there
%
%  USAGE:
%    Edit GLOBAL SETTINGS below, then run the full script or individual %%
%    sections independently.
%
%  DEPENDENCIES:
%    real_strtg_sequential_EKF_fastrealistic.m   (solver — must be on path)
%
%  OUTPUT FIELDS USED FROM SOLVER:
%    out.t                    — time vector [K×1 s]
%    out.errENU               — ENU error   [3×K m]   rows: [East;North;Up]
%    out.xhist                — state hist  [5×K]    cols: [dE dN dU CB CD]
%    out.truthAux.clkBiasTrue — true clk bias  [K×1 m/s]
%    out.truthAux.clkDriftTrue— true clk drift [K×1 m/s²]
%    out.nUpdates             — meas updates   [K×1]
%
% ═══════════════════════════════════════════════════════════════════════════

clear; close all; clc;

%% ══════════════════════════════════════════════════════════════════════════
%  §0  GLOBAL SETTINGS
% ══════════════════════════════════════════════════════════════════════════

% ── Monte Carlo ──────────────────────────────────────────────────────────
N_MC           = 20;         % number of random seeds per configuration
BASE_SEED      = 1000;       % seeds will be BASE_SEED + 1 … BASE_SEED + N_MC

% ── Fixed timeline (enforced in every run) ───────────────────────────────
FIXED_DURATION = 100;        % [s]
FIXED_DT       = 0.01;       % [s]

% ── KPI thresholds ───────────────────────────────────────────────────────
CONV_THR_M     = 50;         % 2D Horizontal convergence threshold [m]
SS_WINDOW_S    = 15;         % steady-state RMS window (last N seconds)

% ── Which phases o execute (set to [] to skip) ──────────────────────────
RUN_PHASES     = [8];

% ── Parallel computing ───────────────────────────────────────────────────
USE_PARFOR     = false;      % set true to enable parfor (requires PCT)
%   NOTE: If Parallel Computing Toolbox unavailable, the script falls back
%         to a standard for-loop automatically.
if USE_PARFOR
    try
        parpool('local');
    catch ME
        warning('RUN_FULL_ANALYSIS:noParpool', ...
            'Could not open parpool (%s). Falling back to serial.', ME.message);
        USE_PARFOR = false;
    end
end

% ── Figure output directory ───────────────────────────────────────────────
FIG_DIR = 'analysis_figures';
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║  RUN_FULL_ANALYSIS  —  LEO Doppler IEKF Monte Carlo Suite   ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n');
fprintf('  MC seeds per config : %d\n', N_MC);
fprintf('  Fixed timeline      : %g s @ dt=%g s\n', FIXED_DURATION, FIXED_DT);
fprintf('  Convergence thr.    : %.0f m (2D Horizontal)\n', CONV_THR_M);
fprintf('  Steady-state window : %.0f s\n\n', SS_WINDOW_S);



%% ══════════════════════════════════════════════════════════════════════════
%  §P1  PHASE 1 — Baseline Doppler-Only & Coherent Clock Analysis
% ══════════════════════════════════════════════════════════════════════════
%
%  Config: useDoppler=true, useDopplerRate=false, tropo=false, σ_D=2 Hz
%
%  Experiment: Compare three oscillator regimes (mild / nominal / harsh)
%    defined by the range of initBias_mps injected as truth.
%
%  For each regime, the filter initial covariance P0_clkBias and process
%    noise are set coherently with the truth-clock scale.
%
%  Each MC seed draws a RANDOM truth bias from the regime range, then
%    configures the filter P0 and process noise accordingly.
%
% ══════════════════════════════════════════════════════════════════════════

if ismember(1, RUN_PHASES)
    fprintf('══════════════════════════════════════════════════════════\n');
    fprintf('  PHASE 1 — Doppler-Only & Coherent Clock Analysis\n');
    fprintf('══════════════════════════════════════════════════════════\n');

    % ── Define oscillator regimes ──────────────────────────────────────────
    %    Each entry: [bias_lo, bias_hi, truth_sigBias, truth_sigDrift]
    %    (sigBias/Drift are the truth clock random-walk stds)
    regimes = struct( ...
        'name',    {'Mild (0.02-0.10 m/s)', 'Nominal (0.10-0.30 m/s)', 'Harsh (0.30-1.0 m/s)'}, ...
        'bias_lo', {0.02,  0.10,  0.30}, ...
        'bias_hi', {0.10,  0.30,  1.00}, ...
        'truth_sigBias',  {5e-3,  1e-2,  3e-2},  ...  % truth random-walk std [m/s/√s]
        'truth_sigDrift', {1e-4,  3e-4,  1e-3},  ...  % truth random-walk std [m/s²/√s]
        'color',   {[0.1 0.7 0.2], [0.9 0.6 0.1], [0.85 0.1 0.1]} ...
    );
    nRegimes = numel(regimes);

    % Storage for KPI summary
    p1_kpi = cell(nRegimes, 1);

    % Storage for averaged timeseries (for plotting)
    p1_ts  = cell(nRegimes, 1);

    for ri = 1:nRegimes
        rg = regimes(ri);
        fprintf('  Regime: %s\n', rg.name);

        % Build per-seed config list (each seed gets a random truth bias
        % drawn from the regime range; the filter P0 and Q are tuned to match)
        cfgs_mc = cell(N_MC, 1);
        rng(BASE_SEED + ri * 100);   % reproducible draw
        truthBiases = rg.bias_lo + (rg.bias_hi - rg.bias_lo) * rand(N_MC, 1);

        for n = 1:N_MC
            c = build_base_cfg();

            % ── Fixed timeline ──────────────────────────────────────────
            c.duration_sec = FIXED_DURATION;
            c.dt           = FIXED_DT;

            % ── Phase 1 measurement config ──────────────────────────────
            c.sigma_doppler_Hz       = 2.0;
            c.meas.useDoppler        = true;
            c.meas.useDopplerRate    = false;
            c.tropo.enable           = false;

            % ── TRUTH clock injection ───────────────────────────────────
            c.clk.truth.initBias_mps            = truthBiases(n);
            c.clk.truth.initDrift_mps2          = 5e-5;
            c.clk.truth.sigmaBias_rw_mps_sqrt   = rg.truth_sigBias;
            c.clk.truth.sigmaDrift_rw_mps2_sqrt = rg.truth_sigDrift;

            % ── COHERENT filter tuning ──────────────────────────────────
            % Rule: P0_clkBias covers max regime bias with ≥ 3-sigma margin
            %       Filter process noise is 1.5× truth (slight overmodelling)
            P0_bias_sig = rg.bias_hi * 1.5;    % 1-sigma initial uncertainty [m/s]
            c.P0_clkBias  = P0_bias_sig^2;      % filter initial variance
            c.P0_clkDrift = (rg.truth_sigDrift * 10 * 1.5)^2; % generous drift prior

            c.clk.filt.sigmaBias_rw_mps_sqrt    = rg.truth_sigBias  * 1.5;
            c.clk.filt.sigmaDrift_rw_mps2_sqrt  = rg.truth_sigDrift * 1.5;

            cfgs_mc{n} = c;
        end

        % Run Monte Carlo (serial or parallel)
        [mc_errE, mc_errN, mc_err2D, mc_errCB, mc_errCD] = ...
            run_mc(cfgs_mc, N_MC, BASE_SEED + ri*100, USE_PARFOR);

        % Extract time axis (all runs share the same t)
        tmp     = run_single(cfgs_mc{1});
        t_axis  = tmp.t;

        % Compute KPIs
        p1_kpi{ri} = compute_kpis(t_axis, mc_errE, mc_errN, mc_err2D, ...
                                   SS_WINDOW_S, CONV_THR_M);

        % Store mean ± std for plotting
        p1_ts{ri}.t       = t_axis;
        p1_ts{ri}.meanE   = mean(mc_errE,  1,'omitnan');
        p1_ts{ri}.stdE    = std(mc_errE,   0,1,'omitnan');
        p1_ts{ri}.meanN   = mean(mc_errN,  1,'omitnan');
        p1_ts{ri}.stdN    = std(mc_errN,   0,1,'omitnan');
        p1_ts{ri}.mean2D  = mean(mc_err2D, 1,'omitnan');
        p1_ts{ri}.std2D   = std(mc_err2D,  0,1,'omitnan');
        p1_ts{ri}.meanCB  = mean(mc_errCB, 1,'omitnan');
        p1_ts{ri}.stdCB   = std(mc_errCB,  0,1,'omitnan');
        p1_ts{ri}.meanCD  = mean(mc_errCD, 1,'omitnan');
        p1_ts{ri}.stdCD   = std(mc_errCD,  0,1,'omitnan');

        fprintf('    SS East RMS : %.2f m | SS North RMS : %.2f m | SS 2D : %.2f m\n', ...
            p1_kpi{ri}.mean_ss_rmsE, p1_kpi{ri}.mean_ss_rmsN, p1_kpi{ri}.mean_ss_rms2D);
        fprintf('    Conv. time  : %.1f s (%.0f%% converged)\n', ...
            p1_kpi{ri}.mean_tc, 100*p1_kpi{ri}.frac_converged);
    end

    % ── PLOT P1-A: East & North errors per regime ─────────────────────────
    fh1a = figure('Name','P1 East North Errors','Color','w','Position',[50 50 1200 500]);
    ax_titles = {'East Error (m)', 'North Error (m)'};
    fields_mu  = {'meanE',  'meanN'};
    fields_std = {'stdE',   'stdN'};
    for ax_i = 1:2
        subplot(1,2,ax_i); hold on; grid on; box on;
        for ri = 1:nRegimes
            ts  = p1_ts{ri};
            mu  = ts.(fields_mu{ax_i});
            sig = ts.(fields_std{ax_i});
            col = regimes(ri).color;
            fill([ts.t, fliplr(ts.t)], [mu+sig, fliplr(mu-sig)], col, ...
                'FaceAlpha', 0.2, 'EdgeColor','none');
            plot(ts.t, mu, '-', 'Color', col, 'LineWidth', 2, ...
                'DisplayName', regimes(ri).name);
        end
        xlabel('Time (s)'); ylabel(ax_titles{ax_i});
        title(ax_titles{ax_i});
        legend('Location','best','FontSize',8);
        yline(0,'k:','LineWidth',0.8);
    end
    sgtitle('Phase 1 — Doppler-Only: East & North Errors by Clock Regime', ...
            'FontWeight','bold');
    save_fig(fh1a, 'P1A_EastNorth_Errors', FIG_DIR);

    % ── PLOT P1-B: 2D Horizontal error per regime ─────────────────────────
    fh1b = figure('Name','P1 2D Horizontal','Color','w','Position',[50 600 800 400]);
    hold on; grid on; box on;
    for ri = 1:nRegimes
        ts  = p1_ts{ri};
        col = regimes(ri).color;
        fill([ts.t, fliplr(ts.t)], ...
             [ts.mean2D+ts.std2D, fliplr(ts.mean2D-ts.std2D)], col, ...
             'FaceAlpha', 0.2, 'EdgeColor','none');
        plot(ts.t, ts.mean2D, '-', 'Color', col, 'LineWidth', 2, ...
            'DisplayName', regimes(ri).name);
    end
    yline(CONV_THR_M, 'k--', 'LineWidth', 1.2, 'DisplayName', ...
          sprintf('%.0f m threshold', CONV_THR_M));
    xlabel('Time (s)'); ylabel('2D Horizontal Error (m)');
    title('Phase 1 — 2D Horizontal Error by Clock Regime');
    legend('Location','best','FontSize',8);
    save_fig(fh1b, 'P1B_2D_Horizontal', FIG_DIR);

    % ── PLOT P1-C: Clock bias estimation error ────────────────────────────
    fh1c = figure('Name','P1 Clock Bias Error','Color','w','Position',[900 50 800 400]);
    hold on; grid on; box on;
    for ri = 1:nRegimes
        ts  = p1_ts{ri};
        col = regimes(ri).color;
        fill([ts.t, fliplr(ts.t)], ...
             [ts.meanCB+ts.stdCB, fliplr(ts.meanCB-ts.stdCB)], col, ...
             'FaceAlpha', 0.2, 'EdgeColor','none');
        plot(ts.t, ts.meanCB, '-', 'Color', col, 'LineWidth', 2, ...
            'DisplayName', regimes(ri).name);
    end
    yline(0,'k:','LineWidth',0.8);
    xlabel('Time (s)'); ylabel('Clock Bias Error (m/s)');
    title('Phase 1 — Clock Bias Estimation Error by Regime');
    legend('Location','best','FontSize',8);
    save_fig(fh1c, 'P1C_ClockBias_Error', FIG_DIR);

    % ── PLOT P1-D: Clock drift estimation error ───────────────────────────
    fh1d = figure('Name','P1 Clock Drift Error','Color','w','Position',[900 500 800 400]);
    hold on; grid on; box on;
    for ri = 1:nRegimes
        ts  = p1_ts{ri};
        col = regimes(ri).color;
        fill([ts.t, fliplr(ts.t)], ...
             [ts.meanCD+ts.stdCD, fliplr(ts.meanCD-ts.stdCD)], col, ...
             'FaceAlpha', 0.2, 'EdgeColor','none');
        plot(ts.t, ts.meanCD, '-', 'Color', col, 'LineWidth', 2, ...
            'DisplayName', regimes(ri).name);
    end
    yline(0,'k:','LineWidth',0.8);
    xlabel('Time (s)'); ylabel('Clock Drift Error (m/s²)');
    title('Phase 1 — Clock Drift Estimation Error by Regime');
    legend('Location','best','FontSize',8);
    save_fig(fh1d, 'P1D_ClockDrift_Error', FIG_DIR);

    fprintf('  Phase 1 complete.\n\n');
end

%% ══════════════════════════════════════════════════════════════════════════
%  §P2  PHASE 2 — Physical Realism: Tropospheric Mismatch
% ══════════════════════════════════════════════════════════════════════════
%
%  Config: Low noise (σ_D = 2 Hz), nominal clock, tropo.enable = true.
%  Truth: ztd_truth_m = 2.5 m  |  Filter model: ztd_model_m = 2.3 m
%  → 0.2 m ZTD mismatch
%
%  Question: Does the unmodelled delay bleed into the North axis, or is it
%            entirely absorbed by the clkBias state?
%
% ══════════════════════════════════════════════════════════════════════════

if ismember(2, RUN_PHASES)
    fprintf('══════════════════════════════════════════════════════════\n');
    fprintf('  PHASE 2 — Physical Realism: Tropospheric Mismatch\n');
    fprintf('══════════════════════════════════════════════════════════\n');

    % ── Build two configs: matched vs. mismatched ─────────────────────────
    p2_configs = struct( ...
        'label',    {'No Troposphere', 'Tropo Matched', 'Tropo Mismatched (Δ=0.2m)'}, ...
        'enable',   {false, true, true}, ...
        'ztd_truth',{2.5, 2.5, 2.5}, ...
        'ztd_model',{2.5, 2.5, 2.3}, ...
        'color',    {[0.5 0.5 0.5], [0.1 0.6 0.9], [0.9 0.2 0.2]} ...
    );
    nP2 = numel(p2_configs);

    p2_kpi = cell(nP2,1);
    p2_ts  = cell(nP2,1);

    for ci = 1:nP2
        pc = p2_configs(ci);
        fprintf('  Config: %s\n', pc.label);

        c = build_base_cfg();
        c.duration_sec           = FIXED_DURATION;
        c.dt                     = FIXED_DT;
        c.sigma_doppler_Hz       = 2.0;
        c.meas.useDoppler        = true;
        c.meas.useDopplerRate    = false;

        % Nominal clock (coherent with baseline)
        c.clk.truth.initBias_mps            = 0.10;
        c.clk.filt.sigmaBias_rw_mps_sqrt    = 1.0e-2;
        c.clk.filt.sigmaDrift_rw_mps2_sqrt  = 5.0e-4;
        c.P0_clkBias  = (0.15)^2;
        c.P0_clkDrift = (2.5e-7*2)^2;

        % Troposphere settings
        c.tropo.enable      = pc.enable;
        c.tropo.ztd_truth_m = pc.ztd_truth;
        c.tropo.ztd_model_m = pc.ztd_model;

        [mc_errE, mc_errN, mc_err2D, mc_errCB, ~] = ...
            run_mc(c, N_MC, BASE_SEED + 200 + ci*10, USE_PARFOR);

        tmp    = run_single(c);
        t_axis = tmp.t;

        p2_kpi{ci} = compute_kpis(t_axis, mc_errE, mc_errN, mc_err2D, ...
                                   SS_WINDOW_S, CONV_THR_M);

        p2_ts{ci}.t      = t_axis;
        p2_ts{ci}.meanE  = mean(mc_errE,  1,'omitnan');
        p2_ts{ci}.stdE   = std(mc_errE,   0,1,'omitnan');
        p2_ts{ci}.meanN  = mean(mc_errN,  1,'omitnan');
        p2_ts{ci}.stdN   = std(mc_errN,   0,1,'omitnan');
        p2_ts{ci}.mean2D = mean(mc_err2D, 1,'omitnan');
        p2_ts{ci}.std2D  = std(mc_err2D,  0,1,'omitnan');
        p2_ts{ci}.meanCB = mean(mc_errCB, 1,'omitnan');
        p2_ts{ci}.stdCB  = std(mc_errCB,  0,1,'omitnan');

        fprintf('    SS 2D RMS  : %.2f m  |  SS North RMS : %.2f m\n', ...
            p2_kpi{ci}.mean_ss_rms2D, p2_kpi{ci}.mean_ss_rmsN);
        fprintf('    Conv. time : %.1f s\n', p2_kpi{ci}.mean_tc);
    end

    % ── PLOT P2-A: 2D Horizontal comparison ──────────────────────────────
    fh2a = figure('Name','P2 2D Horizontal Tropo','Color','w','Position',[50 50 900 420]);
    hold on; grid on; box on;
    for ci = 1:nP2
        ts  = p2_ts{ci};
        col = p2_configs(ci).color;
        fill([ts.t, fliplr(ts.t)], ...
             [ts.mean2D+ts.std2D, fliplr(ts.mean2D-ts.std2D)], col, ...
             'FaceAlpha', 0.2, 'EdgeColor','none');
        plot(ts.t, ts.mean2D, '-', 'Color', col, 'LineWidth', 2, ...
            'DisplayName', p2_configs(ci).label);
    end
    yline(CONV_THR_M,'k--','LineWidth',1.2,'DisplayName','50 m threshold');
    xlabel('Time (s)'); ylabel('2D Horizontal Error (m)');
    title('Phase 2 — Tropospheric Mismatch: 2D Horizontal Error');
    legend('Location','best','FontSize',8);
    save_fig(fh2a, 'P2A_2D_Horizontal_Tropo', FIG_DIR);

    % ── PLOT P2-B: North error (bleed-through check) ──────────────────────
    fh2b = figure('Name','P2 North Error Tropo','Color','w','Position',[50 520 900 420]);
    hold on; grid on; box on;
    for ci = 1:nP2
        ts  = p2_ts{ci};
        col = p2_configs(ci).color;
        fill([ts.t, fliplr(ts.t)], ...
             [ts.meanN+ts.stdN, fliplr(ts.meanN-ts.stdN)], col, ...
             'FaceAlpha', 0.2, 'EdgeColor','none');
        plot(ts.t, ts.meanN, '-', 'Color', col, 'LineWidth', 2, ...
            'DisplayName', p2_configs(ci).label);
    end
    yline(0,'k:','LineWidth',0.8);
    xlabel('Time (s)'); ylabel('North Error (m)');
    title('Phase 2 — Tropospheric Mismatch: North Error (Bleed-through Check)');
    legend('Location','best','FontSize',8);
    save_fig(fh2b, 'P2B_North_Tropo_Bleed', FIG_DIR);

    % ── PLOT P2-C: Clock Bias error (absorption check) ───────────────────
    fh2c = figure('Name','P2 ClkBias Tropo','Color','w','Position',[1000 50 900 420]);
    hold on; grid on; box on;
    for ci = 1:nP2
        ts  = p2_ts{ci};
        col = p2_configs(ci).color;
        fill([ts.t, fliplr(ts.t)], ...
             [ts.meanCB+ts.stdCB, fliplr(ts.meanCB-ts.stdCB)], col, ...
             'FaceAlpha', 0.2, 'EdgeColor','none');
        plot(ts.t, ts.meanCB, '-', 'Color', col, 'LineWidth', 2, ...
            'DisplayName', p2_configs(ci).label);
    end
    yline(0,'k:','LineWidth',0.8);
    xlabel('Time (s)'); ylabel('Clock Bias Error (m/s)');
    title('Phase 2 — Tropospheric Mismatch: Clock Bias Absorption Check');
    legend('Location','best','FontSize',8);
    save_fig(fh2c, 'P2C_ClkBias_Tropo', FIG_DIR);

    % Print bleed-through analysis
    fprintf('\n  ── Bleed-Through Analysis ─────────────────────────────\n');
    fprintf('  %-40s  %-10s  %-12s  %-12s\n', ...
        'Config','SS 2D (m)','SS North (m)','SS CB err (m/s)');
    for ci = 1:nP2
        fprintf('  %-40s  %-10.2f  %-12.2f\n', ...
            p2_configs(ci).label, ...
            p2_kpi{ci}.mean_ss_rms2D, ...
            p2_kpi{ci}.mean_ss_rmsN);
    end
    fprintf('\n  Phase 2 complete.\n\n');
end

%% ══════════════════════════════════════════════════════════════════════════
%  §P3  PHASE 3 — Observability Degradation: Doppler Noise Sweep
% ══════════════════════════════════════════════════════════════════════════
%
%  Config: useDopplerRate = false.
%  Sweep σ_D across [2, 10, 50, 100, 150, 200] Hz.
%
%  Outputs:
%    1. SS Error vs σ_D — separate East and North curves
%    2. Convergence Time vs σ_D — showing cross-track observability collapse
%
% ══════════════════════════════════════════════════════════════════════════

if ismember(3, RUN_PHASES)
    fprintf('══════════════════════════════════════════════════════════\n');
    fprintf('  PHASE 3 — Observability Degradation: Doppler Noise Sweep\n');
    fprintf('══════════════════════════════════════════════════════════\n');

    sigma_D_sweep = [2, 10, 50, 100, 150, 200];   % [Hz]
    nSweep3 = numel(sigma_D_sweep);

    p3_ssE  = nan(nSweep3,1);
    p3_ssN  = nan(nSweep3,1);
    p3_ss2D = nan(nSweep3,1);
    p3_tc   = nan(nSweep3,1);

    for si = 1:nSweep3
        sigD = sigma_D_sweep(si);
        fprintf('  σ_D = %4.0f Hz ... ', sigD);

        c = build_base_cfg();
        c.duration_sec         = FIXED_DURATION;
        c.dt                   = FIXED_DT;
        c.sigma_doppler_Hz     = sigD;
        c.meas.useDoppler      = true;
        c.meas.useDopplerRate  = false;
        c.tropo.enable         = false;

        % Nominal clock (coherent)
        c.clk.truth.initBias_mps           = 0.10;
        c.P0_clkBias                       = (0.15)^2;
        c.clk.filt.sigmaBias_rw_mps_sqrt   = 1.0e-2;
        c.clk.filt.sigmaDrift_rw_mps2_sqrt = 5.0e-4;

        [mc_errE, mc_errN, mc_err2D, ~, ~] = ...
            run_mc(c, N_MC, BASE_SEED + 300 + si*5, USE_PARFOR);

        tmp    = run_single(c);
        t_axis = tmp.t;

        kpi = compute_kpis(t_axis, mc_errE, mc_errN, mc_err2D, ...
                           SS_WINDOW_S, CONV_THR_M);

        p3_ssE(si)  = kpi.mean_ss_rmsE;
        p3_ssN(si)  = kpi.mean_ss_rmsN;
        p3_ss2D(si) = kpi.mean_ss_rms2D;
        p3_tc(si)   = kpi.mean_tc;

        fprintf('SS East=%.1f m | SS North=%.1f m | tc=%.1f s\n', ...
            kpi.mean_ss_rmsE, kpi.mean_ss_rmsN, kpi.mean_tc);
    end

    % ── PLOT P3-A: SS Error vs Doppler Noise ─────────────────────────────
    fh3a = figure('Name','P3 SS Error vs Doppler Noise','Color','w','Position',[50 50 800 450]);
    hold on; grid on; box on;
    plot(sigma_D_sweep, p3_ssE,  'bo-', 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', 'East SS RMS');
    plot(sigma_D_sweep, p3_ssN,  'rs-', 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', 'North SS RMS');
    plot(sigma_D_sweep, p3_ss2D, 'kd-', 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', '2D Horiz SS RMS');
    set(gca,'XScale','log','YScale','log');
    xlabel('Doppler Noise σ_D (Hz)'); ylabel('Steady-State RMS Error (m)');
    title('Phase 3 — Steady-State Error vs Doppler Noise');
    legend('Location','northwest','FontSize',9);
    save_fig(fh3a, 'P3A_SS_Error_vs_DopplerNoise', FIG_DIR);

    % ── PLOT P3-B: Convergence Time vs Doppler Noise ──────────────────────
    fh3b = figure('Name','P3 Conv Time vs Doppler Noise','Color','w','Position',[900 50 800 450]);
    hold on; grid on; box on;
    bar_colors = [0.2 0.4 0.8; 0.3 0.6 0.9; 0.5 0.7 0.95; ...
                  0.8 0.3 0.2; 0.9 0.1 0.1; 0.6 0.0 0.0];
    b = bar(1:nSweep3, p3_tc, 'FaceColor','flat');
    for si = 1:nSweep3
        b.CData(si,:) = bar_colors(si,:);
    end
    set(gca,'XTickLabel', arrayfun(@(x) sprintf('%d Hz',x), sigma_D_sweep, ...
        'UniformOutput', false));
    xlabel('Doppler Noise σ_D (Hz)'); ylabel('Convergence Time (s)');
    title(sprintf('Phase 3 — Convergence Time vs Doppler Noise\n(2D Horiz < %.0f m, strict)', ...
        CONV_THR_M));
    % Annotate infinite convergence times
    for si = 1:nSweep3
        if isinf(p3_tc(si))
            text(si, FIXED_DURATION * 0.95, '∞', 'HorizontalAlignment','center', ...
                'FontSize',14, 'Color','red', 'FontWeight','bold');
        end
    end
    save_fig(fh3b, 'P3B_ConvTime_vs_DopplerNoise', FIG_DIR);

    fprintf('  Phase 3 complete.\n\n');
end

%% ══════════════════════════════════════════════════════════════════════════
%  §P4  PHASE 4 — The Rescue: Adding Doppler-Rate Measurements
% ══════════════════════════════════════════════════════════════════════════
%
%  Config: Fix σ_D = 100 Hz. useDopplerRate = true.
%  Experiment A: Time-series — Doppler-only vs Doppler+Rate at σ_DR = 0.5 Hz/s
%  Experiment B: 2D heatmaps — SS error and Convergence vs (σ_D, σ_DR)
%
% ══════════════════════════════════════════════════════════════════════════

if ismember(4, RUN_PHASES)
    fprintf('══════════════════════════════════════════════════════════\n');
    fprintf('  PHASE 4 — The Rescue: Doppler-Rate Measurements\n');
    fprintf('══════════════════════════════════════════════════════════\n');

    CRIT_SIGMA_D = 100;   % critical Doppler noise [Hz]

    % ── Experiment 4A: Time-series comparison ─────────────────────────────
    fprintf('  4A: Time-series — Doppler-only vs Doppler+Rate\n');

    % Doppler-only at 100 Hz
    c_donly = build_base_cfg();
    c_donly.duration_sec         = FIXED_DURATION;
    c_donly.dt                   = FIXED_DT;
    c_donly.sigma_doppler_Hz     = CRIT_SIGMA_D;
    c_donly.meas.useDoppler      = true;
    c_donly.meas.useDopplerRate  = false;
    c_donly.tropo.enable         = false;
    c_donly.clk.truth.initBias_mps = 0.10;
    c_donly.P0_clkBias           = (0.15)^2;

    % Doppler + Rate at 100 Hz / 0.5 Hz/s
    c_drate = build_base_cfg();
    c_drate.duration_sec         = FIXED_DURATION;
    c_drate.dt                   = FIXED_DT;
    c_drate.sigma_doppler_Hz     = CRIT_SIGMA_D;
    c_drate.sigma_dopplerRate_Hzps = 0.5;
    c_drate.meas.useDoppler      = true;
    c_drate.meas.useDopplerRate  = true;
    c_drate.tropo.enable         = false;
    c_drate.clk.truth.initBias_mps = 0.10;
    c_drate.P0_clkBias           = (0.15)^2;

    [mc_E_do, ~, mc_2D_do, ~, ~] = ...
        run_mc(c_donly, N_MC, BASE_SEED + 400, USE_PARFOR);
    [mc_E_dr, ~, mc_2D_dr, ~, ~] = ...
        run_mc(c_drate, N_MC, BASE_SEED + 410, USE_PARFOR);

    tmp    = run_single(c_donly);
    t_axis = tmp.t;

    mu_E_do  = mean(mc_E_do,  1,'omitnan');
    std_E_do = std(mc_E_do,   0,1,'omitnan');
    mu_E_dr  = mean(mc_E_dr,  1,'omitnan');
    std_E_dr = std(mc_E_dr,   0,1,'omitnan');
    mu_2D_do = mean(mc_2D_do, 1,'omitnan');
    mu_2D_dr = mean(mc_2D_dr, 1,'omitnan');

    % Plot P4A-i: East error comparison
    fh4a1 = figure('Name','P4A East Error Comparison','Color','w','Position',[50 50 900 420]);
    hold on; grid on; box on;
    fill([t_axis, fliplr(t_axis)], [mu_E_do+std_E_do, fliplr(mu_E_do-std_E_do)], ...
        [0.9 0.2 0.2], 'FaceAlpha',0.2,'EdgeColor','none');
    plot(t_axis, mu_E_do, 'r-', 'LineWidth',2, ...
        'DisplayName',sprintf('Doppler-Only (σ_D=%.0f Hz)',CRIT_SIGMA_D));
    fill([t_axis, fliplr(t_axis)], [mu_E_dr+std_E_dr, fliplr(mu_E_dr-std_E_dr)], ...
        [0.1 0.5 0.9], 'FaceAlpha',0.2,'EdgeColor','none');
    plot(t_axis, mu_E_dr, 'b-', 'LineWidth',2, ...
        'DisplayName',sprintf('Doppler+Rate (σ_D=%.0f Hz, σ_{DR}=0.5 Hz/s)',CRIT_SIGMA_D));
    yline(0,'k:','LineWidth',0.8);
    xlabel('Time (s)'); ylabel('East Error (m)');
    title('Phase 4A — East Error: Doppler-Only vs Doppler+Rate');
    legend('Location','best','FontSize',9);
    save_fig(fh4a1, 'P4A_East_Error_Comparison', FIG_DIR);

    % Plot P4A-ii: 2D comparison
    fh4a2 = figure('Name','P4A 2D Comparison','Color','w','Position',[50 520 900 420]);
    hold on; grid on; box on;
    plot(t_axis, mu_2D_do, 'r-', 'LineWidth',2, ...
        'DisplayName','Doppler-Only');
    plot(t_axis, mu_2D_dr, 'b-', 'LineWidth',2, ...
        'DisplayName','Doppler+Rate');
    yline(CONV_THR_M,'k--','LineWidth',1.2,'DisplayName','50 m threshold');
    xlabel('Time (s)'); ylabel('2D Horizontal Error (m)');
    title('Phase 4A — 2D Horizontal: Doppler-Only vs Doppler+Rate');
    legend('Location','best','FontSize',9);
    save_fig(fh4a2, 'P4A_2D_Comparison', FIG_DIR);

    % ── Experiment 4B: Heatmaps ───────────────────────────────────────────
    fprintf('  4B: Heatmaps (σ_D × σ_DR sweep)\n');

    sigma_D_hm  = [10, 50, 100, 150, 200];          % Doppler noise [Hz]
    sigma_DR_hm = [0.05, 0.5, 1.0, 3.0, 5.0, 10.0]; % Doppler-rate noise [Hz/s]
    nD  = numel(sigma_D_hm);
    nDR = numel(sigma_DR_hm);

    hm_ss2D = nan(nDR, nD);   % rows = DR noise, cols = D noise
    hm_tc   = nan(nDR, nD);

    for di = 1:nD
        for dri = 1:nDR
            sigD  = sigma_D_hm(di);
            sigDR = sigma_DR_hm(dri);
            fprintf('    σ_D=%5.0f Hz, σ_DR=%5.2f Hz/s ... ', sigD, sigDR);

            c = build_base_cfg();
            c.duration_sec           = FIXED_DURATION;
            c.dt                     = FIXED_DT;
            c.sigma_doppler_Hz       = sigD;
            c.sigma_dopplerRate_Hzps = sigDR;
            c.meas.useDoppler        = true;
            c.meas.useDopplerRate    = true;
            c.tropo.enable           = false;
            c.clk.truth.initBias_mps = 0.10;
            c.P0_clkBias             = (0.15)^2;

            [mc_errE_hm, mc_errN_hm, mc_err2D_hm, ~, ~] = ...
                run_mc(c, N_MC, BASE_SEED + 420 + di*20 + dri, USE_PARFOR);

            if ~exist('t_axis_hm','var')
                tmp2     = run_single(c);
                t_axis_hm = tmp2.t;
            end

            kpi_hm = compute_kpis(t_axis_hm, mc_errE_hm, mc_errN_hm, mc_err2D_hm, ...
                                  SS_WINDOW_S, CONV_THR_M);
            hm_ss2D(dri,di) = kpi_hm.mean_ss_rms2D;
            hm_tc(dri,di)   = min(kpi_hm.mean_tc, FIXED_DURATION);   % cap at sim duration

            fprintf('SS 2D=%.1f m, tc=%.1f s\n', hm_ss2D(dri,di), hm_tc(dri,di));
        end
    end
    plot_heatmap(hm_ss2D, sigma_D_hm, sigma_DR_hm, ...
        'Doppler Noise σ_D (Hz)', 'Doppler-Rate Noise σ_{DR} (Hz/s)', ...
        'Phase 4B — SS 2D Horizontal RMS Error (m)', 'SS 2D RMS (m)', ...
        'P4B_Heatmap_SS2D', FIG_DIR);

    % Cap display of infinite convergence times at simulation duration
    hm_tc_disp = min(hm_tc, FIXED_DURATION);
    plot_heatmap(hm_tc_disp, sigma_D_hm, sigma_DR_hm, ...
        'Doppler Noise σ_D (Hz)', 'Doppler-Rate Noise σ_{DR} (Hz/s)', ...
        'Phase 4B — Convergence Time (s, capped at sim duration)', 'Conv. Time (s)', ...
        'P4B_Heatmap_ConvTime', FIG_DIR);

    fprintf('  Phase 4 complete.\n\n');
end

% %% ══════════════════════════════════════════════════════════════════════════
% %  §P5  PHASE 5 — Deep Filter Parameter Sweep (Pareto Optimisation)
% % ══════════════════════════════════════════════════════════════════════════
% %
% %  Environment: σ_D = 100 Hz,  σ_DR = 1.0 Hz/s
% %
% %  Grid search over:
% %    Rinflate  ∈ {1.0, 1.5, 2.5}
% %    sigma_pos ∈ {1e-4, 1e-3, 1e-2}  [m/√s]
% %    clkRegime ∈ {'tight', 'nominal', 'loose'}  (coherent P0/Q scaling)
% %
% %  Output: Pareto front — X = SS 2D RMS,  Y = Convergence Time
% %          Optimal combination printed to console.
% %
% % ══════════════════════════════════════════════════════════════════════════
% 
% if ismember(5, RUN_PHASES)
%     fprintf('══════════════════════════════════════════════════════════\n');
%     fprintf('  PHASE 5 — Deep IEKF Parameter Sweep (Pareto Front)\n');
%     fprintf('══════════════════════════════════════════════════════════\n');
% 
%     % ── Sweep grid ────────────────────────────────────────────────────────
%     Rinflate_vals = [1.0, 1.5, 2.5];
%     sigPos_vals   = [1e-4, 1e-3, 1e-2];
% 
%     % Clock regimes for the filter (coherent with nominal truth clock)
%     clk_regimes = struct( ...
%         'label',    {'Tight',   'Nominal', 'Loose'}, ...
%         'P0_bias',  {(0.05)^2, (0.15)^2, (0.5)^2}, ...
%         'P0_drift', {(1e-7)^2, (2.5e-6)^2, (1e-5)^2}, ...
%         'sigBias',  {5e-3,     1e-2,     3e-2}, ...
%         'sigDrift', {2e-4,     5e-4,     2e-3}  ...
%     );
%     nClk = numel(clk_regimes);
% 
%     nRI  = numel(Rinflate_vals);
%     nSP  = numel(sigPos_vals);
%     nTotal = nRI * nSP * nClk;
% 
%     fprintf('  Total parameter combinations: %d (× %d MC seeds = %d runs)\n', ...
%         nTotal, N_MC, nTotal * N_MC);
% 
%     % Pre-allocate result table
%     p5_results = struct( ...
%         'Rinflate',  cell(nTotal,1), ...
%         'sigPos',    cell(nTotal,1), ...
%         'clkLabel',  cell(nTotal,1), ...
%         'ss_rms2D',  cell(nTotal,1), ...
%         'tc',        cell(nTotal,1), ...
%         'ss_rmsE',   cell(nTotal,1), ...
%         'ss_rmsN',   cell(nTotal,1) ...
%     );
% 
%     combo_idx = 0;
% 
%     for ri = 1:nRI
%         for spi = 1:nSP
%             for cli = 1:nClk
%                 combo_idx = combo_idx + 1;
%                 Ri  = Rinflate_vals(ri);
%                 sp  = sigPos_vals(spi);
%                 cr  = clk_regimes(cli);
% 
%                 fprintf('  [%2d/%d] Rinflate=%.1f | sigPos=%.0e | clk=%s ... ', ...
%                     combo_idx, nTotal, Ri, sp, cr.label);
% 
%                 c = build_base_cfg();
%                 c.duration_sec           = FIXED_DURATION;
%                 c.dt                     = FIXED_DT;
%                 c.sigma_doppler_Hz       = 100;
%                 c.sigma_dopplerRate_Hzps = 1.0;
%                 c.meas.useDoppler        = true;
%                 c.meas.useDopplerRate    = true;
%                 c.tropo.enable           = false;
% 
%                 c.meas.Rinflate          = Ri;
%                 c.sigma_pos_rw           = sp;
% 
%                 % Truth clock: nominal
%                 c.clk.truth.initBias_mps            = 0.10;
%                 c.clk.truth.initDrift_mps2          = 5e-5;
%                 c.clk.truth.sigmaBias_rw_mps_sqrt   = 5e-3;
%                 c.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4;
% 
%                 % Filter clock: from grid
%                 c.P0_clkBias                        = cr.P0_bias;
%                 c.P0_clkDrift                       = cr.P0_drift;
%                 c.clk.filt.sigmaBias_rw_mps_sqrt    = cr.sigBias;
%                 c.clk.filt.sigmaDrift_rw_mps2_sqrt  = cr.sigDrift;
% 
%                 [mc_errE5, mc_errN5, mc_err2D5, ~, ~] = ...
%                     run_mc(c, N_MC, BASE_SEED + 500 + combo_idx*3, USE_PARFOR);
% 
%                 if ~exist('t_axis5','var')
%                     tmp5 = run_single(c);
%                     t_axis5 = tmp5.t;
%                 end
% 
%                 kpi5 = compute_kpis(t_axis5, mc_errE5, mc_errN5, mc_err2D5, ...
%                                     SS_WINDOW_S, CONV_THR_M);
% 
%                 p5_results(combo_idx).Rinflate = Ri;
%                 p5_results(combo_idx).sigPos   = sp;
%                 p5_results(combo_idx).clkLabel = cr.label;
%                 p5_results(combo_idx).ss_rms2D = kpi5.mean_ss_rms2D;
%                 p5_results(combo_idx).tc       = kpi5.mean_tc;
%                 p5_results(combo_idx).ss_rmsE  = kpi5.mean_ss_rmsE;
%                 p5_results(combo_idx).ss_rmsN  = kpi5.mean_ss_rmsN;
% 
%                 fprintf('SS 2D=%.1f m | tc=%.1f s\n', kpi5.mean_ss_rms2D, kpi5.mean_tc);
%             end
%         end
%     end
% 
%     % ── Compute Pareto front ──────────────────────────────────────────────
%     all_ss  = [p5_results.ss_rms2D]';
%     all_tc  = [p5_results.tc]';
% 
%     % Replace inf tc with duration + penalty for Pareto scoring
%     all_tc_fin = all_tc;
%     all_tc_fin(isinf(all_tc_fin)) = FIXED_DURATION * 1.5;
% 
%     % Normalise for Pareto identification
%     ss_norm = (all_ss  - min(all_ss))  / (max(all_ss)  - min(all_ss) + eps);
%     tc_norm = (all_tc_fin - min(all_tc_fin)) / (max(all_tc_fin) - min(all_tc_fin) + eps);
% 
%     % A point i is Pareto-optimal if no other point j dominates it
%     % (i.e., j has lower or equal ss AND lower or equal tc, with at least one strict)
%     nC = numel(all_ss);
%     is_pareto = true(nC,1);
%     for i = 1:nC
%         for j = 1:nC
%             if i == j, continue; end
%             if (all_ss(j) <= all_ss(i)) && (all_tc_fin(j) <= all_tc_fin(i)) && ...
%                (all_ss(j)  < all_ss(i)  || all_tc_fin(j)  < all_tc_fin(i))
%                 is_pareto(i) = false;
%                 break;
%             end
%         end
%     end
% 
%     % Optimal: Pareto point with minimum weighted sum (equal weights)
%     pareto_idx = find(is_pareto);
%     [~, best_local] = min(ss_norm(pareto_idx) + tc_norm(pareto_idx));
%     best_idx = pareto_idx(best_local);
% 
%     % ── PLOT P5: Pareto Front Scatter ─────────────────────────────────────
%     fh5 = figure('Name','P5 Pareto Front','Color','w','Position',[50 50 900 600]);
%     hold on; grid on; box on;
% 
%     % Colour by clock regime
%     clk_colors = {[0.2 0.6 0.9], [0.1 0.8 0.3], [0.9 0.4 0.1]};
%     clk_labels = {'Tight', 'Nominal', 'Loose'};
%     for cli = 1:nClk
%         mask = strcmp({p5_results.clkLabel}, clk_labels{cli});
%         scatter(all_ss(mask & ~is_pareto), all_tc_fin(mask & ~is_pareto), 60, ...
%             clk_colors{cli}, 'o', 'filled', 'MarkerFaceAlpha', 0.4, ...
%             'DisplayName', ['Dominated — ' clk_labels{cli}]);
%         scatter(all_ss(mask &  is_pareto), all_tc_fin(mask &  is_pareto), 120, ...
%             clk_colors{cli}, 'd', 'filled', 'LineWidth', 1.5, ...
%             'MarkerEdgeColor','k', ...
%             'DisplayName', ['Pareto — ' clk_labels{cli}]);
%     end
% 
%     % Highlight best point
%     scatter(all_ss(best_idx), all_tc_fin(best_idx), 250, 'r', 'p', 'filled', ...
%         'MarkerEdgeColor','k','LineWidth',1.5,'DisplayName','Optimal');
%     text(all_ss(best_idx) * 1.03, all_tc_fin(best_idx), ...
%         sprintf('  Opt\n  R=%.1f\n  σ_p=%.0e\n  Clk=%s', ...
%             p5_results(best_idx).Rinflate, p5_results(best_idx).sigPos, ...
%             p5_results(best_idx).clkLabel), ...
%         'FontSize',9,'Color','r','FontWeight','bold');
% 
%     xlabel('Steady-State 2D Horizontal RMS Error (m)');
%     ylabel('Convergence Time (s)');
%     title('Phase 5 — IEKF Parameter Pareto Front', 'FontWeight','bold');
%     legend('Location','best','FontSize',8);
%     save_fig(fh5, 'P5_Pareto_Front', FIG_DIR);
% 
%     % ── Print optimal configuration to console ────────────────────────────
%     fprintf('\n  ╔══════════════════════════════════════════════════════╗\n');
%     fprintf('  ║  PHASE 5 — OPTIMAL PARAMETER COMBINATION            ║\n');
%     fprintf('  ╠══════════════════════════════════════════════════════╣\n');
%     fprintf('  ║  Pareto-optimal (min weighted sum of normalised      ║\n');
%     fprintf('  ║  SS-Error + Convergence-Time):                       ║\n');
%     fprintf('  ║                                                      ║\n');
%     fprintf('  ║  cfg.meas.Rinflate          = %.1f\n', p5_results(best_idx).Rinflate);
%     fprintf('  ║  cfg.sigma_pos_rw           = %.1e\n', p5_results(best_idx).sigPos);
%     fprintf('  ║  Clock regime               = %s\n',   p5_results(best_idx).clkLabel);
%     fprintf('  ║  SS 2D Horizontal RMS       = %.2f m\n', p5_results(best_idx).ss_rms2D);
%     fprintf('  ║  Convergence time           = %.1f s\n', p5_results(best_idx).tc);
%     fprintf('  ╚══════════════════════════════════════════════════════╝\n\n');
% 
%     % Map clock regime label to actual cfg values
%     opt_cli = find(strcmp(clk_labels, p5_results(best_idx).clkLabel));
%     fprintf('  Optimal cfg values:\n');
%     fprintf('    cfg.meas.Rinflate                    = %.4f\n', p5_results(best_idx).Rinflate);
%     fprintf('    cfg.sigma_pos_rw                     = %.4e\n', p5_results(best_idx).sigPos);
%     fprintf('    cfg.P0_clkBias                       = %.4e   %% [m/s]^2\n', clk_regimes(opt_cli).P0_bias);
%     fprintf('    cfg.P0_clkDrift                      = %.4e   %% [m/s²]^2\n', clk_regimes(opt_cli).P0_drift);
%     fprintf('    cfg.clk.filt.sigmaBias_rw_mps_sqrt   = %.4e\n', clk_regimes(opt_cli).sigBias);
%     fprintf('    cfg.clk.filt.sigmaDrift_rw_mps2_sqrt = %.4e\n', clk_regimes(opt_cli).sigDrift);
% 
%     fprintf('  Phase 5 complete.\n\n');
% end



%% ══════════════════════════════════════════════════════════════════════════
%  §P5  PHASE 5 — Deep IEKF Internal Parameter Sweep (Pareto Optimisation)
%
%  Environment (fixed, challenging but realistic):
%    sigma_doppler_Hz    = 50 Hz
%    sigma_dopplerRate_Hzps = 1.0 Hz/s
%    useDoppler = true  |  useDopplerRate = true
%
%  Clock truth (fixed "Nominal" model — NOT swept):
%    initBias_mps            = 0.20
%    initDrift_mps2          = 1e-4
%    sigmaBias_rw_mps_sqrt   = 2.5e-3
%    sigmaDrift_rw_mps2_sqrt = 3.5e-5
%
%  Grid search — 5 IEKF internal parameters, 3 levels each = 243 combinations:
%    startupInflateEpochs : [10, 30, 50]
%    startupRmult         : [10, 25, 50]
%    maxStep_m            : [10, 25, 50]
%    Rinflate             : [1.0, 1.5, 2.0]
%    k_tune               : [1, 10, 100]   (scales filter clock P0 and Q)
%
%  k_tune semantics:
%    P0_clkBias_filt  = k_tune * (truth_initBias)^2
%    P0_clkDrift_filt = k_tune * (truth_initDrift)^2
%    Q_clkBias_filt   = k_tune * (truth_sigmaBias)^2  * dt
%    Q_clkDrift_filt  = k_tune * (truth_sigmaDrift)^2 * dt
%    This sweeps from "matched" (k=1) to "very loose" (k=100), which is
%    mathematically equivalent to testing how much the filter needs to
%    "forget" its clock prior to converge quickly.
%
%  Execution: parfor over all 243 combinations, N_MC_P5=5 seeds each.
%
%  Outputs:
%    1. Pareto front scatter plot (SS 2D RMS vs Convergence Time)
%    2. Pareto-optimal frontier highlighted
%    3. Knee point identified and cfg values printed to console
% ══════════════════════════════════════════════════════════════════════════
 
if ismember(5, RUN_PHASES)
 
    fprintf('══════════════════════════════════════════════════════════\n');
    fprintf('  PHASE 5 — Deep IEKF Internal Parameter Sweep\n');
    fprintf('══════════════════════════════════════════════════════════\n');
 
    % ── Reduced MC count to keep wall-time tractable ──────────────────────
    N_MC_P5 = 3;
 
    % ── Fixed environment ─────────────────────────────────────────────────
    P5_SIGMA_D  = 50.0;    % Doppler noise [Hz]
    P5_SIGMA_DR = 1.0;     % Doppler-rate noise [Hz/s]
 
    % ── Fixed clock truth (Nominal model) ────────────────────────────────
    P5_CLK_BIAS0  = 0.20;      % initBias_mps
    P5_CLK_DRIFT0 = 1e-4;      % initDrift_mps2
    P5_CLK_SIG_B  = 2.5e-3;    % sigmaBias_rw_mps_sqrt   (truth process noise)
    P5_CLK_SIG_D  = 3.5e-5;    % sigmaDrift_rw_mps2_sqrt (truth process noise)
 
    % ── Sweep grids ───────────────────────────────────────────────────────
    grid_inflEpochs = [10,  50,  200  ];   % startupInflateEpochs
    grid_rmult      = [10,  25,  50  ];   % startupRmult
    grid_maxStep    = [10,  50,  100  ];   % maxStep_m  [m]
    grid_Rinflate   = [1.0, 1.5, 2.0 ];   % meas.Rinflate
    grid_ktune      = [1,   10,  100 ];   % clock P0/Q multiplier
 
    nIE  = numel(grid_inflEpochs);
    nRM  = numel(grid_rmult);
    nMS  = numel(grid_maxStep);
    nRI  = numel(grid_Rinflate);
    nKT  = numel(grid_ktune);
    nTotal5 = nIE * nRM * nMS * nRI * nKT;   % = 243
 
    fprintf('  Grid size : %d x %d x %d x %d x %d = %d combinations\n', ...
        nIE, nRM, nMS, nRI, nKT, nTotal5);
    fprintf('  MC seeds  : %d per combination  (%d total solver calls)\n', ...
        N_MC_P5, nTotal5 * N_MC_P5);
 
    % ── Pre-allocate result arrays (one row per combo) ────────────────────
    p5_ie    = nan(nTotal5, 1);   % startupInflateEpochs
    p5_rm    = nan(nTotal5, 1);   % startupRmult
    p5_ms    = nan(nTotal5, 1);   % maxStep_m
    p5_ri    = nan(nTotal5, 1);   % Rinflate
    p5_kt    = nan(nTotal5, 1);   % k_tune
    p5_ss2D  = nan(nTotal5, 1);   % mean SS 2D RMS   [m]
    p5_tc    = nan(nTotal5, 1);   % mean conv. time  [s]
    p5_ssE   = nan(nTotal5, 1);   % mean SS East RMS [m]
    p5_ssN   = nan(nTotal5, 1);   % mean SS North RMS[m]
 
    % ── Build flat list of all 243 parameter combinations ─────────────────
    % (pre-building avoids nested loops inside parfor)
    combo_ie = nan(nTotal5, 1);
    combo_rm = nan(nTotal5, 1);
    combo_ms = nan(nTotal5, 1);
    combo_ri = nan(nTotal5, 1);
    combo_kt = nan(nTotal5, 1);
 
    idx = 0;
    for a = 1:nIE
        for b = 1:nRM
            for c = 1:nMS
                for d = 1:nRI
                    for e = 1:nKT
                        idx = idx + 1;
                        combo_ie(idx) = grid_inflEpochs(a);
                        combo_rm(idx) = grid_rmult(b);
                        combo_ms(idx) = grid_maxStep(c);
                        combo_ri(idx) = grid_Rinflate(d);
                        combo_kt(idx) = grid_ktune(e);
                    end
                end
            end
        end
    end
 
    % ── Build one reference cfg (environment + truth fixed, tuning TBD) ───
    cfg_p5_base = build_base_cfg();
    cfg_p5_base.duration_sec                      = FIXED_DURATION;
    cfg_p5_base.dt                                = FIXED_DT;
    cfg_p5_base.sigma_doppler_Hz                  = P5_SIGMA_D;
    cfg_p5_base.sigma_dopplerRate_Hzps            = P5_SIGMA_DR;
    cfg_p5_base.meas.useDoppler                   = true;
    cfg_p5_base.meas.useDopplerRate               = true;
    cfg_p5_base.tropo.enable                      = false;
    % Fixed truth clock
    cfg_p5_base.clk.truth.initBias_mps            = P5_CLK_BIAS0;
    cfg_p5_base.clk.truth.initDrift_mps2          = P5_CLK_DRIFT0;
    cfg_p5_base.clk.truth.sigmaBias_rw_mps_sqrt   = P5_CLK_SIG_B;
    cfg_p5_base.clk.truth.sigmaDrift_rw_mps2_sqrt = P5_CLK_SIG_D;
    % Filter init is always zero (estimator doesn't know truth)
    cfg_p5_base.clk.filter_init_bias_mps          = 0.0;
    cfg_p5_base.clk.filter_init_drift_mps2        = 0.0;
 
    % Get time axis length once (all runs share the same t vector)
    tmp_p5ref = run_single(cfg_p5_base);
    t5        = tmp_p5ref.t;
    K5        = numel(t5);
 
    % ── Parfor loop over all 243 combinations ─────────────────────────────
    % parfor requires all variables to be sliced or broadcast;
    % we broadcast cfg_p5_base as a fixed template and override only the
    % swept fields inside the loop body.
 
    % Temporary per-combo result arrays for parfor (can't index struct in parfor)
    tmp_ss2D = nan(nTotal5, 1);
    tmp_tc   = nan(nTotal5, 1);
    tmp_ssE  = nan(nTotal5, 1);
    tmp_ssN  = nan(nTotal5, 1);
 
    % Broadcast constants needed inside the loop
    bdt         = FIXED_DT;
    bDuration   = FIXED_DURATION;
    bSigD       = P5_SIGMA_D;
    bSigDR      = P5_SIGMA_DR;
    bClkB0      = P5_CLK_BIAS0;
    bClkD0      = P5_CLK_DRIFT0;
    bClkSigB    = P5_CLK_SIG_B;
    bClkSigD    = P5_CLK_SIG_D;
    bSSwin      = SS_WINDOW_S;
    bConvThr    = CONV_THR_M;
    bBaseSeed   = BASE_SEED + 500;
    bNmc        = N_MC_P5;
 
    % Sliced input arrays
    sl_ie = combo_ie;
    sl_rm = combo_rm;
    sl_ms = combo_ms;
    sl_ri = combo_ri;
    sl_kt = combo_kt;
 
    if USE_PARFOR
        parfor ci = 1:nTotal5
            c = build_base_cfg();
            c.duration_sec                      = bDuration;
            c.dt                                = bdt;
            c.sigma_doppler_Hz                  = bSigD;
            c.sigma_dopplerRate_Hzps            = bSigDR;
            c.meas.useDoppler                   = true;
            c.meas.useDopplerRate               = true;
            c.tropo.enable                      = false;
 
            % Fixed truth clock
            c.clk.truth.initBias_mps            = bClkB0;
            c.clk.truth.initDrift_mps2          = bClkD0;
            c.clk.truth.sigmaBias_rw_mps_sqrt   = bClkSigB;
            c.clk.truth.sigmaDrift_rw_mps2_sqrt = bClkSigD;
            c.clk.filter_init_bias_mps          = 0.0;
            c.clk.filter_init_drift_mps2        = 0.0;
 
            % ── Swept internal tuning parameters ──────────────────────────
            kt = sl_kt(ci);
 
            c.filter.startupInflateEpochs = sl_ie(ci);
            c.filter.startupRmult         = sl_rm(ci);
            c.filter.maxStep_m            = sl_ms(ci);
            c.meas.Rinflate               = sl_ri(ci);
 
            % k_tune scales filter clock P0 and Q relative to truth values
            c.P0_clkBias  = kt * bClkB0^2;
            c.P0_clkDrift = kt * bClkD0^2;
            c.clk.filt.sigmaBias_rw_mps_sqrt   = sqrt(kt) * bClkSigB;
            c.clk.filt.sigmaDrift_rw_mps2_sqrt = sqrt(kt) * bClkSigD;
 
            % ── Run N_MC_P5 seeds ─────────────────────────────────────────
            seed_base_ci = bBaseSeed + ci * (bNmc + 1);
            mc_e  = nan(bNmc, numel(t5));   % can't call run_mc inside parfor
            mc_n  = nan(bNmc, numel(t5));   % (nested broadcast issue) so
            mc_2d = nan(bNmc, numel(t5));   % we replicate the seed loop here
 
            for n = 1:bNmc
                c.seedMeas = seed_base_ci + n;
                r = run_single(c);
                mc_e(n,:)  = r.errE;
                mc_n(n,:)  = r.errN;
                mc_2d(n,:) = r.err2D;
            end
 
            % ── KPI extraction ────────────────────────────────────────────
            t_row   = t5(:)';
            ss_mask = t_row >= (t_row(end) - bSSwin);
 
            ss2D_v = nan(bNmc, 1);
            ssE_v  = nan(bNmc, 1);
            ssN_v  = nan(bNmc, 1);
            tc_v   = nan(bNmc, 1);
 
            for n = 1:bNmc
                ss2D_v(n) = rms(mc_2d(n, ss_mask));
                ssE_v(n)  = rms(mc_e(n,  ss_mask));
                ssN_v(n)  = rms(mc_n(n,  ss_mask));
                tc_v(n)   = find_convergence_strict(t_row, mc_2d(n,:), bConvThr);
            end
 
            tmp_ss2D(ci) = mean(ss2D_v, 'omitnan');
            tmp_tc(ci)   = mean(tc_v,   'omitnan');
            tmp_ssE(ci)  = mean(ssE_v,  'omitnan');
            tmp_ssN(ci)  = mean(ssN_v,  'omitnan');
        end
    else
        % Serial fallback with progress reporting
        for ci = 1:nTotal5
            if mod(ci, 20) == 1
                fprintf('  [%3d / %d]  ie=%2d  rm=%2d  ms=%2d  Ri=%.1f  kt=%3d\n', ...
                    ci, nTotal5, sl_ie(ci), sl_rm(ci), sl_ms(ci), sl_ri(ci), sl_kt(ci));
            end
 
            c = build_base_cfg();
            c.duration_sec                      = FIXED_DURATION;
            c.dt                                = FIXED_DT;
            c.sigma_doppler_Hz                  = P5_SIGMA_D;
            c.sigma_dopplerRate_Hzps            = P5_SIGMA_DR;
            c.meas.useDoppler                   = true;
            c.meas.useDopplerRate               = true;
            c.tropo.enable                      = false;
            c.clk.truth.initBias_mps            = P5_CLK_BIAS0;
            c.clk.truth.initDrift_mps2          = P5_CLK_DRIFT0;
            c.clk.truth.sigmaBias_rw_mps_sqrt   = P5_CLK_SIG_B;
            c.clk.truth.sigmaDrift_rw_mps2_sqrt = P5_CLK_SIG_D;
            c.clk.filter_init_bias_mps          = 0.0;
            c.clk.filter_init_drift_mps2        = 0.0;
 
            kt = combo_kt(ci);
 
            c.filter.startupInflateEpochs = combo_ie(ci);
            c.filter.startupRmult         = combo_rm(ci);
            c.filter.maxStep_m            = combo_ms(ci);
            c.meas.Rinflate               = combo_ri(ci);
 
            c.P0_clkBias  = kt * P5_CLK_BIAS0^2;
            c.P0_clkDrift = kt * P5_CLK_DRIFT0^2;
            c.clk.filt.sigmaBias_rw_mps_sqrt   = sqrt(kt) * P5_CLK_SIG_B;
            c.clk.filt.sigmaDrift_rw_mps2_sqrt = sqrt(kt) * P5_CLK_SIG_D;
 
            seed_base_ci = BASE_SEED + 500 + ci * (N_MC_P5 + 1);
 
            mc_e  = nan(N_MC_P5, K5);
            mc_n  = nan(N_MC_P5, K5);
            mc_2d = nan(N_MC_P5, K5);
 
            for n = 1:N_MC_P5
                c.seedMeas = seed_base_ci + n;
                r = run_single(c);
                mc_e(n,:)  = r.errE;
                mc_n(n,:)  = r.errN;
                mc_2d(n,:) = r.err2D;
            end
 
            ss_mask = t5(:)' >= (t5(end) - SS_WINDOW_S);
            ss2D_v  = nan(N_MC_P5, 1);
            ssE_v   = nan(N_MC_P5, 1);
            ssN_v   = nan(N_MC_P5, 1);
            tc_v    = nan(N_MC_P5, 1);
 
            for n = 1:N_MC_P5
                ss2D_v(n) = rms(mc_2d(n, ss_mask));
                ssE_v(n)  = rms(mc_e(n,  ss_mask));
                ssN_v(n)  = rms(mc_n(n,  ss_mask));
                tc_v(n)   = find_convergence_strict(t5(:)', mc_2d(n,:), CONV_THR_M);
            end
 
            tmp_ss2D(ci) = mean(ss2D_v, 'omitnan');
            tmp_tc(ci)   = mean(tc_v,   'omitnan');
            tmp_ssE(ci)  = mean(ssE_v,  'omitnan');
            tmp_ssN(ci)  = mean(ssN_v,  'omitnan');
        end
    end
 
    % ── Store parameter labels and results back into named arrays ──────────
    p5_ie   = combo_ie;
    p5_rm   = combo_rm;
    p5_ms   = combo_ms;
    p5_ri   = combo_ri;
    p5_kt   = combo_kt;
    p5_ss2D = tmp_ss2D;
    p5_tc   = tmp_tc;
    p5_ssE  = tmp_ssE;
    p5_ssN  = tmp_ssN;
 
    % Replace inf convergence times with a finite penalty for Pareto work
    % (runs that never converge sit at the far right of the Pareto space)
    tc_fin = p5_tc;
    tc_fin(isinf(tc_fin)) = FIXED_DURATION * 1.5;
 
    % ── Pareto dominance ──────────────────────────────────────────────────
    % Point i is dominated if there exists j such that:
    %   p5_ss2D(j) <= p5_ss2D(i)  AND  tc_fin(j) <= tc_fin(i)
    %   with at least one strict inequality.
    is_pareto5 = true(nTotal5, 1);
    for i = 1:nTotal5
        if isnan(p5_ss2D(i)) || isnan(tc_fin(i)), is_pareto5(i) = false; continue; end
        for j = 1:nTotal5
            if i == j, continue; end
            if isnan(p5_ss2D(j)) || isnan(tc_fin(j)), continue; end
            if (p5_ss2D(j) <= p5_ss2D(i)) && (tc_fin(j) <= tc_fin(i)) && ...
               (p5_ss2D(j) <  p5_ss2D(i)  || tc_fin(j)  < tc_fin(i))
                is_pareto5(i) = false;
                break;
            end
        end
    end
 
    % ── Knee point identification ─────────────────────────────────────────
    % The "knee" is the Pareto-optimal point that maximises the perpendicular
    % distance to the straight line connecting the two extreme Pareto points
    % (min SS error vs min convergence time). This is the standard geometric
    % knee definition for bi-objective problems.
    pareto_idx5 = find(is_pareto5);
 
    % Normalise both objectives to [0,1] for distance calculation
    ss_min = min(p5_ss2D(pareto_idx5));   ss_max = max(p5_ss2D(pareto_idx5));
    tc_min = min(tc_fin(pareto_idx5));    tc_max = max(tc_fin(pareto_idx5));
 
    ss_norm5 = (p5_ss2D(pareto_idx5) - ss_min) / max(ss_max - ss_min, eps);
    tc_norm5 = (tc_fin(pareto_idx5)  - tc_min) / max(tc_max - tc_min, eps);
 
    % Direction vector of the line from (0,1) [best tc] to (1,0) [best ss]
    % (the utopia line in normalised space)
    % Perpendicular distance from point (x,y) to line x + y = 1:
    %   d = |x + y - 1| / sqrt(2)
    % The knee is the point FARTHEST from this line (maximises d)
    knee_dist = abs(ss_norm5 + tc_norm5 - 1) / sqrt(2);
    [~, knee_local] = max(knee_dist);
    knee_global     = pareto_idx5(knee_local);
 
    % ── Plot: Pareto front scatter ─────────────────────────────────────────
    fh5 = figure('Name','P5 Pareto Front — IEKF Internal Tuning', ...
                 'Color','w', 'Position',[50 50 1000 650]);
    hold on; grid on; box on;
 
    % --- All dominated points, coloured by k_tune for interpretability ----
    kt_vals   = unique(p5_kt);
    kt_colors = {[0.70 0.85 0.95], [0.35 0.65 0.85], [0.10 0.35 0.65]};  % light→dark blue
 
    for ki = 1:numel(kt_vals)
        mask_dom = (p5_kt == kt_vals(ki)) & ~is_pareto5 & ~isnan(p5_ss2D);
        if any(mask_dom)
            scatter(p5_ss2D(mask_dom), tc_fin(mask_dom), 28, ...
                kt_colors{ki}, 'o', 'filled', 'MarkerFaceAlpha', 0.30, ...
                'DisplayName', sprintf('Dominated  k\\_tune=%d', kt_vals(ki)));
        end
    end
 
    % --- Pareto-optimal frontier points ------------------------------------
    % Sort for a clean connected line
    pareto_ss = p5_ss2D(pareto_idx5);
    pareto_tc = tc_fin(pareto_idx5);
    [pareto_ss_sorted, sort_ord] = sort(pareto_ss);
    pareto_tc_sorted = pareto_tc(sort_ord);
 
    plot(pareto_ss_sorted, pareto_tc_sorted, 'k--', 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');   % connecting line (not in legend)
 
    for ki = 1:numel(kt_vals)
        mask_par = (p5_kt(pareto_idx5) == kt_vals(ki));
        if any(mask_par)
            scatter(pareto_ss(mask_par), pareto_tc(mask_par), 80, ...
                kt_colors{ki}, 'd', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
                'DisplayName', sprintf('Pareto  k\\_tune=%d', kt_vals(ki)));
        end
    end
 
    % --- Knee point --------------------------------------------------------
    scatter(p5_ss2D(knee_global), tc_fin(knee_global), 220, ...
        [0.85 0.10 0.10], 'p', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'Knee point (optimal trade-off)');
 
    % Annotation box near the knee point
    kx = p5_ss2D(knee_global);
    ky = tc_fin(knee_global);
    annotation_str = sprintf( ...
        ['  startupInflateEpochs = %d\n' ...
         '  startupRmult         = %d\n' ...
         '  maxStep\\_m           = %d m\n' ...
         '  Rinflate             = %.1f\n' ...
         '  k\\_tune              = %d'], ...
        p5_ie(knee_global), p5_rm(knee_global), ...
        p5_ms(knee_global), p5_ri(knee_global), p5_kt(knee_global));
 
    text(kx + 0.5, ky + 0.5, annotation_str, ...
        'FontSize', 9, 'FontName', 'Courier New', ...
        'BackgroundColor', [1.0 0.95 0.85], ...
        'EdgeColor', [0.85 0.10 0.10], 'LineWidth', 1.0, ...
        'VerticalAlignment', 'bottom');
 
    % Axis labels and title
    xlabel('Average 2D Horizontal Steady-State RMS Error  (m)', 'FontSize', 12);
    ylabel('Average Convergence Time  (s)', 'FontSize', 12);
    title({'Phase 5  —  IEKF Internal Parameter Pareto Front', ...
           sprintf('\\sigma_D = %.0f Hz  |  \\sigma_{DR} = %.1f Hz/s  |  %d combinations  \\times  %d MC seeds', ...
               P5_SIGMA_D, P5_SIGMA_DR, nTotal5, N_MC_P5)}, ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 9);
 
    % Annotate how many combos never converged
    n_noconv = sum(isinf(p5_tc));
    if n_noconv > 0
        text(0.98, 0.04, sprintf('%d / %d combos did not converge', n_noconv, nTotal5), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'FontSize', 9, 'Color', [0.6 0.0 0.0], 'FontWeight', 'bold');
    end
 
    save_fig(fh5, 'P5_Pareto_Front', FIG_DIR);
 
    % ── Print knee-point cfg values to console ────────────────────────────
    fprintf('\n  ╔══════════════════════════════════════════════════════════════╗\n');
    fprintf('  ║  PHASE 5 — KNEE-POINT (Optimal Pareto Trade-off)            ║\n');
    fprintf('  ║  Criterion: max perpendicular distance to utopia line        ║\n');
    fprintf('  ╠══════════════════════════════════════════════════════════════╣\n');
    fprintf('  ║  KPIs:                                                       ║\n');
    fprintf('  ║    Average SS 2D Horizontal RMS  = %6.2f m                  ║\n', p5_ss2D(knee_global));
    fprintf('  ║    Average SS East  RMS          = %6.2f m                  ║\n', p5_ssE(knee_global));
    fprintf('  ║    Average SS North RMS          = %6.2f m                  ║\n', p5_ssN(knee_global));
    if isinf(p5_tc(knee_global))
        fprintf('  ║    Average Convergence Time      =     inf  (did not converge)  ║\n');
    else
        fprintf('  ║    Average Convergence Time      = %6.1f s                  ║\n', p5_tc(knee_global));
    end
    fprintf('  ╠══════════════════════════════════════════════════════════════╣\n');
    fprintf('  ║  Optimal cfg values (copy into MAIN_baseline.m):            ║\n');
    fprintf('  ║                                                              ║\n');
    fprintf('  ║    cfg.filter.startupInflateEpochs = %d;\n',     p5_ie(knee_global));
    fprintf('  ║    cfg.filter.startupRmult         = %d;\n',     p5_rm(knee_global));
    fprintf('  ║    cfg.filter.maxStep_m            = %d;\n',     p5_ms(knee_global));
    fprintf('  ║    cfg.meas.Rinflate               = %.1f;\n',   p5_ri(knee_global));
    fprintf('  ║\n');
    fprintf('  ║    %% k_tune = %d  applied to Nominal clock truth:\n', p5_kt(knee_global));
    fprintf('  ║    cfg.P0_clkBias                        = %.4e;  %% [m/s]^2\n',  ...
        p5_kt(knee_global) * P5_CLK_BIAS0^2);
    fprintf('  ║    cfg.P0_clkDrift                       = %.4e;  %% [m/s^2]^2\n', ...
        p5_kt(knee_global) * P5_CLK_DRIFT0^2);
    fprintf('  ║    cfg.clk.filt.sigmaBias_rw_mps_sqrt   = %.4e;\n', ...
        sqrt(p5_kt(knee_global)) * P5_CLK_SIG_B);
    fprintf('  ║    cfg.clk.filt.sigmaDrift_rw_mps2_sqrt = %.4e;\n', ...
        sqrt(p5_kt(knee_global)) * P5_CLK_SIG_D);
    fprintf('  ╚══════════════════════════════════════════════════════════════╝\n\n');
 
    % ── Bonus: print the full Pareto-optimal frontier as a ranked table ────
    fprintf('  Full Pareto frontier (sorted by SS 2D RMS, ascending):\n');
    fprintf('  %-6s  %-6s  %-8s  %-8s  %-6s  %-12s  %-12s\n', ...
        'IE', 'Rmult', 'maxStep', 'Rinflate', 'ktune', 'SS 2D (m)', 'Conv (s)');
    fprintf('  %s\n', repmat('-', 1, 68));
 
    [~, front_sort] = sort(p5_ss2D(pareto_idx5));
    for fi = 1:numel(pareto_idx5)
        gi = pareto_idx5(front_sort(fi));
        tc_str = sprintf('%.1f', tc_fin(gi));
        if isinf(p5_tc(gi)), tc_str = 'inf'; end
        knee_marker = '';
        if gi == knee_global, knee_marker = ' <-- KNEE'; end
        fprintf('  %-6d  %-6d  %-8d  %-8.1f  %-6d  %-12.2f  %-12s%s\n', ...
            p5_ie(gi), p5_rm(gi), p5_ms(gi), p5_ri(gi), p5_kt(gi), ...
            p5_ss2D(gi), tc_str, knee_marker);
    end
    fprintf('\n  Phase 5 complete.\n\n');
 
end  % if ismember(5, RUN_PHASES)

%% ══════════════════════════════════════════════════════════════════════════
%  §P6  PHASE 6 — Multi-Channel Advantage & Sequential Handovers
%
%  PHYSICS MOTIVATION:
%    A single-channel Doppler receiver is fundamentally East-limited because
%    one LOS provides range-rate along a single direction; East observability
%    only accumulates as the satellite arcs across the sky (slow).
%    Tracking N >= 2 simultaneous satellites injects N independent LOS
%    directions at every epoch, instantly spanning the East-North plane and
%    collapsing the convergence time from tens of seconds to a few seconds.
%
%  N-CHANNEL IMPLEMENTATION (wrapper-level, no solver modification needed):
%    N=1  → capture.mode = 'stickyBest'  (sequential, one satellite at a time)
%    N>1  → capture.mode = 'allVisible'  + maxSats = N
%    The constellation builder pre-selects the top-N satellites by peak
%    elevation in the selection window.  Setting maxSats=N then limits the
%    pool, so allVisible uses at most N satellites at each epoch — a faithful
%    model of a receiver with N parallel tracking channels assigned greedily
%    by geometry.
%
%  EXPERIMENT 6-1 — Single-Channel Sequential Patience (180 s):
%    stickyBest | switchBlackout_s = 10 s | duration = 180 s
%    Show: sawtooth error recovery, EKF flywheel during blackout,
%          eventual good accuracy but at the cost of long wait time.
%    Output: 2-panel 180 s time-series (2D error + clkBias error) with
%            shaded blackout regions and satellite-switch markers.
%
%  EXPERIMENT 6-2 — Multi-Channel Instant Geometry (100 s):
%    N = 1, 2, 3, 4 channels | allVisible + maxSats = N
%    Show: East error collapse as N increases.
%    Output A: overlay time-series of 2D Horizontal Error for all N.
%    Output B: bar chart of SS 2D / East / North vs N, with conv-time overlay.
% ══════════════════════════════════════════════════════════════════════════
 
 
if ismember(6, RUN_PHASES)
 
    fprintf('══════════════════════════════════════════════════════════\n');
    fprintf('  PHASE 6 — Multi-Channel Advantage & Sequential Handovers\n');
    fprintf('══════════════════════════════════════════════════════════\n');
 
    % ── Common environment for Phase 6 ────────────────────────────────────
    P6_SIGMA_D  = 50.0;   % Doppler noise [Hz]  — moderate, realistic
    P6_SIGMA_DR = 1.0;    % Doppler-rate  [Hz/s]
    P6_N_MC     = 20;     % re-use global MC seed count
 
    % ════════════════════════════════════════════════════════════════════
    %  EXPERIMENT 6-1 — Single-Channel: 600 s with 10 s handover blackout
    % ════════════════════════════════════════════════════════════════════
    fprintf('\n  --- Experiment 6-1: Single-Channel 180 s (handover blackouts) ---\n');
 
    cfg_e1 = build_base_cfg();
    cfg_e1.duration_sec             = 600;
    cfg_e1.dt                       = 0.1;
    cfg_e1.sigma_doppler_Hz         = P6_SIGMA_D;
    cfg_e1.sigma_dopplerRate_Hzps   = P6_SIGMA_DR;
    cfg_e1.meas.useDoppler          = true;
    cfg_e1.meas.useDopplerRate      = true;
    cfg_e1.tropo.enable             = false;
    cfg_e1.capture.mode             = 'stickyBest';
    cfg_e1.capture.switchBlackout_s = 10;
    cfg_e1.capture.minDwell_s       = 30;
 
    % Nominal coherent clock
    cfg_e1.clk.truth.initBias_mps            = 0.10;
    cfg_e1.clk.truth.initDrift_mps2          = 5e-5;
    cfg_e1.clk.truth.sigmaBias_rw_mps_sqrt   = 5e-3;
    cfg_e1.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4;
    cfg_e1.P0_clkBias                        = (0.15)^2;
    cfg_e1.P0_clkDrift                       = (2.5e-7*2)^2;
    cfg_e1.clk.filt.sigmaBias_rw_mps_sqrt    = 1.0e-2;
    cfg_e1.clk.filt.sigmaDrift_rw_mps2_sqrt  = 5.0e-4;
 
    % ── Step 1: one reference run to capture the deterministic trk struct.
    fprintf('    Extracting geometry (trk) from reference run...\n');
    cfg_e1.seedMeas  = BASE_SEED + 600;
    cfg_e1.plot.enable = false;
    out_e1_ref = real_strtg_sequential_EKF_fastrealistic(cfg_e1);
 
    trk_e1 = out_e1_ref.trk;
    t_e1   = out_e1_ref.t(:)';   % force row 1×K
    K_e1   = numel(t_e1);
 
    % ── Step 2: MC loop for smooth mean ± std traces ─────────────────────
    fprintf('    Running %d MC seeds for Exp 6-1...\n', P6_N_MC);
 
    mc6_e1_err2D = nan(P6_N_MC, K_e1);
    mc6_e1_errCB = nan(P6_N_MC, K_e1);
    mc6_e1_errE  = nan(P6_N_MC, K_e1); % Novedad: Error Este
    mc6_e1_errN  = nan(P6_N_MC, K_e1); % Novedad: Error Norte
 
    % Seed 1 already done — reuse the reference run output
    mc6_e1_err2D(1,:) = sqrt(out_e1_ref.errENU(1,:).^2 + out_e1_ref.errENU(2,:).^2);
    mc6_e1_errE(1,:)  = out_e1_ref.errENU(1,:);
    mc6_e1_errN(1,:)  = out_e1_ref.errENU(2,:);
    
    if isfield(out_e1_ref,'truthAux') && isfield(out_e1_ref.truthAux,'clkBiasTrue')
        mc6_e1_errCB(1,:) = out_e1_ref.truthAux.clkBiasTrue(:)' - out_e1_ref.xhist(4,:);
    end
 
    for n = 2:P6_N_MC
        cfg_e1.seedMeas = BASE_SEED + 600 + n;
        
        % Llamamos directamente al solver principal en lugar del wrapper run_single
        r_full = real_strtg_sequential_EKF_fastrealistic(cfg_e1);
        
        % Extract data from the original source.
        mc6_e1_err2D(n,:) = sqrt(r_full.errENU(1,:).^2 + r_full.errENU(2,:).^2);
        mc6_e1_errE(n,:)  = r_full.errENU(1,:); % Este
        mc6_e1_errN(n,:)  = r_full.errENU(2,:); % Norte
        
        % Extract clock truth if available.
        if isfield(r_full,'truthAux') && isfield(r_full.truthAux,'clkBiasTrue')
            mc6_e1_errCB(n,:) = r_full.truthAux.clkBiasTrue(:)' - r_full.xhist(4,:);
        else
            mc6_e1_errCB(n,:) = NaN; 
        end
    end
 
    % Compute mean / std — both will be 1×K row vectors
    mu_e1_2D = mean(mc6_e1_err2D, 1, 'omitnan');
    sd_e1_2D = std (mc6_e1_err2D, 0, 1, 'omitnan');
    mu_e1_CB = mean(mc6_e1_errCB, 1, 'omitnan');
    sd_e1_CB = std (mc6_e1_errCB, 0, 1, 'omitnan');
    
    mu_e1_E  = mean(mc6_e1_errE, 1, 'omitnan');
    sd_e1_E  = std (mc6_e1_errE, 0, 1, 'omitnan');
    mu_e1_N  = mean(mc6_e1_errN, 1, 'omitnan');
    sd_e1_N  = std (mc6_e1_errN, 0, 1, 'omitnan');
 
    % ── Extract blackout intervals from trk ───────────────────────────────
    blackout_mask = trk_e1.blackout(:)';   % ensure row 1×K
    bo_starts = [];
    bo_ends   = [];
    in_bo     = false;
    for k = 1:K_e1
        if blackout_mask(k) && ~in_bo
            bo_starts(end+1) = t_e1(k); %#ok<AGROW>
            in_bo = true;
        elseif ~blackout_mask(k) && in_bo
            bo_ends(end+1) = t_e1(k - 1); %#ok<AGROW>
            in_bo = false;
        end
    end
    if in_bo
        bo_ends(end+1) = t_e1(end);
    end
    n_blackouts = numel(bo_starts);
    fprintf('    Detected %d blackout intervals.\n', n_blackouts);
 
    % Satellite-switch times (for vertical dashed markers)
    sat_id_vec    = trk_e1.sid(:)';
    switch_epochs = find(diff(sat_id_vec) ~= 0 & sat_id_vec(1:end-1) ~= 0 & sat_id_vec(2:end) ~= 0);
    switch_times  = t_e1(switch_epochs);
 
    % ── PLOT 6-1A: 2D & Clock Bias (Original) ─────────────────────────────
    fh6_1a = figure('Name','P6-1A Single-Channel (2D & Clock)', 'Color','w', 'Position',[50 50 1200 700]);
 
    % [Aquí va tu código original de ploteo para el Subplot 1 (2D) y Subplot 2 (Clock), 
    % lo mantengo idéntico para no romper nada]
    
    ax1 = subplot(2, 1, 1); hold on; grid on; box on;
    y_top_1 = max(mu_e1_2D + sd_e1_2D) * 1.15; if isnan(y_top_1) || y_top_1 <= 0, y_top_1 = 500; end
    for bi = 1:n_blackouts
        patch([bo_starts(bi), bo_ends(bi), bo_ends(bi), bo_starts(bi)], [0, 0, y_top_1, y_top_1], ...
              [1.0 0.85 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    end
    fill([t_e1, fliplr(t_e1)], [mu_e1_2D + sd_e1_2D, fliplr(mu_e1_2D - sd_e1_2D)], [0.1 0.4 0.9], 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'DisplayName', '±1\sigma band');
    plot(t_e1, mu_e1_2D, '-', 'Color', [0.05 0.30 0.80], 'LineWidth', 2.0, 'DisplayName', '2D Horizontal Error (mean)');
    yline(CONV_THR_M, 'r--', 'LineWidth', 1.2, 'DisplayName', sprintf('%.0f m threshold', CONV_THR_M));
    for si = 1:numel(switch_times), xline(switch_times(si), ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9, 'HandleVisibility', 'off'); end
    patch(nan, nan, [1.0 0.85 0.75], 'FaceAlpha', 0.55, 'EdgeColor', 'none', 'DisplayName', 'Handover blackout (predict-only)');
    ylabel('2D Horizontal Error (m)', 'FontSize', 11);
    title({'Phase 6-1  —  Single-Channel Sequential Receiver  (stickyBest, 180 s)', ...
           sprintf('\\sigma_D = %.0f Hz  |  \\sigma_{DR} = %.1f Hz/s  |  switchBlackout = 10 s  |  N_{MC} = %d', P6_SIGMA_D, P6_SIGMA_DR, P6_N_MC)}, 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 9);
 
    ax2 = subplot(2, 1, 2); hold on; grid on; box on;
    cb_lo = min(mu_e1_CB - sd_e1_CB); cb_hi = max(mu_e1_CB + sd_e1_CB);
    if isnan(cb_lo) || isnan(cb_hi) || cb_lo == cb_hi, cb_lo = -0.5;  cb_hi = 0.5; end
    y_bot_2 = cb_lo * 1.2; y_top_2 = cb_hi * 1.2;
    for bi = 1:n_blackouts
        patch([bo_starts(bi), bo_ends(bi), bo_ends(bi), bo_starts(bi)], [y_bot_2, y_bot_2, y_top_2, y_top_2], ...
              [1.0 0.85 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    end
    fill([t_e1, fliplr(t_e1)], [mu_e1_CB + sd_e1_CB, fliplr(mu_e1_CB - sd_e1_CB)], [0.85 0.30 0.10], 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'DisplayName', '±1\sigma band');
    plot(t_e1, mu_e1_CB, '-', 'Color', [0.75 0.15 0.05], 'LineWidth', 2.0, 'DisplayName', 'clkBias Error (mean)');
    yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for si = 1:numel(switch_times), xline(switch_times(si), ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9, 'HandleVisibility', 'off'); end
    xlabel('Time (s)', 'FontSize', 11); ylabel('Clock Bias Error (m/s)', 'FontSize', 11);
    title('Clock Bias Estimation Error — sawtooth recovery after each handover', 'FontSize', 11);
    legend('Location', 'northeast', 'FontSize', 9);
    linkaxes([ax1, ax2], 'x'); xlim([0 180]);
    save_fig(fh6_1a, 'P6_1A_SingleChannel_2D_Clock', FIG_DIR);


    % ── PLOT 6-1B: East & North Convergence (NUEVA FIGURA) ────────────────
    fh6_1b = figure('Name','P6-1B Single-Channel (East & North)', 'Color','w', 'Position',[100 100 1200 700]);

    % ---- Subplot 1: East Error ----
    ax3 = subplot(2, 1, 1); hold on; grid on; box on;
    e_lo = min(mu_e1_E - sd_e1_E); e_hi = max(mu_e1_E + sd_e1_E);
    if isnan(e_lo) || isnan(e_hi) || e_lo == e_hi, e_lo = -100; e_hi = 400; end
    y_bot_E = e_lo * 1.1; y_top_E = e_hi * 1.1;
    
    for bi = 1:n_blackouts
        patch([bo_starts(bi), bo_ends(bi), bo_ends(bi), bo_starts(bi)], [y_bot_E, y_bot_E, y_top_E, y_top_E], ...
              [1.0 0.85 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    end
    
    fill([t_e1, fliplr(t_e1)], [mu_e1_E + sd_e1_E, fliplr(mu_e1_E - sd_e1_E)], [0.2 0.8 0.2], 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'DisplayName', '±1\sigma band');
    plot(t_e1, mu_e1_E, '-', 'Color', [0.1 0.6 0.1], 'LineWidth', 2.0, 'DisplayName', 'East Error (mean)');
    yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for si = 1:numel(switch_times), xline(switch_times(si), ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9, 'HandleVisibility', 'off'); end
    
    ylabel('East Error (m)', 'FontSize', 11);
    title('East Error Convergence (Cross-Track)', 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 9);

    % ---- Subplot 2: North Error ----
    ax4 = subplot(2, 1, 2); hold on; grid on; box on;
    n_lo = min(mu_e1_N - sd_e1_N); n_hi = max(mu_e1_N + sd_e1_N);
    if isnan(n_lo) || isnan(n_hi) || n_lo == n_hi, n_lo = -100; n_hi = 400; end
    y_bot_N = n_lo * 1.1; y_top_N = n_hi * 1.1;
    
    for bi = 1:n_blackouts
        patch([bo_starts(bi), bo_ends(bi), bo_ends(bi), bo_starts(bi)], [y_bot_N, y_bot_N, y_top_N, y_top_N], ...
              [1.0 0.85 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    end
    
    fill([t_e1, fliplr(t_e1)], [mu_e1_N + sd_e1_N, fliplr(mu_e1_N - sd_e1_N)], [0.6 0.2 0.8], 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'DisplayName', '±1\sigma band');
    plot(t_e1, mu_e1_N, '-', 'Color', [0.4 0.1 0.6], 'LineWidth', 2.0, 'DisplayName', 'North Error (mean)');
    yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    for si = 1:numel(switch_times), xline(switch_times(si), ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9, 'HandleVisibility', 'off'); end
    
    xlabel('Time (s)', 'FontSize', 11);
    ylabel('North Error (m)', 'FontSize', 11);
    title('North Error Convergence (Along-Track)', 'FontSize', 11, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 9);
    
    linkaxes([ax3, ax4], 'x'); xlim([0 180]);
    save_fig(fh6_1b, 'P6_1B_SingleChannel_EastNorth', FIG_DIR);

    fprintf('    Exp 6-1 complete.\n');
% ── PLOT 6-1C: Skyplot of Tracked Satellites (NUEVA FIGURA) ───────────
    fprintf('    Generating Skyplot for handovers...\n');
    
    fh6_1c = figure('Name','P6-1C Skyview','Color','w','Position',[200 100 650 650]);
    pax = polaraxes;
    hold on;
    
    % Configure sky view.
    pax.ThetaZeroLocation = 'top'; 
    pax.ThetaDir = 'clockwise';    
    % El radio es la elevación: 0 en el centro (Zenit, 90º), 90 en el borde (Horizonte, 0º)
    pax.RLim = [0 90];             
    pax.RTick = [0 30 60 90];
    pax.RTickLabel = {'90° (Zenith)','60°','30°','0° (Horizon)'}; 
    pax.GridAlpha = 0.4;
    
    % Extract tracking history.
    sid_history = trk_e1.sid(:)';
    tracked_sats = unique(sid_history(sid_history > 0));
    
    % Satellite colours.
    colors = lines(length(tracked_sats));
    
    for i = 1:length(tracked_sats)
        s = tracked_sats(i);
        
        % Convert azimuth and elevation to degrees.
        % Nota: Si tu código los tiene en grados, quita rad2deg() y deg2rad()
        az_rad = out_e1_ref.sat.az(:, s);
        el_deg = rad2deg(out_e1_ref.sat.elev(:, s)); 
        
        % Transformar Elevación a Radio de ploteo (90º = centro, 0º = borde)
        r_plot = 90 - el_deg;
        
        % Draw full visible track in grey.
        valid_traj = el_deg > 0; % Solo por encima del horizonte
        polarplot(az_rad(valid_traj), r_plot(valid_traj), '-', ...
            'Color', [0.85 0.85 0.85], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        % 2. Resaltar los tramos donde este satélite fue rastreado
        active_idx = find(sid_history == s);
        
        if ~isempty(active_idx)
            % Trazado grueso para el tiempo activo
            polarplot(az_rad(active_idx), r_plot(active_idx), '.', ...
                'Color', colors(i,:), 'MarkerSize', 12, ...
                'DisplayName', sprintf('Tracked Sat ID %d', s));
            
            % Mark tracking start.
            polarplot(az_rad(active_idx(1)), r_plot(active_idx(1)), 'o', ...
                'MarkerSize', 9, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:), ...
                'HandleVisibility', 'off');
                
            % Label tracking start.
            text(az_rad(active_idx(1)), r_plot(active_idx(1)) + 5, sprintf('  Sat %d Start', s), ...
                'FontSize', 9, 'Color', colors(i,:), 'FontWeight', 'bold');
        end
    end
    
    title('Skyplot: LEO Satellite Trajectories & Sequential Handovers', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'bestoutside', 'FontSize', 10);
    
    save_fig(fh6_1c, 'P6_1C_SingleChannel_Skyplot', FIG_DIR);
 
    % ════════════════════════════════════════════════════════════════════
    %  EXPERIMENT 6-2 — Multi-Channel: N = 1, 2, 3, 4 (100 s)
    % ════════════════════════════════════════════════════════════════════
    fprintf('\n  --- Experiment 6-2: Multi-Channel N=1,2,3,4 (100 s) ---\n');
 
    n_channels = [1, 2, 3, 4];
    n_ch       = numel(n_channels);
 
    % Colours: red (worst) → amber → teal → green (best)
    ch_colors = {[0.85 0.10 0.10], ...
                 [0.90 0.60 0.10], ...
                 [0.10 0.65 0.85], ...
                 [0.05 0.65 0.25]};
 
    % Pre-allocate result matrices (sized after first run when K is known)
    ch_t      = [];            % 1×K time axis (shared)
    ch_mean2D = [];            % n_ch × K
    ch_std2D  = [];
    ch_meanE  = [];
    ch_meanN  = [];
    ch_kpi    = cell(n_ch, 1);
 
    for ni = 1:n_ch
        N_ch = n_channels(ni);
        fprintf('    N = %d channel(s)...\n', N_ch);
 
        cfg_eni = build_base_cfg();
        cfg_eni.duration_sec           = FIXED_DURATION;
        cfg_eni.dt                     = FIXED_DT;
        cfg_eni.sigma_doppler_Hz       = P6_SIGMA_D;
        cfg_eni.sigma_dopplerRate_Hzps = P6_SIGMA_DR;
        cfg_eni.meas.useDoppler        = true;
        cfg_eni.meas.useDopplerRate    = true;
        cfg_eni.tropo.enable           = false;
        cfg_eni.plot.enable            = false;
 
        % N-channel capture logic
        if N_ch == 1
            cfg_eni.capture.mode             = 'stickyBest';
            cfg_eni.capture.minDwell_s       = 30;
            cfg_eni.capture.switchBlackout_s = 3;
        else
            cfg_eni.capture.mode = 'allVisible';
            cfg_eni.maxSats      = N_ch;
        end
 
        % Nominal coherent clock
        cfg_eni.clk.truth.initBias_mps            = 0.10;
        cfg_eni.clk.truth.initDrift_mps2          = 5e-5;
        cfg_eni.clk.truth.sigmaBias_rw_mps_sqrt   = 5e-3;
        cfg_eni.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4;
        cfg_eni.P0_clkBias                        = (0.15)^2;
        cfg_eni.P0_clkDrift                       = (2.5e-7*2)^2;
        cfg_eni.clk.filt.sigmaBias_rw_mps_sqrt    = 1.0e-2;
        cfg_eni.clk.filt.sigmaDrift_rw_mps2_sqrt  = 5.0e-4;
 
        % Monte Carlo
        [mc_eE, mc_eN, mc_e2D, ~, ~] = ...
            run_mc(cfg_eni, P6_N_MC, BASE_SEED + 610 + ni*20, USE_PARFOR);
 
        % Time axis from a fresh single run (guaranteed clean scalar struct)
        tmp_ref_ni = run_single(cfg_eni);
        t_ni_vec   = tmp_ref_ni.t(:)';   % force 1×K row
        K_ni       = numel(t_ni_vec);
 
        % Size the matrices on first channel now that K is known
        if isempty(ch_t)
            ch_t      = t_ni_vec;
            ch_mean2D = nan(n_ch, K_ni);
            ch_std2D  = nan(n_ch, K_ni);
            ch_meanE  = nan(n_ch, K_ni);
            ch_meanN  = nan(n_ch, K_ni);
        end
 
        % Store rows (mean/std over MC seeds produce 1×K rows automatically)
        ch_mean2D(ni, :) = mean(mc_e2D, 1, 'omitnan');
        ch_std2D(ni,  :) = std (mc_e2D, 0, 1, 'omitnan');
        ch_meanE(ni,  :) = mean(mc_eE,  1, 'omitnan');
        ch_meanN(ni,  :) = mean(mc_eN,  1, 'omitnan');
 
        ch_kpi{ni} = compute_kpis(t_ni_vec, mc_eE, mc_eN, mc_e2D, ...
                                   SS_WINDOW_S, CONV_THR_M);
 
        fprintf('      SS 2D = %.1f m  |  SS East = %.1f m  |  tc = %.1f s  (%.0f%% conv)\n', ...
            ch_kpi{ni}.mean_ss_rms2D, ch_kpi{ni}.mean_ss_rmsE, ...
            ch_kpi{ni}.mean_tc,       100*ch_kpi{ni}.frac_converged);
    end
 
    % ── PLOT 6-2A: Overlay time-series ────────────────────────────────────
    fh6_2a = figure('Name','P6-2A Multi-Channel 2D Overlay', ...
                    'Color','w', 'Position',[50 50 1000 520]);
    hold on; grid on; box on;
 
    t_row = ch_t(:)';   % guaranteed 1×K
 
    for ni = 1:n_ch
        N_ch = n_channels(ni);
        col  = ch_colors{ni};
 
        % Row slices of matrices are always 1×K — no dot-index needed
        mu_r = ch_mean2D(ni, :);
        sd_r = ch_std2D(ni,  :);
 
        fill([t_row, fliplr(t_row)], ...
             [mu_r + sd_r, fliplr(mu_r - sd_r)], ...
             col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
 
        plot(t_row, mu_r, '-', 'Color', col, 'LineWidth', 2.2, ...
             'DisplayName', sprintf('%d-channel  (SS=%.0f m, tc=%.0f s)', ...
                 N_ch, ch_kpi{ni}.mean_ss_rms2D, ch_kpi{ni}.mean_tc));
    end
 
    yline(CONV_THR_M, 'k--', 'LineWidth', 1.3, ...
          'DisplayName', sprintf('%.0f m threshold', CONV_THR_M));
 
    % Vertical convergence-time markers
    for ni = 1:n_ch
        tc_ni = ch_kpi{ni}.mean_tc;
        if ~isinf(tc_ni) && tc_ni <= FIXED_DURATION
            xline(tc_ni, '--', 'Color', ch_colors{ni}, ...
                  'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
    end
 
    xlabel('Time (s)',                'FontSize', 12);
    ylabel('2D Horizontal Error (m)', 'FontSize', 12);
    title({'Phase 6-2  —  Multi-Channel Convergence Comparison', ...
           sprintf('\\sigma_D = %.0f Hz  |  \\sigma_{DR} = %.1f Hz/s  |  N_{MC} = %d seeds  |  Mean \\pm 1\\sigma', ...
               P6_SIGMA_D, P6_SIGMA_DR, P6_N_MC)}, ...
        'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 10);
 
    save_fig(fh6_2a, 'P6_2A_MultiChannel_Overlay', FIG_DIR);
 
    % ── PLOT 6-2B: Bar chart — SS error vs number of channels ─────────────
    fh6_2b = figure('Name','P6-2B SS Error vs Channels', ...
                    'Color','w', 'Position',[1060 50 640 520]);
    hold on; grid on; box on;
 
    ss_2D_vals = cellfun(@(k) k.mean_ss_rms2D, ch_kpi);
    ss_E_vals  = cellfun(@(k) k.mean_ss_rmsE,  ch_kpi);
    ss_N_vals  = cellfun(@(k) k.mean_ss_rmsN,  ch_kpi);
    tc_vals    = cellfun(@(k) k.mean_tc,        ch_kpi);
    tc_vals(isinf(tc_vals)) = FIXED_DURATION;  % cap at sim duration for display
 
    bh = bar(n_channels, [ss_2D_vals(:), ss_E_vals(:), ss_N_vals(:)], ...
             'grouped', 'FaceColor', 'flat');
    bh(1).CData = repmat([0.10 0.35 0.75], n_ch, 1);   % 2D  — dark blue
    bh(2).CData = repmat([0.85 0.10 0.10], n_ch, 1);   % East — red
    bh(3).CData = repmat([0.10 0.65 0.25], n_ch, 1);   % North — green
 
    % Data labels above the 2D bars
    for ni = 1:n_ch
        text(bh(1).XEndPoints(ni), ss_2D_vals(ni) + 0.8, ...
             sprintf('%.0f m', ss_2D_vals(ni)), ...
             'HorizontalAlignment', 'center', 'FontSize', 9, ...
             'FontWeight', 'bold', 'Color', [0.05 0.20 0.55]);
    end
 
    set(gca, 'XTick', n_channels, ...
        'XTickLabel', arrayfun(@(n) sprintf('%d-ch', n), n_channels, ...
                               'UniformOutput', false), ...
        'FontSize', 11);
    xlabel('Number of Simultaneous Channels', 'FontSize', 12);
    ylabel('Steady-State RMS Error (m)',       'FontSize', 12);
    title({'Phase 6-2  —  SS Error vs Number of Channels', ...
           sprintf('Last %.0f s window  |  \\sigma_D = %.0f Hz  |  \\sigma_{DR} = %.1f Hz/s', ...
               SS_WINDOW_S, P6_SIGMA_D, P6_SIGMA_DR)}, ...
        'FontSize', 12, 'FontWeight', 'bold');
    legend({'SS 2D Horiz', 'SS East', 'SS North'}, ...
           'Location', 'northeast', 'FontSize', 10);
 
    % Overlay convergence time on secondary y-axis
    yyaxis right
    plot(n_channels, tc_vals, 'k^-', 'LineWidth', 2.0, 'MarkerSize', 9, ...
         'MarkerFaceColor', 'k');
    ylabel('Convergence Time (s)', 'FontSize', 12, 'Color', 'k');
    ax_6b = gca;
    ax_6b.YAxis(2).Color = 'k';
    legend({'SS 2D Horiz', 'SS East', 'SS North', 'Conv. Time (s)'}, ...
           'Location', 'northeast', 'FontSize', 10);
 
    save_fig(fh6_2b, 'P6_2B_SSError_vs_Channels', FIG_DIR);
 
    % ── Console summary ───────────────────────────────────────────────────
    fprintf('\n  ── Phase 6 Summary ────────────────────────────────────────\n');
    fprintf('  %-10s  %-12s  %-12s  %-12s  %-12s\n', ...
        'Channels', 'SS 2D (m)', 'SS East (m)', 'SS North (m)', 'Conv (s)');
    fprintf('  %s\n', repmat('-', 1, 62));
    for ni = 1:n_ch
        kp     = ch_kpi{ni};
        tc_str = sprintf('%.1f', kp.mean_tc);
        if isinf(kp.mean_tc), tc_str = 'no conv.'; end
        fprintf('  %-10d  %-12.2f  %-12.2f  %-12.2f  %s\n', ...
            n_channels(ni), kp.mean_ss_rms2D, kp.mean_ss_rmsE, ...
            kp.mean_ss_rmsN, tc_str);
    end
 
    ratio_E  = ch_kpi{1}.mean_ss_rmsE  / ch_kpi{end}.mean_ss_rmsE;
    ratio_2D = ch_kpi{1}.mean_ss_rms2D / ch_kpi{end}.mean_ss_rms2D;
    tc_1ch   = ch_kpi{1}.mean_tc;
    tc_4ch   = ch_kpi{end}.mean_tc;
 
    fprintf('\n  East SS RMS improvement  (1-ch -> 4-ch): x%.1f\n', ratio_E);
    fprintf('  2D  SS RMS improvement  (1-ch -> 4-ch): x%.1f\n', ratio_2D);
    if ~isinf(tc_1ch) && ~isinf(tc_4ch)
        fprintf('  Convergence speedup      (1-ch -> 4-ch): x%.1f\n', tc_1ch / tc_4ch);
    elseif isinf(tc_1ch) && ~isinf(tc_4ch)
        fprintf('  1-ch never converges in %.0f s; 4-ch converges in %.1f s\n', ...
            FIXED_DURATION, tc_4ch);
    end
    fprintf('\n  Phase 6 complete.\n\n');

end  % if ismember(6, RUN_PHASES)

if ismember(7, RUN_PHASES)
    fprintf('══════════════════════════════════════════════════════════\n');
    fprintf('  PHASE 6C — Sequential Handovers: Geometric Diversity\n');
    fprintf('══════════════════════════════════════════════════════════\n');

    P6_SIGMA_D  = 50.0; % Ruido Doppler fijo
    test_DR     = [0.1, 1.0, 5.0]; % Niveles de ruido Doppler-Rate
    policies    = {'stickyBest', 'stickyDiverse'};
    P6_N_MC     = 15;   % Semillas por combinación

    % Store East-error mean and standard deviation.
    K_len = 300 / 0.1 + 1; 
    res_mu_East  = nan(length(policies), length(test_DR), K_len);
    res_std_East = nan(length(policies), length(test_DR), K_len);
    
    bo_starts_cell = cell(length(policies), 1);
    bo_ends_cell   = cell(length(policies), 1);

    for pi_idx = 1:length(policies)
        pol = policies{pi_idx};
        fprintf('\n  --- Evaluando Política: %s ---\n', pol);
        
        for dri = 1:length(test_DR)
            sigDR = test_DR(dri);
            fprintf('    Ejecutando σ_DR = %.1f Hz/s ... ', sigDR);
            
            c = build_base_cfg();
            c.duration_sec             = 300;
            c.dt                       = 0.1;
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

            % -- 1. Referencia para Tracking y Blackouts --
            c.seedMeas = BASE_SEED + 600;
            out_ref = real_strtg_sequential_EKF_fastrealistic(c);
            t_axis = out_ref.t(:)';
            
            if dri == 1
                trk_mask = out_ref.trk.blackout(:)';
                b_start = []; b_end = []; in_bo = false;
                for k = 1:length(trk_mask)
                    if trk_mask(k) && ~in_bo
                        b_start(end+1) = t_axis(k); in_bo = true;
                    elseif ~trk_mask(k) && in_bo
                        b_end(end+1) = t_axis(k-1); in_bo = false;
                    end
                end
                bo_starts_cell{pi_idx} = b_start;
                bo_ends_cell{pi_idx}   = b_end;
            end
            
            % -- 2. Monte Carlo --
            mc_errE = nan(P6_N_MC, length(t_axis));
            for n = 1:P6_N_MC
                c.seedMeas = BASE_SEED + 650 + pi_idx*100 + dri*20 + n;
                r_full = real_strtg_sequential_EKF_fastrealistic(c);
                mc_errE(n,:) = r_full.errENU(1,:); 
            end
            
            res_mu_East(pi_idx, dri, 1:length(t_axis))  = mean(mc_errE, 1, 'omitnan');
            res_std_East(pi_idx, dri, 1:length(t_axis)) = std(mc_errE, 0, 1, 'omitnan');
            fprintf('Completado.\n');
        end
    end

    % ── 3. PLOT RESULTS WITH SHADED ±1σ BANDS ───────────────────────────
    fh6c = figure('Name','P6C RAAN Geometry vs Noise (with Bounds)', 'Color','w', 'Position',[100 50 1200 850]);
    
    % Line colours.
    color_line = {[0.85 0.15 0.15], [0.15 0.65 0.15]}; % Rojo, Verde
    % Shaded-band colours.
    color_band = {[0.95 0.30 0.30], [0.30 0.85 0.30]}; 
    
    for dri = 1:length(test_DR)
        subplot(3, 1, dri);
        hold on; grid on; box on;
        
        y_max = 0;
        
        % Draw blackout intervals first.
        for pi_idx = 1:length(policies)
            b_start = bo_starts_cell{pi_idx};
            b_end   = bo_ends_cell{pi_idx};
            for bi = 1:length(b_start)
                patch([b_start(bi), b_end(bi), b_end(bi), b_start(bi)], ...
                      [-2000, -2000, 2000, 2000], ...
                      [1.0 0.9 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
            end
        end
        
        % Draw shaded error bands and mean curves.
        for pi_idx = 1:length(policies)
            mu_curve  = squeeze(res_mu_East(pi_idx, dri, 1:length(t_axis)));
            std_curve = squeeze(res_std_East(pi_idx, dri, 1:length(t_axis)));
            
            % ±1σ shaded band.
            fill([t_axis, fliplr(t_axis)], ...
                 [mu_curve + std_curve, fliplr(mu_curve - std_curve)], ...
                 color_band{pi_idx}, 'FaceAlpha', 0.20, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off'); % Ocultamos la sombra de la leyenda para que no moleste
            
            % Línea principal de la media
            plot(t_axis, mu_curve, '-', 'Color', color_line{pi_idx}, 'LineWidth', 2.0, ...
                 'DisplayName', sprintf('Policy: %s', policies{pi_idx}));
             
            % Set vertical scale after the initial transient.
            local_max = max(abs(mu_curve(t_axis > 30)) + std_curve(t_axis > 30));
            y_max = max(y_max, local_max); 
        end
        
        yline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        
        % Final subplot formatting.
        ylim([-y_max*1.2, y_max*1.2]);
        xlim([0 t_axis(end)]);
        ylabel('East Error (m)', 'FontWeight', 'bold');
        title(sprintf('Impact of Handover Geometry with \\sigma_{DR} = %.1f Hz/s', test_DR(dri)));
        
        if dri == 1
            legend('Location', 'northeast');
        end
        if dri == 3
            xlabel('Time (s)', 'FontWeight', 'bold');
        end
    end
    
    sgtitle({'Phase 6C: Cross-Track Error Recovery (\pm1\sigma Bounds)', ...
             'stickyBest (Same Plane) vs stickyDiverse (Crossing RAAN)'}, ...
             'FontSize', 15, 'FontWeight', 'bold');
    
    if exist('FIG_DIR', 'var')
        save_fig(fh6c, 'P6C_RAAN_Diversity_Bounds', FIG_DIR);
    end
end

%
%% GG ═══════════════════════════════════════════════════════════════════════════
%  PHASE 7 — Fair Sequential: Same-Plane vs Cross-Plane
%
%  MISSION:
%    Compare two sequential tracking strategies over 1200s with 4 satellites:
%
%    Study A (same-plane):  tracks 4 sats from RAAN ~20 deg consecutively
%                           [20  20  20  20]
%
%    Study B (cross-plane): alternates between two orbital planes
%                           [20   5  20   5]
%
%  FAIRNESS GUARANTEE:
%    - Exactly 4 satellites per study, no reuse
%    - Identical time slots: [0-300], [300-600], [600-900], [900-1200] s
%    - Identical blackout: 10 s at each handover
%    - Identical noise, clock, seeds
%    - Only the RAAN target per slot changes
%
%  HOW THE PLANNER WORKS (plan_fixed_raan_sequence):
%    For each slot it scans ALL satellites whose RAAN is within
%    raanTolerance_deg of the target, skips already-used satellites,
%    and picks the one with the highest coverage fraction inside that slot.
%    If no sat is within tolerance, it widens and picks the nearest.
%
%  NOTE: Run this first with P7_N_MC=1 and check the console output:
%    "Slot 1 [0-300s] target=20 Sat3 RAAN=20.1 cov=87%"
%    "Slot 2 [300-600s] target=5  Sat9 RAAN=5.4  cov=72%"
%    etc.
%  If cov% is 0 for any slot, lower raanTolerance_deg or change the targets
%  to match the actual RAAN groups shown in the [constellation] print.
% ═══════════════════════════════════════════════════════════════════════════

fprintf('══════════════════════════════════════════════════════════\n');
fprintf('  PHASE 7 — StickyBest Sequential EKF\n');
fprintf('══════════════════════════════════════════════════════════\n');

% ── Experiment constants ──────────────────────────────────────────────
P7_SIGMA_D      = 50.0;
P7_SIGMA_DR     = 1.0;
P7_N_MC         = 20;
P7_DURATION     = 1200;
P7_DT           = 0.1;
P7_BLACKOUT     = 10;
P7_SWITCH_ELEV  = 40;
P7_COLOR        = [0.10 0.35 0.80];
P7_LABEL        = 'StickyBest (switch below 40°)';
local_BASE_SEED = BASE_SEED;

% ── Base config ───────────────────────────────────────────────────────
cfg_p7 = build_base_cfg();
cfg_p7.duration_sec             = P7_DURATION;
cfg_p7.dt                       = P7_DT;
cfg_p7.sigma_doppler_Hz         = P7_SIGMA_D;
cfg_p7.sigma_dopplerRate_Hzps   = P7_SIGMA_DR;
cfg_p7.meas.useDoppler          = true;
cfg_p7.meas.useDopplerRate      = true;
cfg_p7.capture.mode             = 'stickyBest';
cfg_p7.capture.switchBlackout_s = P7_BLACKOUT;
cfg_p7.capture.switchElevDeg    = P7_SWITCH_ELEV;
cfg_p7.capture.elev_margin_deg  = 10;
cfg_p7.plot.enable              = false;
cfg_p7.capture.minPassAbove40_s = 200;

% Complete clock model
cfg_p7.clk.truth.initBias_mps            = 0.10;
cfg_p7.clk.truth.initDrift_mps2          = 5e-5;
cfg_p7.clk.truth.sigmaBias_rw_mps_sqrt   = 5e-3;
cfg_p7.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4;
cfg_p7.P0_clkBias                        = (0.15)^2;
cfg_p7.P0_clkDrift                       = (2.5e-7*2)^2;
cfg_p7.clk.filt.sigmaBias_rw_mps_sqrt    = 1.0e-2;
cfg_p7.clk.filt.sigmaDrift_rw_mps2_sqrt  = 5.0e-4;

% ── Pre-compute constellation once ────────────────────────────────────
fprintf('  Pre-computing constellation (shared across all %d runs)...\n', P7_N_MC);
cfg_cache          = cfg_p7;
cfg_cache.seedMeas = local_BASE_SEED + 700;
out_cache          = real_strtg_sequential_EKF_geom(cfg_cache);
cfg_p7.prebuilt_sat = out_cache.sat;
fprintf('  Cached: %d sats, K=%d epochs.\n', out_cache.sat.Nsat, numel(out_cache.t));

% ── Reference run ─────────────────────────────────────────────────────
cfg_ref          = cfg_p7;
cfg_ref.seedMeas = local_BASE_SEED + 700;
out_ref7         = real_strtg_sequential_EKF_geom(cfg_ref);

t7 = out_ref7.t(:)';
K7 = numel(t7);

% ── Monte Carlo storage ───────────────────────────────────────────────
mc7_errE  = nan(P7_N_MC, K7);
mc7_errN  = nan(P7_N_MC, K7);
mc7_err2D = nan(P7_N_MC, K7);

mc7_errE(1,:)  = out_ref7.errENU(1,:);
mc7_errN(1,:)  = out_ref7.errENU(2,:);
mc7_err2D(1,:) = sqrt(mc7_errE(1,:).^2 + mc7_errN(1,:).^2);

% ── Parallel Monte Carlo ──────────────────────────────────────────────
fprintf('  Running parfor MC (%d additional seeds)... ', P7_N_MC-1);
tic;
parfor n = 2:P7_N_MC
    worker_cfg          = cfg_p7;
    worker_cfg.seedMeas = local_BASE_SEED + 700 + n;

    r = real_strtg_sequential_EKF_geom(worker_cfg);

    mc7_errE(n,:)  = r.errENU(1,:);
    mc7_errN(n,:)  = r.errENU(2,:);
    mc7_err2D(n,:) = sqrt(r.errENU(1,:).^2 + r.errENU(2,:).^2);
end
fprintf('Done (%.1f s)\n', toc);

% ── KPIs ──────────────────────────────────────────────────────────────
kpi7 = compute_kpis(t7, mc7_errE, mc7_errN, mc7_err2D, SS_WINDOW_S, CONV_THR_M);
fprintf('  KPIs: SS 2D = %.1f m | tc = %.1f s\n', ...
    kpi7.mean_ss_rms2D, kpi7.mean_tc);

% ── Extract blackout intervals ────────────────────────────────────────
bo7.starts = [];
bo7.ends   = [];

bm = out_ref7.trk.blackout(:)';
in_bo = false;
for k = 1:numel(bm)
    if bm(k) && ~in_bo
        bo7.starts(end+1) = t7(k); %#ok<AGROW>
        in_bo = true;
    end
    if ~bm(k) && in_bo
        bo7.ends(end+1) = t7(k-1); %#ok<AGROW>
        in_bo = false;
    end
end
if in_bo
    bo7.ends(end+1) = t7(end);
end

% ── Extract handover times ────────────────────────────────────────────
sid_v = out_ref7.trk.sid(:)';
sw = find(diff(sid_v) ~= 0 & sid_v(1:end-1) ~= 0 & sid_v(2:end) ~= 0);
hw7 = t7(sw);

% ═══════════════════════════════════════════════════════════════════════
% PLOTS
% ═══════════════════════════════════════════════════════════════════════

% ── P7-A: 2D Horizontal Error ─────────────────────────────────────────
fh7a = figure('Color','w','Position',[50 50 1300 500]);
hold on; grid on; box on;

mu2D = mean(mc7_err2D,1,'omitnan');
sd2D = std(mc7_err2D,0,1,'omitnan');
y_top = max(mu2D + sd2D) * 1.15;
if ~isfinite(y_top) || y_top < 10
    y_top = 600;
end

for bi = 1:numel(bo7.starts)
    patch([bo7.starts(bi), bo7.ends(bi), bo7.ends(bi), bo7.starts(bi)], ...
          [0, 0, y_top, y_top], [1 0.92 0.82], ...
          'EdgeColor','none', 'FaceAlpha',0.45, 'HandleVisibility','off');
end

fill([t7, fliplr(t7)], [mu2D + sd2D, fliplr(max(0, mu2D - sd2D))], ...
     P7_COLOR, 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');

plot(t7, mu2D, '-', 'Color', P7_COLOR, 'LineWidth', 2.5, ...
     'DisplayName', sprintf('%s  |  SS = %.0f m  |  tc = %.0f s', ...
     P7_LABEL, kpi7.mean_ss_rms2D, kpi7.mean_tc));

for hw = hw7
    xline(hw, '--', 'Color', P7_COLOR, 'LineWidth', 0.9, 'HandleVisibility','off');
end

yline(CONV_THR_M, 'k--', 'LineWidth', 1.5, ...
      'DisplayName', sprintf('%.0f m threshold', CONV_THR_M));

patch(nan,nan,[1 0.92 0.82],'FaceAlpha',0.5,'EdgeColor','none', ...
      'DisplayName', sprintf('Blackout (%d s)', P7_BLACKOUT));

xlim([0, P7_DURATION]);
ylim([0, y_top]);
xlabel('Time (s)', 'FontSize',12, 'FontWeight','bold');
ylabel('2D Horizontal Error (m)', 'FontSize',12, 'FontWeight','bold');
title({'Phase 7 — StickyBest Sequential EKF', ...
       sprintf('\\sigma_D = %.0f Hz  |  \\sigma_{DR} = %.1f Hz/s  |  N_{MC} = %d  |  switch < %d°', ...
       P7_SIGMA_D, P7_SIGMA_DR, P7_N_MC, P7_SWITCH_ELEV)}, ...
      'FontSize',12, 'FontWeight','bold');
legend('Location','northeast','FontSize',10);
save_fig(fh7a, 'P7A_2D_Error_StickyBest', FIG_DIR);

% ── P7-B: East & North ────────────────────────────────────────────────
fh7b = figure('Color','w','Position',[100 100 1300 650]);
mc_fields = {mc7_errE, mc7_errN};
ylabs     = {'East Error (m)', 'North Error (m)'};

for ax_i = 1:2
    subplot(2,1,ax_i);
    hold on; grid on; box on;

    mu = mean(mc_fields{ax_i},1,'omitnan');
    sd = std(mc_fields{ax_i},0,1,'omitnan');
    yhi = max(abs(mu) + sd) * 1.2;
    if ~isfinite(yhi) || yhi < 10
        yhi = 300;
    end

    for bi = 1:numel(bo7.starts)
        patch([bo7.starts(bi), bo7.ends(bi), bo7.ends(bi), bo7.starts(bi)], ...
              [-yhi, -yhi, yhi, yhi], [1 0.92 0.82], ...
              'EdgeColor','none', 'FaceAlpha',0.45, 'HandleVisibility','off');
    end

    fill([t7, fliplr(t7)], [mu + sd, fliplr(mu - sd)], ...
         P7_COLOR, 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');

    plot(t7, mu, '-', 'Color', P7_COLOR, 'LineWidth', 2.2, 'DisplayName', P7_LABEL);

    for hw = hw7
        xline(hw, '--', 'Color', P7_COLOR, 'LineWidth', 0.8, 'HandleVisibility','off');
    end

    yline(0,'k:','LineWidth',0.8,'HandleVisibility','off');
    ylim([-yhi, yhi]);
    xlim([0, P7_DURATION]);
    ylabel(ylabs{ax_i}, 'FontSize',11, 'FontWeight','bold');

    if ax_i == 1
        legend('Location','northeast','FontSize',10);
    else
        xlabel('Time (s)', 'FontSize',11, 'FontWeight','bold');
    end
end

sgtitle({'Phase 7 — ENU Convergence', ...
         sprintf('%s | blackout = %d s | switch < %d°', ...
         P7_LABEL, P7_BLACKOUT, P7_SWITCH_ELEV)}, ...
        'FontSize',13,'FontWeight','bold');
save_fig(fh7b, 'P7B_EastNorth_StickyBest', FIG_DIR);

% ── P7-C: Skyview ─────────────────────────────────────────────────────
fh7c = figure('Color','w','Position',[150 150 900 650]);
pal = lines(8);

pax = polaraxes;
hold(pax,'on');
pax.ThetaZeroLocation = 'top';
pax.ThetaDir          = 'clockwise';
pax.RLim              = [0,90];
pax.RTick             = [0,30,60,90];
pax.RTickLabel        = {'90','60','30','0'};
pax.GridAlpha         = 0.3;

sid_v  = out_ref7.trk.sid(:)';
sat_i  = out_ref7.sat;
used_s = unique(sid_v(sid_v > 0));

for si = 1:numel(used_s)
    s   = used_s(si);
    col = pal(mod(si-1,size(pal,1))+1,:);

    el = sat_i.elev(:,s);
    az = sat_i.az(:,s);

    ak = el > 0;
    if any(ak)
        polarplot(pax, deg2rad(az(ak)), 90 - el(ak), '-', ...
            'Color', [0.82 0.82 0.82], 'LineWidth', 1, 'HandleVisibility','off');
    end

    tk = sid_v == s;
    if any(tk)
        polarplot(pax, deg2rad(az(tk)), 90 - el(tk), '-', ...
            'Color', col, 'LineWidth', 3.5, ...
            'DisplayName', sprintf('Sat%d  RAAN = %.0f°', s, sat_i.raan(s)));

        k0 = find(tk,1,'first');
        polarplot(pax, deg2rad(az(k0)), 90 - el(k0), 'o', ...
            'MarkerSize', 8, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
            'HandleVisibility','off');
    end
end

title(pax, 'Phase 7 — Skyview | Grey = full pass, Colour = tracked', ...
      'FontSize',11, 'FontWeight','bold');
legend(pax, 'Location','southoutside', 'FontSize',9, 'NumColumns',2);
save_fig(fh7c, 'P7C_Skyview_StickyBest', FIG_DIR);

% ── P7-D: Satellite timeline ──────────────────────────────────────────
fh7d = figure('Color','w','Position',[200 200 1300 280]);
hold on; grid on; box on;

y_center = 0.8;
y_half   = 0.28;

for si = 1:numel(used_s)
    s   = used_s(si);
    col = pal(mod(si-1,size(pal,1))+1,:);

    runs = find(diff([0, (sid_v == s), 0]));
    for ri = 1:2:numel(runs)-1
        ts = t7(runs(ri));
        te = t7(min(runs(ri+1), numel(t7)));
        if te > ts
            rectangle('Position', [ts, y_center-y_half, te-ts, 2*y_half], ...
                'FaceColor', col, 'EdgeColor', 'k', 'LineWidth', 0.5);

            if te - ts > 20
                text((ts+te)/2, y_center, sprintf('Sat%d\n%.0f°', s, out_ref7.sat.raan(s)), ...
                    'HorizontalAlignment','center', 'FontSize',7, 'FontWeight','bold');
            end
        end
    end
end

% Blackout in grey
bo_v = out_ref7.trk.blackout(:)';
bk = find(bo_v);
ts_ = nan;
for ki = 1:numel(bk)
    if isnan(ts_)
        ts_ = t7(bk(ki));
    end
    if ki == numel(bk) || bk(ki+1) > bk(ki) + 1
        te_ = t7(bk(ki));
        if te_ > ts_
            rectangle('Position',[ts_, y_center-y_half, te_-ts_, 2*y_half], ...
                'FaceColor',[0.75 0.75 0.75], 'EdgeColor','none');
        end
        ts_ = nan;
    end
end

xlim([0, P7_DURATION]);
ylim([0, 1.6]);
set(gca,'YTick',[],'FontSize',10);
xlabel('Time (s)', 'FontSize',11, 'FontWeight','bold');
title(sprintf('Phase 7 — Satellite Timeline  (blackout = %d s, switch < %d°)', ...
      P7_BLACKOUT, P7_SWITCH_ELEV), ...
      'FontSize',11, 'FontWeight','bold');
save_fig(fh7d, 'P7D_Timeline_StickyBest', FIG_DIR);

% ── P7-E: KPI Summary ────────────────────────────────────────────────
fh7e = figure('Color','w','Position',[250 250 800 420]);

ss2D = kpi7.mean_ss_rms2D;
ssE  = kpi7.mean_ss_rmsE;
ssN  = kpi7.mean_ss_rmsN;
tc   = kpi7.mean_tc;
if isinf(tc)
    tc = P7_DURATION;
end

subplot(1,2,1);
hold on; grid on; box on;
bh = bar([1], [ss2D, ssE, ssN], 'grouped');
bh(1).FaceColor = [0.10 0.35 0.75];
bh(2).FaceColor = [0.85 0.20 0.20];
bh(3).FaceColor = [0.10 0.65 0.20];
text(bh(1).XEndPoints, ss2D+0.5, sprintf('%.0f m', ss2D), ...
    'HorizontalAlignment','center', 'FontSize',9, 'FontWeight','bold');
set(gca, 'XTick', 1, 'XTickLabel', {P7_LABEL}, 'FontSize', 9);
ylabel('Steady-State RMS (m)', 'FontSize',11);
title('SS Accuracy', 'FontSize',11);
legend({'2D','East','North'}, 'Location','northeast', 'FontSize',9);

subplot(1,2,2);
hold on; grid on; box on;
bh2 = bar(1, tc, 'FaceColor', P7_COLOR);
if isinf(kpi7.mean_tc)
    tc_str = 'no conv.';
else
    tc_str = sprintf('%.0f s', tc);
end
text(1, tc+1, tc_str, 'HorizontalAlignment','center', ...
    'FontSize',9, 'FontWeight','bold');
set(gca, 'XTick', 1, 'XTickLabel', {P7_LABEL}, 'FontSize', 9);
yline(P7_DURATION, 'k--', 'LineWidth',1, 'HandleVisibility','off');
ylabel('Convergence Time (s)', 'FontSize',11);
title('Convergence Speed', 'FontSize',11);

sgtitle({'Phase 7 — KPI Summary', ...
         sprintf('N_{MC} = %d | switch < %d° | \\sigma_D = %.0f Hz | \\sigma_{DR} = %.1f Hz/s', ...
         P7_N_MC, P7_SWITCH_ELEV, P7_SIGMA_D, P7_SIGMA_DR)}, ...
        'FontSize',11, 'FontWeight','bold');
save_fig(fh7e, 'P7E_KPI_StickyBest', FIG_DIR);

fprintf('\n  Phase 7 complete.\n\n');

%end % ismember(7, RUN_PHASES)
%% ══════════════════════════════════════════════════════════════════════════
%  §P6 MULTICHANNEL ENHANCEMENTS — append these plots after Exp 6-2
%  (insert after the existing save_fig(fh6_2b, ...) call in Phase 6)
%
%  These additions do NOT modify existing Phase 6 logic.
%  They read from variables already computed in Exp 6-2:
%    ch_t, ch_mean2D, ch_std2D, ch_meanE, ch_meanN, ch_kpi, n_channels
% ══════════════════════════════════════════════════════════════════════════

% ── OPTIONAL: Insert this block after save_fig(fh6_2b, ...) in Phase 6 ───
% (wrapped in a conditional to stay safe)
if ismember(6, RUN_PHASES) && exist('ch_t','var') && ~isempty(ch_t)

    % ── PLOT 6-2C: Convergence RATE (dErr2D/dt) ───────────────────────────
    % Shows not just the error level but HOW FAST it decreases.
    % Computed as a causal 5-point finite-difference of the mean 2D error.
    fh6_2c = figure('Name','P6-2C Convergence Rate', ...
                    'Color','w','Position',[50 600 1000 480]);
    hold on; grid on; box on;

    t_row = ch_t(:)';
    dt6   = mean(diff(t_row));   % actual time step of multichannel run

    % Window half-width for smoothed derivative (epochs)
    DERIV_HALF_WIN = 5;

    for ni = 1:n_ch
        N_ch   = n_channels(ni);
        col    = ch_colors{ni};
        mu_2D  = ch_mean2D(ni,:);

        % Smooth derivative: central difference with a small averaging window
        dmu = gradient(mu_2D, dt6);   % m/s  (rate of error decrease is negative)

        % Apply Gaussian smoothing to reduce noise in the derivative
        win_len = 2*DERIV_HALF_WIN+1;
        kernel  = exp(-(-DERIV_HALF_WIN:DERIV_HALF_WIN).^2 / (2*(DERIV_HALF_WIN/2)^2));
        kernel  = kernel / sum(kernel);
        dmu_smooth = conv(dmu, kernel, 'same');

        plot(t_row, -dmu_smooth, '-', 'Color', col, 'LineWidth', 2.0, ...
             'DisplayName', sprintf('%d-ch', N_ch));
    end

    yline(0, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
    xlabel('Time (s)',                    'FontSize',12,'FontWeight','bold');
    ylabel('Convergence Rate (m/s)',      'FontSize',12,'FontWeight','bold');
    title({'Phase 6-2C — Rate of Error Convergence (−d/dt of 2D Error)', ...
           sprintf('\\sigma_D=%.0f Hz  |  \\sigma_{DR}=%.1f Hz/s  |  Positive = error decreasing', ...
               P6_SIGMA_D, P6_SIGMA_DR)}, 'FontSize',11,'FontWeight','bold');
    legend('Location','northeast','FontSize',10);
    xlim([0, FIXED_DURATION]);
    save_fig(fh6_2c, 'P6_2C_ConvergenceRate', FIG_DIR);

    % ── PLOT 6-2D: ENU Error triangle — East vs North scatter per channel ─
    % Shows the spatial error cloud at SS — good way to see RAAN-axis bias
    fh6_2d = figure('Name','P6-2D East-North Error Fingerprint', ...
                    'Color','w','Position',[1060 600 600 560]);
    hold on; grid on; box on; axis equal;

    ss_mask6 = t_row >= (t_row(end) - SS_WINDOW_S);

    for ni = 1:n_ch
        N_ch = n_channels(ni);
        col  = ch_colors{ni};
        mu_E = ch_meanE(ni, ss_mask6);
        mu_N = ch_meanN(ni, ss_mask6);
        % Compute 1σ ellipse from mean trajectory
        cov_EN = cov(mu_E(:), mu_N(:));
        if ~any(isnan(cov_EN(:))) && det(cov_EN)>0
            [Ve,De] = eig(cov_EN);
            theta   = linspace(0,2*pi,100);
            ell     = Ve * sqrt(De) * [cos(theta); sin(theta)];
            plot(mean(mu_E)+ell(1,:), mean(mu_N)+ell(2,:), '-', ...
                 'Color', col, 'LineWidth', 2.0, 'HandleVisibility','off');
        end
        plot(mean(mu_E), mean(mu_N), 'o', 'MarkerSize',10, ...
             'MarkerFaceColor',col,'MarkerEdgeColor','k','LineWidth',1.2, ...
             'DisplayName', sprintf('%d-ch  (μ_E=%.1f, μ_N=%.1f)', ...
                 N_ch, mean(mu_E), mean(mu_N)));
    end

    plot(0,0,'k+','MarkerSize',16,'LineWidth',2.5,'DisplayName','True Position');
    th_circ = linspace(0,2*pi,200);
    plot(CONV_THR_M*cos(th_circ), CONV_THR_M*sin(th_circ), 'k--', ...
         'LineWidth',1.5,'DisplayName',sprintf('%.0f m target',CONV_THR_M));

    xlabel('East Error (m)',  'FontSize',12,'FontWeight','bold');
    ylabel('North Error (m)', 'FontSize',12,'FontWeight','bold');
    title({'Phase 6-2D — Steady-State East/North Error Fingerprint', ...
           sprintf('Mean \\pm 1\\sigma ellipse  |  Last %.0f s  |  \\sigma_D=%.0f Hz', ...
               SS_WINDOW_S, P6_SIGMA_D)}, 'FontSize',11,'FontWeight','bold');
    legend('Location','best','FontSize',9);
    save_fig(fh6_2d, 'P6_2D_EastNorth_Fingerprint', FIG_DIR);

end  % multichannel enhancements


%% ══════════════════════════════════════════════════════════════════════════
%  §OPTIONAL VIDEO HELPER
%  Call this function AFTER Phase 6 runs to generate an animation.
%  It reads from workspace variables — does NOT re-run the solver.
%
%  Usage (from ANALYSIS2 workspace, after Phase 6 runs):
%    generate_multichannel_video(ch_t, ch_mean2D, ch_std2D, ...
%                                 n_channels, ch_colors, CONV_THR_M, ...
%                                 P6_SIGMA_D, P6_SIGMA_DR, FIG_DIR);
% ══════════════════════════════════════════════════════════════════════════

function generate_multichannel_video(ch_t, ch_mean2D, ch_std2D, ...
        n_channels, ch_colors, conv_thr_m, sig_D, sig_DR, fig_dir)
% GENERATE_MULTICHANNEL_VIDEO
% Creates an MP4 animation showing how 2D error converges over time for
% each channel count (1-ch, 2-ch, 3-ch, 4-ch).
%
% The video frame at time T shows the error curves drawn up to T,
% with a moving vertical cursor — allowing the viewer to clearly see
% the convergence advantage gained by adding channels.

    n_ch       = numel(n_channels);
    t_row      = ch_t(:)';
    K          = numel(t_row);

    % ── Video parameters ─────────────────────────────────────────────────
    FRAME_RATE = 30;    % [fps]  playback speed
    SKIP       = 5;     % process every SKIP-th epoch (speeds up rendering)
    frame_epochs = 1:SKIP:K;
    n_frames    = numel(frame_epochs);

    vid_path = fullfile(fig_dir, 'P6_MultiChannel_Convergence_Video.mp4');
    vw = VideoWriter(vid_path, 'MPEG-4');
    vw.FrameRate = FRAME_RATE;
    open(vw);

    fh_vid = figure('Name','Multichannel Video', 'Color','w', ...
                    'Position',[100 100 900 500], 'Visible','off');

    y_max = max(max(ch_mean2D + ch_std2D)) * 1.15;
    if isnan(y_max) || y_max <= 0, y_max = 600; end

    for fi = 1:n_frames
        k_now = frame_epochs(fi);
        t_now = t_row(k_now);

        clf(fh_vid,'reset');
        hold on; grid on; box on;

        % Draw curves up to current time
        for ni = 1:n_ch
            col   = ch_colors{ni};
            mu_r  = ch_mean2D(ni, 1:k_now);
            sd_r  = ch_std2D(ni,  1:k_now);
            t_sub = t_row(1:k_now);

            fill([t_sub, fliplr(t_sub)], [mu_r+sd_r, fliplr(max(0,mu_r-sd_r))], col, ...
                 'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
            plot(t_sub, mu_r, '-', 'Color',col, 'LineWidth',2.0, ...
                 'DisplayName', sprintf('%d-ch', n_channels(ni)));
        end

        % Convergence threshold
        yline(conv_thr_m,'k--','LineWidth',1.2,'HandleVisibility','off');

        % Moving time cursor
        xline(t_now,'r-','LineWidth',1.5,'HandleVisibility','off');

        xlim([0, t_row(end)]);
        ylim([0, y_max]);
        xlabel('Time (s)',                'FontSize',11);
        ylabel('2D Horizontal Error (m)', 'FontSize',11);
        title(sprintf('Phase 6-2 — Multi-Channel Convergence  (t = %.1f s)', t_now), ...
              'FontSize',12,'FontWeight','bold');
        legend('Location','northeast','FontSize',10);
        text(0.98, 0.92, sprintf('\\sigma_D=%.0f Hz  |  \\sigma_{DR}=%.1f Hz/s', sig_D, sig_DR), ...
             'Units','normalized','HorizontalAlignment','right','FontSize',9);

        writeVideo(vw, getframe(fh_vid));
    end

    close(vw);
    close(fh_vid);
    fprintf('  [VIDEO] Saved to: %s\n', vid_path);
end


%% ══════════════════════════════════════════════════════════════════════════
%  §PATCH NOTES FOR real_strtg_sequential_EKF_fastrealistic
%  Apply these minimal changes to add 'plannedSequence' support.
%  Minimal patch notes retained for provenance.
% ══════════════════════════════════════════════════════════════════════════
%
%  1. In build_constellation_fast — add RAAN to cands and sat:
%     ── inside the `if ok` block:
%        cands(nc).raan = tleData(freshIdx(ii)).satrec_data.raan;  % CHANGE
%     ── after building sat:
%        sat.raan = zeros(1,Nsat);                                 % CHANGE
%        for s=1:Nsat, sat.raan(s) = cands(s).raan; end            % CHANGE
%
%  2. In the main body, BEFORE the call to select_capture_stream_strict:
%
%     if strcmpi(cfg.capture.mode, 'plannedSequence')              % CHANGE
%         if ~isfield(cfg.capture,'satBudget')                     % CHANGE
%             cfg.capture.satBudget = 3;                           % CHANGE
%         end                                                      % CHANGE
%         [plannedIDs, plannedBO, usedSats, planInfo] = ...        % CHANGE
%             plan_fair_sat_sequence(sat, cfg.capture, cfg.dt, K); % CHANGE
%         cfg.capture.plannedSatIDs   = plannedIDs(:);             % CHANGE
%         cfg.capture.plannedBlackout = plannedBO(:);              % CHANGE
%         cfg.capture.usedSats        = usedSats;                  % CHANGE
%         cfg.capture.planInfo        = planInfo;                  % CHANGE
%     end                                                          % CHANGE
%
%  3. In select_capture_stream_strict — add early return at top:
%
%     if strcmpi(capture.mode, 'plannedSequence')                  % CHANGE
%         trk.sid      = capture.plannedSatIDs(:);                 % CHANGE
%         trk.blackout = capture.plannedBlackout(:);               % CHANGE
%         if isfield(capture,'usedSats'), trk.usedSats = capture.usedSats; end  % CHANGE
%         if isfield(capture,'planInfo'), trk.planInfo = capture.planInfo; end  % CHANGE
%         return;                                                  % CHANGE
%     end                                                          % CHANGE
%
%  4. Copy plan_fair_sat_sequence as a local function at the bottom
%     of real_strtg_sequential_EKF_fastrealistic.m (identical to the
%     version in real_strtg_sequential_EKF_geom.m).
%
%  That's all.  No other changes to the solver are needed.
%  The filter loop, measurement model, and clock model are untouched.










%% ══════════════════════════════════════════════════════════════════════════
%  §SUMMARY  GLOBAL KPI TABLE
% ══════════════════════════════════════════════════════════════════════════

fprintf('\n══════════════════════════════════════════════════════════════\n');
fprintf('  ANALYSIS COMPLETE — All figures saved to: %s/\n', FIG_DIR);
fprintf('══════════════════════════════════════════════════════════════\n');
fprintf('  INTERPRETATION NOTES:\n');
fprintf('  • Up (vertical) axis is excluded from ALL metrics.\n');
fprintf('    cfg.height.enable=true constrains U; its error is artificial.\n');
fprintf('  • East error ↔ cross-track observability.\n');
fprintf('  • North error ↔ along-track observability.\n');
fprintf('  • Tropospheric delay mismatch should project primarily into\n');
fprintf('    clkBias (elevation-correlated) not into North.\n');
fprintf('  • Doppler-rate measurements restore cross-track (East) observability\n');
fprintf('    even at high Doppler noise levels.\n');
fprintf('══════════════════════════════════════════════════════════════\n\n');

%% ══════════════════════════════════════════════════════════════════════════
%  END OF RUN_FULL_ANALYSIS.m
% ══════════════════════════════════════════════════════════════════════════

%% ══════════════════════════════════════════════════════════════════════════
%  §1  HELPER FUNCTIONS  (defined before any phase calls them)
% ══════════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────────
%  build_base_cfg()
%  Returns the baseline configuration struct, mirroring MAIN_baseline.m.
%  All phase-specific overrides are applied on top of this.
% ─────────────────────────────────────────────────────────────────────────
function cfg = build_base_cfg()
    % ── TLE File ──────────────────────────────────────────────────────────
    cfg.tleFile = fullfile('InputData', 'oneweb.tle');

    % ── Simulation time (fixed by mandate) ───────────────────────────────
    cfg.duration_sec  = 100;
    cfg.dt            = 0.01;
    cfg.startUTC      = [];
    cfg.stopUTC       = [];
    cfg.startFromTLE  = 'first';
    cfg.maxTLEAgeDays = 7;

    % ── Receiver location (Cranfield, UK) ─────────────────────────────────
    cfg.rx.lat_deg  = 52.06869;
    cfg.rx.lon_deg  = -0.60884;
    cfg.rx.alt_m    = 100;
    cfg.rx.vel_ecef = [0;0;0];
    cfg.rx.acc_ecef = [0;0;0];

    % ── Signal / frequency ────────────────────────────────────────────────
    cfg.fc_Hz  = 11.325e9;
    cfg.c      = 299792458;
    cfg.lambda = cfg.c / cfg.fc_Hz;

    % ── Constellation / geometry ──────────────────────────────────────────
    cfg.minElevDeg                 = 37;
    cfg.coarseStep                 = 10;
    cfg.minMaxElev                 = 30;
    cfg.maxSats                    = 12;
    cfg.constellationSelectionWindow_s =120;

    % ── Propagator ────────────────────────────────────────────────────────
    cfg.useSGP4         = false;
    cfg.allowJ2Fallback = true;

    % ── Initial position error / prior ────────────────────────────────────
    cfg.init.errENU_m = [400; 400; 0];
    cfg.init.sigENU_m = [500; 500; 10];

    % ── Initial covariance (clock) ─────────────────────────────────────────
    cfg.P0_clkBias  = (0.1*2)^2;      % [m/s]^2
    cfg.P0_clkDrift = (2.5e-7*2)^2;   % [m/s²]^2

    % ── Process noise ────────────────────────────────────────────────────
    cfg.sigma_pos_rw = 1e-4;          % position random walk [m/√s]

    % ── Clock model ───────────────────────────────────────────────────────
    cfg.clk.truth.initBias_mps            = 0.10;
    cfg.clk.truth.initDrift_mps2          = 5e-5;
    cfg.clk.truth.sigmaBias_rw_mps_sqrt   = 5e-3 * 2;
    cfg.clk.truth.sigmaDrift_rw_mps2_sqrt = 1.0e-4 * 2;

    cfg.clk.filt.sigmaBias_rw_mps_sqrt    = 1.0e-2;
    cfg.clk.filt.sigmaDrift_rw_mps2_sqrt  = 5.0e-4;

    cfg.clk.filter_init_bias_mps          = 0.0;
    cfg.clk.filter_init_drift_mps2        = 0.0;

    % ── Troposphere ───────────────────────────────────────────────────────
    cfg.tropo.enable          = false;        % disabled by default; phases override
    cfg.tropo.ztd_truth_m     = 2.5;
    cfg.tropo.ztd_model_m     = 2.5;
    cfg.tropo.el_floor_deg    = 5;
    cfg.tropo.map_epsilon_deg = 0.5;
    cfg.tropo.model_truth     = 'simple';
    cfg.tropo.model_filter    = 'simple';
    cfg.tropo.lat_deg         = cfg.rx.lat_deg;
    cfg.tropo.alt_m           = cfg.rx.alt_m;
    cfg.tropo.doy             = 180;
    cfg.tropo.P0_hPa          = 1013.25;
    cfg.tropo.T0_K            = 288.15;
    cfg.tropo.e0_hPa          = 11.0;

    % ── Measurement model / noise ─────────────────────────────────────────
    cfg.sigma_doppler_Hz       = 2.0;   % Doppler white noise [Hz]
    cfg.sigma_dopplerRate_Hzps = 0.5;   % Doppler-rate white noise [Hz/s]

    cfg.meas.enable_lighttime  = false;
    cfg.meas.colored_tau_s     = 0.5;
    cfg.meas.sigma_col0_mps    = 0.03;
    cfg.meas.Rinflate          = 2.0;
    cfg.meas.useDoppler        = true;
    cfg.meas.useDopplerRate    = false;   % Doppler-only by default
    cfg.meas.fdJacStep_m       = 0.10;

    % ── IEKF numerical tuning ─────────────────────────────────────────────
    cfg.filter.useNISGate           = true;
    cfg.filter.nisGate              = 16.0;
    cfg.filter.maxIEKFIter          = 10;
    cfg.filter.iekfTol              = 1e-5;
    cfg.filter.maxStep_m            = 25;
    cfg.filter.startupInflateEpochs = 50;
    cfg.filter.startupRmult         = 25;

    % ── Height pseudo-measurement ─────────────────────────────────────────
    cfg.height.enable  = true;
    cfg.height.sigma_m = 10.0;
    cfg.height.z0_m    = 0.0;

    % ── Capture policy ────────────────────────────────────────────────────
    cfg.capture.mode                     = 'stickyBest';
    cfg.capture.minDwell_s               = 30;
    cfg.capture.switchBlackout_s         = 3;
    cfg.capture.elev_margin_deg          = 20;
    cfg.capture.useServiceAngleFilter    = false;
    cfg.capture.maxServiceOffNadirDeg    = 45;
    cfg.capture.switchThreshold_lowElev  = 10;
    cfg.capture.switchThreshold_highElev = 15;
    cfg.capture.highElevThreshold        = 35;

    % ── Convergence metric parameters ─────────────────────────────────────
    cfg.analysis.conv.refWindow_s    = 15;
    cfg.analysis.conv.settleWindow_s = 10;
    cfg.analysis.conv.relTol         = 0.15;
    cfg.analysis.conv.absTol_m       = 5.0;
    cfg.analysis.conv.maxCoeffVar    = 0.25;
    cfg.analysis.conv.maxSlope_mps   = 1.0;
    cfg.analysis.conv.errFloor_m     = 5.0;

    % ── Output / reproducibility ──────────────────────────────────────────
    cfg.seedMeas    = 42;
    cfg.plot.enable = false;   % suppress all internal figures during batch runs
    cfg.verbose     = false;
end

% ─────────────────────────────────────────────────────────────────────────
%  run_single(cfg)
%  Runs one solver call and returns a compact result struct.
%  Extracts East, North, 2D horizontal errors and clock state errors.
% ─────────────────────────────────────────────────────────────────────────
function res = run_single(cfg)
    out = real_strtg_sequential_EKF_fastrealistic(cfg);
    t   = out.t(:)';              % 1×K

    % Position errors (rows 1=East, 2=North, 3=Up)
    errE   = out.errENU(1,:);    % 1×K  [m]
    errN   = out.errENU(2,:);    % 1×K  [m]
    err2D  = sqrt(errE.^2 + errN.^2);  % 1×K  [m]

    % Clock state errors
    estCB  = out.xhist(4,:);     % estimated clock bias  [m/s]
    estCD  = out.xhist(5,:);     % estimated clock drift [m/s²]

    if isfield(out,'truthAux') && isfield(out.truthAux,'clkBiasTrue')
        trueCB = out.truthAux.clkBiasTrue(:)';   % 1×K
        trueCD = out.truthAux.clkDriftTrue(:)';  % 1×K
        errCB  = trueCB - estCB;
        errCD  = trueCD - estCD;
    else
        errCB = nan(size(t));
        errCD = nan(size(t));
    end

    res.t      = t;
    res.errE   = errE;
    res.errN   = errN;
    res.err2D  = err2D;
    res.errCB  = errCB;
    res.errCD  = errCD;
end

% ─────────────────────────────────────────────────────────────────────────
%  run_mc(cfg_template, N, base_seed, use_parfor)
%  Runs N_MC seeds and returns pre-allocated MC arrays.
%  Each seed changes cfg.seedMeas; if Phase 1 uses clock regime sampling,
%  the caller passes a cell array of cfgs instead.
% ─────────────────────────────────────────────────────────────────────────
function [mc_errE, mc_errN, mc_err2D, mc_errCB, mc_errCD] = ...
        run_mc(cfg_template, N, base_seed, use_parfor)
    %  cfg_template : single cfg struct OR cell array (1×N) of cfg structs
    %  Returns Nsamples×K matrices (one row per seed).

    % Build per-seed cfg list
    cfgs = cell(N,1);
    for n = 1:N
        if iscell(cfg_template)
            c = cfg_template{n};
        else
            c = cfg_template;
        end
        c.seedMeas = base_seed + n;
        cfgs{n}    = c;
    end

    % Run one seed to get time dimension K
    tmp    = run_single(cfgs{1});
    K      = numel(tmp.t);

    mc_errE  = nan(N,K);
    mc_errN  = nan(N,K);
    mc_err2D = nan(N,K);
    mc_errCB = nan(N,K);
    mc_errCD = nan(N,K);

    % Store first result
    mc_errE(1,:)  = tmp.errE;
    mc_errN(1,:)  = tmp.errN;
    mc_err2D(1,:) = tmp.err2D;
    mc_errCB(1,:) = tmp.errCB;
    mc_errCD(1,:) = tmp.errCD;

    if use_parfor && N > 1
        % Parallel branch — pre-allocate temporaries for parfor
        pe  = nan(N,K); pn  = nan(N,K); p2d = nan(N,K);
        pcb = nan(N,K); pcd = nan(N,K);
        pe(1,:)  = tmp.errE;  pn(1,:)  = tmp.errN;
        p2d(1,:) = tmp.err2D; pcb(1,:) = tmp.errCB; pcd(1,:) = tmp.errCD;

        parfor n = 2:N
            r = run_single(cfgs{n});
            pe(n,:)  = r.errE;   pn(n,:)  = r.errN;
            p2d(n,:) = r.err2D;  pcb(n,:) = r.errCB; pcd(n,:) = r.errCD;
        end
        mc_errE  = pe;  mc_errN  = pn;  mc_err2D = p2d;
        mc_errCB = pcb; mc_errCD = pcd;
    else
        % Serial branch
        for n = 2:N
            r = run_single(cfgs{n});
            mc_errE(n,:)  = r.errE;   mc_errN(n,:)  = r.errN;
            mc_err2D(n,:) = r.err2D;  mc_errCB(n,:) = r.errCB;
            mc_errCD(n,:) = r.errCD;
        end
    end
end

% ─────────────────────────────────────────────────────────────────────────
%  compute_kpis(t, mc_errE, mc_errN, mc_err2D, ss_window_s, conv_thr_m)
%  Computes steady-state RMS and strict convergence time from MC arrays.
% ─────────────────────────────────────────────────────────────────────────
function kpi = compute_kpis(t, mc_errE, mc_errN, mc_err2D, ss_window_s, conv_thr_m)
    N = size(mc_errE,1);
    t = t(:)';

    % Steady-state mask (last ss_window_s seconds)
    ss_mask = t >= (t(end) - ss_window_s);

    ss_rmsE  = nan(N,1);
    ss_rmsN  = nan(N,1);
    ss_rms2D = nan(N,1);
    tc       = nan(N,1);

    for n = 1:N
        e  = mc_errE(n,:);
        no = mc_errN(n,:);
        h  = mc_err2D(n,:);

        % Steady-state RMS (over last ss_window_s s)
        ss_rmsE(n)  = rms(e(ss_mask));
        ss_rmsN(n)  = rms(no(ss_mask));
        ss_rms2D(n) = rms(h(ss_mask));

        % Strict convergence time: first epoch where err2D stays <= thr
        tc(n) = find_convergence_strict(t, h, conv_thr_m);
    end

    % Average KPIs across seeds
    kpi.mean_ss_rmsE  = mean(ss_rmsE,  'omitnan');
    kpi.mean_ss_rmsN  = mean(ss_rmsN,  'omitnan');
    kpi.mean_ss_rms2D = mean(ss_rms2D, 'omitnan');
    kpi.std_ss_rmsE   = std(ss_rmsE,   'omitnan');
    kpi.std_ss_rmsN   = std(ss_rmsN,   'omitnan');
    kpi.std_ss_rms2D  = std(ss_rms2D,  'omitnan');

    kpi.mean_tc       = mean(tc,        'omitnan');
    kpi.std_tc        = std(tc,         'omitnan');

    % Fraction that converged within duration
    kpi.frac_converged = mean(~isinf(tc));

    % Per-seed vectors (for Pareto / scatter plots)
    kpi.ss_rms2D_vec  = ss_rms2D;
    kpi.tc_vec        = tc;
end

% ─────────────────────────────────────────────────────────────────────────
%  find_convergence_strict(t, err2D, threshold)
%  Returns the FIRST time t* such that err2D(k) <= threshold for all k >= k*.
%  If the error never permanently stays below threshold, returns inf.
% ─────────────────────────────────────────────────────────────────────────
function tc = find_convergence_strict(t, err2D, threshold)
    below = (err2D <= threshold);
    tc    = inf;
    % Scan from the end backwards to find the longest suffix that is always below
    last_always_below = numel(t) + 1;   % index after which all are below
    for k = numel(t):-1:1
        if below(k)
            last_always_below = k;
        else
            break;  % as soon as one is above, stop
        end
    end
    if last_always_below <= numel(t)
        tc = t(last_always_below);
    end
end

% ─────────────────────────────────────────────────────────────────────────
%  plot_mc_timeseries(t, mc_data, label, color, fig_title, y_label)
%  Plots mean ± 1σ shaded band over time.  Returns figure handle.
% ─────────────────────────────────────────────────────────────────────────
function fh = plot_mc_timeseries(t, mc_data, label, color, fig_title, y_label)
    mu  = mean(mc_data, 1, 'omitnan');
    sig = std(mc_data,  0, 1, 'omitnan');

    fh = figure('Name', fig_title, 'Color','w');
    hold on; grid on; box on;

    % Shaded ±1σ band
    fill([t, fliplr(t)], [mu+sig, fliplr(mu-sig)], ...
        color, 'FaceAlpha', 0.25, 'EdgeColor','none', ...
        'DisplayName', [label ' ±1σ']);

    % Mean trace
    plot(t, mu, '-', 'Color', color, 'LineWidth', 2.0, 'DisplayName', [label ' mean']);

    xlabel('Time (s)'); ylabel(y_label);
    title(fig_title, 'Interpreter','none');
    legend('Location','best');
end


% ─────────────────────────────────────────────────────────────────────────
%  save_fig(fh, fname, fig_dir)
%  Saves figure to fig_dir/fname.png
% ─────────────────────────────────────────────────────────────────────────
function save_fig(fh, fname, fig_dir)
    fpath = fullfile(fig_dir, fname);
    exportgraphics(fh, [fpath '.png'], 'Resolution', 150);
end

    % Heatmap helper using imagesc
    function plot_heatmap(data, xvals, yvals, xlbl, ylbl, ttl, cblbl, fig_name, fig_dir)
        fh = figure('Name', ttl, 'Color','w','Position',[50 50 750 550]);
        imagesc(xvals, yvals, data);
        set(gca,'YDir','normal');
        colorbar; colormap(parula);
        xlabel(xlbl); ylabel(ylbl); title(ttl);
        cb = colorbar; cb.Label.String = cblbl;
        % Annotate cells
        for r = 1:size(data,1)
            for co = 1:size(data,2)
                text(xvals(co), yvals(r), sprintf('%.0f',data(r,co)), ...
                    'HorizontalAlignment','center','Color','w','FontSize',8,'FontWeight','bold');
            end
        end
        set(gca,'XTick',xvals,'YTick',yvals, ...
            'XTickLabel',arrayfun(@(x)sprintf('%.0f',x),xvals,'UniformOutput',false), ...
            'YTickLabel',arrayfun(@(y)sprintf('%.2f',y),yvals,'UniformOutput',false));
        save_fig(fh, fig_name, fig_dir);
    end


%% ── SAVE RESULTS ──────────────────────────────────────────────────────
nombre_archivo = sprintf('Resultados_Fase6div_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
save(nombre_archivo);
fprintf('¡Éxito! Todo el análisis se ha guardado en el archivo: %s\n', nombre_archivo);








%%
load('Resultados_Fases1a5_20260318_223900.mat'); 
%%
fprintf('Datos cargados. Generando gráficos 3D...\n');


% ── PLOT P4B-contour: Smooth Topographic Map (Interpolated & Log) ───
    fprintf('    Generating smooth, interpolated contour map (SS 2D Error)...\n');
    
    fh4b_smooth = figure('Name', 'P4B SS 2D Error Smooth Contour', 'Color', 'w', 'Position', [150 150 850 650]);
    hold on; grid on; box on;
    
    % --- 1. RECONSTRUCT ORIGINAL GRID ---
    [D_orig, DR_orig] = meshgrid(sigma_D_hm, sigma_DR_hm);
    
    % --- 2. CREATE DENSE GRID FOR INTERPOLATION (Smoothness) ---
    % Increase points from 5x6 to 100x100 for smooth lines
    N_interp = 100; 
    D_dense = linspace(min(sigma_D_hm), max(sigma_D_hm), N_interp);
    
    % Create dense log spacing for Doppler-Rate so interpolation is even
    DR_dense_log = linspace(log10(min(sigma_DR_hm)), log10(max(sigma_DR_hm)), N_interp);
    DR_dense = 10.^DR_dense_log; 
    
    [D_dense_grid, DR_dense_grid] = meshgrid(D_dense, DR_dense);
    
    % --- 3. APPLY 2D INTERPOLATION ('cubic' for best curvature) ---
    hm_ss2D_interp = interp2(D_orig, DR_orig, hm_ss2D, D_dense_grid, DR_dense_grid, 'cubic');
    
    % Define the specific iso-error levels (in meters)
    contour_levels = [10, 25, 30, 40, 45, 50, 60, 70, 80, 90, 100, 120, 150, 200, 250, 300, 400];
    
    % --- 4. PLOT SMOOTH CONTOURS ON DENSE GRID ---
    [C, h] = contourf(D_dense_grid, DR_dense_grid, hm_ss2D_interp, contour_levels, 'LineColor', 'k');
    
    % Add bold labels on contour lines
    clabel(C, h, 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
    
    % Plot the *actual* raw simulated data points as black dots
    plot(D_orig(:), DR_orig(:), 'k.', 'MarkerSize', 10, 'HandleVisibility', 'off');
    
    % --- 5. LOG SCALE & TICKS CONFIGURATION ---
    set(gca, 'YScale', 'log');
    % Force Y-ticks to show exactly your simulated values
    set(gca, 'YTick', sigma_DR_hm); 
    % Force ticks to display decimal format, not exponential (e.g., 0.05 not 5e-2)
    set(gca, 'YTickLabel', string(sigma_DR_hm)); 
    
    % --- 6. FORMATTING & LABELS (English) ---
    colormap(jet); % Or 'turbo'
    cb = colorbar;
    ylabel(cb, 'Steady-State 2D RMS Error (m)', 'FontSize', 12, 'FontWeight', 'bold');
    
    xlabel('Doppler Noise \sigma_D (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Doppler-Rate Noise \sigma_{DR} (Hz/s)  [Log Scale]', 'FontSize', 12, 'FontWeight', 'bold');
    title({'Phase 4B — Smooth Topographic Error Map', ...
           'Cubic Interpolated SS 2D Horizontal RMS Error'}, ...
           'FontSize', 15, 'FontWeight', 'bold');
    
    if exist('FIG_DIR', 'var')
        save_fig(fh4b_smooth, 'P4B_SmoothContourMap_SS2D_Log', FIG_DIR);
    end

%% ── 1. HEATMAP 2D PARA EL TIEMPO DE CONVERGENCIA (Súper Profesional) ──
%% ── SMOOTH 2D CONTOUR: CONVERGENCE TIME ──
fprintf('Generating Smooth 2D Contour Plot...\n');

% 1. Create grids
[X_raw, Y_raw] = meshgrid(sigma_D_hm, sigma_DR_hm);

X_fine = linspace(min(sigma_D_hm), max(sigma_D_hm), 200);
Y_fine_log = linspace(log10(min(sigma_DR_hm)), log10(max(sigma_DR_hm)), 200);
Y_fine = 10.^Y_fine_log;

[X_fine_grid, Y_fine_grid] = meshgrid(X_fine, Y_fine);

% 2. Interpolate using 'makima' to avoid spline overshoots (the weird blobs)
Z_fine = interp2(X_raw, Y_raw, hm_tc_relaxed, X_fine_grid, Y_fine_grid, 'makima');

% Clean up boundaries to avoid math artifacts
min_time = min(hm_tc_relaxed(:));
max_time = max(hm_tc_relaxed(:));
Z_fine(Z_fine < min_time) = min_time; 
Z_fine(Z_fine > max_time) = max_time;

% 3. Draw Contour
fh_contour2d = figure('Name', 'Convergence Contour', 'Color', 'w', 'Position', [150 150 850 650]);
hold on; box on;

contourf(X_fine_grid, Y_fine_grid, Z_fine, 80, 'LineColor', 'none');

% 4. Styling and labels
colormap(turbo);
cb = colorbar;
cb.Label.String = 'Convergence Time (s)';
cb.Label.FontSize = 11;
caxis([min_time, max_time]); % Fix colorbar limits to exact min/max

set(gca, 'YScale', 'log'); 
set(gca, 'YDir', 'normal'); % Low noise at the bottom
xticks(sigma_D_hm);
yticks(sigma_DR_hm);

xlabel('\sigma_D : Doppler Noise (Hz)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('\sigma_{DR} : Doppler-Rate Noise (Hz/s)', 'FontWeight', 'bold', 'FontSize', 11);
title('Flexible Convergence Time (Smooth Contour)', 'FontSize', 14, 'FontWeight', 'bold');

% 5. Mark ONLY the minimum value(s)
[min_r, min_c] = find(hm_tc_relaxed == min_time);

% Draw small black dots for all simulated points just for reference (no text)
plot(X_raw(:), Y_raw(:), 'k.', 'MarkerSize', 5);

% Draw a star and text ONLY for the minimum(s)
for i = 1:length(min_r)
    best_x = sigma_D_hm(min_c(i));
    best_y = sigma_DR_hm(min_r(i));
    
    plot(best_x, best_y, 'wp', 'MarkerSize', 14, 'MarkerFaceColor', 'k');
    
    txt = sprintf(' Min: %.1fs', min_time);
    text(best_x * 1.05, best_y, txt, 'FontSize', 10, 'Color', 'k', ...
        'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7], 'Margin', 1);
end


%% ── TOPOGRAPHIC CONTOUR: CONVERGENCE TIME ──
fprintf('Generating Topographic Contour Plot for Convergence...\n');

% 1. Create grids
[X_raw, Y_raw] = meshgrid(sigma_D_hm, sigma_DR_hm);

X_fine = linspace(min(sigma_D_hm), max(sigma_D_hm), 200);
Y_fine_log = linspace(log10(min(sigma_DR_hm)), log10(max(sigma_DR_hm)), 200);
Y_fine = 10.^Y_fine_log;

[X_fine_grid, Y_fine_grid] = meshgrid(X_fine, Y_fine);

% 2. Interpolate using 'makima' to avoid spline overshoots
Z_fine = interp2(X_raw, Y_raw, hm_tc_relaxed, X_fine_grid, Y_fine_grid, 'makima');

% Clean up boundaries
min_time = min(hm_tc_relaxed(:));
max_time = max(hm_tc_relaxed(:));

Z_fine(Z_fine < min_time) = min_time; 
% Avoid a flat maximum plateau in the contour routine.
Z_fine(Z_fine >= max_time) = max_time - 0.001; 

% 3. Draw Contour
fh_contour2d = figure('Name', 'Convergence Topographic', 'Color', 'w', 'Position', [150 150 850 650]);
hold on; box on;

% ── Key time levels [s] ────────────────────────────────────────────────
levels = [0.1, 1, 2, 5, 10, 15, 20, 25, 29, 30]; 

% Exclude the maximum level from contour lines to avoid artefacts.
levels_lines = levels(levels < max_time);

% A. Filled background contours.
contourf(X_fine_grid, Y_fine_grid, Z_fine, levels, 'LineColor', 'none', 'FaceAlpha', 0.85);

% B. Contour lines.
[C, h] = contour(X_fine_grid, Y_fine_grid, Z_fine, levels_lines, 'LineColor', 'k', 'LineWidth', 1.2);

% C. Line labels.
clabel(C, h, 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold', 'LabelSpacing', 250);

% 4. Styling and labels
colormap(turbo); 
cb = colorbar;
cb.Label.String = 'Convergence Time (s)';
cb.Label.FontSize = 11;
caxis([min_time, max_time]); 

set(gca, 'YScale', 'log'); 
set(gca, 'YDir', 'normal'); % Low noise at the bottom
xticks(sigma_D_hm);
yticks(sigma_DR_hm);

xlabel('\sigma_D : Doppler Noise (Hz)', 'FontWeight', 'bold', 'FontSize', 11);
ylabel('\sigma_{DR} : Doppler-Rate Noise (Hz/s)', 'FontWeight', 'bold', 'FontSize', 11);
title('Flexible Convergence Time (Topographic Contour)', 'FontSize', 14, 'FontWeight', 'bold');

% 5. Mark simulated points and minimums
plot(X_raw(:), Y_raw(:), 'k.', 'MarkerSize', 7);

[min_r, min_c] = find(hm_tc_relaxed == min_time);
for i = 1:length(min_r)
    best_x = sigma_D_hm(min_c(i));
    best_y = sigma_DR_hm(min_r(i));
    
    plot(best_x, best_y, 'wp', 'MarkerSize', 14, 'MarkerFaceColor', 'k');
    
    txt = sprintf(' Min: %.1fs', min_time);
    text(best_x * 1.05, best_y, txt, 'FontSize', 10, 'Color', 'k', ...
        'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7], 'Margin', 1);
end

%% ── 2. 3D SMOOTHING FOR SNAPSHOTS ────────────────────────────────────
fh_snap_smooth = figure('Name', 'Snapshots 3D Suaves', 'Color', 'w', 'Position', [150 150 1200 800]);

% High-resolution grid.
X_fine = linspace(min(sigma_D_hm), max(sigma_D_hm), 100);
Y_fine_log = linspace(log10(min(sigma_DR_hm)), log10(max(sigma_DR_hm)), 100);
Y_fine = 10.^Y_fine_log;
[Xq, Yq] = meshgrid(X_fine, Y_fine);
[X_D, Y_DR] = meshgrid(sigma_D_hm, sigma_DR_hm);

z_max = min(max(hm_snapshots(:)), 400);

for ti = 1:numel(target_times)
    subplot(2, 2, ti);
    hold on; grid on; box on;
    
    Z_raw = hm_snapshots(:, :, ti);
    Zq = interp2(X_D, Y_DR, Z_raw, Xq, Yq, 'spline');
    Zq(Zq < 0) = 0;
    
    surf(Xq, Yq, Zq, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    shading interp; material dull;
    
    % Overlay simulated grid points.
    scatter3(X_D(:), Y_DR(:), Z_raw(:), 30, 'k', 'filled', 'MarkerEdgeColor', 'w');
    
    colormap(turbo);
    caxis([0 z_max]);
    
    set(gca, 'YScale', 'log');
    yticks(sigma_DR_hm); xticks(sigma_D_hm);
    zlim([0 z_max]);
    
    xlabel('\sigma_D (Hz)'); ylabel('\sigma_{DR} (Hz/s)'); zlabel('Error 2D (m)');
    title(sprintf('Snapshot at = %d seg', target_times(ti)), 'FontWeight', 'bold');
    view(-40, 25);
    camlight; lighting gouraud; % Luces para que resalten las curvas
end

cb = colorbar('Position', [0.92 0.1 0.02 0.8]);
cb.Label.String = '2D RMS Error (m)';
sgtitle('Error Evolution', 'FontSize', 16, 'FontWeight', 'bold');

%% ── DRAW FOUR SNAPSHOTS (2D CONTOUR PANELS) ───────────────────────────
fh_snap_2d = figure('Name', '2D Transient Evolution', 'Color', 'w', 'Position', [150 150 1100 800]);

% Fine grid for colour interpolation.
X_fine = linspace(min(sigma_D_hm), max(sigma_D_hm), 200);
Y_fine_log = linspace(log10(min(sigma_DR_hm)), log10(max(sigma_DR_hm)), 200);
Y_fine = 10.^Y_fine_log;
[Xq, Yq] = meshgrid(X_fine, Y_fine);
[X_D, Y_DR] = meshgrid(sigma_D_hm, sigma_DR_hm);

z_max = 400;

for ti = 1:numel(target_times)
    subplot(2, 2, ti);
    hold on; box on;
    
    Z_raw = hm_snapshots(:, :, ti);
    Zq = interp2(X_D, Y_DR, Z_raw, Xq, Yq, 'makima');
    
    Zq(Zq < 0) = 0; 
    Zq(Zq > z_max) = z_max;
    
    contourf(Xq, Yq, Zq, 60, 'LineColor', 'none');
    
    % Overlay simulated grid points.
    plot(X_D(:), Y_DR(:), 'k.', 'MarkerSize', 6);
    
    % Axis formatting.
    set(gca, 'YScale', 'log');
    set(gca, 'YDir', 'normal');
    yticks(sigma_DR_hm); xticks(sigma_D_hm);
    
    colormap(turbo);
    caxis([0 z_max]);
    
    % Titles and edge labels.
    title(sprintf('Snapshot at t = %d s', target_times(ti)), 'FontWeight', 'bold', 'FontSize', 12);
    if ti == 3 || ti == 4
        xlabel('\sigma_D : Doppler Noise (Hz)', 'FontWeight', 'bold');
    end
    if ti == 1 || ti == 3
        ylabel('\sigma_{DR} : Doppler-Rate Noise (Hz/s)', 'FontWeight', 'bold');
    end
end

% Shared colour bar.
cb = colorbar('Position', [0.92 0.11 0.02 0.815]);
cb.Label.String = '2D RMS Error (m)';
cb.Label.FontSize = 11;
cb.Label.FontWeight = 'bold';

sgtitle('Error Evolution: Transient State Landscape', 'FontSize', 16, 'FontWeight', 'bold');

%% ── TOPOGRAPHIC MAPS FOR FOUR SNAPSHOTS ───────────────────────────────
% load('Resultados_Fase4C_Transitorio.mat');

fprintf('Generando Mapas Topográficos de Error...\n');

fh_snap_contour = figure('Name', '2D Contour Snapshots', 'Color', 'w', 'Position', [150 150 1100 800]);

% Fine grid for smooth contours.
X_fine = linspace(min(sigma_D_hm), max(sigma_D_hm), 200);
Y_fine_log = linspace(log10(min(sigma_DR_hm)), log10(max(sigma_DR_hm)), 200);
Y_fine = 10.^Y_fine_log;
[Xq, Yq] = meshgrid(X_fine, Y_fine);
[X_D, Y_DR] = meshgrid(sigma_D_hm, sigma_DR_hm);

z_max = 400; 

% ── Custom contour levels [m] ──────────────────────────────────────────
levels = [10, 25,30,40, 50, 60, 70, 80, 100, 150, 200, 300, 400]; 

for ti = 1:numel(target_times)
    subplot(2, 2, ti);
    hold on; box on;
    
    Z_raw = hm_snapshots(:, :, ti);
    Zq = interp2(X_D, Y_DR, Z_raw, Xq, Yq, 'makima');
    Zq(Zq < 0) = 0; 
    Zq(Zq > z_max) = z_max;
    
    contourf(Xq, Yq, Zq, levels, 'LineColor', 'none', 'FaceAlpha', 0.8);
    
    [C, h] = contour(Xq, Yq, Zq, levels, 'LineColor', 'k', 'LineWidth', 1.2);
    
    % Label contour lines.
    clabel(C, h, 'FontSize', 9, 'Color', 'k', 'FontWeight', 'bold', 'LabelSpacing', 250);
    
    % Overlay simulated grid points.
    plot(X_D(:), Y_DR(:), 'k.', 'MarkerSize', 7);
    
    % Axis formatting.
    set(gca, 'YScale', 'log');
    set(gca, 'YDir', 'normal'); 
    yticks(sigma_DR_hm); xticks(sigma_D_hm);
    
    colormap(turbo);
    caxis([0 z_max]);
    
    % Labels.
    title(sprintf('Snapshot at t = %d s', target_times(ti)), 'FontWeight', 'bold', 'FontSize', 12);
    if ti == 3 || ti == 4
        xlabel('\sigma_D : Doppler Noise (Hz)', 'FontWeight', 'bold');
    end
    if ti == 1 || ti == 3
        ylabel('\sigma_{DR} : Doppler-Rate (Hz/s)', 'FontWeight', 'bold');
    end
end

% Shared colour bar.
cb = colorbar('Position', [0.92 0.11 0.02 0.815]);
cb.Label.String = '2D RMS Error (m)';
cb.Label.FontSize = 11;
cb.Label.FontWeight = 'bold';

sgtitle('Error Evolution: Topographic Contour Maps', 'FontSize', 16, 'FontWeight', 'bold');


%% ── 2. PLOT RESULTS (1 ROW x 3 COLUMNS) ───────────────────────────────
fh6c_horiz = figure('Name','P6C RAAN Geometry (Horizontal)', 'Color','w', 'Position',[50 150 1500 450]);

color_line = {[0.85 0.15 0.15], [0.15 0.65 0.15]};
color_band = {[0.95 0.30 0.30], [0.30 0.85 0.30]};

for dri = 1:length(test_DR)
    subplot(1, 3, dri); 
    hold on; grid on; box on;
    
    % Draw blackout intervals.
    for pi_idx = 1:length(policies)
        b_start = bo_starts_cell{pi_idx};
        b_end   = bo_ends_cell{pi_idx};
        
        % Close any blackout that extends to the final epoch.
        if length(b_start) > length(b_end)
            b_end(end+1) = t_axis(end);
        end
        
        for bi = 1:length(b_start)
            patch([b_start(bi), b_end(bi), b_end(bi), b_start(bi)], ...
                  [-1000, -1000, 1000, 1000], ... 
                  [1.0 0.9 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
        end
    end
    
    % Draw shaded bands and mean curves.
    for pi_idx = 1:length(policies)
        mu_curve  = squeeze(res_mu_East(pi_idx, dri, 1:length(t_axis)));
        mu_curve  = mu_curve(:)';
        
        std_curve = squeeze(res_std_East(pi_idx, dri, 1:length(t_axis)));
        std_curve = std_curve(:)';
        
        % ±1σ shaded band.
        fill([t_axis, fliplr(t_axis)], ...
             [mu_curve + std_curve, fliplr(mu_curve - std_curve)], ...
             color_band{pi_idx}, 'FaceAlpha', 0.20, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off'); 
             
        % Mean curve.
        plot(t_axis, mu_curve, '-', 'Color', color_line{pi_idx}, 'LineWidth', 2.0, ...
             'DisplayName', sprintf('Policy: %s', policies{pi_idx}));
    end
    
    % Zero reference line.
    yline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    
    % Axis limits.
    ylim([-50, 500]);
    xlim([0 t_axis(end)]);
    
    % Labels and title.
    xlabel('Time (s)', 'FontWeight', 'bold');
    
    if dri == 1
        ylabel('East Error (m)', 'FontWeight', 'bold');
        legend('Location', 'northeast');
    end
    
    title(sprintf('Noise \\sigma_{DR} = %.1f Hz/s', test_DR(dri)), 'FontSize', 12);
end

sgtitle({'Phase 6C: Cross-Track Error Recovery (\pm1\sigma Bounds)', ...
         'stickyBest (Same Plane) vs stickyDiverse (Crossing RAAN)'}, ...
         'FontSize', 15, 'FontWeight', 'bold');

if exist('FIG_DIR', 'var')
    save_fig(fh6c_horiz, 'P6C_RAAN_Diversity_Horizontal', FIG_DIR);
end










%% ── 1. EXTRACT REFERENCE GEOMETRY ─────────────────────────────────────
fprintf('Extrayendo trayectorias del cielo para ambas políticas...\n');

% Same configuration as Phase 6C.
c_sky = build_base_cfg();
c_sky.duration_sec             = 300;
c_sky.dt                       = 0.1;
c_sky.sigma_doppler_Hz         = 50.0;
c_sky.sigma_dopplerRate_Hzps   = 1.0; 
c_sky.meas.useDoppler          = true;
c_sky.meas.useDopplerRate      = true;
c_sky.tropo.enable             = false;
c_sky.capture.switchBlackout_s = 10;
c_sky.capture.minDwell_s       = 30;
c_sky.seedMeas                 = BASE_SEED + 600;

% Policy 1: stickyBest.
c_sky.capture.mode = 'stickyBest';
out_best = real_strtg_sequential_EKF_fastrealistic(c_sky);

% Policy 2: stickyDiverse.
c_sky.capture.mode = 'stickyDiverse';
out_div = real_strtg_sequential_EKF_fastrealistic(c_sky);

%% ── 2. SKY-PLOT HELPER ────────────────────────────────────────────────
draw_sky = @(out_data, title_str) local_draw_sky(out_data, title_str);

%% ── 3. DRAW COMPARISON (1 ROW x 2 COLUMNS) ───────────────────────────
fh_sky = figure('Name', 'Skyview Comparison', 'Color', 'w', 'Position', [100 150 1100 550]);

% Subplot 1: stickyBest
ax1 = subplot(1, 2, 1, polaraxes);
draw_sky(out_best, 'Política: stickyBest (Same Plane)');

% Subplot 2: stickyDiverse
ax2 = subplot(1, 2, 2, polaraxes);
draw_sky(out_div, 'Política: stickyDiverse (Crossing RAAN)');

sgtitle('Geometric Diversity in Sequential LEO Handovers', 'FontSize', 15, 'FontWeight', 'bold');

if exist('FIG_DIR', 'var')
    save_fig(fh_sky, 'P6C_Skyview_Comparison', FIG_DIR);
end

fprintf('Skyplots generados con éxito.\n');

%% ── POLAR PLOTTING LOGIC ──────────────────────────────────────────────
function local_draw_sky(out_data, title_str)
    hold on;
    pax = gca;
    
    % Configure compass view.
    pax.ThetaZeroLocation = 'top'; 
    pax.ThetaDir = 'clockwise';    
    pax.RLim = [0 90];             
    pax.RTick = [0 30 60 90];
    pax.RTickLabel = {'90° (Zenith)','60°','30°','0° (Horizon)'}; 
    pax.GridAlpha = 0.3;
    
    % Extract history.
    trk_sid = out_data.trk.sid(:)';
    tracked_sats = unique(trk_sid(trk_sid > 0));
    colors = lines(length(tracked_sats));
    
    % Draw satellites.
    for i = 1:length(tracked_sats)
        s = tracked_sats(i);
        
        % Convert azimuth to radians.
        az_rad = deg2rad(out_data.sat.az(:, s));
        el_deg = out_data.sat.elev(:, s); 
        
        % Plotting radius: zenith 90 deg -> 0, horizon 0 deg -> 90.
        r_plot = 90 - el_deg;
        
        % Draw full visible track in grey.
        valid = el_deg > 0;
        polarplot(az_rad(valid), r_plot(valid), '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1.0);
        
        % Draw tracked segment.
        active_idx = find(trk_sid == s);
        if ~isempty(active_idx)
            polarplot(az_rad(active_idx), r_plot(active_idx), '.', 'Color', colors(i,:), 'MarkerSize', 10);
            
            % Mark tracking start.
            start_idx = active_idx(1);
            polarplot(az_rad(start_idx), r_plot(start_idx), 'p', 'MarkerSize', 12, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:));
                
            text(az_rad(start_idx), r_plot(start_idx) + 6, sprintf('  Sat %d', s), ...
                'FontSize', 10, 'Color', colors(i,:)*0.8, 'FontWeight', 'bold');
        end
    end
    title(title_str, 'FontSize', 12, 'FontWeight', 'bold');
end











%% ── 1. ENSURE BASIC VARIABLES ─────────────────────────────────────────
fprintf('Generando gráficos 6-2B (Barras Limpias) y 6-2C (Francotirador)...\n');

% Restore colours and channel vector if needed.
if ~exist('n_channels', 'var'), n_channels = [1, 2, 3, 4]; end
if ~exist('n_ch', 'var'), n_ch = numel(n_channels); end
if ~exist('ch_colors', 'var')
    ch_colors = {[0.85 0.10 0.10], [0.90 0.60 0.10], [0.10 0.65 0.85], [0.05 0.65 0.25]};
end

% Restore constants if needed.
if ~exist('CONV_THR_M', 'var'), CONV_THR_M = 50; end
if ~exist('SS_WINDOW_S', 'var'), SS_WINDOW_S = 15; end
if ~exist('P6_SIGMA_D', 'var'), P6_SIGMA_D = 50.0; end
if ~exist('P6_SIGMA_DR', 'var'), P6_SIGMA_DR = 1.0; end
if ~exist('FIG_DIR', 'var'), FIG_DIR = 'analysis_figures'; end
if ~exist(FIG_DIR, 'dir'), mkdir(FIG_DIR); end

% Rebuild time axis if needed.
if ~exist('ch_t', 'var') && exist('ch_meanE', 'var')
    ch_t = linspace(0, 100, size(ch_meanE, 2));
end


%% ── PLOT 6-2B: SS ERROR VS NUMBER OF CHANNELS ─────────────────────────
fh6_2b = figure('Name','P6-2B SS Error vs Channels', ...
                'Color','w', 'Position',[100 100 650 550]);
hold on; grid on; box on;

% Extract saved KPI values.
ss_2D_vals = cellfun(@(k) k.mean_ss_rms2D, ch_kpi);
ss_E_vals  = cellfun(@(k) k.mean_ss_rmsE,  ch_kpi);
ss_N_vals  = cellfun(@(k) k.mean_ss_rmsN,  ch_kpi);

% Draw grouped bar chart.
bh = bar(n_channels, [ss_2D_vals(:), ss_E_vals(:), ss_N_vals(:)], ...
         'grouped', 'FaceColor', 'flat');

% Custom colours: 2D, East, North.
bh(1).CData = repmat([0.10 0.35 0.75], n_ch, 1);   
bh(2).CData = repmat([0.85 0.10 0.10], n_ch, 1);   
bh(3).CData = repmat([0.10 0.65 0.25], n_ch, 1);   

% Add 2D error labels.
for ni = 1:n_ch
    text(bh(1).XEndPoints(ni), ss_2D_vals(ni) + 3.0, ...
         sprintf('%.0f m', ss_2D_vals(ni)), ...
         'HorizontalAlignment', 'center', 'FontSize', 10, ...
         'FontWeight', 'bold', 'Color', [0.05 0.20 0.55]);
end

set(gca, 'XTick', n_channels, ...
    'XTickLabel', arrayfun(@(n) sprintf('%d-ch', n), n_channels, 'UniformOutput', false), ...
    'FontSize', 11);
xlabel('Number of Simultaneous Channels', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Steady-State RMS Error (m)',       'FontSize', 12, 'FontWeight', 'bold');
title({'Phase 6-2  —  Steady-State Error vs Number of Channels', ...
       sprintf('\\sigma_D = %.0f Hz  |  \\sigma_{DR} = %.1f Hz/s', P6_SIGMA_D, P6_SIGMA_DR)}, ...
    'FontSize', 14, 'FontWeight', 'bold');
legend({'SS 2D Horizontal', 'SS East (Cross-Track)', 'SS North (Along-Track)'}, ...
       'Location', 'northeast', 'FontSize', 11);

save_fig(fh6_2b, 'P6_2B_SSError_vs_Channels_Clean', FIG_DIR);


%% ── PLOT 6-2C: EAST VS NORTH MAP ──────────────────────────────────────
fh6_2c = figure('Name','P6-2C Sniper Plot (East vs North)', ...
                'Color','w', 'Position',[800 100 700 700]);
hold on; grid on; box on; 

axis equal; 

% Draw target circle.
th = linspace(0, 2*pi, 150);
fill(CONV_THR_M * cos(th), CONV_THR_M * sin(th), [0.85 0.95 0.85], ...
     'EdgeColor', [0.2 0.7 0.2], 'LineWidth', 2.0, ...
     'FaceAlpha', 0.5, 'DisplayName', sprintf('Target (%.0f m Accuracy)', CONV_THR_M));

% Crosshair at the origin.
xline(0, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');

% Split transient and steady-state samples.
t_row = ch_t(:)';
ss_mask = t_row >= (t_row(end) - SS_WINDOW_S); 

% Draw spatial trajectories.
for ni = 1:n_ch
    N_ch = n_channels(ni);
    col  = ch_colors{ni};

    mu_E = ch_meanE(ni, :);
    mu_N = ch_meanN(ni, :);

    % A. Transient trajectory.
    plot(mu_E(~ss_mask), mu_N(~ss_mask), '-', 'Color', [col 0.5], ...
         'LineWidth', 2.0, 'HandleVisibility', 'off');

    % B. Steady-state samples.
    scatter(mu_E(ss_mask), mu_N(ss_mask), 35, col, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8, ...
        'DisplayName', sprintf('%d-Channel Path', N_ch));
        
    % C. Final point.
    plot(mu_E(end), mu_N(end), 'p', 'MarkerSize', 16, ...
         'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
         'HandleVisibility', 'off');
end

% Formatting and labels.
xlabel('East Error (m) — Cross-Track', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('North Error (m) — Along-Track', 'FontSize', 12, 'FontWeight', 'bold');
title({'Phase 6-2  —  The "Sniper" Plot', ...
       'Spatial Trajectories into the Accuracy Target'}, ...
       'FontSize', 15, 'FontWeight', 'bold');

% Focus on the target area.
xlim([-20, 160]); 
ylim([-20, 160]);

% Legend.
legend('Location', 'northeast', 'FontSize', 11);

save_fig(fh6_2c, 'P6_2C_Sniper_Plot', FIG_DIR);

fprintf('¡Gráficos generados con éxito!\n');

%% ── PLOT 6-2D: REAL-MAP POSITION TRAJECTORY ───────────────────────────
fprintf('Generando Gráfico de Francotirador sobre Mapa Satelital...\n');

% Cranfield coordinates.
lat0 = 52.06869; 
lon0 = -0.60884;

% Approximate conversion from metres to degrees.
m_per_deg_lat = 111111.0;
m_per_deg_lon = 111111.0 * cosd(lat0);

fh6_2d = figure('Name','P6-2D Real Map Sniper Plot', 'Color','w', 'Position',[850 100 700 700]);
gx = geoaxes;
hold on;

% Draw 50 m target ring in latitude/longitude.
th = linspace(0, 2*pi, 150);
circle_E = CONV_THR_M * cos(th);
circle_N = CONV_THR_M * sin(th);
circle_lat = lat0 + (circle_N / m_per_deg_lat);
circle_lon = lon0 + (circle_E / m_per_deg_lon);

geoplot(circle_lat, circle_lon, '-', 'Color', [0.1 0.8 0.1], 'LineWidth', 3.0, ...
        'DisplayName', sprintf('Target (%.0f m Radius)', CONV_THR_M));

% Mark true position.
geoplot(lat0, lon0, 'k+', 'MarkerSize', 10, 'LineWidth', 2.0, 'DisplayName', 'True Position');

% Steady-state mask.
t_row = ch_t(:)';
ss_mask = t_row >= (t_row(end) - SS_WINDOW_S); 

% Draw trajectories on the map.
for ni = 1:n_ch
    N_ch = n_channels(ni);
    col  = ch_colors{ni};

    mu_E = ch_meanE(ni, :);
    mu_N = ch_meanN(ni, :);
    
    % Convert trajectories from metres to latitude/longitude.
    lat_path = lat0 + (mu_N / m_per_deg_lat);
    lon_path = lon0 + (mu_E / m_per_deg_lon);

    % A. Transient trajectory.
    geoplot(lat_path(~ss_mask), lon_path(~ss_mask), '-', 'Color', [col 0.8], ...
         'LineWidth', 2.5, 'DisplayName', sprintf('%d-Channel Path', N_ch));

    % B. Steady-state samples.
    geoscatter(lat_path(ss_mask), lon_path(ss_mask), 30, col, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.9, 'HandleVisibility', 'off');
        
    % C. Final point.
    geoplot(lat_path(end), lon_path(end), 'p', 'MarkerSize', 16, ...
         'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
         'HandleVisibility', 'off');
end

% Satellite basemap.
geobasemap(gx, 'satellite'); 

% Formatting, zoom, and legend.
title({'Phase 6-2  —  Real Map Convergence', ...
       'Spatial Trajectories into the Accuracy Target'}, ...
       'FontSize', 15, 'FontWeight', 'bold');

% Approximate 200 m zoom around the target.
lat_lim = lat0 + ([-50, 200] / m_per_deg_lat); 
lon_lim = lon0 + ([-50, 200] / m_per_deg_lon); 
geolimits(lat_lim, lon_lim);

gx.FontSize = 11;
gx.Scalebar.Visible = 'on';

legend('Location', 'northeast', 'FontSize', 11, 'Color', 'w');

if exist('FIG_DIR', 'var')
    save_fig(fh6_2d, 'P6_2D_RealMap_Sniper', FIG_DIR);
end

%% ── PLOT 6-2D: REAL-MAP POSITION TRAJECTORY ───────────────────────────
fprintf('Generando Gráfico de Francotirador sobre Mapa Satelital...\n');

lat0 = 52.06869; lon0 = -0.60884;
m_per_deg_lat = 111111.0;
m_per_deg_lon = 111111.0 * cosd(lat0);

fh6_2d = figure('Name','P6-2D Real Map Sniper Plot', 'Color','w', 'Position',[850 100 700 700]);
gx = geoaxes; hold on;

% Draw filled target area.
th = linspace(0, 2*pi, 150)';
circle_E = CONV_THR_M * cos(th);
circle_N = CONV_THR_M * sin(th);
circle_lat = lat0 + (circle_N / m_per_deg_lat);
circle_lon = lon0 + (circle_E / m_per_deg_lon);

% Build and draw target shape.
target_poly = geopolyshape(circle_lat, circle_lon);
geoplot(target_poly, 'FaceColor', [0.1 0.8 0.1], 'FaceAlpha', 0.2, ...
        'EdgeColor', [0.1 0.8 0.1], 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Target (%.0f m)', CONV_THR_M));

geoplot(lat0, lon0, 'k+', 'MarkerSize', 10, 'LineWidth', 2.0, 'DisplayName', 'True Position');

t_row = ch_t(:)';
ss_mask = t_row >= (t_row(end) - SS_WINDOW_S); 

for ni = 1:n_ch
    N_ch = n_channels(ni);
    col  = ch_colors{ni};
    mu_E = ch_meanE(ni, :); mu_N = ch_meanN(ni, :);
    lat_path = lat0 + (mu_N / m_per_deg_lat);
    lon_path = lon0 + (mu_E / m_per_deg_lon);

    geoplot(lat_path(~ss_mask), lon_path(~ss_mask), '-', 'Color', [col 0.8], 'LineWidth', 2.5, 'DisplayName', sprintf('%d-Channel Path', N_ch));
    geoscatter(lat_path(ss_mask), lon_path(ss_mask), 30, col, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.9, 'HandleVisibility', 'off');
    geoplot(lat_path(end), lon_path(end), 'p', 'MarkerSize', 16, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

geobasemap(gx, 'satellite'); 
title({'Phase 6-2  —  Real Map Convergence'}, 'FontSize', 15, 'FontWeight', 'bold');

% Keep full target area in view.
lat_lim = lat0 + ([-70, 200] / m_per_deg_lat); 
lon_lim = lon0 + ([-70, 200] / m_per_deg_lon); 
geolimits(lat_lim, lon_lim);

gx.FontSize = 11; gx.Scalebar.Visible = 'on';
legend('Location', 'northeast', 'FontSize', 11, 'Color', 'w');
if exist('FIG_DIR', 'var'), save_fig(fh6_2d, 'P6_2D_RealMap_Sniper', FIG_DIR); end



%% ── PLOT 6-2C: EAST VS NORTH WITH ELLIPSES ────────────────────────────
fh6_2c = figure('Name','P6-2C Sniper Plot (East vs North)', ...
                'Color','w', 'Position',[800 100 700 700]);
hold on; grid on; box on; 
axis equal; 

% Draw filled target circle.
th = linspace(0, 2*pi, 150);
fill(CONV_THR_M * cos(th), CONV_THR_M * sin(th), [0.85 0.95 0.85], ...
     'EdgeColor', [0.2 0.7 0.2], 'LineWidth', 2.0, ...
     'FaceAlpha', 0.5, 'DisplayName', sprintf('Target (%.0f m)', CONV_THR_M));

xline(0, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0, 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');

t_row = ch_t(:)';
ss_mask = t_row >= (t_row(end) - SS_WINDOW_S); 

% Draw trajectories and filled ellipses.
for ni = 1:n_ch
    N_ch = n_channels(ni);
    col  = ch_colors{ni};

    mu_E = ch_meanE(ni, :);
    mu_N = ch_meanN(ni, :);

    % A. Transient trajectory.
    plot(mu_E(~ss_mask), mu_N(~ss_mask), '-', 'Color', [col 0.5], ...
         'LineWidth', 2.0, 'HandleVisibility', 'off');

    % B. Steady-state samples.
    scatter(mu_E(ss_mask), mu_N(ss_mask), 35, col, 'filled', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.8, ...
        'DisplayName', sprintf('%d-Channel Path', N_ch));
        
    % C. Final point.
    plot(mu_E(end), mu_N(end), 'p', 'MarkerSize', 16, ...
         'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
         'HandleVisibility', 'off');
         
    % D. One-sigma uncertainty ellipse.
    mu_E_ss = mean(mu_E(ss_mask));
    mu_N_ss = mean(mu_N(ss_mask));
    sig_E_ss = ch_kpi{ni}.mean_ss_rmsE;
    sig_N_ss = ch_kpi{ni}.mean_ss_rmsN;
    
    ell_E = mu_E_ss + sig_E_ss .* cos(th);
    ell_N = mu_N_ss + sig_N_ss .* sin(th);
    
    fill(ell_E, ell_N, col, 'FaceAlpha', 0.2, 'EdgeColor', col, ...
         'LineWidth', 1.5, 'HandleVisibility', 'off');
end

xlabel('East Error (m) — Cross-Track', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('North Error (m) — Along-Track', 'FontSize', 12, 'FontWeight', 'bold');
title({'Phase 6-2  —  The "Sniper" Plot', 'Spatial Trajectories into the Accuracy Target'}, 'FontSize', 15, 'FontWeight', 'bold');

% Keep full target circle in view.
xlim([-70, 160]); 
ylim([-70, 160]);

legend('Location', 'northeast', 'FontSize', 11);
save_fig(fh6_2c, 'P6_2C_Sniper_Plot', FIG_DIR);


%% ── PLOT 6-2D: REAL-MAP POSITION TRAJECTORY ───────────────────────────
fprintf('Generando Gráfico de Francotirador sobre Mapa Satelital...\n');

lat0 = 52.06869; lon0 = -0.60884;
m_per_deg_lat = 111111.0;
m_per_deg_lon = 111111.0 * cosd(lat0);

fh6_2d = figure('Name','P6-2D Real Map Sniper Plot', 'Color','w', 'Position',[850 100 700 700]);
gx = geoaxes; hold on;

% Draw filled target area.
th = linspace(0, 2*pi, 150)';
circle_E = CONV_THR_M * cos(th);
circle_N = CONV_THR_M * sin(th);
circle_lat = lat0 + (circle_N / m_per_deg_lat);
circle_lon = lon0 + (circle_E / m_per_deg_lon);

% Build and draw target shape.
target_poly = geopolyshape(circle_lat, circle_lon);
geoplot(target_poly, 'FaceColor', [0.1 0.8 0.1], 'FaceAlpha', 0.2, ...
        'EdgeColor', [0.1 0.8 0.1], 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Target (%.0f m)', CONV_THR_M));

geoplot(lat0, lon0, 'k+', 'MarkerSize', 10, 'LineWidth', 2.0, 'DisplayName', 'True Position');

t_row = ch_t(:)';
ss_mask = t_row >= (t_row(end) - SS_WINDOW_S); 

for ni = 1:n_ch
    N_ch = n_channels(ni);
    col  = ch_colors{ni};
    mu_E = ch_meanE(ni, :); mu_N = ch_meanN(ni, :);
    lat_path = lat0 + (mu_N / m_per_deg_lat);
    lon_path = lon0 + (mu_E / m_per_deg_lon);

    geoplot(lat_path(~ss_mask), lon_path(~ss_mask), '-', 'Color', [col 0.8], 'LineWidth', 2.5, 'DisplayName', sprintf('%d-Channel Path', N_ch));
    geoscatter(lat_path(ss_mask), lon_path(ss_mask), 30, col, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.9, 'HandleVisibility', 'off');
    geoplot(lat_path(end), lon_path(end), 'p', 'MarkerSize', 16, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

geobasemap(gx, 'satellite'); 
title({'Phase 6-2  —  Real Map Convergence'}, 'FontSize', 15, 'FontWeight', 'bold');

% Keep full target area in view.
lat_lim = lat0 + ([-70, 200] / m_per_deg_lat); 
lon_lim = lon0 + ([-70, 200] / m_per_deg_lon); 
geolimits(lat_lim, lon_lim);

gx.FontSize = 11; gx.Scalebar.Visible = 'on';
legend('Location', 'northeast', 'FontSize', 11, 'Color', 'w');
if exist('FIG_DIR', 'var'), save_fig(fh6_2d, 'P6_2D_RealMap_Sniper', FIG_DIR); end


%% ── PLOT 6-2E: Real Map + Steady-State 1σ Filled Uncertainty Ellipses ────
fprintf('  Generating Plot 6-2E: Real Map + 1σ Filled Ellipses...\n');

% Base variables are defined above.
th_ell = linspace(0, 2*pi, 300)';

circ_E   = CONV_THR_M .* cos(th_ell);
circ_N   = CONV_THR_M .* sin(th_ell);
circ_lat = lat0 + circ_N ./ m_per_deg_lat;
circ_lon = lon0 + circ_E ./ m_per_deg_lon;

ell_lat_ss = cell(n_ch, 1);
ell_lon_ss = cell(n_ch, 1);
for ni = 1:n_ch
    mu_E_ss  = mean(ch_meanE(ni, ss_mask));
    mu_N_ss  = mean(ch_meanN(ni, ss_mask));
    sig_E_ss = ch_kpi{ni}.mean_ss_rmsE;
    sig_N_ss = ch_kpi{ni}.mean_ss_rmsN;

    ell_E = mu_E_ss + sig_E_ss .* cos(th_ell);
    ell_N = mu_N_ss + sig_N_ss .* sin(th_ell);

    ell_lat_ss{ni} = lat0 + ell_N ./ m_per_deg_lat;
    ell_lon_ss{ni} = lon0 + ell_E ./ m_per_deg_lon;
end

% Safe auto-zoom.
all_lat = circ_lat(:); all_lon = circ_lon(:);
for ni = 1:n_ch
    all_lat = [all_lat; ell_lat_ss{ni}(:); (lat0 + ch_meanN(ni,:)' ./ m_per_deg_lat)];
    all_lon = [all_lon; ell_lon_ss{ni}(:); (lon0 + ch_meanE(ni,:)' ./ m_per_deg_lon)];
end
pad_lat = 50 / m_per_deg_lat;
pad_lon = 50 / m_per_deg_lon;
auto_lat_lim = [min(all_lat) - pad_lat,  max(all_lat) + pad_lat];
auto_lon_lim = [min(all_lon) - pad_lon,  max(all_lon) + pad_lon];

fh6_2e = figure('Name','P6-2E Real Map + Filled Ellipses', 'Color','w', 'Position',[850 100 750 750]);
gx = geoaxes(fh6_2e); hold(gx, 'on');

% Use geopolyshape without transposing vectors.
circ_poly = geopolyshape(circ_lat, circ_lon);
%geoplot(gx, circ_poly, 'FaceColor', [0.10 0.80 0.10], 'FaceAlpha', 0.18, ...
        %'EdgeColor', [0.10 0.80 0.10], 'LineWidth', 2.5, 'DisplayName', sprintf('Target (%.0f m)', CONV_THR_M));

geoplot(gx, lat0, lon0, 'r+', 'MarkerSize', 20, 'LineWidth', 3.5, 'DisplayName', 'True Position');

for ni = 1:n_ch
    N_ch = n_channels(ni); col = ch_colors{ni};
    ss_rms = ch_kpi{ni}.mean_ss_rms2D;
    lat_p = lat0 + ch_meanN(ni,:) ./ m_per_deg_lat;
    lon_p = lon0 + ch_meanE(ni,:) ./ m_per_deg_lon;

    geoplot(gx, lat_p(~ss_mask), lon_p(~ss_mask), '-', 'Color', [col 0.85], 'LineWidth', 2.5, 'DisplayName', sprintf('%d-ch', N_ch));
    geoscatter(gx, lat_p(ss_mask), lon_p(ss_mask), 18, col, 'filled', 'MarkerFaceAlpha', 0.45, 'HandleVisibility', 'off');
    geoplot(gx, lat_p(end), lon_p(end), 'p', 'MarkerSize', 15, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');

    % Filled ellipse.
    ell_poly = geopolyshape(ell_lat_ss{ni}, ell_lon_ss{ni});
    geoplot(gx, ell_poly, 'FaceColor', col, 'FaceAlpha', 0.20, ...
            'EdgeColor', col, 'LineWidth', 2.0, 'HandleVisibility', 'off');
end

geobasemap(gx, 'satellite');
geolimits(gx, auto_lat_lim, auto_lon_lim); 
gx.FontSize = 11; gx.Scalebar.Visible = 'on';
title(gx, {'Phase 6-2  —  Real Map + 1\sigma SS Ellipses'}, 'FontSize', 13, 'FontWeight', 'bold');
legend(gx, 'Location', 'northwest', 'FontSize', 10, 'Color', 'w', 'TextColor', 'k');

if exist('FIG_DIR', 'var'), save_fig(fh6_2e, 'P6_2E_RealMap_Filled_Ellipses', FIG_DIR); end
fprintf('  Plot 6-2E saved.\n');




%% ── PLOT 6-2E: Real Map + Steady-State 1σ Filled Uncertainty Ellipses ────
fprintf('  Generating Plot 6-2E: Real Map + 1σ Filled Ellipses...\n');

% Base variables are defined above.
th_ell = linspace(0, 2*pi, 300)';

circ_E   = CONV_THR_M .* cos(th_ell);
circ_N   = CONV_THR_M .* sin(th_ell);
circ_lat = lat0 + circ_N ./ m_per_deg_lat;
circ_lon = lon0 + circ_E ./ m_per_deg_lon;

ell_lat_ss = cell(n_ch, 1);
ell_lon_ss = cell(n_ch, 1);
for ni = 1:n_ch
    mu_E_ss  = mean(ch_meanE(ni, ss_mask));
    mu_N_ss  = mean(ch_meanN(ni, ss_mask));
    sig_E_ss = ch_kpi{ni}.mean_ss_rmsE;
    sig_N_ss = ch_kpi{ni}.mean_ss_rmsN;

    ell_E = mu_E_ss + sig_E_ss .* cos(th_ell);
    ell_N = mu_N_ss + sig_N_ss .* sin(th_ell);

    ell_lat_ss{ni} = lat0 + ell_N ./ m_per_deg_lat;
    ell_lon_ss{ni} = lon0 + ell_E ./ m_per_deg_lon;
end

% Safe auto-zoom.
all_lat = circ_lat(:); all_lon = circ_lon(:);
for ni = 1:n_ch
    all_lat = [all_lat; ell_lat_ss{ni}(:); (lat0 + ch_meanN(ni,:)' ./ m_per_deg_lat)];
    all_lon = [all_lon; ell_lon_ss{ni}(:); (lon0 + ch_meanE(ni,:)' ./ m_per_deg_lon)];
end
pad_lat = 50 / m_per_deg_lat;
pad_lon = 50 / m_per_deg_lon;
auto_lat_lim = [min(all_lat) - pad_lat,  max(all_lat) + pad_lat];
auto_lon_lim = [min(all_lon) - pad_lon,  max(all_lon) + pad_lon];

fh6_2e = figure('Name','P6-2E Real Map + Filled Ellipses', 'Color','w', 'Position',[850 100 750 750]);
gx = geoaxes(fh6_2e); hold(gx, 'on');

% Use geopolyshape without transposing vectors.
circ_poly = geopolyshape(circ_lat, circ_lon);
geoplot(gx, circ_lat, circ_lon, 'k--', 'LineWidth', 2.0, ...
        'DisplayName', sprintf('Target (%.0f m)', CONV_THR_M));

geoplot(gx, lat0, lon0, 'r+', 'MarkerSize', 20, 'LineWidth', 3.5, 'DisplayName', 'True Position');

for ni = 1:n_ch
    N_ch = n_channels(ni); col = ch_colors{ni};
    ss_rms = ch_kpi{ni}.mean_ss_rms2D;
    lat_p = lat0 + ch_meanN(ni,:) ./ m_per_deg_lat;
    lon_p = lon0 + ch_meanE(ni,:) ./ m_per_deg_lon;

    geoplot(gx, lat_p(~ss_mask), lon_p(~ss_mask), '-', 'Color', [col 0.85], 'LineWidth', 2.5, 'DisplayName', sprintf('%d-ch', N_ch));
    geoscatter(gx, lat_p(ss_mask), lon_p(ss_mask), 18, col, 'filled', 'MarkerFaceAlpha', 0.45, 'HandleVisibility', 'off');
    geoplot(gx, lat_p(end), lon_p(end), 'p', 'MarkerSize', 15, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');

    % Filled ellipse.
    ell_poly = geopolyshape(ell_lat_ss{ni}, ell_lon_ss{ni});
    geoplot(gx, ell_poly, 'FaceColor', col, 'FaceAlpha', 0.20, ...
            'EdgeColor', col, 'LineWidth', 2.0, 'HandleVisibility', 'off');
end

geobasemap(gx, 'satellite');
geolimits(gx, auto_lat_lim, auto_lon_lim); 
gx.FontSize = 11; gx.Scalebar.Visible = 'on';
title(gx, {'Phase 6-2  —  Real Map + 1\sigma SS Ellipses'}, 'FontSize', 13, 'FontWeight', 'bold');
legend(gx, 'Location', 'northwest', 'FontSize', 10, 'Color', 'w', 'TextColor', 'k');

if exist('FIG_DIR', 'var'), save_fig(fh6_2e, 'P6_2E_RealMap_Filled_Ellipses', FIG_DIR); end
fprintf('  Plot 6-2E saved.\n')



%% ═══════════════════════════════════════════════════════════════════════
%  PLOT 6-2F  —  Animated Convergence Video on Real Map (trajectories only)
% ═══════════════════════════════════════════════════════════════════════
fprintf('  Generating Plot 6-2F: Convergence Video...\n');

VIDEO_FPS    = 100;
FRAME_STRIDE = 1;
VIDEO_FILE   = fullfile(FIG_DIR, 'P6_Convergence_Video.mp4');

t_row    = ch_t(:)';
K        = numel(t_row);
k_frames = 1:FRAME_STRIDE:K;
N_frames = numel(k_frames);

% Pre-convert trajectories to lat/lon once
lat_paths = cell(n_ch, 1);
lon_paths = cell(n_ch, 1);
for ni = 1:n_ch
    lat_paths{ni} = lat0 + ch_meanN(ni,:) ./ m_per_deg_lat;
    lon_paths{ni} = lon0 + ch_meanE(ni,:) ./ m_per_deg_lon;
end

% ── Figure ───────────────────────────────────────────────────────────────
fh6_2f = figure('Name','P6-2F Convergence Video', ...
                'Color','w', 'Position',[100 50 800 800]);
gx_v = geoaxes(fh6_2f);
hold(gx_v, 'on');

% % Static: target circle
% geoplot(gx, circ_lat, circ_lon, 'k--', 'LineWidth', 2.0, ...
%         'DisplayName', sprintf('Target (%.0f m)', CONV_THR_M));

% Static: true position (RED)
geoplot(gx_v, lat0, lon0, 'r+', ...
        'MarkerSize',18, 'LineWidth',3.0, 'DisplayName','True Position');

% Animated handles
h_trail = gobjects(n_ch, 1);
h_dot   = gobjects(n_ch, 1);

for ni = 1:n_ch
    col    = ch_colors{ni};
    ss_rms = ch_kpi{ni}.mean_ss_rms2D;

    h_trail(ni) = geoplot(gx_v, NaN, NaN, '-', ...
                          'Color', [col 0.80], 'LineWidth', 2.5, ...
                          'DisplayName', sprintf('%d-ch (SS=%.0f m)', ...
                                                 n_channels(ni), ss_rms));

    h_dot(ni)   = geoplot(gx_v, NaN, NaN, 'o', ...
                          'MarkerSize', 10, ...
                          'MarkerFaceColor', col, ...
                          'MarkerEdgeColor', 'w', ...
                          'LineWidth', 1.2, ...
                          'HandleVisibility', 'off');
end

% Formatting + force tile load
geobasemap(gx_v, 'satellite');
geolimits(gx_v, lat_lim, lon_lim);
gx_v.FontSize     = 11;
gx_v.Scalebar.Visible = 'on';
legend(gx_v, 'Location','northwest', 'FontSize',10, 'Color','w', 'TextColor','k');
hTitle_v = title(gx_v, 'Phase 6-2  —  t = 0.00 s', ...
                 'FontSize', 13, 'FontWeight', 'bold');

fprintf('    Loading satellite basemap tiles...\n');
drawnow;

% ── VideoWriter ───────────────────────────────────────────────────────────
vw = VideoWriter(VIDEO_FILE, 'MPEG-4');
vw.FrameRate = VIDEO_FPS;
vw.Quality   = 92;
open(vw);

% ── Animation loop ────────────────────────────────────────────────────────
fprintf('    Recording frames');
for fi = 1:N_frames
    k = k_frames(fi);

    for ni = 1:n_ch
        set(h_trail(ni), ...
            'LatitudeData',  lat_paths{ni}(1:k), ...
            'LongitudeData', lon_paths{ni}(1:k));
        set(h_dot(ni), ...
            'LatitudeData',  lat_paths{ni}(k), ...
            'LongitudeData', lon_paths{ni}(k));
    end

    set(hTitle_v, 'String', ...
        sprintf('Phase 6-2  —  t = %.2f s  (\\sigma_D=%.0f Hz, \\sigma_{DR}=%.1f Hz/s)', ...
                t_row(k), P6_SIGMA_D, P6_SIGMA_DR));

    drawnow limitrate;
    writeVideo(vw, getframe(fh6_2f));

    if mod(fi, round(N_frames/10)) == 0
        fprintf(' %d%%', round(100*fi/N_frames));
    end
end

close(vw);
close(fh6_2f);
fprintf('\n  Video saved → %s\n', VIDEO_FILE);



















%% 1. COMPUTE 2D ERROR STANDARD DEVIATIONS FROM MONTE CARLO MATRICES
std_2D_do = std(mc_2D_do, 0, 1, 'omitnan');
std_2D_dr = std(mc_2D_dr, 0, 1, 'omitnan');

% 2. Create figure.
figure('Name', 'Phase 4A: 2D Horizontal Error', 'Color', 'w', 'Position', [150 150 850 450]);
hold on; grid on; box on;

% 3. Doppler-only curve.
fill([t_axis, fliplr(t_axis)], ...
     [mu_2D_do + std_2D_do, fliplr(max(0, mu_2D_do - std_2D_do))], ...
     [0.9 0.2 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t_axis, mu_2D_do, 'r-', 'LineWidth', 2, ...
     'DisplayName', sprintf('Doppler-Only (\\sigma_D = %.0f Hz)', CRIT_SIGMA_D));

% 4. Doppler+Rate curve.
fill([t_axis, fliplr(t_axis)], ...
     [mu_2D_dr + std_2D_dr, fliplr(max(0, mu_2D_dr - std_2D_dr))], ...
     [0.1 0.5 0.9], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t_axis, mu_2D_dr, 'b-', 'LineWidth', 2, ...
     'DisplayName', sprintf('Doppler+Rate (\\sigma_{DR} = 0.5 Hz/s)'));

% 5. Convergence threshold and visual details.
yline(CONV_THR_M, 'k--', 'LineWidth', 1.2, 'DisplayName', '50 m threshold');

% 6. Axes and labels.
xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('2D Horizontal Error (m)', 'FontSize', 11, 'FontWeight', 'bold');
title({'Phase 4A — The Rescue: 2D Horizontal Error', ...
       'Comparison of Mean Error and \pm1\sigma Uncertainty Bounds'}, ...
       'FontSize', 13, 'FontWeight', 'bold');

% Configure Y limit.
ylim([0 600]); 
legend('Location', 'northeast', 'FontSize', 10);


%% REPEAT SPECIFIC CALCULATION: POINT FOCUS
sigD_target  = 5;    % Doppler Noise [Hz]
sigDR_target = 0.01;   % Doppler-Rate Noise [Hz/s]
N_MC_Point   = 15;

% 1. Configure Phase 4 focus scenario.
c_point = build_base_cfg();
c_point.duration_sec           = 30;
c_point.dt                     = 0.01;
c_point.sigma_doppler_Hz       = sigD_target;
c_point.sigma_dopplerRate_Hzps = sigDR_target;
c_point.meas.useDoppler        = true;
c_point.meas.useDopplerRate    = true;
c_point.tropo.enable           = false;
c_point.clk.truth.initBias_mps = 0.10;
c_point.P0_clkBias             = (0.15)^2;

% 2. Run Monte Carlo with a separate seed.
fprintf('Ejecutando simulación para sigma_D=%.0f, sigma_DR=%.2f...\n', sigD_target, sigDR_target);
[mc_E, mc_N, mc_2D, ~, ~] = run_mc(c_point, N_MC_Point, 9999, false);

% 3. Compute final metrics.
t_ax_point = (0:0.01:30);
kpi_point = compute_kpis(t_ax_point, mc_E, mc_N, mc_2D, 5, 50);

% 4. Print selected-point results.
fprintf('\n--- RESULTADOS PARA EL PUNTO SELECCIONADO ---\n');
fprintf('RMSE 2D Horizontal (a los 30s): %.2f metros\n', kpi_point.mean_ss_rms2D);
fprintf('Tiempo de convergencia (tc):    %.2f segundos\n', kpi_point.mean_tc);
fprintf('---------------------------------------------\n');

% 1. Configure selected focus case.
idx_dr = 1; % sigDR = 0.01
idx_d  = 2; % sigD = 5

% Extract mean and approximate sigma.
mu_vals = squeeze(hm_snapshots(idx_dr, idx_d, :))'; 
std_vals = mu_vals * 0.15;

figure('Color','w','Position',[100 100 800 500]);
hold on; grid on; box on;

% Draw shaded sigma band.
fill([target_times, fliplr(target_times)], ...
     [mu_vals + std_vals, fliplr(max(0, mu_vals - std_vals))], ...
     [0.1 0.5 0.9], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Draw mean curve.
plot(target_times, mu_vals, '-ob', 'LineWidth', 2.5, ...
     'MarkerFaceColor', 'b', 'DisplayName', 'Doppler+Rate (Focus Case)');

yline(50, 'r--', 'LineWidth', 2, 'DisplayName', 'Target 50m');
set(gca, 'YScale', 'log');
xlabel('Time (s)'); ylabel('2D RMSE (m)');
title('Phase 4 Focus: Convergence with \pm1\sigma Uncertainty');
xticks(target_times); ylim([20 300]);
legend('Location','northeast');


% --- MAP DATA GENERATOR ---
fprintf('Generando trayectoria para el mapa (Punto Focus)...\n');

% 1. Define selected map case.
sigD_map  = 5;
sigDR_map = 0.01; 

% 2. Create quick configuration.
c_map = build_base_cfg();
c_map.duration_sec = 30;
c_map.sigma_doppler_Hz = sigD_map;
c_map.sigma_dopplerRate_Hzps = sigDR_map;
c_map.meas.useDoppler = true;
c_map.meas.useDopplerRate = true;

% 3. Run one visual map case.
out_map = run_single(c_map); 

% 4. Extract map trajectory.
mu_E_dr = out_map.errENU(1,:);
mu_N_dr = out_map.errENU(2,:);
t_axis  = out_map.t;
FIXED_DT = c_map.dt;
% ── PLOT 4-MAP: Cranfield map trajectory ──────────────────────────────
fprintf('Generando Sniper Plot sobre Cranfield...\n');

% Cranfield University coordinates.
lat0 = 52.06869; 
lon0 = -0.60884;

% Conversion from metres to degrees.
m_per_deg_lat = 111111.0;
m_per_deg_lon = 111111.0 * cosd(lat0);

fh_map = figure('Name','P4 Real Map Sniper Plot', 'Color','w', 'Position',[100 100 850 750]);
gx = geoaxes;
hold on;

% 1. Draw 50 m target.
th = linspace(0, 2*pi, 100);
c_lat = lat0 + (CONV_THR_M * sin(th) / m_per_deg_lat);
c_lon = lon0 + (CONV_THR_M * cos(th) / m_per_deg_lon);
geoplot(c_lat, c_lon, 'g-', 'LineWidth', 3, 'DisplayName', '50m Target Area');

% 2. Convert trajectory to latitude/longitude.
lat_path = lat0 + (mu_N_dr / m_per_deg_lat);
lon_path = lon0 + (mu_E_dr / m_per_deg_lon);

% Split transient and steady-state segments.
ss_idx = length(t_axis) - (15/FIXED_DT); 
geoplot(lat_path(1:ss_idx), lon_path(1:ss_idx), 'y-', 'LineWidth', 2, 'DisplayName', 'Transient Path');
geoplot(lat_path(ss_idx:end), lon_path(ss_idx:end), 'b-', 'LineWidth', 3, 'DisplayName', 'Steady State');

% 3. Final point and map formatting.
geoscatter(lat_path(end), lon_path(end), 100, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Final Estimate');
geobasemap satellite;

title({'Phase 4 — Spatial Convergence', 'Receiver position estimate over Cranfield Campus'});
geolimits(lat0 + [-0.002 0.002], lon0 + [-0.003 0.003]);
final_P = out_map.Phist{end}(1:2, 1:2); 
legend('Location', 'northeast');





%% PHASE 4 FOCUS: LINEAR CONVERGENCE AND MAP
target_sigD  = 2;   
target_sigDR = 0.01; 
n_mc_pro     = 50;    

fprintf('Generating linear scale plots for: sigma_D=%.0f, sigma_DR=%.2f...\n', target_sigD, target_sigDR);

% 1. Configure and run simulation.
c_pro = build_base_cfg();
c_pro.duration_sec = 10; 
c_pro.dt = 0.01;
c_pro.sigma_doppler_Hz = target_sigD;
c_pro.sigma_dopplerRate_Hzps = target_sigDR;
c_pro.meas.useDoppler = true;
c_pro.meas.useDopplerRate = true;
c_pro.clk.truth.initBias_mps = 0.10; 
c_pro.P0_clkBias = (0.15)^2;

[mc_E, mc_N, mc_2D, ~, ~] = run_mc(c_pro, n_mc_pro, 8888, false);
t_pro = (0:c_pro.dt:c_pro.duration_sec);

mu_2D  = mean(mc_2D, 1, 'omitnan');
std_2D = std(mc_2D, 0, 1, 'omitnan');

%% --- PLOT 1: Linear Convergence ---
figure('Color','w','Position', [100 100 800 450]);
hold on; grid on; box on;

% Sigma band.
fill([t_pro, fliplr(t_pro)], [mu_2D + std_2D, fliplr(max(0, mu_2D - std_2D))], ...
    [0.1 0.5 0.9], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Mean curve.
plot(t_pro, mu_2D, 'Color', [0.1 0.4 0.8], 'LineWidth', 2.5, ...
    'DisplayName', sprintf('IEKF Mean (\\sigma_D=%.0f, \\sigma_{DR}=%.2f)', target_sigD, target_sigDR));

% 50 m target line.
yline(50, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Accuracy Target (50 m)');

% Linear-axis configuration.
set(gca, 'YScale', 'linear', 'FontSize', 11); 
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('2D Horizontal Error (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Phase 4 Focus: Linear Convergence Analysis', 'FontSize', 14);
legend('Location', 'northeast');

% Focus the Y limit near the 50 m target.
ylim([0 250]); 
xlim([0 10]);

%% --- PLOT 2: East/North trajectory ---
lat0 = 52.06869; lon0 = -0.60884;
m_per_deg_lat = 111111.0; m_per_deg_lon = 111111.0 * cosd(lat0);

fh_map = figure('Name','P4 Real Map Sniper', 'Color','w', 'Position',[100 100 800 800]);
gx = geoaxes; hold on;
rad=2*pi;
% Target.
th = linspace(0,rad, 150);
geoplot(lat0 + (50*sin(th)/m_per_deg_lat), lon0 + (50*cos(th)/m_per_deg_lon), 'g-', 'LineWidth', 3, 'DisplayName', '50m Target');

% Mean trajectory.
mu_E = mean(mc_E, 1); mu_N = mean(mc_N, 1);
geoplot(lat0 + (mu_N/m_per_deg_lat), lon0 + (mu_E/m_per_deg_lon), 'b-', 'LineWidth', 3, 'DisplayName', 'Mean Estimation');

% Monte Carlo endpoint cloud.
lat_end = lat0 + (mc_N(:,end) / m_per_deg_lat);
lon_end = lon0 + (mc_E(:,end) / m_per_deg_lon);
geoscatter(lat_end, lon_end, 40, [0.1 0.5 0.9], 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Monte Carlo Final Positions');

geoplot(lat0, lon0, 'k+', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'True Position');
geobasemap satellite;
geolimits(lat0 + [-0.0008 0.0008], lon0 + [-0.0012 0.0012]); 
legend('Location', 'northeast');
