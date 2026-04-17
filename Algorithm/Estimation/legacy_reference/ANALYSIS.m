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

% ── Which phases to execute (set to [] to skip) ──────────────────────────
RUN_PHASES     = [1 2 3 4 5];

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

%% ══════════════════════════════════════════════════════════════════════════
%  §P5  PHASE 5 — Deep Filter Parameter Sweep (Pareto Optimisation)
% ══════════════════════════════════════════════════════════════════════════
%
%  Environment: σ_D = 100 Hz,  σ_DR = 1.0 Hz/s
%
%  Grid search over:
%    Rinflate  ∈ {1.0, 1.5, 2.5}
%    sigma_pos ∈ {1e-4, 1e-3, 1e-2}  [m/√s]
%    clkRegime ∈ {'tight', 'nominal', 'loose'}  (coherent P0/Q scaling)
%
%  Output: Pareto front — X = SS 2D RMS,  Y = Convergence Time
%          Optimal combination printed to console.
%
% ══════════════════════════════════════════════════════════════════════════

if ismember(5, RUN_PHASES)
    fprintf('══════════════════════════════════════════════════════════\n');
    fprintf('  PHASE 5 — Deep IEKF Parameter Sweep (Pareto Front)\n');
    fprintf('══════════════════════════════════════════════════════════\n');

    % ── Sweep grid ────────────────────────────────────────────────────────
    Rinflate_vals = [1.0, 1.5, 2.5];
    sigPos_vals   = [1e-4, 1e-3, 1e-2];

    % Clock regimes for the filter (coherent with nominal truth clock)
    clk_regimes = struct( ...
        'label',    {'Tight',   'Nominal', 'Loose'}, ...
        'P0_bias',  {(0.05)^2, (0.15)^2, (0.5)^2}, ...
        'P0_drift', {(1e-7)^2, (2.5e-6)^2, (1e-5)^2}, ...
        'sigBias',  {5e-3,     1e-2,     3e-2}, ...
        'sigDrift', {2e-4,     5e-4,     2e-3}  ...
    );
    nClk = numel(clk_regimes);

    nRI  = numel(Rinflate_vals);
    nSP  = numel(sigPos_vals);
    nTotal = nRI * nSP * nClk;

    fprintf('  Total parameter combinations: %d (× %d MC seeds = %d runs)\n', ...
        nTotal, N_MC, nTotal * N_MC);

    % Pre-allocate result table
    p5_results = struct( ...
        'Rinflate',  cell(nTotal,1), ...
        'sigPos',    cell(nTotal,1), ...
        'clkLabel',  cell(nTotal,1), ...
        'ss_rms2D',  cell(nTotal,1), ...
        'tc',        cell(nTotal,1), ...
        'ss_rmsE',   cell(nTotal,1), ...
        'ss_rmsN',   cell(nTotal,1) ...
    );

    combo_idx = 0;

    for ri = 1:nRI
        for spi = 1:nSP
            for cli = 1:nClk
                combo_idx = combo_idx + 1;
                Ri  = Rinflate_vals(ri);
                sp  = sigPos_vals(spi);
                cr  = clk_regimes(cli);

                fprintf('  [%2d/%d] Rinflate=%.1f | sigPos=%.0e | clk=%s ... ', ...
                    combo_idx, nTotal, Ri, sp, cr.label);

                c = build_base_cfg();
                c.duration_sec           = FIXED_DURATION;
                c.dt                     = FIXED_DT;
                c.sigma_doppler_Hz       = 100;
                c.sigma_dopplerRate_Hzps = 1.0;
                c.meas.useDoppler        = true;
                c.meas.useDopplerRate    = true;
                c.tropo.enable           = false;

                c.meas.Rinflate          = Ri;
                c.sigma_pos_rw           = sp;

                % Truth clock: nominal
                c.clk.truth.initBias_mps            = 0.10;
                c.clk.truth.initDrift_mps2          = 5e-5;
                c.clk.truth.sigmaBias_rw_mps_sqrt   = 5e-3;
                c.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4;

                % Filter clock: from grid
                c.P0_clkBias                        = cr.P0_bias;
                c.P0_clkDrift                       = cr.P0_drift;
                c.clk.filt.sigmaBias_rw_mps_sqrt    = cr.sigBias;
                c.clk.filt.sigmaDrift_rw_mps2_sqrt  = cr.sigDrift;

                [mc_errE5, mc_errN5, mc_err2D5, ~, ~] = ...
                    run_mc(c, N_MC, BASE_SEED + 500 + combo_idx*3, USE_PARFOR);

                if ~exist('t_axis5','var')
                    tmp5 = run_single(c);
                    t_axis5 = tmp5.t;
                end

                kpi5 = compute_kpis(t_axis5, mc_errE5, mc_errN5, mc_err2D5, ...
                                    SS_WINDOW_S, CONV_THR_M);

                p5_results(combo_idx).Rinflate = Ri;
                p5_results(combo_idx).sigPos   = sp;
                p5_results(combo_idx).clkLabel = cr.label;
                p5_results(combo_idx).ss_rms2D = kpi5.mean_ss_rms2D;
                p5_results(combo_idx).tc       = kpi5.mean_tc;
                p5_results(combo_idx).ss_rmsE  = kpi5.mean_ss_rmsE;
                p5_results(combo_idx).ss_rmsN  = kpi5.mean_ss_rmsN;

                fprintf('SS 2D=%.1f m | tc=%.1f s\n', kpi5.mean_ss_rms2D, kpi5.mean_tc);
            end
        end
    end

    % ── Compute Pareto front ──────────────────────────────────────────────
    all_ss  = [p5_results.ss_rms2D]';
    all_tc  = [p5_results.tc]';

    % Replace inf tc with duration + penalty for Pareto scoring
    all_tc_fin = all_tc;
    all_tc_fin(isinf(all_tc_fin)) = FIXED_DURATION * 1.5;

    % Normalise for Pareto identification
    ss_norm = (all_ss  - min(all_ss))  / (max(all_ss)  - min(all_ss) + eps);
    tc_norm = (all_tc_fin - min(all_tc_fin)) / (max(all_tc_fin) - min(all_tc_fin) + eps);

    % A point i is Pareto-optimal if no other point j dominates it
    % (i.e., j has lower or equal ss AND lower or equal tc, with at least one strict)
    nC = numel(all_ss);
    is_pareto = true(nC,1);
    for i = 1:nC
        for j = 1:nC
            if i == j, continue; end
            if (all_ss(j) <= all_ss(i)) && (all_tc_fin(j) <= all_tc_fin(i)) && ...
               (all_ss(j)  < all_ss(i)  || all_tc_fin(j)  < all_tc_fin(i))
                is_pareto(i) = false;
                break;
            end
        end
    end

    % Optimal: Pareto point with minimum weighted sum (equal weights)
    pareto_idx = find(is_pareto);
    [~, best_local] = min(ss_norm(pareto_idx) + tc_norm(pareto_idx));
    best_idx = pareto_idx(best_local);

    % ── PLOT P5: Pareto Front Scatter ─────────────────────────────────────
    fh5 = figure('Name','P5 Pareto Front','Color','w','Position',[50 50 900 600]);
    hold on; grid on; box on;

    % Colour by clock regime
    clk_colors = {[0.2 0.6 0.9], [0.1 0.8 0.3], [0.9 0.4 0.1]};
    clk_labels = {'Tight', 'Nominal', 'Loose'};
    for cli = 1:nClk
        mask = strcmp({p5_results.clkLabel}, clk_labels{cli});
        scatter(all_ss(mask & ~is_pareto), all_tc_fin(mask & ~is_pareto), 60, ...
            clk_colors{cli}, 'o', 'filled', 'MarkerFaceAlpha', 0.4, ...
            'DisplayName', ['Dominated — ' clk_labels{cli}]);
        scatter(all_ss(mask &  is_pareto), all_tc_fin(mask &  is_pareto), 120, ...
            clk_colors{cli}, 'd', 'filled', 'LineWidth', 1.5, ...
            'MarkerEdgeColor','k', ...
            'DisplayName', ['Pareto — ' clk_labels{cli}]);
    end

    % Highlight best point
    scatter(all_ss(best_idx), all_tc_fin(best_idx), 250, 'r', 'p', 'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.5,'DisplayName','Optimal');
    text(all_ss(best_idx) * 1.03, all_tc_fin(best_idx), ...
        sprintf('  Opt\n  R=%.1f\n  σ_p=%.0e\n  Clk=%s', ...
            p5_results(best_idx).Rinflate, p5_results(best_idx).sigPos, ...
            p5_results(best_idx).clkLabel), ...
        'FontSize',9,'Color','r','FontWeight','bold');

    xlabel('Steady-State 2D Horizontal RMS Error (m)');
    ylabel('Convergence Time (s)');
    title('Phase 5 — IEKF Parameter Pareto Front', 'FontWeight','bold');
    legend('Location','best','FontSize',8);
    save_fig(fh5, 'P5_Pareto_Front', FIG_DIR);

    % ── Print optimal configuration to console ────────────────────────────
    fprintf('\n  ╔══════════════════════════════════════════════════════╗\n');
    fprintf('  ║  PHASE 5 — OPTIMAL PARAMETER COMBINATION            ║\n');
    fprintf('  ╠══════════════════════════════════════════════════════╣\n');
    fprintf('  ║  Pareto-optimal (min weighted sum of normalised      ║\n');
    fprintf('  ║  SS-Error + Convergence-Time):                       ║\n');
    fprintf('  ║                                                      ║\n');
    fprintf('  ║  cfg.meas.Rinflate          = %.1f\n', p5_results(best_idx).Rinflate);
    fprintf('  ║  cfg.sigma_pos_rw           = %.1e\n', p5_results(best_idx).sigPos);
    fprintf('  ║  Clock regime               = %s\n',   p5_results(best_idx).clkLabel);
    fprintf('  ║  SS 2D Horizontal RMS       = %.2f m\n', p5_results(best_idx).ss_rms2D);
    fprintf('  ║  Convergence time           = %.1f s\n', p5_results(best_idx).tc);
    fprintf('  ╚══════════════════════════════════════════════════════╝\n\n');

    % Map clock regime label to actual cfg values
    opt_cli = find(strcmp(clk_labels, p5_results(best_idx).clkLabel));
    fprintf('  Optimal cfg values:\n');
    fprintf('    cfg.meas.Rinflate                    = %.4f\n', p5_results(best_idx).Rinflate);
    fprintf('    cfg.sigma_pos_rw                     = %.4e\n', p5_results(best_idx).sigPos);
    fprintf('    cfg.P0_clkBias                       = %.4e   %% [m/s]^2\n', clk_regimes(opt_cli).P0_bias);
    fprintf('    cfg.P0_clkDrift                      = %.4e   %% [m/s²]^2\n', clk_regimes(opt_cli).P0_drift);
    fprintf('    cfg.clk.filt.sigmaBias_rw_mps_sqrt   = %.4e\n', clk_regimes(opt_cli).sigBias);
    fprintf('    cfg.clk.filt.sigmaDrift_rw_mps2_sqrt = %.4e\n', clk_regimes(opt_cli).sigDrift);

    fprintf('  Phase 5 complete.\n\n');
end

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
    cfg.minElevDeg                 = 25;
    cfg.coarseStep                 = 10;
    cfg.minMaxElev                 = 30;
    cfg.maxSats                    = 12;
    cfg.constellationSelectionWindow_s = 60;

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
nombre_archivo = sprintf('Resultados_Fases1a5_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));

save(nombre_archivo);

fprintf('¡Éxito! Todo el análisis se ha guardado en el archivo: %s\n', nombre_archivo);

%% ── 1. LOAD SAVED DATA ────────────────────────────────────────────────
load('Resultados_Fases1a5_20260318_223900.mat'); 


%% ── 2. SAFE EXTRACTION ────────────────────────────────────────────────
N = numel(p5_results);

all_ss     = nan(N, 1);
Rinflate_v = nan(N, 1);
sigPos_v   = nan(N, 1);
clkLabel_v = cell(N, 1);

for i = 1:N
    if ~isempty(p5_results(i).ss_rms2D)
        all_ss(i) = p5_results(i).ss_rms2D;
    end
    if ~isempty(p5_results(i).Rinflate)
        Rinflate_v(i) = p5_results(i).Rinflate;
    end
    if ~isempty(p5_results(i).sigPos)
        sigPos_v(i) = p5_results(i).sigPos;
    end
    if ~isempty(p5_results(i).clkLabel)
        clkLabel_v{i} = p5_results(i).clkLabel;
    else
        clkLabel_v{i} = 'Unknown';
    end
end

valid_idx = isfinite(all_ss);

ss_clean    = all_ss(valid_idx);
Rinflate_c  = Rinflate_v(valid_idx);
sigPos_c    = sigPos_v(valid_idx);

clkLabel_c  = categorical(clkLabel_v(valid_idx), {'Tight', 'Nominal', 'Loose'});

fprintf('Generando Boxplots para %d simulaciones válidas...\n', numel(ss_clean));

%% ── 3. PLOT SUBPANELS WITHOUT STATISTICS TOOLBOX ──────────────────────
fh_sens = figure('Name','P5 Sensitivity Analysis','Color','w','Position',[100 200 1200 400]);

% --- Subplot 1: Rinflate ---
subplot(1,3,1); hold on; grid on; box on;
    u_R = unique(Rinflate_c);
    for i = 1:length(u_R)
        idx = Rinflate_c == u_R(i);
        vals = ss_clean(idx);
        scatter(i + (rand(sum(idx),1)-0.5)*0.2, vals, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
        plot([i-0.3 i+0.3], [median(vals) median(vals)], 'r-', 'LineWidth', 2.5);
end
xticks(1:length(u_R)); xticklabels(string(u_R));
xlabel('R-inflate Multiplier');
ylabel('Steady-State 2D Error (m)');
title('Impact of Measurement Inflation');

% --- Subplot 2: sigma_pos ---
subplot(1,3,2); hold on; grid on; box on;
u_S = unique(sigPos_c);
for i = 1:length(u_S)
    idx = sigPos_c == u_S(i);
    vals = ss_clean(idx);
    scatter(i + (rand(sum(idx),1)-0.5)*0.2, vals, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
    plot([i-0.3 i+0.3], [median(vals) median(vals)], 'r-', 'LineWidth', 2.5);
end
xticks(1:length(u_S)); xticklabels(string(u_S));
xlabel('Pos Random Walk (\sigma_{pos})');
ylabel('Steady-State 2D Error (m)');
title('Impact of Position Process Noise');

% --- Subplot 3: Clock Regime ---
subplot(1,3,3); hold on; grid on; box on;
cats = {'Tight', 'Nominal', 'Loose'};
for i = 1:length(cats)
    idx = clkLabel_c == cats{i};
    vals = ss_clean(idx);
    if ~isempty(vals)
        scatter(i + (rand(sum(idx),1)-0.5)*0.2, vals, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
        plot([i-0.3 i+0.3], [median(vals) median(vals)], 'r-', 'LineWidth', 2.5);
    end
end
xticks(1:length(cats)); xticklabels(cats);
xlabel('Filter Clock Tuning Regime');
ylabel('Steady-State 2D Error (m)');
title('Impact of Clock Assumptions');

% Título general
sgtitle('Phase 5: Steady-State 2D Error Sensitivity (Median & Scatter)', 'FontWeight', 'bold', 'FontSize', 14);
