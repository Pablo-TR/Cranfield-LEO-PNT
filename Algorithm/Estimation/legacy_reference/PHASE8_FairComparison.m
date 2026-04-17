%% ═══════════════════════════════════════════════════════════════════════════
%  PHASE 8 — Fair Sequential Comparison: stickyBest vs stickyDiverse
%             *** PARFOR-OPTIMISED EDITION ***
%
%  PURPOSE:
%    Apples-to-apples comparison of the two handover policies using
%    real_strtg_sequential_EKF_fastrealistic as the solver.
%
%  PARALLEL STRATEGY:
%    All N_MC × nPol × nDR jobs are FLATTENED into one linear job list
%    and dispatched in a SINGLE parfor.  There is NO nested parallelism,
%    no parfor-inside-parfor, no broadcast of large objects mid-loop.
%    A DataQueue streams live progress to the Command Window without
%    stalling workers.
%
%  FAIRNESS GUARANTEES:
%    • Identical constellation, identical noise, identical seeds
%    • Both policies use capture.maxSwitches = N_SATS-1  (hard cap)
%    • Identical blackout duration at every handover
%    • Duration = N_SATS × DWELL_S + (N_SATS-1) × BLACKOUT_S
%
%  OUTPUTS:
%    Fig 1 — Skyplot comparison  (dots, clean style, no star markers)
%    Fig 2 — East Error + 2D convergence with ±1σ shaded bands
%    .mat   — Full results workspace saved automatically
%
%  SOLVER:  real_strtg_sequential_EKF_geom  (must be on MATLAB path)
%  EDIT:    §0 SETTINGS only.
% ═══════════════════════════════════════════════════════════════════════════

clear; close all; clc;

%% ══════════════════════════════════════════════════════════════════════════
%  §0  SETTINGS — edit here only
% ══════════════════════════════════════════════════════════════════════════

% ── Fairness ─────────────────────────────────────────────────────────────
N_SATS      = 5;          % both policies track this many sats (hard cap)
DWELL_S     = 160;         % [s] minimum dwell per satellite
BLACKOUT_S  = 10;         % [s] handover blackout (identical for both)
DURATION_S  = N_SATS * DWELL_S + (N_SATS - 1) * BLACKOUT_S;  % auto
DT          = 0.1;        % [s] time step

% ── Noise ────────────────────────────────────────────────────────────────
SIGMA_D     = 50.0;       % Doppler noise [Hz]   — fixed across sweep
TEST_DR     = [0.1, 1.0, 5.0];  % Doppler-rate noise [Hz/s] — sweep levels

% ── Monte Carlo ──────────────────────────────────────────────────────────
N_MC        = 20;         % seeds per (policy × noise) combination
BASE_SEED   = 1000;       % deterministic seed base

% ── Parallel pool ────────────────────────────────────────────────────────
%    Set to 0 to auto-detect (uses all physical cores).
%    Set to a positive integer to pin the pool size.
POOL_SIZE   = 0;

% ── Output ───────────────────────────────────────────────────────────────
FIG_DIR     = 'analysis_figures_p8';
SAVE_MAT    = true;       % save full workspace to .mat when done

% ── Visual style ─────────────────────────────────────────────────────────
COLOR_LINE  = {[0.85 0.15 0.15]};  % red / blue / green
COLOR_BAND  = {[0.97 0.70 0.70]};  % pastel versions

%% ══════════════════════════════════════════════════════════════════════════
%  §1  INITIALISE
% ══════════════════════════════════════════════════════════════════════════

if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

policies   = {'stickyBest'};
POL_LABELS = {'stickyBest (Same Plane)'};
nPol       = numel(policies);
nDR        = numel(TEST_DR);
K_len      = round(DURATION_S / DT) + 1;

fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║  PHASE 8 — Fair Sequential EKF  [parfor edition]            ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n');
fprintf('  Solver      : real_strtg_sequential_EKF_geom\n');
fprintf('  Policies    : %s\n', strjoin(policies, ' | '));
fprintf('  N sats each : %d  (maxSwitches cap = %d)\n', N_SATS, N_SATS-1);
fprintf('  Dwell / sat : %d s  |  Blackout : %d s  |  Duration : %d s\n',...
        DWELL_S, BLACKOUT_S, DURATION_S);
fprintf('  sigma_D     : %.1f Hz\n', SIGMA_D);
fprintf('  sigma_DR    : [%s] Hz/s\n', num2str(TEST_DR));
fprintf('  N_MC / combo: %d   |  Total jobs : %d\n\n', ...
        N_MC, nPol * nDR * N_MC);

%% ══════════════════════════════════════════════════════════════════════════
%  §2  OPEN PARALLEL POOL (once, before any computation)
% ══════════════════════════════════════════════════════════════════════════

pool = gcp('nocreate');
if isempty(pool)
    if POOL_SIZE > 0
        pool = parpool('local', POOL_SIZE);
    else
        pool = parpool('local');          % auto-detect physical cores
    end
    fprintf('  Parallel pool opened: %d workers\n\n', pool.NumWorkers);
else
    fprintf('  Reusing existing pool: %d workers\n\n', pool.NumWorkers);
end
N_WORKERS = pool.NumWorkers;

%% ══════════════════════════════════════════════════════════════════════════
%  §3  SERIAL REFERENCE RUNS
%      One run per (policy × noise) for:
%        • canonical t_axis  (all MC jobs must share the same grid)
%        • blackout window extraction
%        • sky-plot data (median noise level = index 2)
%      These are fast and must stay serial (state-dependent extraction).
% ══════════════════════════════════════════════════════════════════════════

fprintf('  [Serial] Running %d reference seeds...\n', nPol * nDR);
t_serial = tic;

out_ref        = cell(nPol, 1);          % sky-plot source (di == 2)
bo_starts_cell = cell(nPol, 1);
bo_ends_cell   = cell(nPol, 1);
t_axis         = [];                     % will be filled on first ref run

for pi = 1:nPol
    for di = 1:nDR
        c = p8_build_cfg(DURATION_S, DT, SIGMA_D, TEST_DR(di), ...
                         policies{pi}, DWELL_S, BLACKOUT_S, N_SATS);
        c.seedMeas = BASE_SEED + 600 + pi*10 + di;
        out_r      = real_strtg_sequential_EKF_geom(c);

        % Lock t_axis on very first run
        if isempty(t_axis)
            t_axis = out_r.t(:)';
        end

        % Sky-plot reference: use middle noise level
        if di == 2
            out_ref{pi} = out_r;
        end

        % Blackout windows: extract once per policy (noise-level agnostic)
        if di == 1
            bm = logical(out_r.trk.blackout(:)');
            K  = numel(t_axis);
            bs = []; be = []; in_bo = false;
            for k = 1:K
                if bm(k) && ~in_bo,  bs(end+1) = t_axis(k);   in_bo = true;  end %#ok<AGROW>
                if ~bm(k) && in_bo,  be(end+1) = t_axis(k-1); in_bo = false; end %#ok<AGROW>
            end
            if in_bo, be(end+1) = t_axis(end); end
            bo_starts_cell{pi} = bs;
            bo_ends_cell{pi}   = be;
        end

        n_used = numel(unique(out_r.trk.sid(out_r.trk.sid > 0)));
        fprintf('    Ref [pol=%d di=%d]: N_sats_used=%d\n', pi, di, n_used);
    end
end
K_actual = numel(t_axis);
fprintf('  Reference runs done  (%.1f s)\n\n', toc(t_serial));

%% ══════════════════════════════════════════════════════════════════════════
%  §4  BUILD FLAT JOB LIST  (all configs pre-computed before parfor)
%
%      parfor requirement: workers must not call helper functions that
%      capture workspace variables or use closures.  Pre-building every
%      config struct here avoids any closure/broadcast issue.
%
%      Job index j encodes:  pi (policy), di (noise), n (seed)
%      j = (pi-1)*nDR*N_MC + (di-1)*N_MC + n
% ══════════════════════════════════════════════════════════════════════════

N_JOBS   = nPol * nDR * N_MC;
job_cfgs = cell(N_JOBS, 1);   % pre-built config structs
job_pi   = zeros(N_JOBS, 1);  % which policy
job_di   = zeros(N_JOBS, 1);  % which noise level
job_n    = zeros(N_JOBS, 1);  % which MC seed

j = 0;
for pi = 1:nPol
    for di = 1:nDR
        for n = 1:N_MC
            j = j + 1;
            c             = p8_build_cfg(DURATION_S, DT, SIGMA_D, TEST_DR(di), ...
                                         policies{pi}, DWELL_S, BLACKOUT_S, N_SATS);
            c.seedMeas    = BASE_SEED + 700 + pi*5000 + di*500 + n;
            job_cfgs{j}   = c;
            job_pi(j)     = pi;
            job_di(j)     = di;
            job_n(j)      = n;
        end
    end
end
fprintf('  Job list built: %d configs ready for parfor.\n\n', N_JOBS);

%% ══════════════════════════════════════════════════════════════════════════
%  §5  PARALLEL MONTE CARLO  (single flat parfor — no nesting)
%
%      Live progress via DataQueue (non-blocking, no fprintf inside parfor).
%      Each worker returns a thin result struct (no solver internals stored).
% ══════════════════════════════════════════════════════════════════════════

% ── Live counter ─────────────────────────────────────────────────────────
dq          = parallel.pool.DataQueue;
done_count  = 0;
t_par_start = tic;

afterEach(dq, @(~) p8_tick(done_count, N_JOBS, t_par_start));

% ── Pre-allocate result containers (cell array — parfor-safe) ────────────
raw_errE  = cell(N_JOBS, 1);
raw_err2D = cell(N_JOBS, 1);

fprintf('  [Parallel] Dispatching %d jobs across %d workers...\n', ...
        N_JOBS, N_WORKERS);

parfor j = 1:N_JOBS                           %#ok<PFUIX>
    cfg  = job_cfgs{j};                        % broadcast: one struct per job
    r    = real_strtg_sequential_EKF_geom(cfg); %#ok<PFBNS>

    K_r  = numel(r.t);
    eE   = nan(1, K_actual);
    e2D  = nan(1, K_actual);
    Kc   = min(K_r, K_actual);
    eE(1:Kc)  = r.errENU(1, 1:Kc);
    e2D(1:Kc) = sqrt(r.errENU(1,1:Kc).^2 + r.errENU(2,1:Kc).^2);

    raw_errE{j}  = eE;
    raw_err2D{j} = e2D;

    send(dq, j);                               % non-blocking progress tick
end

fprintf('\n  Parallel sweep done  (%.1f s total)\n\n', toc(t_par_start));

%% ══════════════════════════════════════════════════════════════════════════
%  §6  AGGREGATE RESULTS  (serial — fast, negligible cost)
% ══════════════════════════════════════════════════════════════════════════

res_mu_E   = nan(nPol, nDR, K_actual);
res_std_E  = nan(nPol, nDR, K_actual);
res_mu_2D  = nan(nPol, nDR, K_actual);
res_std_2D = nan(nPol, nDR, K_actual);

for pi = 1:nPol
    for di = 1:nDR
        % Collect all N_MC runs for this (pi, di) combination
        mask  = (job_pi == pi) & (job_di == di);
        idx   = find(mask);             % N_MC indices into flat arrays

        mc_E  = cell2mat(raw_errE(idx));   % [N_MC x K_actual]
        mc_2D = cell2mat(raw_err2D(idx));

        res_mu_E(pi,di,:)   = mean(mc_E,  1, 'omitnan');
        res_std_E(pi,di,:)  = std(mc_E,   0, 1, 'omitnan');
        res_mu_2D(pi,di,:)  = mean(mc_2D, 1, 'omitnan');
        res_std_2D(pi,di,:) = std(mc_2D,  0, 1, 'omitnan');
    end
end
fprintf('  Results aggregated.\n\n');

%% ══════════════════════════════════════════════════════════════════════════
%  §7  SAVE WORKSPACE
% ══════════════════════════════════════════════════════════════════════════

if SAVE_MAT
    mat_name = fullfile(FIG_DIR, ...
        sprintf('P8_Results_%s.mat', datestr(now,'yyyymmdd_HHMMSS')));
    save(mat_name, ...
         'res_mu_E','res_std_E','res_mu_2D','res_std_2D', ...
         't_axis','out_ref','bo_starts_cell','bo_ends_cell', ...
         'policies','POL_LABELS','TEST_DR','SIGMA_D', ...
         'N_SATS','DWELL_S','BLACKOUT_S','DURATION_S','N_MC');
    fprintf('  Workspace saved -> %s\n\n', mat_name);
end

%% ══════════════════════════════════════════════════════════════════════════
%  §8  FIGURE 1 — SKYPLOT COMPARISON (clean, dots, no star markers)
% ══════════════════════════════════════════════════════════════════════════

fh_sky = figure('Name','P8 Skyplot', 'Color','w', ...
                'Position',[80 80 380*nPol 540]);

for pi = 1:nPol
    subplot(1, nPol, pi, polaraxes);
    p8_draw_sky(out_ref{pi}, POL_LABELS{pi}, N_SATS);
end

%p8_save(fh_sky, FIG_DIR, 'P8_Skyplot_Fair');

%% ══════════════════════════════════════════════════════════════════════════
%  §9  FIGURE 2 — EAST ERROR CONVERGENCE WITH +-1sigma BANDS (3 subplots)
% ══════════════════════════════════════════════════════════════════════════

fh_east = figure('Name','P8 East Error +-1sigma', 'Color','w', ...
                 'Position',[100 60 1100 870]);

for di = 1:nDR
    ax = subplot(3, 1, di);
    hold on; grid on; box on;
    set(ax, 'GridAlpha', 0.20, 'GridLineStyle', ':', 'FontSize', 10);

    % Blackout patches (background, drawn first)
    p8_draw_blackouts(ax, bo_starts_cell{1}, bo_ends_cell{1});

    % +-1sigma bands then mean lines
    y_abs_max = 10;    % minimum y-scale
    for pi = 1:nPol
        mu  = squeeze(res_mu_E(pi, di, :))';
        sig = squeeze(res_std_E(pi, di, :))';

        % Shaded +-1sigma band
        fill([t_axis, fliplr(t_axis)], ...
             [mu + sig, fliplr(mu - sig)], ...
             COLOR_BAND{pi}, ...
             'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');

        % Mean line (on top of band)
        plot(t_axis, mu, '-', ...
             'Color', COLOR_LINE{pi}, 'LineWidth', 2.2, ...
             'DisplayName', sprintf('\\mu: %s', policies{pi}));

        % Auto y-scale: skip first 15 s transient
        mask      = t_axis > 15;
        y_abs_max = max(y_abs_max, max(abs(mu(mask)) + sig(mask), [], 'omitnan'));
    end

    yline(0, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    y_lim = min(y_abs_max * 1.30, 1500);   % cap at 1500 m
    ylim([-100, y_lim]);
    xlim([0, 500]);
    ylabel('East Error (m)', 'FontSize', 11, 'FontWeight', 'bold');

    title(sprintf('\\sigma_{DR} = %.1f Hz/s', TEST_DR(di)), ...
          'FontSize', 11, 'FontWeight', 'bold');

    % Legend with band proxies — top subplot only, built dynamically for nPol
    if di == 1
        h_lines = findobj(ax, 'Type', 'line', '-not', 'LineStyle', ':');
        h_lines = flipud(h_lines(1:min(nPol, numel(h_lines))));
        hp = gobjects(nPol, 1);
        lbl_lines = cell(nPol, 1);
        lbl_bands = cell(nPol, 1);
        for pi2 = 1:nPol
            hp(pi2)       = patch(nan, nan, COLOR_BAND{pi2}, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
            lbl_lines{pi2} = sprintf('\\mu: %s', policies{pi2});
            lbl_bands{pi2} = sprintf('\\pm1\\sigma %s', policies{pi2});
        end
        legend([h_lines; hp], [lbl_lines; lbl_bands], ...
               'Location', 'northeast', 'FontSize', 9, 'Box', 'on');
    end

    if di == nDR
        xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    end
end

%p8_save(fh_east, FIG_DIR, 'P8_EastError_Fair_Sigma');

%% ══════════════════════════════════════════════════════════════════════════
%  §10  FIGURE 3 — 2D HORIZONTAL ERROR WITH +-1sigma BANDS
% ══════════════════════════════════════════════════════════════════════════

fh_2d = figure('Name','P8 2D Horiz Error +-1sigma', 'Color','w', ...
               'Position',[120 60 1100 870]);

for di = 1:nDR
    ax = subplot(3, 1, di);
    hold on; grid on; box on;
    set(ax, 'GridAlpha', 0.20, 'GridLineStyle', ':', 'FontSize', 10);

    p8_draw_blackouts(ax, bo_starts_cell{1}, bo_ends_cell{1});

    y_max = 10;
    for pi = 1:nPol
        mu  = squeeze(res_mu_2D(pi, di, :))';
        sig = squeeze(res_std_2D(pi, di, :))';

        fill([t_axis, fliplr(t_axis)], ...
             [mu + sig, fliplr(max(0, mu - sig))], ...
             COLOR_BAND{pi}, ...
             'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');

        plot(t_axis, mu, '-', ...
             'Color', COLOR_LINE{pi}, 'LineWidth', 2.2, ...
             'DisplayName', sprintf('\\mu: %s', policies{pi}));

        mask  = t_axis > 15;
        %y_max = max(y_max, max(mu(mask) + sig(mask), [], 'omitnan'));
    end

    yline(0, 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');

    ylim([0, 600]);
    xlim([0, 500]);
    ylabel('2D Horiz Error (m)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('\\sigma_{DR} = %.1f Hz/s', TEST_DR(di)), ...
          'FontSize', 11, 'FontWeight', 'bold');

    if di == 1
        h_lines = findobj(ax, 'Type', 'line', '-not', 'LineStyle', ':');
        h_lines = flipud(h_lines(1:min(nPol, numel(h_lines))));
        hp = gobjects(nPol, 1);
        lbl_lines = cell(nPol, 1);
        lbl_bands = cell(nPol, 1);
        for pi2 = 1:nPol
            hp(pi2)        = patch(nan, nan, COLOR_BAND{pi2}, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
            lbl_lines{pi2} = sprintf('\\mu: %s', policies{pi2});
            lbl_bands{pi2} = sprintf('\\pm1\\sigma %s', policies{pi2});
        end
        legend([h_lines; hp], [lbl_lines; lbl_bands], ...
               'Location', 'northeast', 'FontSize', 9, 'Box', 'on');
    end
    if di == nDR
        xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    end
end


p8_save(fh_2d, FIG_DIR, 'P8_2DError_Fair_Sigma');

fprintf('\n  Phase 8 complete.  All figures saved to: %s\n\n', FIG_DIR);

%% ══════════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS  (must be at end of script / in same file)
% ══════════════════════════════════════════════════════════════════════════

% -------------------------------------------------------------------------
%  p8_build_cfg
%  Assembles a fully self-contained config struct.
%  NO closures, NO globals — safe to broadcast into parfor workers.
% -------------------------------------------------------------------------
function cfg = p8_build_cfg(duration_s, dt, sigma_D, sigma_DR, ...
                             policy, dwell_s, blackout_s, n_sats)
    cfg = build_base_cfg();                  % your existing base builder

    cfg.duration_sec                      = duration_s;
    cfg.dt                                = dt;
    cfg.sigma_doppler_Hz                  = sigma_D;
    cfg.sigma_dopplerRate_Hzps            = sigma_DR;
    cfg.meas.useDoppler                   = true;
    cfg.meas.useDopplerRate               = true;
    cfg.tropo.enable                      = false;

    cfg.capture.mode                      = policy;
    cfg.capture.switchBlackout_s          = blackout_s;
    cfg.capture.minDwell_s                = dwell_s;
    cfg.capture.maxSwitches               = n_sats - 1;  % FAIRNESS hard cap

    cfg.clk.truth.initBias_mps            = 0.10;
    cfg.clk.truth.initDrift_mps2          = 5e-5;
    cfg.clk.truth.sigmaBias_rw_mps_sqrt   = 5e-3;
    cfg.clk.truth.sigmaDrift_rw_mps2_sqrt = 1e-4;
    cfg.P0_clkBias                        = (0.15)^2;
    cfg.P0_clkDrift                       = (2.5e-7 * 2)^2;
    cfg.clk.filt.sigmaBias_rw_mps_sqrt    = 1.0e-2;
    cfg.clk.filt.sigmaDrift_rw_mps2_sqrt  = 5.0e-4;
    cfg.plot.enable                       = false;
end

% -------------------------------------------------------------------------
%  p8_tick
%  DataQueue callback: prints a live progress bar (runs in serial context).
%  Uses a persistent counter — safe because afterEach callbacks are serial.
% -------------------------------------------------------------------------
function p8_tick(~, n_total, t_start)
    persistent cnt
    if isempty(cnt), cnt = 0; end
    cnt  = cnt + 1;
    pct  = 100 * cnt / n_total;
    filled  = floor(pct / 5);
    bar_str = [repmat(char(9608), 1, filled), ...    % filled blocks
               repmat(char(9617), 1, 20 - filled)];  % empty blocks
    el   = toc(t_start);
    eta  = el / cnt * (n_total - cnt);
    fprintf('\r  [%s] %3.0f%%  (%d/%d)  elapsed %.0fs  ETA ~%.0fs   ', ...
            bar_str, pct, cnt, n_total, el, max(0, eta));
    if cnt == n_total
        fprintf('\n');
        cnt = [];   % reset persistent for potential re-runs
    end
end

% -------------------------------------------------------------------------
%  p8_draw_sky
%  Clean polar sky plot. Filled circles for acquisition, dots for tracking.
%  No star/pentagon markers.
% -------------------------------------------------------------------------
function p8_draw_sky(out_data, title_str, n_sats_expected)

    pax = gca;
    hold(pax, 'on');

    % Axis cosmetics
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir          = 'clockwise';
    pax.RLim              = [0, 90];
    pax.RTick             = [0, 30, 60, 90];
    pax.RTickLabel        = {'90°','60°','30°','0°'};
    pax.GridAlpha         = 0.22;
    pax.GridLineStyle     = ':';
    pax.FontSize          = 9;
    pax.Color             = [0.97 0.97 0.99];

    % Cardinal labels (outside ring)
    r_lbl = 105;
    text(pax, deg2rad(0),   r_lbl, 'N', 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.45 0.45 0.45]);
    text(pax, deg2rad(90),  r_lbl, 'E', 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.45 0.45 0.45]);
    text(pax, deg2rad(180), r_lbl, 'S', 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.45 0.45 0.45]);
    text(pax, deg2rad(270), r_lbl, 'W', 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.45 0.45 0.45]);

    % Satellite geometry
    trk_sid   = out_data.trk.sid(:)';
    tracked   = unique(trk_sid(trk_sid > 0));
    n_tracked = numel(tracked);
    pal       = lines(max(n_tracked, 6));

    lh_legend = gobjects(n_tracked, 1);

    for i = 1:n_tracked
        s   = tracked(i);
        col = pal(i, :);

        az_rad = deg2rad(out_data.sat.az(:, s));
        el_deg = out_data.sat.elev(:, s);
        r_plot = 90 - el_deg;

        % Full pass arc: thin light grey
        vis = el_deg > 2;
        if any(vis)
            polarplot(pax, az_rad(vis), r_plot(vis), '-', ...
                'Color', [0.80 0.80 0.80], 'LineWidth', 0.9, ...
                'HandleVisibility', 'off');
        end

        % Tracked segment: coloured dots, subsampled for clarity
        act = find(trk_sid == s);
        if ~isempty(act)
            step = max(1, floor(numel(act) / 40));   % ~40 dots max
            sub  = act(1:step:end);
            polarplot(pax, az_rad(sub), r_plot(sub), '.', ...
                'Color', col, 'MarkerSize', 8, ...
                'HandleVisibility', 'off');

            % Acquisition marker: bold filled circle (NOT star)
            k0 = act(1);
            polarplot(pax, az_rad(k0), r_plot(k0), 'o', ...
                'MarkerSize', 10, ...
                'MarkerFaceColor', col, ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.4, ...
                'HandleVisibility', 'off');

            % Label: satellite number, offset radially outward
            r_txt = min(r_plot(k0) + 7, 87);
            text(pax, az_rad(k0), r_txt, ...
                sprintf(' Sat %d', s), ...
                'FontSize', 9, 'FontWeight', 'bold', ...
                'Color', col .* 0.72, ...
                'HorizontalAlignment', 'left');
        end

        % Legend proxy
        lh_legend(i) = polarplot(pax, nan, nan, 'o', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
            'MarkerSize', 7, 'LineWidth', 1.0, ...
            'DisplayName', sprintf('Sat %d', s));
    end

    legend(pax, lh_legend, ...
           'Location', 'southoutside', ...
           'NumColumns', n_tracked, ...
           'FontSize', 8, 'Box', 'off');

    title(pax,    title_str, 'FontSize', 11, 'FontWeight', 'bold');
    subtitle(pax, sprintf('N_{tracked} = %d  (expected %d)', ...
                          n_tracked, n_sats_expected), ...
             'FontSize', 8, 'Color', [0.45 0.45 0.45]);
end

% -------------------------------------------------------------------------
%  p8_draw_blackouts
%  Draw pale orange patch for every handover blackout window.
% -------------------------------------------------------------------------
function p8_draw_blackouts(ax, b_start, b_end)
    n_bo = min(numel(b_start), numel(b_end));
    for bi = 1:n_bo
        patch(ax, ...
              [b_start(bi), b_end(bi), b_end(bi), b_start(bi)], ...
              [-1e5, -1e5, 1e5, 1e5], ...
              [1.00 0.88 0.72], ...
              'EdgeColor', 'none', 'FaceAlpha', 0.42, ...
              'HandleVisibility', 'off');
    end
end

% -------------------------------------------------------------------------
%  p8_save
%  Export figure to high-resolution PNG.
% -------------------------------------------------------------------------
function p8_save(fh, fig_dir, fname)
    out_path = fullfile(fig_dir, [fname, '.png']);
    exportgraphics(fh, out_path, 'Resolution', 200);
    fprintf('  Saved -> %s\n', out_path);
end

function cfg = build_base_cfg()
    % ── TLE File ──────────────────────────────────────────────────────────
    cfg.tleFile = fullfile('InputData', 'oneweb.tle');

    % ── Simulation time (fixed by mandate) ───────────────────────────────
    cfg.duration_sec  = 1000;
    cfg.dt            = 0.1;
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
    cfg.sigma_doppler_Hz       = 50.0;   % Doppler white noise [Hz]
    cfg.sigma_dopplerRate_Hzps = 1.0;   % Doppler-rate white noise [Hz/s]

    cfg.meas.enable_lighttime  = false;
    cfg.meas.colored_tau_s     = 0.5;
    cfg.meas.sigma_col0_mps    = 0.03;
    cfg.meas.Rinflate          = 2.0;
    cfg.meas.useDoppler        = true;
    cfg.meas.useDopplerRate    = true;   % Doppler-only by default
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
