%% ═══════════════════════════════════════════════════════════════════════
%  MAIN_PHASE0_VALIDATION
%  Phase 0 — Full Filter Validation for LEO Doppler Sequential IEKF
%
%  Produces all validation evidence for the filter chapter of the GDP report.
%  Six independent blocks (B1–B6) selectable via RUN_BLOCK flags.
%
%  Pre-computed targets (adopted baseline run):
%    meanNIS = 0.982,  ACF_D = 2%,  ACF_DR = 4%,
%    kurtosis_D = 5.91,  kurtosis_DR = 0.49
%
%  Dependencies:
%    real_strtg_sequential_EKF_fastrealistic.m  (solver – do NOT modify)
%    ekf_extract_metrics.m
%    InputData/oneweb.tle
% ═══════════════════════════════════════════════════════════════════════

clear; close all; clc;

%% ── Run-block flags ──────────────────────────────────────────────────
RUN_BLOCK = struct( ...
    'B1_progressive', true,  ...   % progressive complexity cascade
    'B2_ekf_vs_iekf', true,  ...   % EKF vs IEKF comparison
    'B3_Rinflate',    true,  ...   % R_inflate sweep
    'B4_P0',          true,  ...   % initial covariance sweep
    'B5_Q',           true,  ...   % process noise sweep
    'B6_NIS_ACF',     true   ...   % final NIS / ACF / kurtosis diagnostics
);

N_MC      = 10;    % Monte Carlo seeds per configuration (use 20 for final run)
BASE_SEED = 42;
SAVE_FIGS = true;
FIG_DIR   = 'phase0_figures';

% Seed block offsets for reproducibility
SEED_OFF = struct('B1',0,'B2',100,'B3',200,'B4',300,'B5',400,'B6',500);

%% ── Output directory ─────────────────────────────────────────────────
if SAVE_FIGS && ~exist(FIG_DIR,'dir')
    mkdir(FIG_DIR);
end

%% ── Initialise results structs ───────────────────────────────────────
B1_results = struct(); B2_results = struct(); B3_results = struct();
B4_results = struct(); B5_results = struct(); B6_results = struct();
baseline_cfg = build_baseline_cfg();

fprintf('╔═══════════════════════════════════════════════════════╗\n');
fprintf('║  MAIN_PHASE0_VALIDATION — LEO Doppler IEKF            ║\n');
fprintf('║  N_MC = %d seeds per config                           ║\n', N_MC);
fprintf('╚═══════════════════════════════════════════════════════╝\n\n');

%% ════════════════════════════════════════════════════════════════════
%%  BLOCK B1 — Progressive Complexity Cascade
%%  PURPOSE: Validate filter correctness step by step by adding one
%%  physical complexity layer at a time to the truth model.
%%  ALL steps use Doppler+Rate. What changes is what truth injects.
%%  Validation criterion: filter consistency (NIS, ACF, innovations),
%%  NOT position accuracy (RMSE is always large for single-satellite).
%% ════════════════════════════════════════════════════════════════════
if RUN_BLOCK.B1_progressive
    fprintf('══════ BLOCK B1: Progressive Complexity Cascade ══════\n');

    %  Step definitions — what changes in TRUTH MODEL per step
    %  Filter always uses: Dopp+Rate, Rinflate=2, startup inflation, etc.
    %
    %  clkZero  : truth clock bias=0, filter init=0, clock RW=0 (known clock)
    %  whiteNoise: sigma_D level to inject
    %  colored  : whether AR(1) coloured noise is active
    %  tropo    : whether tropospheric delay is active
    %  label    : description for plot title and console
    steps = struct( ...
        'name',       {'S1','S2','S3','S4','S5','S6'}, ...
        'label',      {'Geometry + Jacobian check (near-zero noise, known clock)', ...
                       'Small noise, known clock (NIS calibration)', ...
                       'Realistic noise + clock mismatch (cold start)', ...
                       'Add troposphere (matched ZTD)', ...
                       'Add coloured noise (kurtosis check)', ...
                       'Full baseline (target: NIS~0.98, ACF<5%)'}, ...
        'clkZero',    {true,  true,  false, false, false, false}, ...
        'sigma_D_Hz', {0.1,   5.0,   50.0,  50.0,  50.0,  50.0}, ...
        'sigma_DR',   {0.001, 0.05,  1.0,   1.0,   1.0,   1.0 }, ...
        'colored',    {false, false, false, false, true,  true }, ...
        'tropo',      {false, false, false, true,  true,  true }  ...
    );

    nSteps = numel(steps);
    B1_results.metrics    = cell(nSteps,1);
    B1_results.out        = cell(nSteps,1);
    B1_results.steps      = steps;
    B1_results.max_innov  = nan(nSteps,1);   % for S1 Jacobian check

    for si = 1:nSteps
        st    = steps(si);
        cfg_s = build_baseline_cfg();
        cfg_s.duration_sec        = 50;    % 50 s sufficient
        cfg_s.plot.enable         = false;
        cfg_s.verbose             = false;
        cfg_s.seedMeas            = BASE_SEED + SEED_OFF.B1 + si;

        % All steps use Doppler + Doppler-rate
        cfg_s.meas.useDoppler     = true;
        cfg_s.meas.useDopplerRate = true;

        % Measurement noise level (truth injection)
        cfg_s.sigma_doppler_Hz       = st.sigma_D_Hz;
        cfg_s.sigma_dopplerRate_Hzps = st.sigma_DR;

        % Troposphere
        cfg_s.tropo.enable = st.tropo;

        % Coloured noise
        if ~st.colored
            cfg_s.meas.sigma_col0_mps = 0;
        end

        % Clock: zero mismatch (S1, S2) or realistic cold-start (S3-S6)
        if st.clkZero
            % Truth clock: zero bias, zero drift, zero random walk
            cfg_s.clk.truth.initBias_mps            = 0.0;
            cfg_s.clk.truth.initDrift_mps2          = 0.0;
            cfg_s.clk.truth.sigmaBias_rw_mps_sqrt   = 0.0;
            cfg_s.clk.truth.sigmaDrift_rw_mps2_sqrt = 0.0;
            % Filter clock: initialise at truth (zero) with tight prior
            cfg_s.clk.filter_init_bias_mps          = 0.0;
            cfg_s.clk.filter_init_drift_mps2        = 0.0;
            cfg_s.P0_clkBias  = (0.001)^2;   % tight — we know the clock
            cfg_s.P0_clkDrift = (1e-8)^2;
        end
        % S3-S6: keep baseline (truth bias=0.10 m/s, filter init=0) —
        % this is the realistic cold-start with unknown clock bias.

        fprintf('  Running %s: %s ...\n', st.name, st.label);
        try
            out_s = real_strtg_sequential_EKF_fastrealistic(cfg_s);
            m_s   = ekf_extract_metrics(out_s, cfg_s);
            B1_results.metrics{si} = m_s;
            B1_results.out{si}     = out_s;

            % For S1: compute max innovation to validate Jacobian
            if si == 1
                nu_raw = out_s.nu_hist(:);
                nu_valid = nu_raw(~isnan(nu_raw) & nu_raw ~= 0);
                if ~isempty(nu_valid)
                    B1_results.max_innov(1) = max(abs(nu_valid));
                    thresh = 3 * (st.sigma_D_Hz * cfg_s.c / cfg_s.fc_Hz);
                    if B1_results.max_innov(1) < thresh
                        fprintf('  S1 Jacobian check PASSED: max|nu| = %.4f m/s < 3*sigma = %.4f m/s\n', ...
                            B1_results.max_innov(1), thresh);
                    else
                        fprintf('  *** S1 Jacobian check FAILED: max|nu| = %.4f m/s >= 3*sigma = %.4f m/s\n', ...
                            B1_results.max_innov(1), thresh);
                        fprintf('  *** This suggests a bug in the geometric model or Jacobians.\n');
                    end
                else
                    fprintf('  *** S1: No valid innovations — all measurements rejected by NIS gate.\n');
                    fprintf('  *** Check that sigma_D is not too small (>0) and clock is correctly zeroed.\n');
                end
            end

        catch ME
            fprintf('  ERROR in %s: %s\n', st.name, ME.message);
            B1_results.metrics{si} = struct('RMSE_h_m',NaN,'meanNIS',NaN, ...
                'fracSigACF_dopp',NaN,'excessKurtosis_dopp',NaN,'convTime_horiz_s',NaN);
            B1_results.out{si} = [];
        end
    end

    % ── Console table ─────────────────────────────────────────────────
    fprintf('\n  Step | sigma_D  | clkZero | tropo | colored | RMSE_h   | NIS    | ACF%%  | kurtosis\n');
    fprintf('  -----|----------|---------|-------|---------|----------|--------|-------|----------\n');
    for si = 1:nSteps
        m_s  = B1_results.metrics{si};
        st   = steps(si);
        clkStr   = tern_str(st.clkZero,'yes','no');
        tropoStr = tern_str(st.tropo,'yes','no');
        colStr   = tern_str(st.colored,'yes','no');
        fprintf('  %-4s | %8.2f | %-7s | %-5s | %-7s | %7.1f m | %6.3f | %4.0f%% | %6.2f\n', ...
            st.name, st.sigma_D_Hz, clkStr, tropoStr, colStr, ...
            m_s.RMSE_h_m, m_s.meanNIS, ...
            100*m_s.fracSigACF_dopp, m_s.excessKurtosis_dopp);
    end
    fprintf('\n  NOTE: RMSE is large for all steps — this is correct physics for single-satellite.\n');
    fprintf('  The validation criterion is NIS~1.0 and ACF<5%%, not RMSE.\n');
    fprintf('  S1 NIS is not meaningful (near-zero noise); check max|innovation| instead.\n\n');

    % ── Plot B1: 6-panel progressive cascade ──────────────────────────
    fig1 = figure('Name','B1_progressive_cascade','Position',[50 50 1400 900]);
    y_max = 600;
    for si = 1:nSteps
        subplot(2,3,si);
        out_s = B1_results.out{si};
        m_s   = B1_results.metrics{si};
        st    = steps(si);

        if ~isempty(out_s)
            t_s  = out_s.t;
            eE_s = out_s.errENU(1,:);
            eN_s = out_s.errENU(2,:);
            eH_s = sqrt(eE_s.^2 + eN_s.^2);
            plot(t_s, eH_s, 'b-', 'LineWidth', 1.2); hold on;
            plot(t_s, abs(eE_s), 'r--', 'LineWidth', 0.8);
            plot(t_s, abs(eN_s), 'g--', 'LineWidth', 0.8);
            if si == 1
                legend('2D Horiz','|East|','|North|','Location','northeast','FontSize',8);
            end
        end
        yline(50, 'r-', 'LineWidth', 1.5);
        ylim([0, y_max]);
        xlim([0, 50]);
        xlabel('Time [s]', 'FontSize', 10);
        ylabel('Error [m]', 'FontSize', 10);
        % Shorten title to step name + key parameter
        title(sprintf('%s: sigma_D=%.1f Hz, clkZero=%s', ...
            st.name, st.sigma_D_Hz, tern_str(st.clkZero,'T','F')), ...
            'FontSize', 10, 'Interpreter', 'none');

        % Annotation: show NIS and ACF, not RMSE (which is large by design)
        if si == 1 && isfinite(B1_results.max_innov(1))
            txt = sprintf('max|nu| = %.4f m/s\nNIS = n/a (near-zero noise)\nValidation: Jacobian check', ...
                B1_results.max_innov(1));
        else
            txt = sprintf('NIS = %.3f\nACF = %.0f%%\nKurt = %.2f', ...
                m_s.meanNIS, 100*m_s.fracSigACF_dopp, m_s.excessKurtosis_dopp);
        end
        text(0.98, 0.97, txt, 'FontSize', 8, 'Units','normalized', ...
            'VerticalAlignment','top', 'HorizontalAlignment','right', ...
            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k');
        grid on;
        set(gca, 'FontSize', 10);
    end
    sgtitle('B1 — Progressive Complexity Cascade (single run, 50 s)', ...
        'FontSize', 12, 'FontWeight', 'bold');
    if SAVE_FIGS
        print(fig1, fullfile(FIG_DIR, 'B1_progressive_cascade'), '-dpng', '-r150');
        fprintf('  Saved: B1_progressive_cascade.png\n');
    end
end

%% ════════════════════════════════════════════════════════════════════
%%  BLOCK B2 — EKF vs IEKF Comparison
%% ════════════════════════════════════════════════════════════════════
if RUN_BLOCK.B2_ekf_vs_iekf
    fprintf('══════ BLOCK B2: EKF vs IEKF Comparison (N_MC=%d) ══════\n', N_MC);

    configs = struct( ...
        'label', {'EKF','IEKF'}, ...
        'maxIter', {1, 10} ...  % EKF = single linearisation per epoch
    );
    nCfg = numel(configs);

    B2_results.metrics = cell(nCfg, N_MC);
    B2_results.eH_cell = cell(nCfg, N_MC);
    B2_results.nis_cell = cell(nCfg, N_MC);
    B2_results.t_ref   = [];

    for ci = 1:nCfg
        fprintf('  Config: %s (maxIEKFIter=%d)\n', configs(ci).label, configs(ci).maxIter);

        % Hoist loop-invariant scalars so parfor workers can broadcast them
        maxIter_ci   = configs(ci).maxIter;
        seed_base_ci = BASE_SEED + SEED_OFF.B2 + (ci-1)*N_MC;
        K_est_ci     = round(100/0.01) + 1;

        % Sliced output temporaries (parfor writes mc-indexed cells locally)
        metrics_ci = cell(N_MC,1);
        eH_ci      = cell(N_MC,1);
        nis_ci     = cell(N_MC,1);
        tref_ci    = cell(N_MC,1);  % collect t from each run; t_ref set after

        parfor mc = 1:N_MC
            cfg_r = build_baseline_cfg();
            cfg_r.filter.maxIEKFIter = maxIter_ci;  %#ok<PFBNS>
            cfg_r.plot.enable = false;
            cfg_r.verbose     = false;
            cfg_r.seedMeas    = seed_base_ci + mc;   %#ok<PFBNS>
            try
                out_r = real_strtg_sequential_EKF_fastrealistic(cfg_r);
                m_r   = ekf_extract_metrics(out_r, cfg_r);
                metrics_ci{mc} = m_r;
                eH_ci{mc}      = sqrt(out_r.errENU(1,:).^2 + out_r.errENU(2,:).^2);
                nis_ci{mc}     = out_r.nis(:)';
                tref_ci{mc}    = out_r.t(:)';
            catch ME
                fprintf('    [B2 ci=%d mc=%d] failed: %s\n', ci, mc, ME.message);
                metrics_ci{mc} = struct('RMSE_h_m',NaN,'meanNIS',NaN,'convTime_horiz_s',Inf);
                eH_ci{mc}      = nan(1, K_est_ci);  %#ok<PFBNS>
                nis_ci{mc}     = nan(1, K_est_ci);
                tref_ci{mc}    = [];
            end
        end

        % Copy back into B2_results
        B2_results.metrics(ci,:)  = metrics_ci';
        B2_results.eH_cell(ci,:)  = eH_ci';
        B2_results.nis_cell(ci,:) = nis_ci';

        % Pick t_ref from first successful run of this config
        if isempty(B2_results.t_ref)
            for mc = 1:N_MC
                if ~isempty(tref_ci{mc})
                    B2_results.t_ref = tref_ci{mc};
                    break;
                end
            end
        end
        fprintf('    Done %d/%d seeds.\n', N_MC, N_MC);
    end

    % Aggregate
    t_ref = B2_results.t_ref;
    K_ref = numel(t_ref);
    cfg_colours = {'r','b'};

    % Pad/trim all series to K_ref
    eH_mat  = cell(nCfg,1);
    nis_mat = cell(nCfg,1);
    for ci = 1:nCfg
        eH_mat{ci}  = nan(N_MC, K_ref);
        nis_mat{ci} = nan(N_MC, K_ref);
        for mc = 1:N_MC
            v = B2_results.eH_cell{ci,mc};
            n = min(numel(v), K_ref);
            eH_mat{ci}(mc,1:n)  = v(1:n);
            v2 = B2_results.nis_cell{ci,mc};
            n2 = min(numel(v2), K_ref);
            nis_mat{ci}(mc,1:n2) = v2(1:n2);
        end
    end

    % Console output
    fprintf('\n  Config | RMSE_h (mean) | RMSE_h (std) | meanNIS | frac_converged\n');
    fprintf('  -------|---------------|--------------|---------|---------------\n');
    for ci = 1:nCfg
        rmse_vals = cellfun(@(m) m.RMSE_h_m, B2_results.metrics(ci,:));
        nis_vals  = cellfun(@(m) m.meanNIS,   B2_results.metrics(ci,:));
        conv_vals = cellfun(@(m) m.convTime_horiz_s, B2_results.metrics(ci,:));
        frac_conv = sum(isfinite(conv_vals) & conv_vals < Inf);
        fprintf('  %-6s | %13.1f | %12.1f | %7.3f | %d/%d\n', ...
            configs(ci).label, mean(rmse_vals, 'omitnan'), std(rmse_vals, 'omitnan'), ...
            mean(nis_vals, 'omitnan'), frac_conv, N_MC);
    end
    fprintf('\n');

    % ── Plot B2a: 2D horizontal error ────────────────────────────────
    fig2a = figure('Name','B2a_ekf_vs_iekf_error','Position',[50 50 900 500]);
    leg_handles = [];
    leg_labels  = {};
    for ci = 1:nCfg
        mu  = mean(eH_mat{ci}, 1, 'omitnan');
        sg  = std(eH_mat{ci}, 0, 1, 'omitnan');
        col = cfg_colours{ci};
        fh = fill([t_ref, fliplr(t_ref)], [mu+sg, fliplr(max(mu-sg,0))], col, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on;
        ph = plot(t_ref, mu, col, 'LineWidth', 2.0);
        leg_handles(end+1) = ph; %#ok<AGROW>
        leg_labels{end+1}  = sprintf('%s (mean ±1σ)', configs(ci).label); %#ok<AGROW>
    end
    yline(50, 'r--', 'LineWidth', 1.5, 'Label', '50 m');
    xlabel('Time [s]', 'FontSize', 12);
    ylabel('2D Horizontal Error [m]', 'FontSize', 12);
    title(sprintf('EKF vs IEKF: 2D Horizontal Error (N_{MC}=%d)', N_MC), ...
        'FontSize', 13, 'FontWeight', 'bold');
    legend(leg_handles, leg_labels, 'Location', 'northeast', 'FontSize', 11);
    grid on; set(gca,'FontSize',11);
    if SAVE_FIGS
        print(fig2a, fullfile(FIG_DIR,'B2a_ekf_vs_iekf_error'), '-dpng', '-r150');
        fprintf('  Saved: B2a_ekf_vs_iekf_error.png\n');
    end

    % ── Plot B2b: NIS evolution ───────────────────────────────────────
    fig2b = figure('Name','B2b_ekf_vs_iekf_nis','Position',[50 50 900 400]);
    for ci = 1:nCfg
        col = cfg_colours{ci};
        mu_nis = mean(nis_mat{ci}, 1, 'omitnan');
        plot(t_ref, mu_nis, col, 'LineWidth', 2.0, 'DisplayName', configs(ci).label); hold on;
    end
    yline(1.0, 'k--', 'LineWidth', 1.5, 'Label', 'NIS=1');
    xlabel('Time [s]', 'FontSize', 12);
    ylabel('Mean NIS', 'FontSize', 12);
    title('EKF vs IEKF: NIS Evolution', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 11);
    grid on; set(gca,'FontSize',11);
    if SAVE_FIGS
        print(fig2b, fullfile(FIG_DIR,'B2b_ekf_vs_iekf_nis'), '-dpng', '-r150');
        fprintf('  Saved: B2b_ekf_vs_iekf_nis.png\n');
    end
end

%% ════════════════════════════════════════════════════════════════════
%%  BLOCK B3 — R_inflate Sweep
%% ════════════════════════════════════════════════════════════════════
if RUN_BLOCK.B3_Rinflate
    Rinflate_vals = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 5.0];
    nR = numel(Rinflate_vals);
    adopted_Rinflate = 2.0;

    fprintf('══════ BLOCK B3: R_inflate Sweep (N_MC=%d per value) ══════\n', N_MC);

    B3_results.Rinflate_vals = Rinflate_vals;
    B3_results.RMSE_h_mean   = nan(1,nR);
    B3_results.RMSE_h_std    = nan(1,nR);
    B3_results.NIS_mean      = nan(1,nR);
    B3_results.NIS_std       = nan(1,nR);
    B3_results.frac_consistent = nan(1,nR);  % fraction with meanNIS in [0.5, 2.0]

    for vi = 1:nR
        Ri           = Rinflate_vals(vi);
        seed_base_vi = BASE_SEED + SEED_OFF.B3 + (vi-1)*N_MC;  % broadcast scalar
        rmse_v = nan(N_MC,1);
        nis_v  = nan(N_MC,1);
        parfor mc = 1:N_MC
            cfg_r = build_baseline_cfg();
            cfg_r.meas.Rinflate = Ri;           %#ok<PFBNS>
            cfg_r.plot.enable   = false;
            cfg_r.verbose       = false;
            cfg_r.seedMeas      = seed_base_vi + mc;  %#ok<PFBNS>
            try
                out_r = real_strtg_sequential_EKF_fastrealistic(cfg_r);
                m_r   = ekf_extract_metrics(out_r, cfg_r);
                rmse_v(mc) = m_r.RMSE_h_m;
                nis_v(mc)  = m_r.meanNIS;
            catch
            end
        end
        B3_results.RMSE_h_mean(vi)     = mean(rmse_v, 'omitnan');
        B3_results.RMSE_h_std(vi)      = std(rmse_v, 'omitnan');
        B3_results.NIS_mean(vi)        = mean(nis_v, 'omitnan');
        B3_results.NIS_std(vi)         = std(nis_v, 'omitnan');
        B3_results.frac_consistent(vi) = mean(nis_v >= 0.5 & nis_v <= 2.0, 'omitnan');
        fprintf('  Rinflate=%5.2f: RMSE_h=%.1f±%.1f m, NIS=%.3f±%.3f, cons=%.0f%%\n', ...
            Ri, B3_results.RMSE_h_mean(vi), B3_results.RMSE_h_std(vi), ...
            B3_results.NIS_mean(vi), B3_results.NIS_std(vi), ...
            100*B3_results.frac_consistent(vi));
    end

    % Console table
    fprintf('\n  Rinflate | RMSE_h [m] |  meanNIS | cons.\n');
    fprintf('  ---------|------------|----------|------\n');
    for vi = 1:nR
        marker = '';
        if Rinflate_vals(vi) == adopted_Rinflate
            marker = [marker, ' *'];
        end
        if B3_results.NIS_mean(vi) < 0.5 || B3_results.NIS_mean(vi) > 2.0
            marker = [marker, ' !'];
        end
        fprintf('  %8.2f | %10.1f | %8.3f | %4.0f%%%s\n', ...
            Rinflate_vals(vi), B3_results.RMSE_h_mean(vi), ...
            B3_results.NIS_mean(vi), 100*B3_results.frac_consistent(vi), marker);
    end
    fprintf('  (* adopted value,  ! NIS outside [0.5,2.0])\n\n');

    % ── Plot B3 ──────────────────────────────────────────────────────
    adopted_idx = find(Rinflate_vals == adopted_Rinflate, 1);
    fig3 = figure('Name','B3_Rinflate_sweep','Position',[50 50 1100 500]);

    subplot(1,2,1);
    errorbar(Rinflate_vals, B3_results.RMSE_h_mean, B3_results.RMSE_h_std, ...
        'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 7); hold on;
    yline(50, 'r--', 'LineWidth', 1.5, 'Label', '50 m');
    if ~isempty(adopted_idx)
        xline(adopted_Rinflate, 'k-', 'LineWidth', 2, 'Label', 'Adopted');
    end
    xlabel('R_{inflate}', 'FontSize', 12);
    ylabel('RMSE_h [m]', 'FontSize', 12);
    title('RMSE vs R_{inflate}', 'FontSize', 12);
    grid on; set(gca,'FontSize',11);

    subplot(1,2,2);
    plot(Rinflate_vals, B3_results.NIS_mean, 'b-o', 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'b', 'MarkerSize', 7); hold on;
    yline(0.5, 'g--', 'LineWidth', 1.5, 'Label', '0.5');
    yline(2.0, 'g--', 'LineWidth', 1.5, 'Label', '2.0');
    if ~isempty(adopted_idx)
        xline(adopted_Rinflate, 'k-', 'LineWidth', 2, 'Label', 'Adopted');
    end
    xlabel('R_{inflate}', 'FontSize', 12);
    ylabel('Mean NIS', 'FontSize', 12);
    title('Mean NIS vs R_{inflate}', 'FontSize', 12);
    grid on; set(gca,'FontSize',11);

    sgtitle('Effect of Measurement Noise Inflation on RMSE and NIS', ...
        'FontSize', 13, 'FontWeight', 'bold');
    if SAVE_FIGS
        print(fig3, fullfile(FIG_DIR,'B3_Rinflate_sweep'), '-dpng', '-r150');
        fprintf('  Saved: B3_Rinflate_sweep.png\n');
    end
end

%% ════════════════════════════════════════════════════════════════════
%%  BLOCK B4 — Initial Covariance Sweep
%% ════════════════════════════════════════════════════════════════════
if RUN_BLOCK.B4_P0
    fprintf('══════ BLOCK B4: Initial Covariance Sweep (N_MC=%d) ══════\n', N_MC);

    % B4a — Position prior σ_EN sweep
    pos_std_vals = [100, 200, 350, 500, 700, 1000, 1500];
    nP = numel(pos_std_vals);
    adopted_pos_std = 500;

    B4_results.pos_std_vals = pos_std_vals;
    B4_results.a_RMSE_h_mean  = nan(1,nP);
    B4_results.a_RMSE_h_std   = nan(1,nP);
    B4_results.a_convTime_mean = nan(1,nP);

    fprintf('  B4a — Position prior sweep\n');
    for vi = 1:nP
        pstd         = pos_std_vals(vi);
        seed_base_vi = BASE_SEED + SEED_OFF.B4 + (vi-1)*N_MC;
        rmse_v = nan(N_MC,1);
        conv_v = nan(N_MC,1);
        parfor mc = 1:N_MC
            cfg_r = build_baseline_cfg();
            cfg_r.init.sigENU_m = [pstd; pstd; 10];  %#ok<PFBNS>
            cfg_r.plot.enable   = false;
            cfg_r.verbose       = false;
            cfg_r.seedMeas      = seed_base_vi + mc;   %#ok<PFBNS>
            try
                out_r = real_strtg_sequential_EKF_fastrealistic(cfg_r);
                m_r   = ekf_extract_metrics(out_r, cfg_r);
                rmse_v(mc) = m_r.RMSE_h_m;
                conv_v(mc) = m_r.convTime_horiz_s;
            catch
            end
        end
        B4_results.a_RMSE_h_mean(vi)  = mean(rmse_v, 'omitnan');
        B4_results.a_RMSE_h_std(vi)   = std(rmse_v, 'omitnan');
        B4_results.a_convTime_mean(vi) = mean(conv_v(isfinite(conv_v)), 'omitnan');
        fprintf('    σ_EN=%4d m: RMSE_h=%.1f±%.1f m\n', pstd, ...
            B4_results.a_RMSE_h_mean(vi), B4_results.a_RMSE_h_std(vi));
    end

    % B4b — Clock bias prior σ_b sweep
    clk_std_vals = [0.02, 0.05, 0.10, 0.20, 0.50, 1.00, 2.00];
    nC = numel(clk_std_vals);
    adopted_clk_std = 0.20;

    B4_results.clk_std_vals  = clk_std_vals;
    B4_results.b_RMSE_h_mean = nan(1,nC);
    B4_results.b_NIS_mean    = nan(1,nC);

    fprintf('  B4b — Clock bias prior sweep\n');
    for vi = 1:nC
        cstd         = clk_std_vals(vi);
        seed_base_vi = BASE_SEED + SEED_OFF.B4 + 50 + (vi-1)*N_MC;
        rmse_v = nan(N_MC,1);
        nis_v  = nan(N_MC,1);
        parfor mc = 1:N_MC
            cfg_r = build_baseline_cfg();
            cfg_r.P0_clkBias  = cstd^2;            %#ok<PFBNS>
            cfg_r.plot.enable = false;
            cfg_r.verbose     = false;
            cfg_r.seedMeas    = seed_base_vi + mc;  %#ok<PFBNS>
            try
                out_r = real_strtg_sequential_EKF_fastrealistic(cfg_r);
                m_r   = ekf_extract_metrics(out_r, cfg_r);
                rmse_v(mc) = m_r.RMSE_h_m;
                nis_v(mc)  = m_r.meanNIS;
            catch
            end
        end
        B4_results.b_RMSE_h_mean(vi) = mean(rmse_v, 'omitnan');
        B4_results.b_NIS_mean(vi)    = mean(nis_v, 'omitnan');
        fprintf('    σ_clk=%.3f m/s: RMSE_h=%.1f m, NIS=%.3f\n', cstd, ...
            B4_results.b_RMSE_h_mean(vi), B4_results.b_NIS_mean(vi));
    end

    % Check plateau (B4a)
    low_end  = B4_results.a_RMSE_h_mean(pos_std_vals <= 350);
    high_end = B4_results.a_RMSE_h_mean(pos_std_vals >= 500);
    plateau_ok = (std(high_end, 'omitnan') / max(mean(high_end, 'omitnan'), 1)) < 0.20;
    B4_results.plateau_confirmed = plateau_ok;
    fprintf('  P0 position insensitivity plateau: %s\n\n', tern_str(plateau_ok,'YES','NO'));

    % ── Plot B4 ──────────────────────────────────────────────────────
    fig4 = figure('Name','B4_P0_sweep','Position',[50 50 1100 500]);

    subplot(1,2,1);
    errorbar(pos_std_vals, B4_results.a_RMSE_h_mean, B4_results.a_RMSE_h_std, ...
        'b-o','LineWidth',1.5,'MarkerFaceColor','b','MarkerSize',7); hold on;
    xline(adopted_pos_std,'k-','LineWidth',2,'Label','Adopted 500 m');
    set(gca,'XScale','log');
    xlabel('\sigma_{P0,pos} [m]', 'FontSize', 12);
    ylabel('RMSE_h [m]', 'FontSize', 12);
    title('B4a: RMSE vs Initial Position Prior', 'FontSize', 12);
    grid on; set(gca,'FontSize',11);

    subplot(1,2,2);
    yyaxis left;
    plot(clk_std_vals, B4_results.b_RMSE_h_mean, 'b-o', ...
        'LineWidth',1.5,'MarkerFaceColor','b','MarkerSize',7); hold on;
    ylabel('RMSE_h [m]', 'FontSize', 12);
    yyaxis right;
    plot(clk_std_vals, B4_results.b_NIS_mean, 'r-s', ...
        'LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',7);
    ylabel('Mean NIS', 'FontSize', 12);
    xline(adopted_clk_std,'k-','LineWidth',2,'Label','Adopted 0.2 m/s');
    set(gca,'XScale','log');
    xlabel('\sigma_{P0,b_c} [m/s]', 'FontSize', 12);
    title('B4b: RMSE & NIS vs Clock Bias Prior', 'FontSize', 12);
    grid on; set(gca,'FontSize',11);

    sgtitle('Sensitivity to Initial Covariance', 'FontSize', 13, 'FontWeight', 'bold');
    if SAVE_FIGS
        print(fig4, fullfile(FIG_DIR,'B4_P0_sweep'), '-dpng', '-r150');
        fprintf('  Saved: B4_P0_sweep.png\n');
    end
end

%% ════════════════════════════════════════════════════════════════════
%%  BLOCK B5 — Process Noise Sweep
%% ════════════════════════════════════════════════════════════════════
if RUN_BLOCK.B5_Q
    fprintf('══════ BLOCK B5: Process Noise Sweep (N_MC=%d) ══════\n', N_MC);

    % B5a — Position RW
    pos_rw_vals    = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
    nPR = numel(pos_rw_vals);
    adopted_pos_rw = 1e-4;

    B5_results.pos_rw_vals   = pos_rw_vals;
    B5_results.a_RMSE_h_mean = nan(1,nPR);
    B5_results.a_RMSE_h_std  = nan(1,nPR);

    fprintf('  B5a — Position RW sweep\n');
    for vi = 1:nPR
        prw          = pos_rw_vals(vi);
        seed_base_vi = BASE_SEED + SEED_OFF.B5 + (vi-1)*N_MC;
        rmse_v = nan(N_MC,1);
        parfor mc = 1:N_MC
            cfg_r = build_baseline_cfg();
            cfg_r.sigma_pos_rw = prw;              %#ok<PFBNS>
            cfg_r.plot.enable  = false;
            cfg_r.verbose      = false;
            cfg_r.seedMeas     = seed_base_vi + mc; %#ok<PFBNS>
            try
                out_r = real_strtg_sequential_EKF_fastrealistic(cfg_r);
                m_r   = ekf_extract_metrics(out_r, cfg_r);
                rmse_v(mc) = m_r.RMSE_h_m;
            catch
            end
        end
        B5_results.a_RMSE_h_mean(vi) = mean(rmse_v, 'omitnan');
        B5_results.a_RMSE_h_std(vi)  = std(rmse_v, 'omitnan');
        fprintf('    σ_pos=%s m/√s: RMSE_h=%.1f±%.1f m\n', ...
            num2str(prw,'%.2e'), B5_results.a_RMSE_h_mean(vi), B5_results.a_RMSE_h_std(vi));
    end

    % Check plateau
    mid_idx  = pos_rw_vals >= 1e-5 & pos_rw_vals <= 1e-3;
    plateau_ok_q = (std(B5_results.a_RMSE_h_mean(mid_idx), 'omitnan') / ...
                    max(mean(B5_results.a_RMSE_h_mean(mid_idx), 'omitnan'), 1)) < 0.20;
    B5_results.plateau_confirmed = plateau_ok_q;

    % B5b — Clock Q scale
    clk_q_scale = [0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0];
    nQS = numel(clk_q_scale);
    % baseline values
    sigma_b_base  = 1e-2;
    sigma_bd_base = 5e-4;

    B5_results.clk_q_scale   = clk_q_scale;
    B5_results.b_RMSE_h_mean = nan(1,nQS);
    B5_results.b_NIS_mean    = nan(1,nQS);
    B5_results.b_clkBiasRMSE = nan(1,nQS);
    B5_results.b_clkDriftRMSE= nan(1,nQS);

    fprintf('  B5b — Clock Q scale sweep\n');
    for vi = 1:nQS
        sc              = clk_q_scale(vi);
        seed_base_vi    = BASE_SEED + SEED_OFF.B5 + 50 + (vi-1)*N_MC;
        sb_scaled       = sigma_b_base  * sc;   % broadcast scalars
        sbd_scaled      = sigma_bd_base * sc;
        rmse_v = nan(N_MC,1);
        nis_v  = nan(N_MC,1);
        cb_v   = nan(N_MC,1);
        cd_v   = nan(N_MC,1);
        parfor mc = 1:N_MC
            cfg_r = build_baseline_cfg();
            cfg_r.clk.filt.sigmaBias_rw_mps_sqrt   = sb_scaled;   %#ok<PFBNS>
            cfg_r.clk.filt.sigmaDrift_rw_mps2_sqrt = sbd_scaled;  %#ok<PFBNS>
            cfg_r.plot.enable = false;
            cfg_r.verbose     = false;
            cfg_r.seedMeas    = seed_base_vi + mc;                 %#ok<PFBNS>
            try
                out_r = real_strtg_sequential_EKF_fastrealistic(cfg_r);
                m_r   = ekf_extract_metrics(out_r, cfg_r);
                rmse_v(mc) = m_r.RMSE_h_m;
                nis_v(mc)  = m_r.meanNIS;
                cb_v(mc)   = m_r.clkBiasRMSE_mps;
                cd_v(mc)   = m_r.clkDriftRMSE_mps2;
            catch
            end
        end
        B5_results.b_RMSE_h_mean(vi)  = mean(rmse_v, 'omitnan');
        B5_results.b_NIS_mean(vi)     = mean(nis_v, 'omitnan');
        B5_results.b_clkBiasRMSE(vi)  = mean(cb_v, 'omitnan');
        B5_results.b_clkDriftRMSE(vi) = mean(cd_v, 'omitnan');
        fprintf('    scale=%5.2f: RMSE_h=%.1f m, NIS=%.3f\n', sc, ...
            B5_results.b_RMSE_h_mean(vi), B5_results.b_NIS_mean(vi));
    end
    fprintf('  Q_pos plateau: %s\n\n', tern_str(plateau_ok_q,'YES','NO'));

    % ── Plot B5 ──────────────────────────────────────────────────────
    fig5 = figure('Name','B5_process_noise_sweep','Position',[50 50 1100 500]);

    subplot(1,2,1);
    % Use finite values only for log-x plot
    valid_pr = pos_rw_vals > 0;
    errorbar(pos_rw_vals(valid_pr), B5_results.a_RMSE_h_mean(valid_pr), ...
        B5_results.a_RMSE_h_std(valid_pr), 'b-o','LineWidth',1.5, ...
        'MarkerFaceColor','b','MarkerSize',7); hold on;
    if ~isnan(B5_results.a_RMSE_h_mean(1))
        plot(1e-7, B5_results.a_RMSE_h_mean(1), 'bs','MarkerFaceColor','b','MarkerSize',9, ...
            'DisplayName','\sigma_{pos}=0');
    end
    xline(adopted_pos_rw,'k-','LineWidth',2,'Label','Adopted 1e-4');
    set(gca,'XScale','log');
    xlabel('\sigma_{pos} [m/\surd{}s]', 'FontSize', 12);
    ylabel('RMSE_h [m]', 'FontSize', 12);
    title('B5a: RMSE vs Position RW', 'FontSize', 12);
    % Bracket annotation for plateau
    xlim_cur = get(gca,'XLim');
    text(1e-3, 0.95*max(B5_results.a_RMSE_h_mean(valid_pr)), '← flat plateau →', ...
        'FontSize', 10, 'Color', [0 0.5 0]);
    grid on; set(gca,'FontSize',11);

    subplot(1,2,2);
    yyaxis left;
    plot(clk_q_scale, B5_results.b_RMSE_h_mean, 'b-o', ...
        'LineWidth',1.5,'MarkerFaceColor','b','MarkerSize',7);
    ylabel('RMSE_h [m]', 'FontSize', 12);
    yyaxis right;
    plot(clk_q_scale, B5_results.b_NIS_mean, 'r-s', ...
        'LineWidth',1.5,'MarkerFaceColor','r','MarkerSize',7);
    ylabel('Mean NIS', 'FontSize', 12);
    xline(1.0,'k-','LineWidth',2,'Label','Adopted (scale=1)');
    yline(1.0,'k:','LineWidth',1,'Label','NIS=1');
    set(gca,'XScale','log');
    xlabel('Clock Q scale factor', 'FontSize', 12);
    title('B5b: RMSE & NIS vs Clock Q Scale', 'FontSize', 12);
    grid on; set(gca,'FontSize',11);

    sgtitle('B5 — Process Noise Sensitivity', 'FontSize', 13, 'FontWeight', 'bold');
    if SAVE_FIGS
        print(fig5, fullfile(FIG_DIR,'B5_process_noise_sweep'), '-dpng', '-r150');
        fprintf('  Saved: B5_process_noise_sweep.png\n');
    end
end

%% ════════════════════════════════════════════════════════════════════
%%  BLOCK B6 — Final Filter Consistency Diagnostics (Adopted Config)
%% ════════════════════════════════════════════════════════════════════
if RUN_BLOCK.B6_NIS_ACF
    fprintf('══════ BLOCK B6: Final Diagnostics (N_MC=%d, adopted config) ══════\n', N_MC);

    % Sliced output temporaries — parfor writes to these, assigned to
    % B6_results after the loop (parfor cannot write to struct fields)
    metrics_mc  = cell(N_MC,1);
    out_mc_list = cell(N_MC,1);
    seed_base_b6 = BASE_SEED + SEED_OFF.B6;   % broadcast scalar

    parfor mc = 1:N_MC
        cfg_r = build_baseline_cfg();
        cfg_r.plot.enable = false;
        cfg_r.verbose     = false;
        cfg_r.seedMeas    = seed_base_b6 + mc;  %#ok<PFBNS>
        try
            out_r = real_strtg_sequential_EKF_fastrealistic(cfg_r);
            m_r   = ekf_extract_metrics(out_r, cfg_r);
            metrics_mc{mc}  = m_r;
            out_mc_list{mc} = out_r;
        catch ME
            fprintf('  [B6 mc=%d] failed: %s\n', mc, ME.message);
            metrics_mc{mc}  = struct('RMSE_h_m',NaN,'meanNIS',NaN,...
                'fracSigACF_dopp',NaN,'fracSigACF_rate',NaN,...
                'excessKurtosis_dopp',NaN,'excessKurtosis_rate',NaN,...
                'convTime_horiz_s',Inf);
            out_mc_list{mc} = [];
        end
    end

    % Assign to results struct after parfor completes
    B6_results.metrics  = metrics_mc;
    B6_results.out_list = out_mc_list;

    % Per-seed summary (printed after parfor, order is deterministic here)
    for mc = 1:N_MC
        fprintf('  Seed %2d/%2d: RMSE_h=%.1f m, NIS=%.3f\n', mc, N_MC, ...
            B6_results.metrics{mc}.RMSE_h_m, B6_results.metrics{mc}.meanNIS);
    end

    % Find representative seed (closest to mean RMSE_h)
    rmse_all = cellfun(@(m) m.RMSE_h_m, B6_results.metrics);
    valid_runs = find(isfinite(rmse_all));
    if isempty(valid_runs)
        rep_idx = 1;
    else
        [~, ri_local] = min(abs(rmse_all(valid_runs) - mean(rmse_all, 'omitnan')));
        rep_idx = valid_runs(ri_local);
    end
    out_rep = B6_results.out_list{rep_idx};
    fprintf('  Representative seed: %d (RMSE_h=%.1f m)\n', rep_idx, rmse_all(rep_idx));

    % Aggregate statistics
    B6_results.RMSE_h_mean        = mean(rmse_all, 'omitnan');
    B6_results.RMSE_h_std         = std(rmse_all, 'omitnan');
    B6_results.meanNIS_ensemble   = mean(cellfun(@(m) m.meanNIS, B6_results.metrics), 'omitnan');
    B6_results.fracACF_D_mean     = mean(cellfun(@(m) m.fracSigACF_dopp, B6_results.metrics), 'omitnan');
    B6_results.fracACF_DR_mean    = mean(cellfun(@(m) m.fracSigACF_rate, B6_results.metrics), 'omitnan');
    B6_results.kurtosis_D_mean    = mean(cellfun(@(m) m.excessKurtosis_dopp, B6_results.metrics), 'omitnan');
    B6_results.kurtosis_DR_mean   = mean(cellfun(@(m) m.excessKurtosis_rate, B6_results.metrics), 'omitnan');
    conv_vals = cellfun(@(m) m.convTime_horiz_s, B6_results.metrics);
    B6_results.frac_converged     = sum(isfinite(conv_vals) & conv_vals < Inf);

    %% ── Plot B6a: Innovation sequence (representative seed) ─────────
    if ~isempty(out_rep)
        t_rep     = out_rep.t(:);
        nu_D      = out_rep.nu_hist(:);
        if isfield(out_rep, 'nu_rate_hist') && ~isempty(out_rep.nu_rate_hist)
            nu_DR = out_rep.nu_rate_hist(:);
        else
            nu_DR = nan(size(t_rep));
        end
        nis_rep   = out_rep.nis(:);
        K_rep     = numel(t_rep);

        % Running mean NIS (window 401 epochs)
        win = 401;
        nis_run = running_mean(nis_rep, win);

        % Chi-squared 95% CI for chi^2(1)  (χ²(1) 2.5th/97.5th percentiles)
        if exist('chi2inv','file') || exist('chi2inv','builtin')
            chi2_lo = chi2inv(0.025, 1);   % ≈ 0.0010
            chi2_hi = chi2inv(0.975, 1);   % ≈ 5.024
        else
            chi2_lo = 0.001;   % hard-coded χ²(1) 2.5th percentile
            chi2_hi = 5.024;   % hard-coded χ²(1) 97.5th percentile
        end

        fig6a = figure('Name','B6a_innovations','Position',[50 50 1100 800]);
        subplot(4,1,1);
        plot(t_rep, nu_D, 'b.', 'MarkerSize', 3); hold on;
        yline(0,'k:','LineWidth',1);
        xlabel('Time [s]','FontSize',11);
        ylabel('\nu_D [m/s]','FontSize',11);
        title('Doppler Innovation','FontSize',12);
        grid on; set(gca,'FontSize',11);

        subplot(4,1,2);
        plot(t_rep, nu_DR, 'r.', 'MarkerSize', 3); hold on;
        yline(0,'k:','LineWidth',1);
        xlabel('Time [s]','FontSize',11);
        ylabel('\nu_{DR} [m/s^2]','FontSize',11);
        title('Doppler-Rate Innovation','FontSize',12);
        grid on; set(gca,'FontSize',11);

        subplot(4,1,3);
        plot(t_rep, nis_rep, 'b.', 'MarkerSize', 2); hold on;
        yline(chi2_lo,'g--','LineWidth',1.5,'Label',sprintf('\\chi^2(1) 95%% lo=%.2f',chi2_lo));
        yline(chi2_hi,'g--','LineWidth',1.5,'Label',sprintf('\\chi^2(1) 95%% hi=%.2f',chi2_hi));
        yline(1.0,'k:','LineWidth',1.5,'Label','NIS=1');
        xlabel('Time [s]','FontSize',11);
        ylabel('NIS','FontSize',11);
        title('NIS Time Series','FontSize',12);
        ylim([0, min(50, max(chi2_hi*3, 5))]);
        grid on; set(gca,'FontSize',11);

        subplot(4,1,4);
        plot(t_rep, nis_run, 'b-', 'LineWidth', 1.5); hold on;
        yline(0.5,'g--','LineWidth',1.5,'Label','0.5');
        yline(2.0,'g--','LineWidth',1.5,'Label','2.0');
        yline(1.0,'k:','LineWidth',1.5,'Label','1.0');
        xlabel('Time [s]','FontSize',11);
        ylabel('Running Mean NIS','FontSize',11);
        title(sprintf('Running Mean NIS (window=%d epochs)', win),'FontSize',12);
        grid on; set(gca,'FontSize',11);

        sgtitle(sprintf('B6a — Innovations (representative seed %d)', rep_idx), ...
            'FontSize',13,'FontWeight','bold');
        if SAVE_FIGS
            print(fig6a, fullfile(FIG_DIR,'B6a_innovations'), '-dpng', '-r150');
            fprintf('  Saved: B6a_innovations.png\n');
        end

        %% ── Plot B6b: Innovation ACF ─────────────────────────────────
        fig6b = figure('Name','B6b_ACF','Position',[50 50 1000 450]);

        subplot(1,2,1);
        [frac_acf_D, ~, acf_D, lags_D, conf_D] = compute_acf_stats(nu_D);
        bar(lags_D(lags_D>0), acf_D(lags_D>0), 'b', 'FaceAlpha', 0.7); hold on;
        yline( conf_D, 'r--', 'LineWidth', 1.5, 'Label', '+1.96/√N');
        yline(-conf_D, 'r--', 'LineWidth', 1.5, 'Label', '-1.96/√N');
        yline(0,'k:','LineWidth',1);
        xlabel('Lag','FontSize',11); ylabel('ACF','FontSize',11);
        title(sprintf('Doppler Innovation ACF — WHITE: %.0f%% lags sig.', 100*frac_acf_D), ...
            'FontSize',11);
        xlim([0 50]); grid on; set(gca,'FontSize',11);

        subplot(1,2,2);
        [frac_acf_DR, ~, acf_DR, lags_DR, conf_DR] = compute_acf_stats(nu_DR);
        bar(lags_DR(lags_DR>0), acf_DR(lags_DR>0), 'r', 'FaceAlpha', 0.7); hold on;
        yline( conf_DR, 'b--', 'LineWidth', 1.5, 'Label', '+1.96/√N');
        yline(-conf_DR, 'b--', 'LineWidth', 1.5, 'Label', '-1.96/√N');
        yline(0,'k:','LineWidth',1);
        xlabel('Lag','FontSize',11); ylabel('ACF','FontSize',11);
        title(sprintf('Doppler-Rate Innovation ACF — WHITE: %.0f%% lags sig.', 100*frac_acf_DR), ...
            'FontSize',11);
        xlim([0 50]); grid on; set(gca,'FontSize',11);

        sgtitle(sprintf('B6b — Innovation ACF (seed %d)', rep_idx), ...
            'FontSize',13,'FontWeight','bold');
        if SAVE_FIGS
            print(fig6b, fullfile(FIG_DIR,'B6b_ACF'), '-dpng', '-r150');
            fprintf('  Saved: B6b_ACF.png\n');
        end

        %% ── Plot B6c: Innovation histograms ──────────────────────────
        fig6c = figure('Name','B6c_innovation_histograms','Position',[50 50 1000 450]);

        subplot(1,2,1);
        [kurt_D, mean_D, std_D] = plot_innovation_hist(nu_D, 'b', 'Doppler'); hold on;
        x_pdf = linspace(-5,5,200);
        plot(x_pdf, normpdf(x_pdf), 'r-', 'LineWidth', 2, 'DisplayName', 'N(0,1)');
        % Annotation
        kurt_note = '';
        if kurt_D > 1
            kurt_note = sprintf('\n(AR(1) coloured noise\ntau=0.5 s; ACF white)');
        end
        txt_d = sprintf('mean=%.3f\nstd=%.3f\nkurt=%.2f%s', mean_D, std_D, kurt_D, kurt_note);
        text(-4.5, 0.9*normpdf(0), txt_d, 'FontSize', 9, 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k');
        legend('Location','northeast','FontSize',10); grid on; set(gca,'FontSize',11);

        subplot(1,2,2);
        [kurt_DR, mean_DR, std_DR] = plot_innovation_hist(nu_DR, [1 0.5 0], 'Doppler-Rate'); hold on;
        plot(x_pdf, normpdf(x_pdf), 'r-', 'LineWidth', 2, 'DisplayName', 'N(0,1)');
        txt_dr = sprintf('mean=%.3f\nstd=%.3f\nkurt=%.2f', mean_DR, std_DR, kurt_DR);
        text(-4.5, 0.9*normpdf(0), txt_dr, 'FontSize', 9, 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k');
        legend('Location','northeast','FontSize',10); grid on; set(gca,'FontSize',11);

        sgtitle(sprintf('B6c — Normalised Innovation Histograms (seed %d)', rep_idx), ...
            'FontSize',13,'FontWeight','bold');
        if SAVE_FIGS
            print(fig6c, fullfile(FIG_DIR,'B6c_innovation_histograms'), '-dpng', '-r150');
            fprintf('  Saved: B6c_innovation_histograms.png\n');
        end

        %% ── Plot B6d: NIS distribution ───────────────────────────────
        fig6d = figure('Name','B6d_NIS_distribution','Position',[50 50 1000 450]);

        % Collect all NIS values for histogram
        nis_all = [];
        for mc = 1:N_MC
            if ~isempty(B6_results.out_list{mc})
                n_mc = B6_results.out_list{mc}.nis(:);
                nis_all = [nis_all; n_mc(~isnan(n_mc) & n_mc>0)]; %#ok<AGROW>
            end
        end

        subplot(1,2,1);
        if ~isempty(nis_all)
            histogram(nis_all, 'Normalization','pdf','FaceColor',[0.2 0.5 0.8], ...
                'FaceAlpha',0.7,'EdgeColor','none'); hold on;
        end
        x_chi2_max = max(20, prctile(nis_all, 99) * 1.1);
        if isempty(nis_all) || ~isfinite(x_chi2_max), x_chi2_max = 20; end
        x_chi2 = linspace(0, x_chi2_max, 200);
        if exist('chi2pdf','file') || exist('chi2pdf','builtin')
            chi2_curve = chi2pdf(x_chi2, 1);
        else
            chi2_curve = (1./(sqrt(2*pi)*sqrt(x_chi2))) .* exp(-x_chi2/2);  % χ²(1) PDF
            chi2_curve(x_chi2 <= 0) = 0;
        end
        plot(x_chi2, chi2_curve, 'r-', 'LineWidth', 2, 'DisplayName', '\chi^2(1)');
        xlabel('NIS','FontSize',11); ylabel('PDF','FontSize',11);
        title(sprintf('NIS Distribution (all MC epochs)\nmeanNIS=%.3f', mean(nis_all, 'omitnan')), ...
            'FontSize',11);
        if ~isempty(nis_all) && isfinite(prctile(nis_all,99))
            xlim([0, min(20, prctile(nis_all,99))]);
        end
        legend('Location','northeast','FontSize',10); grid on; set(gca,'FontSize',11);

        subplot(1,2,2);
        hold on;
        nis_cumul = zeros(N_MC, K_rep);
        for mc = 1:N_MC
            if ~isempty(B6_results.out_list{mc})
                n_mc = B6_results.out_list{mc}.nis(:);
                n_mc(isnan(n_mc)) = NaN;
                % running mean
                K_mc = min(numel(n_mc), K_rep);
                rm_mc = running_mean(n_mc(1:K_mc), win);
                nis_cumul(mc, 1:K_mc) = rm_mc;
                lh = plot(t_rep(1:K_mc), rm_mc, 'b-', 'LineWidth', 0.8);
                lh.Color(4) = 0.35;   % transparency (R2022a+); silently ignored in older MATLAB
            end
        end
        ens_mean = mean(nis_cumul, 1, 'omitnan');
        plot(t_rep, ens_mean, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Ensemble mean');
        yline(1.0,'k--','LineWidth',1.5,'Label','NIS=1');
        xlabel('Time [s]','FontSize',11); ylabel('Running Mean NIS','FontSize',11);
        title(sprintf('NIS Convergence over %d MC Seeds', N_MC),'FontSize',11);
        grid on; set(gca,'FontSize',11);
        legend('Location','northeast','FontSize',10);

        sgtitle('B6d — NIS Distribution and Convergence', 'FontSize',13,'FontWeight','bold');
        if SAVE_FIGS
            print(fig6d, fullfile(FIG_DIR,'B6d_NIS_distribution'), '-dpng', '-r150');
            fprintf('  Saved: B6d_NIS_distribution.png\n');
        end
    end

    %% ── Plot B6e: MC error bands ─────────────────────────────────────
    if ~isempty(out_rep)
        t_ref_b6 = out_rep.t(:)';
    else
        % Fall back to time from any successful run
        t_ref_b6 = [];
        for mc_fb = 1:N_MC
            if ~isempty(B6_results.out_list{mc_fb})
                t_ref_b6 = B6_results.out_list{mc_fb}.t(:)';
                break;
            end
        end
    end
    if isempty(t_ref_b6)
        warning('B6e: no successful runs found — skipping error band plot.');
        K_b6 = 0;
    else
        K_b6 = numel(t_ref_b6);
    end

    eE_mat  = nan(N_MC, max(K_b6,1));
    eN_mat  = nan(N_MC, max(K_b6,1));
    eH_mat  = nan(N_MC, max(K_b6,1));
    eCB_mat = nan(N_MC, max(K_b6,1));

    for mc = 1:N_MC
        out_mc = B6_results.out_list{mc};
        if isempty(out_mc), continue; end
        K_mc = min(size(out_mc.errENU,2), K_b6);
        eE_mat(mc,1:K_mc)  = out_mc.errENU(1,1:K_mc);
        eN_mat(mc,1:K_mc)  = out_mc.errENU(2,1:K_mc);
        eH_mat(mc,1:K_mc)  = sqrt(out_mc.errENU(1,1:K_mc).^2 + out_mc.errENU(2,1:K_mc).^2);
        % clock bias error: estimated - truth
        xCB  = out_mc.xhist(4,1:K_mc);
        cbT  = out_mc.truthAux.clkBiasTrue(1:K_mc)';
        eCB_mat(mc,1:K_mc) = xCB - cbT;
    end

    if K_b6 > 0
        fig6e = figure('Name','B6e_MC_error_bands','Position',[50 50 1200 800]);

        panel_data = { ...
            eE_mat,  'East Error [m]',          'r',  'East Error (ENU)'; ...
            eN_mat,  'North Error [m]',         'g',  'North Error (ENU)'; ...
            eH_mat,  '2D Horiz. Error [m]',     'b',  '2D Horizontal Error'; ...
            eCB_mat, 'Clock Bias Error [m/s]',  [0.5 0 0.5], 'Clock Bias Error' ...
        };

        for pp = 1:4
            subplot(2,2,pp);
            dat  = panel_data{pp,1};
            ylbl = panel_data{pp,2};
            col  = panel_data{pp,3};
            ttle = panel_data{pp,4};

            mu_p  = mean(dat, 1, 'omitnan');
            sg_p  = std(dat, 0, 1, 'omitnan');
            sg_lo = mu_p - sg_p;   % keep y-range symmetric
            fill([t_ref_b6, fliplr(t_ref_b6)], [mu_p+sg_p, fliplr(sg_lo)], col, ...
                'FaceAlpha',0.2,'EdgeColor','none'); hold on;
            plot(t_ref_b6, mu_p, '-', 'Color', col, 'LineWidth', 2.0);
            if pp == 3
                yline(50,'r--','LineWidth',1.5,'Label','50 m');
            end
            xlabel('Time [s]','FontSize',11); ylabel(ylbl,'FontSize',11);
            title(ttle,'FontSize',12);
            grid on; set(gca,'FontSize',11);
        end
        sgtitle(sprintf('B6e — MC Error Bands (N_{MC}=%d, adopted config)', N_MC), ...
            'FontSize',13,'FontWeight','bold');
        if SAVE_FIGS
            print(fig6e, fullfile(FIG_DIR,'B6e_MC_error_bands'), '-dpng', '-r150');
            fprintf('  Saved: B6e_MC_error_bands.png\n');
        end
    end
end

%% ════════════════════════════════════════════════════════════════════
%%  FINAL CONSOLE SUMMARY
%% ════════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║  PHASE 0 — FILTER VALIDATION SUMMARY                        ║\n');
fprintf('╠═══════════════════════════════════════════════════════════════╣\n');

% B1
if RUN_BLOCK.B1_progressive && ~isempty(B1_results.metrics) && ~isempty(B1_results.metrics{1})
    fprintf('║  B1 S1 Jacobian check: max|nu| = %.4f m/s              ║\n', ...
        B1_results.max_innov(1));
    fprintf('║  B1 S6 Full model NIS: %.3f (target ~0.98)             ║\n', ...
        B1_results.metrics{6}.meanNIS);
else
    fprintf('║  B1: not run                                                  ║\n');
end
fprintf('║                                                               ║\n');

% B2
if RUN_BLOCK.B2_ekf_vs_iekf && isfield(B2_results,'metrics') && ~isempty(B2_results.metrics)
    fprintf('║  B2 EKF vs IEKF:                                              ║\n');
    for ci = 1:2
        rmse_v = cellfun(@(m) m.RMSE_h_m, B2_results.metrics(ci,:));
        conv_v = cellfun(@(m) m.convTime_horiz_s, B2_results.metrics(ci,:));
        fc = sum(isfinite(conv_v) & conv_v < Inf);
        cfg_lbl = tern_str(ci==1,'EKF ','IEKF');
        fprintf('║    %s RMSE_h = %5.1f m,  frac_conv = %d/%d\n', ...
            cfg_lbl, mean(rmse_v, 'omitnan'), fc, N_MC);
    end
else
    fprintf('║  B2: not run                                                  ║\n');
end
fprintf('║                                                               ║\n');

% B3
if RUN_BLOCK.B3_Rinflate
    idx_ad = find(B3_results.Rinflate_vals == 2.0, 1);
    if ~isempty(idx_ad)
        is_cons = (B3_results.NIS_mean(idx_ad) >= 0.5 && B3_results.NIS_mean(idx_ad) <= 2.0);
        fprintf('║  B3 Adopted Rinflate=2.0:                                     ║\n');
        fprintf('║    RMSE_h = %5.1f m,  meanNIS = %.3f  [CONSISTENT: %s]   ║\n', ...
            B3_results.RMSE_h_mean(idx_ad), B3_results.NIS_mean(idx_ad), ...
            tern_str(is_cons,'Y','N'));
    end
else
    fprintf('║  B3: not run                                                  ║\n');
end
fprintf('║                                                               ║\n');

% B4 & B5
if RUN_BLOCK.B4_P0 && isfield(B4_results,'plateau_confirmed')
    fprintf('║  B4 P0 insensitivity: plateau confirmed %s                  ║\n', ...
        tern_str(B4_results.plateau_confirmed,'Y/N — Y','Y/N — N'));
else
    fprintf('║  B4: not run                                                  ║\n');
end
if RUN_BLOCK.B5_Q && isfield(B5_results,'plateau_confirmed')
    fprintf('║  B5 Q_pos insensitivity: plateau confirmed %s               ║\n', ...
        tern_str(B5_results.plateau_confirmed,'Y/N — Y','Y/N — N'));
else
    fprintf('║  B5: not run                                                  ║\n');
end
fprintf('║                                                               ║\n');

% B6
if RUN_BLOCK.B6_NIS_ACF && isfield(B6_results,'metrics') && ...
        ~isempty(B6_results.metrics) && ~isempty(B6_results.metrics{1})
    kD_note  = tern_str(B6_results.kurtosis_D_mean > 1, '  (AR(1) coloured noise)', '');
    kDR_note = tern_str(B6_results.kurtosis_DR_mean > 1, '  (AR(1))', '');
    fprintf('║  B6 Final diagnostics (N_MC=%2d):                              ║\n', N_MC);
    fprintf('║    meanNIS        = %.3f  (target ~1.0)                   ║\n', B6_results.meanNIS_ensemble);
    fprintf('║    fracSigACF_D   = %.1f%%   (target <5%%)                    ║\n', 100*B6_results.fracACF_D_mean);
    fprintf('║    fracSigACF_DR  = %.1f%%   (target <5%%)                    ║\n', 100*B6_results.fracACF_DR_mean);
    fprintf('║    kurtosis_D     = %.2f%s\n', B6_results.kurtosis_D_mean, kD_note);
    fprintf('║    kurtosis_DR    = %.2f%s\n', B6_results.kurtosis_DR_mean, kDR_note);
    fprintf('║    frac_converged = %d/%d\n', B6_results.frac_converged, N_MC);
else
    fprintf('║  B6: not run                                                  ║\n');
end
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

%% ── Save results ─────────────────────────────────────────────────────
save('phase0_results.mat', 'B1_results', 'B2_results', 'B3_results', ...
    'B4_results', 'B5_results', 'B6_results', 'baseline_cfg');
fprintf('Results saved to phase0_results.mat\n\n');

%% ════════════════════════════════════════════════════════════════════
%%  LOCAL HELPER FUNCTIONS
%% ════════════════════════════════════════════════════════════════════

% ─────────────────────────────────────────────────────────────────────
function cfg = build_baseline_cfg()
% BASELINE CONFIGURATION for Phase 0 validation.
% Follows exact parameter values from the validation specification.

cfg.tleFile              = fullfile('InputData','oneweb.tle');
cfg.duration_sec         = 100;
cfg.dt                   = 0.01;
cfg.startUTC             = [];
cfg.stopUTC              = [];
cfg.startFromTLE         = 'first';
cfg.maxTLEAgeDays        = 7;

%% Receiver (Cranfield, UK)
cfg.rx.lat_deg           = 52.069;
cfg.rx.lon_deg           = -0.609;
cfg.rx.alt_m             = 100;
cfg.rx.vel_ecef          = [0;0;0];
cfg.rx.acc_ecef          = [0;0;0];

%% Signal / frequency
cfg.fc_Hz                = 11.325e9;   % 11.325 GHz
cfg.c                    = 299792458;
% lambda derived internally by solver: c/fc = 0.026473... m

%% Constellation
cfg.minElevDeg           = 25;
cfg.coarseStep           = 10;
cfg.minMaxElev           = 30;
cfg.maxSats              = 12;
cfg.constellationSelectionWindow_s = 60;

%% Propagator
cfg.useSGP4              = false;
cfg.allowJ2Fallback      = true;

%% Initial position error / prior
cfg.init.errENU_m        = [400; 400; 0];   % guessed offset [m]
cfg.init.sigENU_m        = [500; 500; 10];  % prior std [dE dN dU] [m]

%% Initial covariance (clock)
cfg.P0_clkBias           = (0.2)^2;        % [m/s]^2
cfg.P0_clkDrift          = (5e-7)^2;       % [m/s^2]^2

%% Process noise
cfg.sigma_pos_rw         = 1e-4;           % position RW [m/sqrt(s)]

%% Clock model — truth
cfg.clk.truth.initBias_mps            =  0.10;     % [m/s]
cfg.clk.truth.initDrift_mps2          =  5e-5;     % [m/s^2]
cfg.clk.truth.sigmaBias_rw_mps_sqrt   =  1e-2;     % [m/s / sqrt(s)]
cfg.clk.truth.sigmaDrift_rw_mps2_sqrt =  3e-4;     % [m/s^2 / sqrt(s)]

%% Clock model — filter
cfg.clk.filt.sigmaBias_rw_mps_sqrt    =  1e-2;     % overmodels truth
cfg.clk.filt.sigmaDrift_rw_mps2_sqrt  =  5e-4;

cfg.clk.filter_init_bias_mps          =  0.0;
cfg.clk.filter_init_drift_mps2        =  0.0;

%% Troposphere
cfg.tropo.enable          = true;
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

%% Measurement model / noise
cfg.sigma_doppler_Hz       = 50.0;         % Doppler white noise [Hz]
cfg.sigma_dopplerRate_Hzps =  1.0;         % Doppler-rate white noise [Hz/s]

cfg.meas.enable_lighttime  = false;
cfg.meas.colored_tau_s     = 0.5;          % AR(1) time constant [s]
cfg.meas.sigma_col0_mps    = 0.03;         % coloured noise amplitude [m/s]
cfg.meas.Rinflate          = 2.0;          % adopted inflation factor

cfg.meas.useDoppler        = true;
cfg.meas.useDopplerRate    = true;
cfg.meas.fdJacStep_m       = 0.10;

%% IEKF numerical tuning
cfg.filter.useNISGate           = true;
cfg.filter.nisGate              = 16.0;
cfg.filter.maxIEKFIter          = 10;
cfg.filter.iekfTol              = 1e-5;
cfg.filter.maxStep_m            = 25;
cfg.filter.startupInflateEpochs = 50;
cfg.filter.startupRmult         = 25;

%% Height pseudo-measurement
cfg.height.enable  = true;
cfg.height.sigma_m = 10.0;
cfg.height.z0_m    = 0.0;

%% Capture policy
cfg.capture.mode                     = 'stickyBest';
cfg.capture.minDwell_s               = 30;
cfg.capture.switchBlackout_s         = 3;
cfg.capture.elev_margin_deg          = 20;
cfg.capture.useServiceAngleFilter    = false;
cfg.capture.maxServiceOffNadirDeg    = 45;
cfg.capture.switchThreshold_lowElev  = 10;
cfg.capture.switchThreshold_highElev = 15;
cfg.capture.highElevThreshold        = 35;

%% Convergence metric parameters
cfg.analysis.conv.refWindow_s    = 15;
cfg.analysis.conv.settleWindow_s = 10;
cfg.analysis.conv.relTol         = 0.15;
cfg.analysis.conv.absTol_m       = 5.0;
cfg.analysis.conv.maxCoeffVar    = 0.25;
cfg.analysis.conv.maxSlope_mps   = 1.0;
cfg.analysis.conv.errFloor_m     = 5.0;

%% Output control
cfg.seedMeas    = 42;
cfg.plot.enable = false;   % suppress per-run plots
cfg.verbose     = false;
cfg.lambda      = cfg.c / cfg.fc_Hz;   % wavelength [m]
end

% ─────────────────────────────────────────────────────────────────────
function s = tern_str(cond, a, b)
% Ternary string: return a if cond is true, else b
if cond
    s = a;
else
    s = b;
end
end

% ─────────────────────────────────────────────────────────────────────
function y = running_mean(x, win)
% Compute causal running mean with window 'win'.  NaN-safe.
x  = x(:)';
K  = numel(x);
y  = nan(1, K);
for k = 1:K
    i0 = max(1, k - floor(win/2));
    i1 = min(K, k + floor(win/2));
    v  = x(i0:i1);
    v  = v(~isnan(v));
    if ~isempty(v)
        y(k) = mean(v);
    end
end
end

% ─────────────────────────────────────────────────────────────────────
function [fracSig, exKurt, acf, lags, conf_bound] = compute_acf_stats(nu_raw)
% Compute ACF whiteness fraction, excess kurtosis, and ACF vector.
nu_raw = nu_raw(:);
valid  = ~isnan(nu_raw) & nu_raw ~= 0;
nu_v   = nu_raw(valid);

if numel(nu_v) < 20
    fracSig = NaN; exKurt = NaN;
    acf = []; lags = []; conf_bound = NaN;
    return;
end

sigma = std(nu_v);
if sigma < eps
    fracSig = NaN; exKurt = NaN;
    acf = []; lags = []; conf_bound = NaN;
    return;
end

nu_norm = (nu_v - mean(nu_v)) / sigma;
maxLag  = min(50, floor(numel(nu_norm)/4));
[acf_full, lags_full] = xcorr(nu_norm, maxLag, 'normalized');

acf  = acf_full;
lags = lags_full;

conf_bound = 1.96 / sqrt(numel(nu_norm));
pos_lags   = lags > 0;
fracSig    = mean(abs(acf(pos_lags)) > conf_bound);

n      = numel(nu_norm);
m4     = sum(nu_norm.^4) / n;
exKurt = m4 - 3;
end

% ─────────────────────────────────────────────────────────────────────
function [exKurt, mu_n, std_n] = plot_innovation_hist(nu_raw, col, chan_name)
% Plot normalised innovation histogram.
% Returns excess kurtosis, mean, std of the normalised sequence.
nu_raw = nu_raw(:);
valid  = ~isnan(nu_raw) & nu_raw ~= 0;
nu_v   = nu_raw(valid);

if numel(nu_v) < 10
    exKurt = NaN; mu_n = NaN; std_n = NaN; return;
end

mu_n  = mean(nu_v);
std_n = std(nu_v);
nu_n  = (nu_v - mu_n) / std_n;

histogram(nu_n, 50, 'Normalization','pdf', 'FaceColor', col, ...
    'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', [chan_name, ' (normalised)']);

n      = numel(nu_n);
m4     = sum(nu_n.^4) / n;
exKurt = m4 - 3;

xlabel('Normalised Innovation', 'FontSize', 11);
ylabel('PDF', 'FontSize', 11);
title(sprintf('%s Innovation vs N(0,1)', chan_name), 'FontSize', 12);
end
