%% ═══════════════════════════════════════════════════════════════════════
%  MAIN_SENSITIVITY   1D sweeps and 2D heatmaps for parameter sensitivity
%
%  Studies how performance (RMSE_h, meanNIS, convTime) changes as each
%  parameter is swept over a physically motivated range.
%
%  Usage:  Run this script (F5). Set sweep_mode to control which sweeps run.
%
%  Dependencies:
%    real_strtg_sequential_EKF_fastrealistic.m
%    ekf_extract_metrics.m
%    ekf_tropo_models.m
% ═══════════════════════════════════════════════════════════════════════

clear; close all; clc;

%% ── Sweep selection ──────────────────────────────────────────────────────
run_1d_sweeps = true;
run_2d_heatmaps = true;
nRepeat = 3;   % runs per parameter value (averaged) — set 1 for speed
saveResults = false;

%% ── Load baseline config ─────────────────────────────────────────────────
cfg_base = baseline_cfg();

%% ═══════════════════════════════════════════════════════════════════════════
%  1D SWEEPS
%% ═══════════════════════════════════════════════════════════════════════════
if run_1d_sweeps
    fprintf('=== 1D SENSITIVITY SWEEPS ===\n\n');

    sweeps = define_sweeps();
    nSweeps = numel(sweeps);
    sweep_results = struct();

    for si = 1:nSweeps
        sw = sweeps(si);
        fprintf('  Sweep [%d/%d]: %s  (%d points x %d repeats)\n', ...
            si, nSweeps, sw.label, numel(sw.values), nRepeat);

        result = run_sweep_1d(cfg_base, sw.param_path, sw.values, nRepeat);
        sweep_results.(sw.tag) = result;
        sweep_results.(sw.tag).label  = sw.label;
        sweep_results.(sw.tag).units  = sw.units;
        sweep_results.(sw.tag).values = sw.values;
    end

    plot_1d_sweeps(sweep_results, sweeps);
    plot_tornado(sweep_results, sweeps);
end

%% ═══════════════════════════════════════════════════════════════════════════
%  2D HEATMAPS
%% ═══════════════════════════════════════════════════════════════════════════
if run_2d_heatmaps
    fprintf('\n=== 2D HEATMAPS ===\n\n');

    % Key pair 1: sigma_pos_rw vs sigma_doppler_Hz
    vals1 = logspace(-5, -2, 8);   % sigma_pos_rw
    vals2 = logspace(0, 2, 8);     % sigma_doppler_Hz
    fprintf('  Heatmap 1: sigma_pos_rw vs sigma_doppler_Hz...\n');
    hm1 = run_sweep_2d(cfg_base, 'sigma_pos_rw', vals1, ...
                                 'sigma_doppler_Hz', vals2);
    plot_heatmap(hm1, 'sigma\_pos\_rw', 'sigma\_doppler\_Hz', ...
        vals1, vals2, 'RMSE_h_m', 'RMSE Horizontal (m)');

    % Key pair 2: P0_clkBias vs height.sigma_m
    vals3 = logspace(-4, 0, 8);    % P0_clkBias (variance: sqrt gives std)
    vals4 = logspace(0, 2.5, 8);   % height.sigma_m
    fprintf('  Heatmap 2: sqrt(P0_clkBias) [std] vs height.sigma_m...\n');
    hm2 = run_sweep_2d(cfg_base, 'P0_clkBias', vals3.^2, ...
                                 'height.sigma_m', vals4);
    plot_heatmap(hm2, 'P0\_clkBias std (m/s)', 'height.sigma\_m (m)', ...
        vals3, vals4, 'RMSE_h_m', 'RMSE Horizontal (m)');
end

%% ── Save ─────────────────────────────────────────────────────────────────
if saveResults
    if run_1d_sweeps
        save('sensitivity_1d.mat', 'sweep_results', 'sweeps', 'cfg_base');
        fprintf('\n  Saved sensitivity_1d.mat\n');
    end
end


%% ═══════════════════════════════════════════════════════════════════════════
%  SWEEP DEFINITIONS
%% ═══════════════════════════════════════════════════════════════════════════
function sweeps = define_sweeps()
% Each entry: param_path (dot-separated), values, label, units, tag
s = struct;
i = 0;

i=i+1;
s(i).param_path = 'sigma_pos_rw';
s(i).values     = logspace(-5,-2,12);
s(i).label      = 'sigma\_pos\_rw';
s(i).units      = 'm/sqrt(s)';
s(i).tag        = 'sigma_pos_rw';

i=i+1;
s(i).param_path = 'sigma_doppler_Hz';
s(i).values     = logspace(0,2,12);
s(i).label      = 'sigma\_doppler\_Hz';
s(i).units      = 'Hz';
s(i).tag        = 'sigma_doppler_Hz';

i=i+1;
s(i).param_path = 'sigma_dopplerRate_Hzps';
s(i).values     = logspace(-2,1,12);
s(i).label      = 'sigma\_dopplerRate\_Hzps';
s(i).units      = 'Hz/s';
s(i).tag        = 'sigma_dopplerRate_Hzps';

i=i+1;
s(i).param_path = 'meas.Rinflate';
s(i).values     = [0.5, 0.75, 1, 1.5, 2, 3, 5, 8, 12];
s(i).label      = 'meas.Rinflate';
s(i).units      = '(-)';
s(i).tag        = 'Rinflate';

i=i+1;
s(i).param_path = 'height.sigma_m';
s(i).values     = [1, 2, 5, 10, 20, 50, 100, 500, 1000];
s(i).label      = 'height.sigma\_m';
s(i).units      = 'm';
s(i).tag        = 'height_sigma';

i=i+1;
s(i).param_path = 'P0_clkBias';
s(i).values     = logspace(-4,0,10).^2;   % variance
s(i).label      = 'sqrt(P0\_clkBias)';
s(i).units      = 'm/s';
s(i).tag        = 'P0_clkBias';

i=i+1;
s(i).param_path = 'P0_clkDrift';
s(i).values     = logspace(-6,-2,10).^2;
s(i).label      = 'sqrt(P0\_clkDrift)';
s(i).units      = 'm/s^2';
s(i).tag        = 'P0_clkDrift';

i=i+1;
s(i).param_path = 'clk.filt.sigmaBias_rw_mps_sqrt';
s(i).values     = logspace(-4,-1,10);
s(i).label      = 'clk.filt.sigmaBias\_rw';
s(i).units      = 'm/s/sqrt(s)';
s(i).tag        = 'clk_bias_rw';

i=i+1;
s(i).param_path = 'clk.filt.sigmaDrift_rw_mps2_sqrt';
s(i).values     = logspace(-6,-2,10);
s(i).label      = 'clk.filt.sigmaDrift\_rw';
s(i).units      = 'm/s^2/sqrt(s)';
s(i).tag        = 'clk_drift_rw';

sweeps = s;
end

%% ═══════════════════════════════════════════════════════════════════════════
%  1D SWEEP RUNNER
%% ═══════════════════════════════════════════════════════════════════════════
function result = run_sweep_1d(cfg_base, param_path, values, nRepeat)
N = numel(values);
fields_to_extract = {'RMSE_h_m','RMSE_E_m','RMSE_N_m','finalHoriz_m', ...
                     'meanNIS','convTime_horiz_s','fractionRejectedByGate'};

nF = numel(fields_to_extract);
data = nan(N, nF);

for vi = 1:N
    vals_rep = nan(nRepeat, nF);
    for rr = 1:nRepeat
        cfg_i = set_cfg_field(cfg_base, param_path, values(vi));
        cfg_i.plot.enable = false;
        cfg_i.verbose     = false;
        cfg_i.seedMeas    = cfg_base.seedMeas + rr - 1;
        try
            out_i = real_strtg_sequential_EKF_fastrealistic(cfg_i);
            mi    = ekf_extract_metrics(out_i, cfg_i);
            for fi = 1:nF
                fn = fields_to_extract{fi};
                if isfield(mi, fn), vals_rep(rr,fi) = mi.(fn); end
            end
        catch
            % leave as NaN
        end
    end
    data(vi,:) = nanmean(vals_rep, 1);
end

result.values = values;
result.data   = data;
result.fields = fields_to_extract;
end

%% ═══════════════════════════════════════════════════════════════════════════
%  2D SWEEP RUNNER
%% ═══════════════════════════════════════════════════════════════════════════
function hm = run_sweep_2d(cfg_base, path1, vals1, path2, vals2)
N1 = numel(vals1);
N2 = numel(vals2);
rmse_grid    = nan(N1, N2);
nis_grid     = nan(N1, N2);

for i = 1:N1
    for j = 1:N2
        cfg_ij = set_cfg_field(cfg_base, path1, vals1(i));
        cfg_ij = set_cfg_field(cfg_ij,   path2, vals2(j));
        cfg_ij.plot.enable = false;
        cfg_ij.verbose     = false;
        try
            out_ij = real_strtg_sequential_EKF_fastrealistic(cfg_ij);
            mi = ekf_extract_metrics(out_ij, cfg_ij);
            rmse_grid(i,j) = mi.RMSE_h_m;
            nis_grid(i,j)  = mi.meanNIS;
        catch
            % leave as NaN
        end
    end
end

hm.path1 = path1;  hm.path2 = path2;
hm.vals1 = vals1;  hm.vals2 = vals2;
hm.rmse_grid = rmse_grid;
hm.nis_grid  = nis_grid;
end

%% ═══════════════════════════════════════════════════════════════════════════
%  PLOTTING
%% ═══════════════════════════════════════════════════════════════════════════
function plot_1d_sweeps(sweep_results, sweeps)
% One figure per sweep: RMSE_h and meanNIS on twin axes.
for si = 1:numel(sweeps)
    sw  = sweeps(si);
    res = sweep_results.(sw.tag);
    v   = res.values;
    fi_rmse = find(strcmp(res.fields, 'RMSE_h_m'),       1);
    fi_nis  = find(strcmp(res.fields, 'meanNIS'),         1);
    fi_conv = find(strcmp(res.fields, 'convTime_horiz_s'),1);

    rmse = res.data(:,fi_rmse);
    nis  = res.data(:,fi_nis);
    conv = res.data(:,fi_conv);

    figure('Name',['Sweep: ' sw.tag],'Position',[60 80 1000 420]);

    subplot(1,2,1);
    semilogx(v, rmse, 'b-o','LineWidth',2,'MarkerSize',5); hold on; grid on;
    xlabel([sw.label ' (' sw.units ')']);
    ylabel('RMSE Horizontal (m)');
    title([sw.label '  vs  RMSE_h']);

    subplot(1,2,2);
    semilogx(v, nis, 'r-s','LineWidth',2,'MarkerSize',5); hold on; grid on;
    yline(1.0,'g--','Ideal','LineWidth',1.5);
    xlabel([sw.label ' (' sw.units ')']);
    ylabel('Mean NIS');
    title([sw.label '  vs  Mean NIS']);

    sgtitle(['Sensitivity: ' strrep(sw.label,'\_','_')], ...
        'FontSize',12,'FontWeight','bold');
end
end

function plot_tornado(sweep_results, sweeps)
% Tornado plot: effect of sweeping each parameter from low to high.
nS = numel(sweeps);
effect_rmse = nan(nS,1);
label_str   = cell(nS,1);

for si = 1:nS
    sw  = sweeps(si);
    res = sweep_results.(sw.tag);
    fi  = find(strcmp(res.fields,'RMSE_h_m'),1);
    rmse = res.data(:,fi);
    valid = isfinite(rmse);
    if sum(valid) < 2, continue; end
    effect_rmse(si) = max(rmse(valid)) - min(rmse(valid));
    label_str{si}   = strrep(sw.label,'\_','_');
end

valid = isfinite(effect_rmse);
if sum(valid) == 0, return; end

[~,ord] = sort(effect_rmse(valid),'descend');
eff_v   = effect_rmse(valid);
lbl_v   = label_str(valid);

figure('Name','Tornado: RMSE_h sensitivity','Position',[60 80 700 500]);
barh(eff_v(ord),'FaceColor',[0.2 0.5 0.85],'FaceAlpha',0.8);
set(gca,'YTick',1:numel(ord),'YTickLabel',lbl_v(ord));
xlabel('\DeltaRMSE Horizontal (m)  [max - min over sweep]');
title('Sensitivity tornado: parameter influence on RMSE_h');
grid on;
end

function plot_heatmap(hm, lbl1, lbl2, v1, v2, field, fieldLabel)
if strcmp(field,'RMSE_h_m')
    grid_data = hm.rmse_grid;
else
    grid_data = hm.nis_grid;
end

figure('Name',['Heatmap: ' lbl1 ' vs ' lbl2],'Position',[60 80 700 520]);
imagesc(log10(v2), log10(v1), grid_data);
colorbar; colormap(jet);
xlabel(lbl2); ylabel(lbl1);
title(['Heatmap: ' fieldLabel '  (' lbl1 '  vs  ' lbl2 ')']);
set(gca,'YDir','normal');

% Overlay NIS=1 contour if plotting RMSE
if strcmp(field,'RMSE_h_m') && any(isfinite(hm.nis_grid(:)))
    hold on;
    contour(log10(v2), log10(v1), hm.nis_grid, [1 1], 'w-','LineWidth',2.5);
    text(0.5, 0.95,'NIS=1 contour','Units','normalized','Color','white', ...
        'HorizontalAlignment','center','FontSize',10);
end
end

%% ═══════════════════════════════════════════════════════════════════════════
%  HELPERS
%% ═══════════════════════════════════════════════════════════════════════════
function cfg = set_cfg_field(cfg, path, value)
% Set a nested cfg field by dot-separated path string.
parts = strsplit(path, '.');
switch numel(parts)
    case 1, cfg.(parts{1}) = value;
    case 2, cfg.(parts{1}).(parts{2}) = value;
    case 3, cfg.(parts{1}).(parts{2}).(parts{3}) = value;
    otherwise
        error('set_cfg_field: unsupported path depth: %s', path);
end
end

function v = nanmean(x, dim)
if nargin < 2, dim = 1; end
mask = ~isnan(x);
s = sum(x .* mask, dim);
n = sum(mask, dim);
v = s ./ max(n, 1);
v(n == 0) = NaN;
end

%% ═══════════════════════════════════════════════════════════════════════════
%  BASELINE CONFIGURATION
%% ═══════════════════════════════════════════════════════════════════════════
function cfg = baseline_cfg()
cfg.tleFile = fullfile('InputData','oneweb.tle');
cfg.duration_sec = 5;  cfg.dt = 0.01;
cfg.startUTC = [];  cfg.stopUTC = [];
cfg.startFromTLE = 'first';  cfg.maxTLEAgeDays = 7;
cfg.rx.lat_deg = 52.06869;  cfg.rx.lon_deg = -0.60884;
cfg.rx.alt_m = 100;  cfg.rx.vel_ecef = [0;0;0];  cfg.rx.acc_ecef = [0;0;0];
cfg.fc_Hz = 11.325e9;  cfg.c = 299792458;
cfg.minElevDeg = 25;  cfg.coarseStep = 10;  cfg.minMaxElev = 30;
cfg.maxSats = 12;  cfg.constellationSelectionWindow_s = 60;
cfg.useSGP4 = false;  cfg.allowJ2Fallback = true;
cfg.init.errENU_m = [500;500;0];  cfg.init.sigENU_m = [550;550;10];
cfg.P0_clkBias = (0.1)^2;  cfg.P0_clkDrift = (5e-5)^2;
cfg.sigma_pos_rw = 1e-4;
cfg.clk.truth.initBias_mps = 0.10;  cfg.clk.truth.initDrift_mps2 = 5e-5;
cfg.clk.truth.sigmaBias_rw_mps_sqrt = 5e-3;
cfg.clk.truth.sigmaDrift_rw_mps2_sqrt = 1.0e-4;
cfg.clk.filt.sigmaBias_rw_mps_sqrt = 1.0e-2;
cfg.clk.filt.sigmaDrift_rw_mps2_sqrt = 1.0e-4;
cfg.clk.filter_init_bias_mps = 0.0;  cfg.clk.filter_init_drift_mps2 = 0.0;
cfg.tropo.enable = true;  cfg.tropo.ztd_truth_m = 2.5;
cfg.tropo.ztd_model_m = 2.3;  cfg.tropo.el_floor_deg = 5;
cfg.tropo.map_epsilon_deg = 0.5;  cfg.tropo.model_truth = 'simple';
cfg.tropo.model_filter = 'simple';  cfg.tropo.lat_deg = 52.06869;
cfg.tropo.alt_m = 100;  cfg.tropo.doy = 180;
cfg.tropo.P0_hPa = 1013.25;  cfg.tropo.T0_K = 288.15;  cfg.tropo.e0_hPa = 11.0;
cfg.sigma_doppler_Hz = 5.0;  cfg.sigma_dopplerRate_Hzps = 0.05;
cfg.meas.enable_lighttime = true;  cfg.meas.colored_tau_s = 30;
cfg.meas.sigma_col0_mps = 0.03;  cfg.meas.Rinflate = 1.0;
cfg.meas.useDoppler = true;  cfg.meas.useDopplerRate = true;
cfg.meas.fdJacStep_m = 0.10;
cfg.filter.useNISGate = true;  cfg.filter.nisGate = 16.0;
cfg.filter.maxIEKFIter = 10;  cfg.filter.iekfTol = 1e-5;
cfg.filter.maxStep_m = 25;  cfg.filter.startupInflateEpochs = 20;
cfg.filter.startupRmult = 25;
cfg.height.enable = true;  cfg.height.sigma_m = 10.0;  cfg.height.z0_m = 0.0;
cfg.capture.mode = 'stickyBest';  cfg.capture.minDwell_s = 30;
cfg.capture.switchBlackout_s = 3;  cfg.capture.elev_margin_deg = 20;
cfg.capture.useServiceAngleFilter = false;  cfg.capture.maxServiceOffNadirDeg = 45;
cfg.capture.switchThreshold_lowElev = 10;  cfg.capture.switchThreshold_highElev = 15;
cfg.capture.highElevThreshold = 35;
cfg.seedMeas = 42;  cfg.plot.enable = false;  cfg.verbose = false;
cfg.analysis.conv.refWindow_s = 15;  cfg.analysis.conv.settleWindow_s = 10;
cfg.analysis.conv.relTol = 0.15;  cfg.analysis.conv.absTol_m = 5.0;
cfg.analysis.conv.maxCoeffVar = 0.25;  cfg.analysis.conv.maxSlope_mps = 1.0;
cfg.analysis.conv.errFloor_m = 5.0;
end
