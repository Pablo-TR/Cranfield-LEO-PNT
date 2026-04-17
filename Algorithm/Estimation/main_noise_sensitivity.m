%% MAIN_NOISE_SENSITIVITY
% Report entry point for measurement/process-noise sensitivity studies.
%
% This is the cleaned version of the repeated sweep logic in MAIN_sensitivity
% and ANALYSIS/ANALYSIS2. It keeps the same sweep values used in
% MAIN_sensitivity.m and centralises the run, aggregation, and plotting.

clear; close all; clc;

root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root_dir, 'functions')));
setup_report_paths(root_dir);

out_dir = fullfile(root_dir, 'outputs', 'noise_sensitivity');
ensure_dir(out_dir);

settings.n_repeat = 3;
settings.base_seed = 4200;
settings.run_2d_heatmaps = true;

cfg_base = cfg_report_baseline('sensitivity');
sweeps = define_noise_sweeps();

fprintf('Running 1D noise/tuning sensitivity sweeps...\n');
results_1d = struct();
for si = 1:numel(sweeps)
    sw = sweeps(si);
    fprintf('  [%d/%d] %s (%d values x %d repeats)\n', ...
        si, numel(sweeps), sw.label, numel(sw.values), settings.n_repeat);
    results_1d.(sw.tag) = run_sweep1d(cfg_base, sw.param_path, sw.values, ...
        settings.n_repeat, settings.base_seed + 100*si);
    results_1d.(sw.tag).label = sw.label;
    results_1d.(sw.tag).units = sw.units;
end

results_2d = struct();
if settings.run_2d_heatmaps
    fprintf('Running 2D sensitivity heatmaps...\n');
    results_2d.pos_rw_vs_doppler = run_sweep2d(cfg_base, ...
        'sigma_pos_rw', logspace(-5, -2, 8), ...
        'sigma_doppler_Hz', logspace(0, 2, 8), ...
        settings.base_seed + 1000);
    results_2d.p0_clock_vs_height = run_sweep2d(cfg_base, ...
        'P0_clkBias', logspace(-4, 0, 8).^2, ...
        'height.sigma_m', logspace(0, 2.5, 8), ...
        settings.base_seed + 2000);
end

plot_noise_sensitivity(results_1d, results_2d, sweeps, out_dir);
save(fullfile(out_dir, 'noise_sensitivity_results.mat'), ...
    'cfg_base', 'settings', 'sweeps', 'results_1d', 'results_2d');
fprintf('Saved noise sensitivity outputs to %s\n', out_dir);
