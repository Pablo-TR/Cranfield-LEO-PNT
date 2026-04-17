%% MAIN_BASELINE
% Report entry point for one adopted LEO Doppler IEKF baseline run.
%
% Outputs:
%   outputs/baseline/baseline_results.mat
%   outputs/baseline/baseline_errors.png

clear; close all; clc;

root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root_dir, 'functions')));
setup_report_paths(root_dir);

out_dir = fullfile(root_dir, 'outputs', 'baseline');
ensure_dir(out_dir);

cfg = cfg_report_baseline('phase0');
cfg.plot.enable = false;
cfg.verbose = false;

fprintf('Running report baseline IEKF...\n');
[out, metrics] = run_case(cfg);
print_metric_summary(metrics);
plot_baseline_report(out, metrics, out_dir);

%save(fullfile(out_dir, 'baseline_results.mat'), 'cfg', 'out', 'metrics');
fprintf('Saved baseline outputs to %s\n', out_dir);
