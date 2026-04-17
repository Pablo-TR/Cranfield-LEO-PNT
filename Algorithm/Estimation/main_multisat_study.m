%% MAIN_MULTISAT_STUDY
% Report entry point for multisatellite / multichannel LEO Doppler studies.
%
% The study follows the multichannel logic developed in ANALYSIS2: one
% channel uses stickyBest tracking; multiple channels use allVisible with a
% capped number of satellites.

clear; close all; clc;

root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root_dir, 'functions')));
setup_report_paths(root_dir);

out_dir = fullfile(root_dir, 'outputs', 'multisat_study');
ensure_dir(out_dir);

settings.channel_counts = [1 2 3 4];
settings.n_mc = 20;
settings.base_seed = 6100;
settings.duration_sec = 100;
settings.dt = 0.01;
settings.sigma_doppler_Hz = 50.0;
settings.sigma_dopplerRate_Hzps = 1.0;

cfg_base = cfg_report_baseline('multisat');
study = run_channel_sweep(cfg_base, settings);

plot_multisat_summary(study, out_dir);
%save(fullfile(out_dir, 'multisat_study_results.mat'), 'cfg_base', 'settings', 'study');
fprintf('Saved multisat study outputs to %s\n', out_dir);
