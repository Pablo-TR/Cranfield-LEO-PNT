%% MAIN_PHASE0_VALIDATION
% Report entry point for the full Phase 0 validation and tuning workflow.
%
% The detailed B1-B6 validation script is preserved in legacy_reference to
% avoid changing numerical behaviour late in the report process. This wrapper
% provides the clean report entry point, paths, and package-local execution.

clear; close all; clc;

root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root_dir, 'functions')));
setup_report_paths(root_dir);

run_phase0_validation(root_dir);
