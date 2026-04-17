function run_phase0_validation(root_dir)
%RUN_PHASE0_VALIDATION Execute the preserved Phase 0 B1-B6 validation script.
%
% The Phase 0 script contains many report-specific plot layouts and console
% tables. It is kept as a single preserved workflow to avoid numerical drift.

if nargin < 1 || isempty(root_dir)
    root_dir = report_root();
end

script_path = fullfile(root_dir, 'legacy_reference', 'MAIN_phase0_validation.m');
if ~exist(script_path, 'file')
    error('Missing Phase 0 validation script: %s', script_path);
end

old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(root_dir);
run(script_path);
end
