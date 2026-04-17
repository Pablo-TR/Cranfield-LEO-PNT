function setup_report_paths(root_dir)
%SETUP_REPORT_PATHS Add package functions and switch to the package root.
%
% Running from the report-code root keeps relative outputs and the copied
% InputData folder self-contained.

if nargin < 1 || isempty(root_dir)
    root_dir = report_root();
end

addpath(genpath(fullfile(root_dir, 'functions')));
cd(root_dir);
end
