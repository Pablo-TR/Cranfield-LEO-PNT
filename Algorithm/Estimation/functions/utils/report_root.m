function root_dir = report_root()
%REPORT_ROOT Return the absolute path to the report-code package root.

this_file = mfilename('fullpath');
root_dir = fileparts(fileparts(fileparts(this_file)));
end
