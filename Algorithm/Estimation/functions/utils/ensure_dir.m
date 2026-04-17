function ensure_dir(path_name)
%ENSURE_DIR Create a directory if it does not already exist.

if ~exist(path_name, 'dir')
    mkdir(path_name);
end
end
