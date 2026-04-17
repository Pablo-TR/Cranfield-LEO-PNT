function cfg = set_field(cfg, dot_path, value)
%SET_FIELD Set a nested struct field using a dot-separated path.
%
% Example:
%   cfg = set_field(cfg, 'meas.Rinflate', 2.0);

parts = strsplit(dot_path, '.');
cfg = set_field_rec(cfg, parts, value);
end

function s = set_field_rec(s, parts, value)
if isscalar(parts)
    s.(parts{1}) = value;
else
    head = parts{1};
    if ~isfield(s, head) || ~isstruct(s.(head))
        s.(head) = struct();
    end
    s.(head) = set_field_rec(s.(head), parts(2:end), value);
end
end
