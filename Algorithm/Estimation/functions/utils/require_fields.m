function require_fields(s, names, prefix)
for i = 1:numel(names)
    if ~isfield(s, names{i})
        error('Missing %s%s', prefix, names{i});
    end
end
end

%% ═══════════════════════ CONSTELLATION BUILDER ═══════════════════════════
