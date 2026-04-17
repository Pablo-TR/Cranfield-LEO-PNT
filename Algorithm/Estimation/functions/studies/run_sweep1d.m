function result = run_sweep1d(cfg_base, param_path, values, n_repeat, base_seed)
%RUN_SWEEP1D Run a one-dimensional parameter sweep with repeated seeds.

fields = {'RMSE_h_m','RMSE_E_m','RMSE_N_m','finalHoriz_m', ...
    'meanNIS','convTime_horiz_s','fractionRejectedByGate'};

n_values = numel(values);
n_fields = numel(fields);
samples = nan(n_values, n_repeat, n_fields);

for vi = 1:n_values
    for ri = 1:n_repeat
        cfg = set_field(cfg_base, param_path, values(vi));
        cfg.seedMeas = base_seed + 100*vi + ri;
        try
            [~, m] = run_case(cfg);
            for fi = 1:n_fields
                if isfield(m, fields{fi})
                    samples(vi, ri, fi) = m.(fields{fi});
                end
            end
        catch ME
            warning('run_sweep1d:runFailed', '%s value %d repeat %d failed: %s', ...
                param_path, vi, ri, ME.message);
        end
    end
end

result.param_path = param_path;
result.values = values;
result.fields = fields;
result.samples = samples;
result.mean = squeeze(mean(samples, 2, 'omitnan'));
result.std = squeeze(std(samples, 0, 2, 'omitnan'));
end
