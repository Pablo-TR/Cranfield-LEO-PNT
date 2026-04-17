function metrics = extract_metrics(out, cfg)
%EXTRACT_METRICS Short-name wrapper around ekf_extract_metrics.
%
% Adds a defensive placeholder for Doppler-rate innovations when a scenario
% intentionally disables Doppler-rate measurements.

if ~isfield(out, 'nu_rate_hist') || isempty(out.nu_rate_hist)
    out.nu_rate_hist = nan(size(out.nu_hist));
end
metrics = ekf_extract_metrics(out, cfg);
end
