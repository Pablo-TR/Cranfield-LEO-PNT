function [out, metrics] = run_case(cfg)
%RUN_CASE Execute one IEKF run and extract report metrics.

out = leo_iekf(cfg);
metrics = extract_metrics(out, cfg);
end
