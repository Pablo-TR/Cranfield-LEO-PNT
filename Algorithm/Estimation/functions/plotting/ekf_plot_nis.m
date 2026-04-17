function ekf_plot_nis(out, cfg)
% EKF_PLOT_NIS   Three-panel NIS diagnostic (time-series / histogram / running mean).
if ~cfg.plot.enable, return; end

nis = out.nis(:);
t   = out.t(:);

valid = ~isnan(nis) & nis > 0;
nis_v = nis(valid);
t_v   = t(valid);

if isempty(nis_v)
    fprintf('  [NIS] No valid samples.\n');
    return;
end

n        = numel(nis_v);
mean_n   = mean(nis_v);
std_n    = std(nis_v);
median_n = median(nis_v);

ci_half  = 1.96 * sqrt(2/n);
ci_lo    = max(0, 1 - ci_half);
ci_hi    = 1 + ci_half;

frac_gated = mean(nis_v > cfg.filter.nisGate);

fprintf('\n╔════════════════════════════════════════════╗\n');
fprintf('║           NIS DIAGNOSTIC REPORT           ║\n');
fprintf('╠════════════════════════════════════════════╣\n');
fprintf('║  Samples:          %-6d                 ║\n', n);
fprintf('║  Mean NIS:         %-8.4f               ║\n', mean_n);
fprintf('║  Median NIS:       %-8.4f               ║\n', median_n);
fprintf('║  Std  NIS:         %-8.4f               ║\n', std_n);
fprintf('║  95%% CI for mean:  [%.4f, %.4f]       ║\n', ci_lo, ci_hi);
fprintf('║  Fraction gated:   %-6.1f%%               ║\n', 100*frac_gated);
fprintf('╠════════════════════════════════════════════╣\n');
if mean_n > ci_hi
    sug = sqrt(mean_n);
    fprintf('║  STATUS: OVERCONFIDENT                    ║\n');
    fprintf('║  → Increase Q or inflate R                ║\n');
    fprintf('║  → Try sigma_pos_rw  × %.2f               ║\n', sug);
    fprintf('║  → Try Rinflate      × %.2f               ║\n', sug);
elseif mean_n < ci_lo
    sug = sqrt(mean_n);
    fprintf('║  STATUS: UNDERCONFIDENT                   ║\n');
    fprintf('║  → Decrease Q or deflate R                ║\n');
    fprintf('║  → Try sigma_pos_rw  × %.2f               ║\n', sug);
else
    fprintf('║  STATUS: WELL-TUNED  (mean NIS in CI)     ║\n');
end
fprintf('╚════════════════════════════════════════════╝\n\n');

tmin = t_v / 60;
figure('Name','NIS Diagnostics','Position',[100 100 1350 450]);

% Panel 1: NIS time-series
subplot(1,3,1);
plot(tmin, nis_v, 'b.', 'MarkerSize',4); hold on; grid on;
fill([tmin; flipud(tmin)], [ci_lo*ones(size(t_v)); ci_hi*ones(size(t_v))], ...
    [0.6 1.0 0.6], 'EdgeColor','none','FaceAlpha',0.25);
yline(1.0,              'g-', 'Ideal=1',                   'LineWidth',2.0);
yline(mean_n,           'k-', sprintf('Mean=%.2f',mean_n), 'LineWidth',1.5);
yline(cfg.filter.nisGate,'r--','Gate',                     'LineWidth',1.5);
xlabel('Time (min)'); ylabel('NIS'); title('NIS time-series');
ylim([0, min(max(nis_v)*1.1+0.5, cfg.filter.nisGate*1.5)]);
legend('NIS','95% CI','Ideal','Mean','Gate','Location','northeast');

% Panel 2: NIS histogram vs chi²(1) PDF
subplot(1,3,2);
clip_val = min(20, prctile_manual(nis_v, 99));
nis_clip = nis_v(nis_v <= clip_val);
edges    = linspace(0, clip_val, 40);
histogram(nis_clip, edges, 'Normalization','pdf', ...
    'FaceColor',[0.25 0.45 0.85], 'EdgeColor','white','FaceAlpha',0.75);
hold on; grid on;
x_chi = linspace(0.01, clip_val, 400);
y_chi = exp(-x_chi/2) ./ sqrt(2*pi*x_chi);
plot(x_chi, y_chi, 'r-', 'LineWidth',2.5);
xline(mean_n, 'k-',  sprintf('Mean=%.2f',mean_n), 'LineWidth',1.5);
xline(ci_lo,  'm--', 'CI lo', 'LineWidth',1.0);
xline(ci_hi,  'm--', 'CI hi', 'LineWidth',1.0);
xlabel('NIS value'); ylabel('PDF');
title(sprintf('NIS distribution  (n=%d)', n));
legend('NIS','\chi^2(1)','Mean','95% CI','Location','northeast');

% Panel 3: running mean
subplot(1,3,3);
running_mean = cumsum(nis_v) ./ (1:n)';
k_seq = (1:n)';
ci_r  = 1.96 * sqrt(2 ./ k_seq);
lo_r  = max(0, 1 - ci_r);
hi_r  = 1 + ci_r;
fill([tmin; flipud(tmin)], [lo_r; flipud(hi_r)], ...
    [0.6 1.0 0.6], 'EdgeColor','none','FaceAlpha',0.25);
hold on; grid on;
plot(tmin, running_mean, 'b-', 'LineWidth',2.0);
yline(1.0, 'g-', 'Target=1', 'LineWidth',2.0);
xline(tmin(end), 'k:');
xlabel('Time (min)'); ylabel('Cumulative mean NIS');
title({'Running mean NIS','(converge to 1.0)'});
legend('95% CI','Running mean','Target','Location','northeast');

sgtitle('EKF Tuning Diagnostics: NIS Analysis','FontSize',13,'FontWeight','bold');
end

% ═══════════════════════════════════════════════════════════════════════════

function p = prctile_manual(x, pct)
% Toolbox-free percentile.
x = sort(x(:));
n = numel(x);
if n == 0, p = NaN; return; end
idx = max(1, min(n, round(pct/100 * n)));
p   = x(idx);
end
