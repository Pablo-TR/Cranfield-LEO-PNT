function plot_baseline_report(out, metrics, out_dir)
%PLOT_BASELINE_REPORT Save a compact baseline diagnostic figure.

ensure_dir(out_dir);
t = out.t(:);
eE = out.errENU(1,:).';
eN = out.errENU(2,:).';
eH = sqrt(eE.^2 + eN.^2);

fig = figure('Name', 'baseline_errors', 'Color', 'w', 'Position', [80 80 1000 650]);

subplot(2,1,1);
plot(t, eE, 'r-', 'LineWidth', 1.1); hold on;
plot(t, eN, 'g-', 'LineWidth', 1.1);
plot(t, eH, 'b-', 'LineWidth', 1.4);
yline(50, 'k--', '50 m');
grid on;
xlabel('Time [s]');
ylabel('Error [m]');
legend('East', 'North', '2D horizontal', 'Location', 'best');
title(sprintf('Baseline position error, RMSE_h = %.1f m', metrics.RMSE_h_m));

subplot(2,1,2);
plot(t, out.nis(:), 'k.', 'MarkerSize', 3); hold on;
yline(1.0, 'r--', 'NIS = 1');
grid on;
xlabel('Time [s]');
ylabel('NIS');
title(sprintf('Consistency, mean NIS = %.3f', metrics.meanNIS));

print(fig, fullfile(out_dir, 'baseline_errors'), '-dpng', '-r150');
end
