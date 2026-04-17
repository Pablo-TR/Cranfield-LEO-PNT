function plot_multisat_summary(study, out_dir)
%PLOT_MULTISAT_SUMMARY Save the multichannel summary figure.

ensure_dir(out_dir);

fig = figure('Name', 'multisat_summary', 'Color', 'w', 'Position', [80 80 900 420]);

yyaxis left;
errorbar(study.channel_counts, study.mean_RMSE_h_m, study.std_RMSE_h_m, ...
    'b-o', 'LineWidth', 1.4, 'MarkerFaceColor', 'b');
ylabel('RMSE_h [m]');

yyaxis right;
plot(study.channel_counts, study.mean_NIS, 'r-s', 'LineWidth', 1.4, 'MarkerFaceColor', 'r');
ylabel('mean NIS');

grid on;
xlabel('Number of tracked channels / satellites');
title('Multichannel sensitivity: position error and consistency');
legend('RMSE_h mean +/- 1 sigma', 'mean NIS', 'Location', 'best');
print(fig, fullfile(out_dir, 'multisat_summary'), '-dpng', '-r150');
end
