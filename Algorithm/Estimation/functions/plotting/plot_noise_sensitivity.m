function plot_noise_sensitivity(results_1d, results_2d, sweeps, out_dir)
%PLOT_NOISE_SENSITIVITY Save summary figures for sensitivity studies.

ensure_dir(out_dir);

fig = figure('Name', 'noise_sensitivity_1d', 'Color', 'w', 'Position', [80 80 1200 800]);
n = numel(sweeps);
for si = 1:n
    sw = sweeps(si);
    res = results_1d.(sw.tag);
    subplot(3, 3, si);
    y = res.mean(:, strcmp(res.fields, 'RMSE_h_m'));
    plot(res.values, y, 'b-o', 'LineWidth', 1.2, 'MarkerFaceColor', 'b');
    set(gca, 'XScale', axis_scale(res.values));
    grid on;
    xlabel(sprintf('%s [%s]', sw.label, sw.units), 'Interpreter', 'none');
    ylabel('RMSE_h [m]');
    title(sw.tag, 'Interpreter', 'none');
end
print(fig, fullfile(out_dir, 'noise_sensitivity_1d'), '-dpng', '-r150');

if ~isempty(fieldnames(results_2d))
    names = fieldnames(results_2d);
    for i = 1:numel(names)
        hm = results_2d.(names{i});
        fig_h = figure('Name', names{i}, 'Color', 'w', 'Position', [120 120 700 560]);
        imagesc(hm.values_x, hm.values_y, hm.RMSE_h_m);
        set(gca, 'YDir', 'normal', 'XScale', axis_scale(hm.values_x), ...
            'YScale', axis_scale(hm.values_y));
        colorbar;
        xlabel(hm.path_x, 'Interpreter', 'none');
        ylabel(hm.path_y, 'Interpreter', 'none');
        title(sprintf('%s: RMSE_h [m]', names{i}), 'Interpreter', 'none');
        print(fig_h, fullfile(out_dir, names{i}), '-dpng', '-r150');
    end
end
end

function scale = axis_scale(v)
if all(v > 0) && max(v) / min(v) > 20
    scale = 'log';
else
    scale = 'linear';
end
end
