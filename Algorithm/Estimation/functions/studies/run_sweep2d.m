function hm = run_sweep2d(cfg_base, path_x, values_x, path_y, values_y, base_seed)
%RUN_SWEEP2D Run a two-parameter heatmap using one seed per grid cell.

nx = numel(values_x);
ny = numel(values_y);
hm.path_x = path_x;
hm.path_y = path_y;
hm.values_x = values_x;
hm.values_y = values_y;
hm.RMSE_h_m = nan(ny, nx);
hm.meanNIS = nan(ny, nx);
hm.convTime_horiz_s = nan(ny, nx);

for ix = 1:nx
    for iy = 1:ny
        cfg = set_field(cfg_base, path_x, values_x(ix));
        cfg = set_field(cfg, path_y, values_y(iy));
        cfg.seedMeas = base_seed + 100*ix + iy;
        try
            [~, m] = run_case(cfg);
            hm.RMSE_h_m(iy, ix) = m.RMSE_h_m;
            hm.meanNIS(iy, ix) = m.meanNIS;
            hm.convTime_horiz_s(iy, ix) = m.convTime_horiz_s;
        catch ME
            warning('run_sweep2d:runFailed', 'Grid (%d,%d) failed: %s', ix, iy, ME.message);
        end
    end
end
end
