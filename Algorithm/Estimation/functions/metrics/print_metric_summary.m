function print_metric_summary(m)
%PRINT_METRIC_SUMMARY Compact console summary for one solver run.

fprintf('\nBaseline metrics\n');
fprintf('  RMSE East/North/Horiz : %.2f / %.2f / %.2f m\n', ...
    m.RMSE_E_m, m.RMSE_N_m, m.RMSE_h_m);
fprintf('  Final East/North/Horiz: %.2f / %.2f / %.2f m\n', ...
    m.finalEast_m, m.finalNorth_m, m.finalHoriz_m);
fprintf('  mean NIS              : %.3f\n', m.meanNIS);
fprintf('  ACF significant lags  : Doppler %.1f%%, Rate %.1f%%\n', ...
    100*m.fracSigACF_dopp, 100*m.fracSigACF_rate);
fprintf('  Excess kurtosis       : Doppler %.2f, Rate %.2f\n', ...
    m.excessKurtosis_dopp, m.excessKurtosis_rate);
fprintf('  Convergence time      : %.2f s\n\n', m.convTime_horiz_s);
end
