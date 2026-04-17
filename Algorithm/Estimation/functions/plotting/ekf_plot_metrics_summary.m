function ekf_plot_metrics_summary(m, cfg) %#ok<INUSD>
% EKF_PLOT_METRICS_SUMMARY   Print a formatted metrics summary to console.

fprintf('\n╔══════════════════════════════════════════════════════╗\n');
fprintf('║              METRICS SUMMARY                        ║\n');
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  POSITION ERRORS                                    ║\n');
fprintf('║  Final East          : %+8.2f m                   ║\n', m.finalEast_m);
fprintf('║  Final North         : %+8.2f m                   ║\n', m.finalNorth_m);
fprintf('║  Final Horizontal    : %8.2f m                   ║\n', m.finalHoriz_m);
fprintf('║  Min  Horizontal     : %8.2f m                   ║\n', m.minHoriz_m);
fprintf('║  RMSE East           : %8.2f m                   ║\n', m.RMSE_E_m);
fprintf('║  RMSE North          : %8.2f m                   ║\n', m.RMSE_N_m);
fprintf('║  RMSE Horizontal     : %8.2f m                   ║\n', m.RMSE_h_m);
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  CLOCK                                              ║\n');
fprintf('║  Bias  RMSE          : %.4e m/s                ║\n', m.clkBiasRMSE_mps);
fprintf('║  Drift RMSE          : %.4e m/s^2              ║\n', m.clkDriftRMSE_mps2);
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  FILTER CONSISTENCY                                 ║\n');
fprintf('║  Mean  NIS           : %8.4f  (ideal=1.0)      ║\n', m.meanNIS);
fprintf('║  Median NIS          : %8.4f                   ║\n', m.medianNIS);
fprintf('║  Std NIS             : %8.4f                   ║\n', m.stdNIS);
fprintf('║  Fraction gated      : %8.2f%%                  ║\n', 100*m.fractionRejectedByGate);
fprintf('║  ACF frac sig (Dopp) : %8.2f%%                  ║\n', 100*m.fracSigACF_dopp);
fprintf('║  Excess kurtosis     : %8.4f  (Gauss→ 0)      ║\n', m.excessKurtosis_dopp);
fprintf('╠══════════════════════════════════════════════════════╣\n');
fprintf('║  CONVERGENCE                                        ║\n');
fprintf('║  Horiz settling time : %8.2f s                  ║\n', m.convTime_horiz_s);
fprintf('║  East  settling time : %8.2f s                  ║\n', m.convTime_east_s);
fprintf('║  North settling time : %8.2f s                  ║\n', m.convTime_north_s);
fprintf('╠══════════════════════════════════════════════════════╣\n');
if m.success
    fprintf('║  STATUS: SUCCESS                                    ║\n');
else
    fprintf('║  STATUS: FAILED — %s\n', m.failMsg);
end
fprintf('╚══════════════════════════════════════════════════════╝\n\n');
end
