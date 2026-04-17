function ekf_plot_clock(out, cfg)
% EKF_PLOT_CLOCK  Clock bias and drift error + running RMSE.
if ~cfg.plot.enable, return; end

t    = out.t(:);
tmin = t / 60;

cBe = out.xhist(4,:).';
cDe = out.xhist(5,:).';
cBt = out.truthAux.clkBiasTrue(:);
cDt = out.truthAux.clkDriftTrue(:);

Lb = min(numel(cBe), numel(cBt));
Ld = min(numel(cDe), numel(cDt));

vB = ~isnan(cBe(1:Lb)) & ~isnan(cBt(1:Lb));
vD = ~isnan(cDe(1:Ld)) & ~isnan(cDt(1:Ld));

biasErr  = cBe(1:Lb) - cBt(1:Lb);
driftErr = cDe(1:Ld) - cDt(1:Ld);

% Running RMSE
biasRMSE_run  = nan(Lb,1);
driftRMSE_run = nan(Ld,1);
bIdx = find(vB);
for i = 1:numel(bIdx)
    k = bIdx(i);
    e = biasErr(vB & ((1:Lb)' <= k));
    biasRMSE_run(k) = sqrt(mean(e.^2));
end
dIdx = find(vD);
for i = 1:numel(dIdx)
    k = dIdx(i);
    e = driftErr(vD & ((1:Ld)' <= k));
    driftRMSE_run(k) = sqrt(mean(e.^2));
end

finalB = sqrt(mean(biasErr(vB).^2));
finalD = sqrt(mean(driftErr(vD).^2));

fprintf('\n=== CLOCK RMSE ===\n');
fprintf('  Final clock bias RMSE : %.6f m/s\n',   finalB);
fprintf('  Final clock drift RMSE: %.6e m/s^2\n', finalD);

figure('Name','Clock bias/drift RMSE','Position',[120 120 1200 700]);

subplot(2,2,1);
plot(tmin(1:Lb), biasErr,'b-','LineWidth',1.2); hold on; grid on;
yline(0,'k:','LineWidth',0.8); xlabel('Time (min)'); ylabel('Bias error (m/s)');
title('Clock bias error');

subplot(2,2,2);
plot(tmin(1:Ld), driftErr,'r-','LineWidth',1.2); hold on; grid on;
yline(0,'k:','LineWidth',0.8); xlabel('Time (min)'); ylabel('Drift error (m/s^2)');
title('Clock drift error');

subplot(2,2,3);
plot(tmin(1:Lb), biasRMSE_run,'b-','LineWidth',1.8); hold on; grid on;
yline(finalB,'k--',sprintf('Final=%.4g',finalB),'LineWidth',1.0);
xlabel('Time (min)'); ylabel('Running RMSE (m/s)'); title('Clock bias running RMSE');

subplot(2,2,4);
plot(tmin(1:Ld), driftRMSE_run,'r-','LineWidth',1.8); hold on; grid on;
yline(finalD,'k--',sprintf('Final=%.4g',finalD),'LineWidth',1.0);
xlabel('Time (min)'); ylabel('Running RMSE (m/s^2)'); title('Clock drift running RMSE');

sgtitle('Clock estimation performance','FontSize',12,'FontWeight','bold');
end

% ═══════════════════════════════════════════════════════════════════════════
