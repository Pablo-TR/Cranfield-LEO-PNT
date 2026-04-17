function ekf_plot_measurements(out, cfg)
% EKF_PLOT_MEASUREMENTS  Doppler, Doppler-rate, tracking timeline, skyview.
if ~cfg.plot.enable, return; end

t    = out.t(:);
tmin = t / 60;
trk  = out.trk;
meas = out.meas;
satU = out.satUsed(:);
sat  = out.sat;

figure('Name','Measurement diagnostics','Position',[60 500 1300 700]);

% Visibility / tracking
subplot(2,2,1); hold on; grid on;
plot(tmin, meas.nVisible,'k-','LineWidth',1.2);
sid_plot = satU;  sid_plot(sid_plot==0) = NaN;
plot(tmin, sid_plot,'r.','MarkerSize',6);
xlabel('Time (min)'); ylabel('Count / sat idx');
title('Visibility and tracked satellite');
legend('# visible','tracked','Location','best');

% Doppler
subplot(2,2,2); hold on; grid on;
dM = nan(size(t));  dT = nan(size(t));  dP = nan(size(t));
for k = 1:numel(t)
    s = satU(k); if s <= 0, continue; end
    dM(k) = -meas.dopp.obs(k,s)   / cfg.lambda;
    dT(k) = -meas.dopp.true(k,s)  / cfg.lambda;
    dP(k) = -out.zPredDoppUsed(k) / cfg.lambda;
end
plot(tmin, dM,'b.','MarkerSize',4);
plot(tmin, dP,'r-','LineWidth',1.2);
plot(tmin, dT,'k--','LineWidth',1.0);
xlabel('Time (min)'); ylabel('Doppler (Hz)');
title('Doppler stream (used satellite)');
legend('measured','predicted','truth','Location','best');

% Doppler-rate
subplot(2,2,3); hold on; grid on;
rM = nan(size(t)); rT = nan(size(t)); rP = nan(size(t));
for k = 1:numel(t)
    s = satU(k); if s <= 0, continue; end
    rM(k) = -meas.rate.obs(k,s)   / cfg.lambda;
    rT(k) = -meas.rate.true(k,s)  / cfg.lambda;
    rP(k) = -out.zPredRateUsed(k) / cfg.lambda;
end
plot(tmin, rM,'b.','MarkerSize',4);
plot(tmin, rP,'r-','LineWidth',1.2);
plot(tmin, rT,'k--','LineWidth',1.0);
xlabel('Time (min)'); ylabel('Doppler rate (Hz/s)');
title('Doppler-rate stream (used satellite)');
legend('measured','predicted','truth','Location','best');

% Tracking timeline
subplot(2,2,4); hold on; grid on;
plot(tmin, sid_plot,'.','MarkerSize',6);
if any(trk.blackout)
    plot(tmin(trk.blackout), zeros(sum(trk.blackout),1),'rx','MarkerSize',6);
end
xlabel('Time (min)'); ylabel('Satellite index');
title('Capture timeline');

% Skyview
figure('Name','Skyview','Position',[1200 100 500 500]);
hold on; axis equal; grid on;
th = linspace(0,2*pi,400);
for elR = [0 30 60 80]
    rho = 90 - elR;
    plot(rho*sin(th), rho*cos(th),'k:','LineWidth',0.6);
end
text(0, 94,'N','HorizontalAlignment','center','FontWeight','bold');
text(94, 0,'E','HorizontalAlignment','center','FontWeight','bold');
text(0,-94,'S','HorizontalAlignment','center','FontWeight','bold');
text(-94,0,'W','HorizontalAlignment','center','FontWeight','bold');
for k = find(satU > 0).'
    s  = satU(k); az = sat.az(k,s); el = sat.elev(k,s);
    rho = 90 - el;
    plot(rho*sind(az), rho*cosd(az),'r.','MarkerSize',3);
end
xlim([-95 95]); ylim([-95 95]);
title('Skyview of tracked satellite');
end

% ═══════════════════════════════════════════════════════════════════════════
