function plot_capture_and_measurements(t, trk, meas, satUsed, zPredDoppUsed, zPredRateUsed, sat, cfg)

tmin = t/60;

figure('Name','Capture summary','Position',[60 500 1200 700]);

subplot(2,2,1); hold on; grid on;
plot(tmin, meas.nVisible, 'k-', 'LineWidth',1.2);
plot(tmin, satUsed, 'r.', 'MarkerSize',8);
xlabel('Time (min)');
ylabel('Count / sat idx');
title('Visibility and used satellite');
legend('# visible','used sat','Location','best');

subplot(2,2,2); hold on; grid on;
dopp_meas = nan(size(t));
dopp_true = nan(size(t));
dopp_pred = nan(size(t));
for k = 1:numel(t)
    s = satUsed(k);
    if s <= 0
        continue;
    end
    dopp_meas(k) = -meas.dopp.obs(k,s)  / cfg.lambda;
    dopp_true(k) = -meas.dopp.true(k,s) / cfg.lambda;
    dopp_pred(k) = -zPredDoppUsed(k)    / cfg.lambda;
end
plot(tmin, dopp_meas, 'b.', 'MarkerSize',5);
plot(tmin, dopp_pred, 'r-', 'LineWidth',1.2);
plot(tmin, dopp_true, 'k--','LineWidth',1.0);
xlabel('Time (min)');
ylabel('Doppler (Hz)');
title('Used Doppler stream');
legend('measured','predicted','truth','Location','best');

subplot(2,2,3); hold on; grid on;
rate_meas = nan(size(t));
rate_true = nan(size(t));
rate_pred = nan(size(t));
for k = 1:numel(t)
    s = satUsed(k);
    if s <= 0
        continue;
    end
    rate_meas(k) = -meas.rate.obs(k,s)  / cfg.lambda;
    rate_true(k) = -meas.rate.true(k,s) / cfg.lambda;
    rate_pred(k) = -zPredRateUsed(k)    / cfg.lambda;
end
plot(tmin, rate_meas, 'b.', 'MarkerSize',5);
plot(tmin, rate_pred, 'r-', 'LineWidth',1.2);
plot(tmin, rate_true, 'k--','LineWidth',1.0);
xlabel('Time (min)');
ylabel('Doppler rate (Hz/s)');
title('Used Doppler-rate stream');
legend('measured','predicted','truth','Location','best');

subplot(2,2,4); hold on; grid on;
sid = trk.sid;
sid(sid==0) = NaN;
plot(tmin, sid, '.', 'MarkerSize',6);
if any(trk.blackout)
    plot(tmin(trk.blackout), zeros(sum(trk.blackout),1), 'rx', 'MarkerSize',6);
end
xlabel('Time (min)');
ylabel('sat idx');
title('Single-channel capture timeline');

figure('Name','Skyview','Position',[1200 120 500 500]);
hold on; axis equal; grid on;
th = linspace(0,2*pi,400);
for elR = [0 30 60]
    rho = 90 - elR;
    plot(rho*sin(th), rho*cos(th), 'k:');
end
text( 0, 94,'N','HorizontalAlignment','center','FontWeight','bold');
text(94,  0,'E','HorizontalAlignment','center','FontWeight','bold');
text( 0,-94,'S','HorizontalAlignment','center','FontWeight','bold');
text(-94, 0,'W','HorizontalAlignment','center','FontWeight','bold');

for k = find(satUsed>0).'
    s  = satUsed(k);
    az = sat.az(k,s);
    el = sat.elev(k,s);
    rho = 90 - el;
    plot(rho*sind(az), rho*cosd(az), 'r.', 'MarkerSize',4);
end
xlim([-95 95]);
ylim([-95 95]);
title('Skyview of used satellite');
end

%% ═══════════════════════ GEOMETRY HELPERS ════════════════════════════════
