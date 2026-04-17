function ekf_plot_position_errors(out, cfg)
% EKF_PLOT_POSITION_ERRORS  ENU error with 2-sigma bounds + horizontal error.
if ~cfg.plot.enable, return; end

t    = out.t_prior(:) / 60;   % minutes
K    = numel(t);
eENU = out.errENU_prior;      % 3 x (K)

% Extract sigma from filter covariance
sigENU = nan(3,K);
for k = 1:K
    Pk = out.Pfull_prior{k};
    if isempty(Pk), continue; end
    for d = 1:3
        sigENU(d,k) = sqrt(max(Pk(d,d),0));
    end
end

labels = {'East','North','Up'};
cols   = {[0.1 0.4 0.9],[0.1 0.7 0.2],[0.85 0.33 0.1]};

figure('Name','ENU position errors','Position',[60 80 1350 400]);
for d = 1:3
    subplot(1,3,d); hold on; grid on;
    fill([t; flipud(t)], ...
         [2*sigENU(d,:)'; flipud(-2*sigENU(d,:)')], ...
         [1 0.7 0.7],'EdgeColor','none','FaceAlpha',0.35);
    plot(t, eENU(d,:)','-','Color',cols{d},'LineWidth',1.5);
    plot(t,  2*sigENU(d,:)','r--','LineWidth',0.9);
    plot(t, -2*sigENU(d,:)','r--','LineWidth',0.9);
    yline(0,'k:','LineWidth',0.8);
    xlabel('Time (min)'); ylabel('Error (m)'); title([labels{d} ' error']);
    legend('2\sigma band','error','+2\sigma','-2\sigma','Location','best');
end
sgtitle('ENU position error with filter uncertainty','FontSize',12,'FontWeight','bold');

% Horizontal error
eH = sqrt(eENU(1,:).^2 + eENU(2,:).^2);
figure('Name','Horizontal error','Position',[60 500 700 350]);
plot(t, eH,'b-','LineWidth',1.5); hold on; grid on;
xlabel('Time (min)'); ylabel('Horizontal error (m)');
title('Horizontal positioning error  \sqrt{E^2+N^2}');
yline(50,'r--','50 m','LineWidth',1.0);
yline(100,'k--','100 m','LineWidth',1.0);
end

% ═══════════════════════════════════════════════════════════════════════════
