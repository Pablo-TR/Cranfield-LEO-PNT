function plot_enu_fastrealistic(t, errENU, Pfull, enuIdx)
tmin = t;
K    = numel(t);

sigENU = nan(3,K);
for k = 1:K
    Pk = Pfull{k};
    if isempty(Pk)
        continue;
    end
    for d = 1:3
        sigENU(d,k) = sqrt(max(Pk(enuIdx(d), enuIdx(d)), 0));
    end
end

labels = {'East','North','Up'};
cols   = {[0.1 0.4 0.9], [0.1 0.7 0.2], [0.85 0.33 0.1]};

for d = 1:3
    figure('Name',[labels{d} ' error'], 'Position', [60+430*(d-1), 120, 410, 330]);
    hold on; grid on;
    fill([tmin; flipud(tmin)], ...
         [2*sigENU(d,:), fliplr(-2*sigENU(d,:))], ...
         [1 0.7 0.7], 'EdgeColor','none','FaceAlpha',0.35);
    plot(tmin, errENU(d,:), '-', 'Color', cols{d}, 'LineWidth', 1.5);
    plot(tmin,  2*sigENU(d,:), 'r--', 'LineWidth', 0.9);
    plot(tmin, -2*sigENU(d,:), 'r--', 'LineWidth', 0.9);
    yline(0,'k:','LineWidth',0.8);
    xlabel('Time (sec)');
    ylabel('Error (m)');
    title([labels{d} ' error']);
    legend('2\sigma band','error','2\sigma','-2\sigma','Location','best');
end
end
