function ekf_plot_innovation_whiteness(out, cfg)
% EKF_PLOT_INNOVATION_WHITENESS  ACF / histogram for Doppler and Doppler-rate.
if ~cfg.plot.enable, return; end

plot_one_channel(out.nu_hist(:),      out.t(:), 'Doppler innovation');
plot_one_channel(out.nu_rate_hist(:), out.t(:), 'Doppler-rate innovation');
end


function plot_one_channel(nu_raw, t, label)
valid = ~isnan(nu_raw) & nu_raw ~= 0;
nu_v  = nu_raw(valid);
t_v   = t(valid);

if numel(nu_v) < 20
    fprintf('  [Whiteness] %s: only %d valid samples. Skipping.\n', label, numel(nu_v));
    return;
end

sigma  = std(nu_v);
if sigma < eps, return; end

nu_norm = (nu_v - mean(nu_v)) / sigma;
maxLag  = min(50, floor(numel(nu_norm)/4));
[acf, lags] = acf_manual(nu_norm, maxLag);
conf_bound  = 1.96 / sqrt(numel(nu_norm));
num_sig     = sum(abs(acf(lags ~= 0)) > conf_bound);
frac_sig    = num_sig / sum(lags ~= 0);

figure('Name', label, 'Position',[100 600 1300 420]);

subplot(1,3,1);
plot(t_v/60, nu_norm, 'b.', 'MarkerSize',3); hold on; grid on;
yline( 2,'r--','+2\sigma','LineWidth',1.2);
yline(-2,'r--','-2\sigma','LineWidth',1.2);
yline( 0,'k-',            'LineWidth',0.8);
xlabel('Time (min)'); ylabel('Normalised innovation'); title('Innovation sequence');
ylim([-5 5]);

subplot(1,3,2);
stem(lags, acf, 'filled', 'Color',[0.2 0.4 0.8], 'MarkerSize',3); hold on; grid on;
yline( conf_bound,'r--','+95% bound','LineWidth',1.5);
yline(-conf_bound,'r--','-95% bound','LineWidth',1.5);
yline(0,'k-','LineWidth',0.8);
xlabel('Lag (epochs)'); ylabel('Autocorrelation');
title({'Innovation ACF','(target: all lags ≈ 0)'});

subplot(1,3,3);
edges = linspace(-4,4,35);
histogram(nu_norm, edges, 'Normalization','pdf', ...
    'FaceColor',[0.25 0.45 0.85], 'EdgeColor','white', 'FaceAlpha',0.75);
hold on; grid on;
x_g = linspace(-4,4,300);
y_g = exp(-0.5*x_g.^2) / sqrt(2*pi);
plot(x_g, y_g, 'r-', 'LineWidth',2.5);
xlabel('Normalised innovation'); ylabel('PDF');
title('vs N(0,1)');
legend('Innovations','N(0,1)','Location','best');

if frac_sig > 0.10
    status_str = sprintf('COLORED: %.0f%% lags significant', 100*frac_sig);
    col = [0.8 0.1 0.1];
else
    status_str = sprintf('WHITE: %.0f%% lags significant',   100*frac_sig);
    col = [0.1 0.6 0.1];
end
sgtitle(sprintf('%s  —  %s', label, status_str), ...
    'FontSize',12,'FontWeight','bold','Color',col);
end

% ═══════════════════════════════════════════════════════════════════════════
