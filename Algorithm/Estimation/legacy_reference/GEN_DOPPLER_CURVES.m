%% ═══════════════════════════════════════════════════════════════════════════
%  GEN_DOPPLER_CURVES.m
%
%  PURPOSE:
%    Stand-alone script that, given a precise recording timestamp and
%    duration, propagates the OneWeb constellation and computes the
%    clean theoretical Doppler curves for the top-5 satellites by
%    peak elevation.
%
%    Output is a ready-to-use struct  doppler_data  saved as a .mat file
%    that the tracking team can load directly for correlation and
%    satellite identification.
%
%  ┌─────────────────────────────────────────────────────────────────────┐
%  │  EDIT THIS SECTION ONLY — everything else is automatic              │
%  └─────────────────────────────────────────────────────────────────────┘
%
%  USAGE:
%    1.  Set the recording timestamp and duration below.
%    2.  Run.  Two files are written to OUTPUT_DIR:
%          doppler_curves_<YYYYMMDD_HHMMSS>.mat   ← data for tracking team
%          doppler_curves_<YYYYMMDD_HHMMSS>.png   ← summary figure (PDF also)
%    3.  Pass the .mat to the correlator.  The key matrices are documented
%        at the end of this file under §OUTPUT FORMAT.
%
%  DEPENDENCIES:
%    InputData/oneweb.tle   (same TLE file used by the main analysis)
%    No additional toolboxes required beyond core MATLAB.
%
% ═══════════════════════════════════════════════════════════════════════════

clear; close all; clc;

%% ══════════════════════════════════════════════════════════════════════════
%  §1  USER SETTINGS  ← EDIT HERE
% ══════════════════════════════════════════════════════════════════════════

% ── Recording timestamp (UTC) ─────────────────────────────────────────────
%    Set the exact moment the signal was captured.
REC_YEAR   = 2025;
REC_MONTH  = 10;
REC_DAY    = 14;
REC_HOUR   = 14;
REC_MIN    = 04;
REC_SEC    = 28.05;    % fractional seconds OK (e.g.  47.5  =  47 s 500 ms)

% ── Recording duration ────────────────────────────────────────────────────
REC_DURATION_S = 3;    % [s]  how many seconds of signal were recorded
REC_DT_S       = 0.01;    % [s]  desired time resolution of Doppler curves

% ── Receiver location ─────────────────────────────────────────────────────
RX_LAT_DEG = 52.0714774;   % Cranfield University (change if needed)
RX_LON_DEG = -0.6274441;
RX_ALT_M   = 100;

% ── Signal parameters ─────────────────────────────────────────────────────
FC_HZ   = 11.325e9;      % carrier frequency [Hz]  — OneWeb Ku-band downlink
C_MPS   = 299792458;     % speed of light [m/s]

% ── How many top satellites to output ─────────────────────────────────────
N_TOP   = 5;

% ── Minimum elevation to include a satellite in the output ────────────────
MIN_ELEV_DEG = 35;       % [deg]

% ── TLE file ──────────────────────────────────────────────────────────────
TLE_FILE = fullfile('InputData', 'oneweb.tle');

% ── Output directory ──────────────────────────────────────────────────────
OUTPUT_DIR = 'doppler_output';

%% ══════════════════════════════════════════════════════════════════════════
%  §2  DERIVED CONSTANTS  (do not edit)
% ══════════════════════════════════════════════════════════════════════════

LAMBDA  = C_MPS / FC_HZ;          % wavelength [m]
if ~exist(OUTPUT_DIR,'dir'), mkdir(OUTPUT_DIR); end

% Compose UTC timestamp
t0_UTC = datetime(REC_YEAR, REC_MONTH, REC_DAY, ...
                  REC_HOUR, REC_MIN, REC_SEC, ...
                  'TimeZone','UTC');

% Build time vector
t_rel   = (0 : REC_DT_S : REC_DURATION_S)';   % relative time [s]  K×1
K       = numel(t_rel);
t_UTC   = t0_UTC + seconds(t_rel);

% Timestamp string for filenames
ts_str  = datestr(t0_UTC, 'yyyymmdd_HHMMSS');

fprintf('╔══════════════════════════════════════════════════════════════╗\n');
fprintf('║  GEN_DOPPLER_CURVES — Tracking Team Sidequest               ║\n');
fprintf('╚══════════════════════════════════════════════════════════════╝\n');
fprintf('  Recording start : %s UTC\n', datestr(t0_UTC, 'yyyy-mm-dd HH:MM:SS.FFF'));
fprintf('  Duration        : %.1f s  (dt = %.3f s → %d epochs)\n', ...
        REC_DURATION_S, REC_DT_S, K);
fprintf('  Receiver        : %.5f° N,  %.5f° E,  %.0f m ASL\n', ...
        RX_LAT_DEG, RX_LON_DEG, RX_ALT_M);
fprintf('  Carrier         : %.6f GHz  (λ = %.4f mm)\n', ...
        FC_HZ/1e9, LAMBDA*1e3);
fprintf('  TLE file        : %s\n', TLE_FILE);
fprintf('\n');

%% ══════════════════════════════════════════════════════════════════════════
%  §3  GEOMETRY PROPAGATION
% ══════════════════════════════════════════════════════════════════════════

% Receiver ECEF position (static receiver)
[rx_x, rx_y, rx_z] = geodetic2ecef_dc(RX_LAT_DEG, RX_LON_DEG, RX_ALT_M);
recECEF = [rx_x; rx_y; rx_z];
R_e2n   = ecef2enu_rotation(deg2rad(RX_LAT_DEG), deg2rad(RX_LON_DEG));

fprintf('  Loading TLE file...\n');
tleData = read_tle_file_dc(TLE_FILE);
nTLE    = numel(tleData);
if nTLE == 0
    error('No TLEs loaded from %s. Check the path.', TLE_FILE);
end
fprintf('  Loaded %d TLE entries.\n', nTLE);

% ── Propagate all satellites ───────────────────────────────────────────────
fprintf('  Propagating %d satellites for %d epochs (may take a moment)...\n', nTLE, K);

% Pre-allocate: store results for every TLE, trim later
all_pos     = nan(3, K, nTLE);   % ECEF position  [3 × K × nTLE]  m
all_vel     = nan(3, K, nTLE);   % ECEF velocity  [3 × K × nTLE]  m/s
all_elev    = nan(K,    nTLE);   % elevation angle [K × nTLE]      deg
all_az      = nan(K,    nTLE);   % azimuth        [K × nTLE]       deg
all_range   = nan(K,    nTLE);   % slant range    [K × nTLE]       m
all_rdot    = nan(K,    nTLE);   % range-rate     [K × nTLE]       m/s
all_fdop    = nan(K,    nTLE);   % Doppler shift  [K × nTLE]       Hz
all_valid   = false(1,  nTLE);   % propagation succeeded?

for i = 1:nTLE
    ok = true;
    for k = 1:K
        [p, v] = propagate_j2_dc(tleData(i), t_UTC(k));
        if isempty(p), ok = false; break; end
        all_pos(:,k,i) = p;
        all_vel(:,k,i) = v;
    end
    if ~ok, continue; end
    all_valid(i) = true;

    for k = 1:K
        p  = all_pos(:,k,i);
        v  = all_vel(:,k,i);
        dr = p - recECEF;
        rho = norm(dr);
        u   = dr / rho;
        drENU = R_e2n * dr;

        all_elev(k,i)  = rad2deg(atan2(drENU(3), hypot(drENU(1),drENU(2))));
        all_az(k,i)    = wrapTo360_dc(atan2d(drENU(1), drENU(2)));
        all_range(k,i) = rho;

        % Range-rate: positive = satellite moving away
        rdot = v.' * u;    % (v_sat - v_rx) · u,  v_rx = 0 (static)
        all_rdot(k,i)  = rdot;

        % Doppler frequency shift:
        %   f_received = f_c * (1 - rdot/c)
        %   f_D = f_received - f_c = -f_c * rdot / c
        % Sign convention: approaching satellite → rdot < 0 → f_D > 0 (blueshift)
        all_fdop(k,i)  = -(FC_HZ / C_MPS) * rdot;
    end
end

nValid = sum(all_valid);
fprintf('  Successfully propagated %d / %d satellites.\n', nValid, nTLE);

%% ══════════════════════════════════════════════════════════════════════════
%  §4  SELECT TOP-N SATELLITES BY PEAK ELEVATION
% ══════════════════════════════════════════════════════════════════════════

% Compute peak elevation for each valid satellite during the recording window
peak_elev = nan(1, nTLE);
for i = 1:nTLE
    if all_valid(i)
        vis_k = all_elev(:,i) >= MIN_ELEV_DEG;
        if any(vis_k)
            peak_elev(i) = max(all_elev(vis_k, i));
        end
    end
end

% Sort and pick top N
[sorted_peak, sorted_idx] = sort(peak_elev, 'descend', 'MissingPlacement','last');
top_idx = sorted_idx(1:min(N_TOP, nTLE));
top_idx = top_idx(~isnan(sorted_peak(1:numel(top_idx))));   % remove NaN entries
if numel(top_idx) < N_TOP
    warning('GEN_DOPPLER:fewSats', ...
        'Only %d satellites visible above %.0f° during the recording window.', ...
        numel(top_idx), MIN_ELEV_DEG);
end
n_out = numel(top_idx);

fprintf('\n  Top %d satellites by peak elevation during recording window:\n', n_out);
fprintf('  %-4s  %-30s  %-10s  %-12s  %-12s\n', ...
        'Rank','TLE Name','Peak Elev','Min Doppler','Max Doppler');
fprintf('  %s\n', repmat('-',1,74));

for ri = 1:n_out
    i     = top_idx(ri);
    fd    = all_fdop(:,i);
    vis_k = all_elev(:,i) >= MIN_ELEV_DEG;
    fprintf('  %-4d  %-30s  %-10.2f  %-12.2f  %-12.2f\n', ri, ...
            tleData(i).name, sorted_peak(ri), ...
            min(fd(vis_k)), max(fd(vis_k)));
end

%% ══════════════════════════════════════════════════════════════════════════
%  §5  BUILD OUTPUT STRUCT
% ══════════════════════════════════════════════════════════════════════════
%
%  doppler_data is the hand-off struct for the tracking / correlation team.
%  See §OUTPUT FORMAT at the bottom for field-by-field documentation.

doppler_data.meta.recording_start_UTC   = t0_UTC;
doppler_data.meta.recording_start_str   = datestr(t0_UTC,'yyyy-mm-dd HH:MM:SS.FFF');
doppler_data.meta.duration_s            = REC_DURATION_S;
doppler_data.meta.dt_s                  = REC_DT_S;
doppler_data.meta.K_epochs              = K;
doppler_data.meta.rx_lat_deg            = RX_LAT_DEG;
doppler_data.meta.rx_lon_deg            = RX_LON_DEG;
doppler_data.meta.rx_alt_m              = RX_ALT_M;
doppler_data.meta.rx_ecef_m             = recECEF;
doppler_data.meta.fc_Hz                 = FC_HZ;
doppler_data.meta.lambda_m              = LAMBDA;
doppler_data.meta.c_mps                 = C_MPS;
doppler_data.meta.min_elev_deg          = MIN_ELEV_DEG;
doppler_data.meta.n_candidates          = n_out;
doppler_data.meta.tle_file              = TLE_FILE;

% Time axes
doppler_data.t_rel_s    = t_rel;           % K×1  relative time from recording start [s]
doppler_data.t_UTC      = t_UTC;           % K×1  datetime array (UTC)
doppler_data.t_unix_s   = posixtime(t_UTC);% K×1  POSIX / Unix timestamp [s]

% Per-candidate matrices  (K × n_out)
doppler_data.elev_deg   = nan(K, n_out);   % elevation angle            [deg]
doppler_data.az_deg     = nan(K, n_out);   % azimuth (N=0, E=90)        [deg]
doppler_data.range_m    = nan(K, n_out);   % slant range                [m]
doppler_data.rdot_mps   = nan(K, n_out);   % range-rate (+ = receding)  [m/s]
doppler_data.fdop_Hz    = nan(K, n_out);   % Doppler shift              [Hz]
doppler_data.fdot_Hzps  = nan(K, n_out);   % Doppler rate (d/dt fdop)   [Hz/s]
doppler_data.vis_mask   = false(K, n_out); % true when elev >= MIN_ELEV [logical]

% Per-candidate scalars (1 × n_out)
doppler_data.sat_names      = cell(1, n_out);
doppler_data.sat_raan_deg   = nan(1, n_out);
doppler_data.sat_incl_deg   = nan(1, n_out);
doppler_data.peak_elev_deg  = nan(1, n_out);
doppler_data.tle_index      = nan(1, n_out);   % index into tleData array

for ri = 1:n_out
    i  = top_idx(ri);
    fd = all_fdop(:,i);

    doppler_data.elev_deg(:,ri)  = all_elev(:,i);
    doppler_data.az_deg(:,ri)    = all_az(:,i);
    doppler_data.range_m(:,ri)   = all_range(:,i);
    doppler_data.rdot_mps(:,ri)  = all_rdot(:,i);
    doppler_data.fdop_Hz(:,ri)   = fd;
    doppler_data.vis_mask(:,ri)  = (all_elev(:,i) >= MIN_ELEV_DEG);

    % Doppler rate: numerical derivative (central difference)
    fdot = gradient(fd, REC_DT_S);
    doppler_data.fdot_Hzps(:,ri) = fdot;

    doppler_data.sat_names{ri}     = strtrim(tleData(i).name);
    doppler_data.sat_raan_deg(ri)  = tleData(i).satrec.raan;
    doppler_data.sat_incl_deg(ri)  = tleData(i).satrec.incl;
    doppler_data.peak_elev_deg(ri) = sorted_peak(ri);
    doppler_data.tle_index(ri)     = i;
end

% Save .mat
mat_fname = fullfile(OUTPUT_DIR, sprintf('doppler_curves_%s.mat', ts_str));
save(mat_fname, 'doppler_data');
fprintf('\n  Data saved → %s\n', mat_fname);

%% ══════════════════════════════════════════════════════════════════════════
%  §6  PLOTS
% ══════════════════════════════════════════════════════════════════════════

% Consistent colour palette for all plots (one colour per candidate)
PAL = [0.08 0.47 0.71;   % blue
       0.90 0.27 0.15;   % red
       0.16 0.63 0.30;   % green
       0.87 0.55 0.06;   % amber
       0.52 0.18 0.65];  % purple

% Short labels for legend (trim TLE name)
short_names = cell(1, n_out);
for ri = 1:n_out
    nm = strtrim(doppler_data.sat_names{ri});
    if numel(nm) > 22, nm = [nm(1:20) '...']; end
    short_names{ri} = sprintf('#%d  %s  (%.0f°)', ri, nm, ...
                               doppler_data.peak_elev_deg(ri));
end

t_axis = doppler_data.t_rel_s;   % convenience alias

% ── FIGURE 1: Doppler Frequency Curves ────────────────────────────────────
fh1 = figure('Name','Doppler Frequency Curves', 'Color','w', ...
             'Position',[50 50 1150 560]);
hold on; grid on; box on;

for ri = 1:n_out
    fd   = doppler_data.fdop_Hz(:,ri);
    vm   = doppler_data.vis_mask(:,ri);
    col  = PAL(ri,:);

    % Light dashed line for invisible portion
    invisible = ~vm;
    if any(invisible)
        plot(t_axis(invisible), fd(invisible)/1e3, ':', ...
             'Color', col*0.55 + 0.45, 'LineWidth', 0.9, ...
             'HandleVisibility','off');
    end
    % Solid thick line when satellite is visible
    plot(t_axis(vm), fd(vm)/1e3, '-', 'Color', col, 'LineWidth', 2.4, ...
         'DisplayName', short_names{ri});

    % Mark pass start and end
    if any(vm)
        k_s = find(vm,1,'first'); k_e = find(vm,1,'last');
        plot(t_axis(k_s), fd(k_s)/1e3, 'o', 'MarkerSize', 8, ...
             'MarkerFaceColor', col, 'MarkerEdgeColor','k', 'LineWidth',1.2, ...
             'HandleVisibility','off');
        plot(t_axis(k_e), fd(k_e)/1e3, 's', 'MarkerSize', 8, ...
             'MarkerFaceColor', col, 'MarkerEdgeColor','k', 'LineWidth',1.2, ...
             'HandleVisibility','off');
    end
end

yline(0,'k-','LineWidth',0.8,'HandleVisibility','off');
xlabel('Time from Recording Start (s)','FontSize',12,'FontWeight','bold');
ylabel('Doppler Shift  f_D  (kHz)',    'FontSize',12,'FontWeight','bold');
title({'Theoretical Doppler Frequency Curves — Top-5 Satellites by Elevation', ...
       sprintf('Recording: %s UTC  |  f_c = %.4f GHz  |  Rx: %.5f° N, %.5f° E', ...
           datestr(t0_UTC,'dd-mmm-yyyy HH:MM:SS'), FC_HZ/1e9, RX_LAT_DEG, RX_LON_DEG)}, ...
      'FontSize',12,'FontWeight','bold');
legend('Location','eastoutside','FontSize',9);
xlim([0, REC_DURATION_S]);

% Secondary x-axis showing actual UTC time
ax1 = gca;
ax_top = axes('Position',ax1.Position,'XAxisLocation','top','Color','none', ...
              'XLim',[0, REC_DURATION_S],'YLim',ax1.YLim,'YTick',[],'Box','off');
t_tick_s  = linspace(0, REC_DURATION_S, 7);
tick_lbls = cell(1,7);
for ti=1:7
    tick_lbls{ti} = datestr(t0_UTC + seconds(t_tick_s(ti)),'HH:MM:SS');
end
set(ax_top,'XTick',t_tick_s,'XTickLabel',tick_lbls,'FontSize',8);
xlabel(ax_top,'UTC','FontSize',9);
axes(ax1);   % restore current axes

save_fig_dc(fh1, fullfile(OUTPUT_DIR, sprintf('P_Doppler_Freq_%s', ts_str)));

% ── FIGURE 2: Doppler Rate Curves ─────────────────────────────────────────
fh2 = figure('Name','Doppler Rate Curves', 'Color','w', ...
             'Position',[100 100 1150 480]);
hold on; grid on; box on;

for ri = 1:n_out
    fd_dot = doppler_data.fdot_Hzps(:,ri);
    vm     = doppler_data.vis_mask(:,ri);
    col    = PAL(ri,:);
    if any(vm)
        plot(t_axis(vm), fd_dot(vm), '-', 'Color', col, 'LineWidth', 2.2, ...
             'DisplayName', short_names{ri});
    end
end

yline(0,'k-','LineWidth',0.8,'HandleVisibility','off');
xlabel('Time from Recording Start (s)','FontSize',12,'FontWeight','bold');
ylabel('Doppler Rate  \dotf_D  (Hz/s)','FontSize',12,'FontWeight','bold');
title({'Doppler Rate (d/dt) — Discriminator for Satellite Identification', ...
       sprintf('Recording: %s UTC', datestr(t0_UTC,'dd-mmm-yyyy HH:MM:SS'))}, ...
      'FontSize',12,'FontWeight','bold');
legend('Location','eastoutside','FontSize',9);
xlim([0, REC_DURATION_S]);
save_fig_dc(fh2, fullfile(OUTPUT_DIR, sprintf('P_Doppler_Rate_%s', ts_str)));

% ── FIGURE 3: Elevation & Azimuth Profiles ────────────────────────────────
fh3 = figure('Name','Elevation & Azimuth Profiles', 'Color','w', ...
             'Position',[150 150 1150 650]);

ax_el = subplot(2,1,1); hold on; grid on; box on;
for ri = 1:n_out
    el  = doppler_data.elev_deg(:,ri);
    vm  = doppler_data.vis_mask(:,ri);
    col = PAL(ri,:);
    plot(t_axis, el, '-', 'Color', col*0.4+0.6, 'LineWidth', 1.0, 'HandleVisibility','off');
    if any(vm)
        plot(t_axis(vm), el(vm), '-', 'Color', col, 'LineWidth', 2.5, ...
             'DisplayName', short_names{ri});
    end
end
yline(MIN_ELEV_DEG,'k--','LineWidth',1.2,'DisplayName',sprintf('%.0f° min elev.',MIN_ELEV_DEG));
ylabel('Elevation (deg)','FontSize',11,'FontWeight','bold');
title('Elevation Profiles During Recording Window','FontSize',11,'FontWeight','bold');
legend('Location','eastoutside','FontSize',8);
xlim([0, REC_DURATION_S]); ylim([0, 95]);

ax_az = subplot(2,1,2); hold on; grid on; box on;
for ri = 1:n_out
    az  = doppler_data.az_deg(:,ri);
    vm  = doppler_data.vis_mask(:,ri);
    col = PAL(ri,:);
    if any(vm)
        plot(t_axis(vm), az(vm), '-', 'Color', col, 'LineWidth', 2.5, ...
             'DisplayName', short_names{ri});
    end
end
set(gca,'YTick',0:45:360,'YTickLabel',{'N','NE','E','SE','S','SW','W','NW','N'});
ylabel('Azimuth (deg)','FontSize',11,'FontWeight','bold');
xlabel('Time from Recording Start (s)','FontSize',11,'FontWeight','bold');
title('Azimuth Profiles (N=0°, E=90°, clockwise)','FontSize',11,'FontWeight','bold');
legend('Location','eastoutside','FontSize',8);
xlim([0, REC_DURATION_S]);
linkaxes([ax_el, ax_az],'x');
sgtitle({'Geometry Profiles — Top-5 Satellites', ...
         sprintf('Recording: %s UTC', datestr(t0_UTC,'dd-mmm-yyyy HH:MM:SS'))}, ...
        'FontSize',12,'FontWeight','bold');
save_fig_dc(fh3, fullfile(OUTPUT_DIR, sprintf('P_Elevation_Azimuth_%s', ts_str)));

% ── FIGURE 4: Skyview (Polar) ──────────────────────────────────────────────
fh4 = figure('Name','Skyview — Polar Plot', 'Color','w', ...
             'Position',[200 200 680 680]);
pax = polaraxes;
hold on;
pax.ThetaZeroLocation = 'top';    % North at top
pax.ThetaDir          = 'clockwise';
pax.RLim              = [0 90];
pax.RTick             = [0 30 60 90];
pax.RTickLabel        = {'90°  ','60°  ','30°  ','0°  '};
pax.GridAlpha         = 0.35;
pax.FontSize          = 9;

for ri = 1:n_out
    el_r  = doppler_data.elev_deg(:,ri);
    az_r  = doppler_data.az_deg(:,ri);
    vm    = doppler_data.vis_mask(:,ri);
    col   = PAL(ri,:);

    % Full arc in pale colour
    all_above_horiz = el_r > 0;
    if any(all_above_horiz)
        r_full = 90 - el_r(all_above_horiz);
        polarplot(deg2rad(az_r(all_above_horiz)), r_full, '-', ...
                  'Color', col*0.3+0.7, 'LineWidth', 1.0, 'HandleVisibility','off');
    end

    % Visible segment (thick)
    if any(vm)
        r_vis = 90 - el_r(vm);
        polarplot(deg2rad(az_r(vm)), r_vis, '-', 'Color', col, 'LineWidth', 3.0, ...
                  'DisplayName', short_names{ri});

        % Start marker
        k_s = find(vm,1,'first');
        polarplot(deg2rad(az_r(k_s)), 90-el_r(k_s), 'o', ...
                  'MarkerSize',9,'MarkerFaceColor',col,'MarkerEdgeColor','k', ...
                  'LineWidth',1.2,'HandleVisibility','off');
        % Direction arrow at midpoint
        k_m = find(vm, floor(sum(vm)/2), 'first'); k_m = k_m(end);
        if k_m < find(vm,1,'last')
            polarplot(deg2rad(az_r(k_m)), 90-el_r(k_m), '>', ...
                      'MarkerSize',7,'MarkerFaceColor',col,'MarkerEdgeColor','k', ...
                      'LineWidth',1.2,'HandleVisibility','off');
        end
    end
end

legend('Location','southoutside','FontSize',8,'NumColumns',2);
title({'Skyview from Rx Station — Top-5 Satellites', ...
       sprintf('Recording: %s UTC  |  Rx: %.5f°N, %.5f°E', ...
           datestr(t0_UTC,'dd-mmm-yyyy HH:MM:SS'), RX_LAT_DEG, RX_LON_DEG), ...
       'Thick arc = visible during recording  |  Circle = pass start'}, ...
      'FontSize',10,'FontWeight','bold');
save_fig_dc(fh4, fullfile(OUTPUT_DIR, sprintf('P_Skyview_%s', ts_str)));

% ── FIGURE 5: Doppler + Elevation — Combined Panel (Correlation Sheet) ────
% This is the "cheat sheet" the tracking team uses for manual correlation:
% top row = Doppler, bottom row = elevation, one column per satellite.
fh5 = figure('Name','Correlation Sheet', 'Color','w', ...
             'Position',[250 50 1400 700]);

for ri = 1:n_out
    col   = PAL(ri,:);
    fd    = doppler_data.fdop_Hz(:,ri);
    el    = doppler_data.elev_deg(:,ri);
    vm    = doppler_data.vis_mask(:,ri);
    fdot  = doppler_data.fdot_Hzps(:,ri);

    % ── Top panel: Doppler frequency ──────────────────────────────────────
    ax_d = subplot(2, n_out, ri);
    hold on; grid on; box on;
    plot(t_axis, fd/1e3, ':', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.9, ...
         'HandleVisibility','off');
    if any(vm)
        plot(t_axis(vm), fd(vm)/1e3, '-', 'Color', col, 'LineWidth', 2.5);
        % Annotate start and end frequency
        k_s = find(vm,1,'first');
        k_e = find(vm,1,'last');
        text(t_axis(k_s), fd(k_s)/1e3, sprintf(' %.2f kHz', fd(k_s)/1e3), ...
             'FontSize',7,'Color',col,'VerticalAlignment','bottom');
        text(t_axis(k_e), fd(k_e)/1e3, sprintf(' %.2f kHz', fd(k_e)/1e3), ...
             'FontSize',7,'Color',col,'VerticalAlignment','top');
    end
    yline(0,'k-','LineWidth',0.7,'HandleVisibility','off');

    % Add Doppler rate as thin dotted second line on same axes (scaled)
    yyaxis right;
    if any(vm)
        plot(t_axis(vm), fdot(vm), ':', 'Color', col*0.65, 'LineWidth', 1.2, ...
             'HandleVisibility','off');
    end
    set(gca,'YColor', col*0.65);
    ylabel('f_D rate (Hz/s)','FontSize',7);
    yyaxis left;

    nm   = strtrim(doppler_data.sat_names{ri});
    if numel(nm)>18, nm=[nm(1:16) '...']; end
    title(sprintf('%s\nRAAN=%.0f°  i=%.1f°', nm, ...
                  doppler_data.sat_raan_deg(ri), doppler_data.sat_incl_deg(ri)), ...
          'FontSize',9,'FontWeight','bold','Color',col);
    if ri==1, ylabel('Doppler f_D (kHz)','FontSize',9,'FontWeight','bold'); end
    set(gca,'FontSize',8); xlim([0, REC_DURATION_S]);

    % ── Bottom panel: Elevation ────────────────────────────────────────────
    ax_e = subplot(2, n_out, n_out + ri);
    hold on; grid on; box on;
    % Full arc pale
    plot(t_axis, el, '-', 'Color', col*0.3+0.7, 'LineWidth', 1.0, ...
         'HandleVisibility','off');
    % Visible portion thick
    if any(vm)
        plot(t_axis(vm), el(vm), '-', 'Color', col, 'LineWidth', 2.5);
        % Peak annotation
        [el_pk, k_pk] = max(el);
        text(t_axis(k_pk), el_pk+2, sprintf('%.0f°', el_pk), ...
             'FontSize',8,'Color',col,'HorizontalAlignment','center', ...
             'FontWeight','bold');
    end
    yline(MIN_ELEV_DEG,'k--','LineWidth',0.9,'HandleVisibility','off');

    xlim([0, REC_DURATION_S]); ylim([0 95]);
    set(gca,'FontSize',8);
    if ri==1, ylabel('Elevation (deg)','FontSize',9,'FontWeight','bold'); end
    xlabel('t (s)','FontSize',8);
end

sgtitle({'Doppler Correlation Sheet — Recording ' datestr(t0_UTC,'dd-mmm-yyyy HH:MM:SS') ' UTC', ...
         sprintf('f_c = %.4f GHz  |  Rx: %.5f°N, %.5f°E  |  Top %d satellites by peak elevation', ...
             FC_HZ/1e9, RX_LAT_DEG, RX_LON_DEG, n_out)}, ...
        'FontSize',11,'FontWeight','bold');
save_fig_dc(fh5, fullfile(OUTPUT_DIR, sprintf('P_Correlation_Sheet_%s', ts_str)));

%% ══════════════════════════════════════════════════════════════════════════
%  §7  CONSOLE SUMMARY
% ══════════════════════════════════════════════════════════════════════════

fprintf('\n╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║  DOPPLER CURVE GENERATION COMPLETE                              ║\n');
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Recording: %s UTC                         ║\n', ...
        datestr(t0_UTC,'dd-mmm-yyyy HH:MM:SS'));
fprintf('║  Duration : %.1f s  |  dt = %.3f s  |  %d epochs             ║\n', ...
        REC_DURATION_S, REC_DT_S, K);
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Candidate Satellites:                                          ║\n');
for ri = 1:n_out
    fd_vis = doppler_data.fdop_Hz(doppler_data.vis_mask(:,ri), ri);
    if isempty(fd_vis), fd_vis = [0;0]; end
    vis_s  = sum(doppler_data.vis_mask(:,ri)) * REC_DT_S;
    fprintf('║  [%d] %-22s  PeakEl=%5.1f°  VisTime=%5.0f s      ║\n', ...
            ri, strtrim(doppler_data.sat_names{ri}), ...
            doppler_data.peak_elev_deg(ri), vis_s);
    fprintf('║      f_D range: [%+8.2f, %+8.2f] kHz  RAAN=%.0f°  i=%.1f°      ║\n', ...
            min(fd_vis)/1e3, max(fd_vis)/1e3, ...
            doppler_data.sat_raan_deg(ri), doppler_data.sat_incl_deg(ri));
end
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Output files:                                                  ║\n');
fprintf('║    %s  ║\n', mat_fname);
fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');

fprintf('  Hand-off to tracking team:\n');
fprintf('    load(''%s'')\n', mat_fname);
fprintf('    %% Key matrices  (all K×n_out, rows=time, cols=satellite rank)\n');
fprintf('    %%   doppler_data.fdop_Hz    — Doppler shift [Hz]\n');
fprintf('    %%   doppler_data.fdot_Hzps  — Doppler rate  [Hz/s]\n');
fprintf('    %%   doppler_data.rdot_mps   — range-rate    [m/s]\n');
fprintf('    %%   doppler_data.elev_deg   — elevation     [deg]\n');
fprintf('    %%   doppler_data.vis_mask   — visibility    [logical]\n');
fprintf('    %%   doppler_data.t_rel_s    — time vector   [s]\n');
fprintf('    %%   doppler_data.sat_names  — cell, 1×n_out\n\n');


%% ══════════════════════════════════════════════════════════════════════════
%  §8  LOCAL HELPER FUNCTIONS
% ══════════════════════════════════════════════════════════════════════════

function [x,y,z] = geodetic2ecef_dc(lat_deg, lon_deg, alt_m)
    a=6378137.0; e2=0.00669437999014;
    lat=deg2rad(lat_deg); lon=deg2rad(lon_deg);
    N=a/sqrt(1-e2*sin(lat)^2);
    x=(N+alt_m)*cos(lat)*cos(lon);
    y=(N+alt_m)*cos(lat)*sin(lon);
    z=(N*(1-e2)+alt_m)*sin(lat);
end

function R = ecef2enu_rotation(lat, lon)
    sl=sin(lat); cl=cos(lat); sn=sin(lon); cn=cos(lon);
    R=[-sn,cn,0; -sl*cn,-sl*sn,cl; cl*cn,cl*sn,sl];
end

function w = wrapTo360_dc(a), w = mod(a,360); end

function tleData = read_tle_file_dc(fn)
    fid=fopen(fn,'r'); L={}; tleData=struct('name',{},'epochUTC',{},'satrec',{});
    if fid<0, error('Cannot open TLE file: %s', fn); end
    while ~feof(fid)
        s=strtrim(fgetl(fid));
        if ~isempty(s) && ischar(s), L{end+1}=s; end %#ok<AGROW>
    end
    fclose(fid);
    c=0; i=1;
    while i<=numel(L)-2
        if length(L{i+1})>=69 && length(L{i+2})>=69
            c=c+1;
            tleData(c).name = L{i};
            ey=str2double(L{i+1}(19:20)); ed=str2double(L{i+1}(21:32));
            yr=2000+ey; if ey>=57, yr=1900+ey; end
            tleData(c).epochUTC = datetime(yr,1,1,'TimeZone','UTC')+days(ed-1);
            mm=str2double(L{i+2}(53:63)); nr=mm*2*pi/86400; mu=3.986004418e14;
            tleData(c).satrec = struct( ...
                'a',   (mu/nr^2)^(1/3), ...
                'ecc', str2double(['0.' L{i+2}(27:33)]), ...
                'incl',str2double(L{i+2}(9:16)), ...
                'raan',str2double(L{i+2}(18:25)), ...
                'argp',str2double(L{i+2}(35:42)), ...
                'M0',  str2double(L{i+2}(44:51)), ...
                'n_rad_s', nr);
            i=i+3;
        else
            i=i+1;
        end
    end
end

function [p,v] = propagate_j2_dc(tle, tgt_datetime)
    sr=tle.satrec; mu=3.986e14; J2=1.08263e-3; Re=6378137;
    dt_=seconds(tgt_datetime - tle.epochUTC);
    a=sr.a; e=sr.ecc; inc=deg2rad(sr.incl);
    O0=deg2rad(sr.raan); w0=deg2rad(sr.argp); M0=deg2rad(sr.M0); n=sr.n_rad_s;
    pp=a*(1-e^2); ci=cos(inc);
    Od=-1.5*n*J2*(Re/pp)^2*ci;
    wd= 0.75*n*J2*(Re/pp)^2*(5*ci^2-1);
    Md= n + 0.75*n*J2*(Re/pp)^2*sqrt(1-e^2)*(3*ci^2-1);
    O=O0+Od*dt_; w=w0+wd*dt_; M=mod(M0+Md*dt_,2*pi);
    E=M;
    for it=1:50
        dE=(M-E+e*sin(E))/(1-e*cos(E)); E=E+dE;
        if abs(dE)<1e-13, break; end
    end
    nu=2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
    r_=a*(1-e*cos(E));
    rp_=[r_*cos(nu); r_*sin(nu); 0];
    h_=sqrt(mu*pp); vp_=(mu/h_)*[-sin(nu); e+cos(nu); 0];
    cO=cos(O); sO=sin(O); cw=cos(w); sw=sin(w); si=sin(inc);
    Rot=[cO*cw-sO*sw*ci, -cO*sw-sO*cw*ci,  sO*si;
         sO*cw+cO*sw*ci, -sO*sw+cO*cw*ci, -cO*si;
         sw*si,           cw*si,            ci];
    pi_=Rot*rp_; vi_=Rot*vp_;
    we=7.2921159e-5;
    jd=juliandate(tgt_datetime);
    T=(jd-2451545)/36525;
    g=280.4606+360.9856*(jd-2451545)+0.00038*T^2;
    th=deg2rad(mod(g,360));
    ct=cos(th); st=sin(th);
    R_=[ct,st,0; -st,ct,0; 0,0,1];
    dR=we*[-st,ct,0; -ct,-st,0; 0,0,0];
    p=R_*pi_; v=R_*vi_+dR*pi_;
end

function save_fig_dc(fh, base_path)
    % Save as both .png (quick preview) and .pdf (vector, for reports)
    try
        exportgraphics(fh, [base_path '.png'], 'Resolution', 200);
        exportgraphics(fh, [base_path '.pdf'], 'ContentType','vector');
    catch
        % Fallback for older MATLAB versions
        print(fh, [base_path '.png'], '-dpng', '-r200');
        print(fh, [base_path '.pdf'], '-dpdf');
    end
end


%% ══════════════════════════════════════════════════════════════════════════
%  §OUTPUT FORMAT — Reference for the tracking team
% ══════════════════════════════════════════════════════════════════════════
%
%  load('doppler_curves_YYYYMMDD_HHMMSS.mat')   → struct doppler_data
%
%  ┌─────────────────────────────────────────────────────────────────────┐
%  │  MATRICES  (all K × n_out, rows = time epochs, cols = sat rank)     │
%  ├──────────────────────────┬──────────┬──────────────────────────────┤
%  │ Field                    │ Unit     │ Description                  │
%  ├──────────────────────────┼──────────┼──────────────────────────────┤
%  │ fdop_Hz(:,ri)            │ Hz       │ Doppler frequency shift      │
%  │                          │          │ f_D = -(fc/c) * rdot         │
%  │                          │          │ + = blueshift (approaching)  │
%  │ fdot_Hzps(:,ri)          │ Hz/s     │ Doppler rate df_D/dt         │
%  │ rdot_mps(:,ri)           │ m/s      │ Range-rate (+ = receding)    │
%  │ range_m(:,ri)            │ m        │ Slant range receiver↔sat     │
%  │ elev_deg(:,ri)           │ deg      │ Elevation angle              │
%  │ az_deg(:,ri)             │ deg      │ Azimuth (N=0, E=90, CW)      │
%  │ vis_mask(:,ri)           │ logical  │ true when elev ≥ MIN_ELEV    │
%  ├──────────────────────────┼──────────┼──────────────────────────────┤
%  │  TIME VECTORS  (K × 1)                                             │
%  ├──────────────────────────┼──────────┼──────────────────────────────┤
%  │ t_rel_s                  │ s        │ Time from recording start    │
%  │ t_unix_s                 │ s        │ POSIX / Unix epoch time      │
%  │ t_UTC                    │ datetime │ UTC datetime array           │
%  ├──────────────────────────┼──────────┼──────────────────────────────┤
%  │  PER-SATELLITE VECTORS  (1 × n_out)                                │
%  ├──────────────────────────┼──────────┼──────────────────────────────┤
%  │ sat_names{ri}            │ string   │ TLE name                     │
%  │ sat_raan_deg(ri)         │ deg      │ RAAN from TLE epoch          │
%  │ sat_incl_deg(ri)         │ deg      │ Inclination from TLE epoch   │
%  │ peak_elev_deg(ri)        │ deg      │ Peak elevation during window │
%  └──────────────────────────┴──────────┴──────────────────────────────┘
%
%  SIGN CONVENTION (Doppler):
%    Satellite approaching receiver  →  rdot < 0  →  fdop_Hz > 0  (blueshift)
%    Satellite receding  from receiver  →  rdot > 0  →  fdop_Hz < 0  (redshift)
%    Crossover (maximum range-rate → 0) marks the point of closest approach.
%
%  CORRELATION WORKFLOW (suggested):
%    1.  Load the .mat and select the candidate column by rank (col index).
%    2.  Extract fdop_Hz(:,ri) as the reference template.
%    3.  Compute cross-correlation between the received signal's measured
%        Doppler and the template.  Peak lag → time offset.
%    4.  Use fdot_Hzps as a secondary discriminator: the rate of change
%        is unique per orbital plane and can reject same-frequency ambiguities.
%    5.  Use elev_deg + az_deg for geometry-based plausibility checks.
% ══════════════════════════════════════════════════════════════════════════
