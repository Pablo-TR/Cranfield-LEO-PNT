function out = real_strtg_sequential_EKF_geom(cfg)
% REAL_STRTG_SEQUENTIAL_EKF_FASTREALISTIC
% Sequential IEKF for static receiver positioning using:
%   - Doppler
%   - Doppler rate
%   - Optional height pseudo-measurement
%
% State:
%   x = [dE dN dU clkBias clkDrift]^T
%
% Notes:
%   - No bootstrap
%   - No solver-side tuning defaults
%   - Truth is used only in simulation and diagnostics
%   - Doppler-rate is modeled as:
%         rhoddot + tropoAcc + clkDrift
%
% Required cfg fields must be set from MAIN.

if nargin < 1
    error('cfg must be provided from main.');
end

cfg = validate_and_derive_cfg(cfg);


%% ── SETUP ────────────────────────────────────────────────────────────────
[rx,ry,rz] = geodetic2ecef_local(cfg.rx.lat_deg, cfg.rx.lon_deg, cfg.rx.alt_m);
recTrue    = [rx;ry;rz];
recVelTrue = cfg.rx.vel_ecef(:);

R_e2n = compute_R_ecef2enu_local(deg2rad(cfg.rx.lat_deg), deg2rad(cfg.rx.lon_deg));
R_n2e = R_e2n.';

% Initial guessed position used by estimator
r0_guess = recTrue + R_n2e * cfg.init.errENU_m(:);


%% ── TLE & TIME ───────────────────────────────────────────────────────────
[tleData, ~] = read_tle_file_fast(cfg.tleFile, cfg.useSGP4);
nTLE = numel(tleData);
if nTLE == 0
    error('No TLEs read. Check cfg.tleFile.');
end

epochs = NaT(nTLE,1,'TimeZone','UTC');
for i = 1:nTLE
    epochs(i) = tleData(i).epochUTC;
end
validEpoch = ~isnat(epochs);

if isempty(cfg.startUTC)
    vIdx = find(validEpoch);
    if isempty(vIdx)
        error('No valid TLE epochs found.');
    end
    switch lower(cfg.startFromTLE)
        case 'first'
            cfg.startUTC = epochs(vIdx(1));
        otherwise
            cfg.startUTC = median(epochs(vIdx));
    end
end

if ~isempty(cfg.stopUTC)
    cfg.duration_sec = seconds(cfg.stopUTC - cfg.startUTC);
end
if cfg.duration_sec <= 0
    error('cfg.duration_sec must be > 0');
end

ages   = abs(days(cfg.startUTC - epochs));
hasOrb = arrayfun(@(d) ~isempty(d.satrec_data), tleData).';
freshIdx = find(validEpoch & (ages <= cfg.maxTLEAgeDays) & hasOrb);
if isempty(freshIdx)
    warning('No TLEs within %.1f days. Using all valid orbital TLEs.', cfg.maxTLEAgeDays);
    freshIdx = find(validEpoch & hasOrb);
end
if isempty(freshIdx)
    error('No usable TLEs with orbital data.');
end

t    = (0:cfg.dt:cfg.duration_sec)';
K    = numel(t);
tUTC = cfg.startUTC + seconds(t);
tUTC.TimeZone = 'UTC';


%% ── CONSTELLATION / GEOMETRY ─────────────────────────────────────────────
[sat, ~] = build_constellation_fast(tleData, freshIdx, tUTC, recTrue, R_e2n, cfg);

Nsat = sat.Nsat;
if Nsat == 0
    error('No candidate satellites found.');
end

% Satellite acceleration from propagated velocity
sat.acc = nan(size(sat.vel));
for s = 1:sat.Nsat
    V = squeeze(sat.vel(:,1,s,:));      % 3 x K
    Aecef = deriv1_cols(V, cfg.dt);     % 3 x K
    sat.acc(:,1,s,:) = reshape(Aecef, [3,1,1,size(V,2)]);
end


%% ── VISIBILITY / CAPTURE ─────────────────────────────────────────────────
vis = sat.elev >= cfg.minElevDeg;

sat.az = nan(K,Nsat);
for s = 1:Nsat
    for k = find(vis(:,s)).'
        drENU = R_e2n * (sat.pos(:,1,s,k) - recTrue);
        sat.az(k,s) = wrapTo360(atan2d(drENU(1), drENU(2)));
    end
end

trk = select_capture_stream_online(vis, sat.elev, sat.az, sat.raan, cfg.capture, cfg.dt);


%% ── MEASUREMENT SIMULATION ───────────────────────────────────────────────
rng(cfg.seedMeas);
[meas, truthAux] = generate_measurements_fastrealistic(sat, vis, recTrue, recVelTrue, cfg);

% Attach visibility summary to meas struct
meas.vis      = vis;
meas.nVisible = sum(vis,2);
meas.anyVis   = meas.nVisible > 0;

%% ── FILTER DEFINITION ────────────────────────────────────────────────────

nx = 5;
iDE = 1; iDN = 2; iDU = 3;
iCB = 4; iCD = 5;

% x contains correction relative to r0_guess
x = zeros(nx,1);
x(iCB) = cfg.clk.filter_init_bias_mps;
x(iCD) = cfg.clk.filter_init_drift_mps2;

P = diag([ ...
    cfg.init.sigENU_m(1)^2, ...
    cfg.init.sigENU_m(2)^2, ...
    cfg.init.sigENU_m(3)^2, ...
    cfg.P0_clkBias, ...
    cfg.P0_clkDrift]);

A = eye(nx);
A(iCB,iCD) = cfg.dt;

qPos = cfg.sigma_pos_rw^2;
Q = diag([ ...
    qPos * cfg.dt, ...
    qPos * cfg.dt, ...
    qPos * cfg.dt, ...
    (cfg.clk.filt.sigmaBias_rw_mps_sqrt^2)   * cfg.dt, ...
    (cfg.clk.filt.sigmaDrift_rw_mps2_sqrt^2) * cfg.dt]);

Inx = eye(nx);


% Prior state before any measurement update
rEst0    = r0_guess + R_n2e * x(iDE:iDU);
errECEF0 = rEst0 - recTrue;
errENU0  = R_e2n * errECEF0;
err3d0   = norm(errECEF0);


%% ── STORAGE ──────────────────────────────────────────────────────────────
xhist         = nan(nx,K);
Pfull         = cell(K,1);
rEstHist      = nan(3,K);
errECEF       = nan(3,K);
errENU        = nan(3,K);
err3d         = nan(K,1);
nisHist       = nan(K,1);
nUpdates      = zeros(K,1);
satUsed       = zeros(K,1);

zPredDoppUsed = nan(K,1);
zPredRateUsed = nan(K,1);
nu_hist       = nan(K,1);    % Doppler innovation history
nu_rate_hist  = nan(K,1);    % Doppler-rate innovation history

% Prior-augmented histories
t_prior        = [-cfg.dt; t];
xhist_prior    = nan(nx, K+1);
Pfull_prior    = cell(K+1,1);
rEstHist_prior = nan(3,K+1);
errECEF_prior  = nan(3,K+1);
errENU_prior   = nan(3,K+1);
err3d_prior    = nan(K+1,1);

xhist_prior(:,1)    = x;
Pfull_prior{1}      = P;
rEstHist_prior(:,1) = rEst0;
errECEF_prior(:,1)  = errECEF0;
errENU_prior(:,1)   = errENU0;
err3d_prior(1)      = err3d0;
% Track consecutive epochs each satellite has been used
consecutive_tracking = zeros(1, Nsat);
%% ── MAIN LOOP ────────────────────────────────────────────────────────────
for k = 1:K

    % Predict
    x = A * x;
    P = A * P * A' + Q;

    % Determine available real measurements at this epoch
    avail = false(1,Nsat);
    if cfg.meas.useDoppler
        avail = avail | ~isnan(meas.dopp.obs(k,:));
    end
    if cfg.meas.useDopplerRate
        avail = avail | ~isnan(meas.rate.obs(k,:));
    end

    if strcmpi(cfg.capture.mode,'allVisible')
        useList = find(vis(k,:) & avail);
        if ~isempty(useList)
            [~,ord] = sort(sat.elev(k,useList), 'descend');
            useList = useList(ord);
        end
    else
        sSel = trk.sid(k);
        if trk.blackout(k) || sSel == 0 || ~vis(k,sSel) || ~avail(sSel)
            useList = [];
        else
            useList = sSel;
        end
    end
    % Reset tracking counters for satellites not in use or during blackout
    for ss = 1:Nsat
        if ~ismember(ss, useList) || (exist('trk','var') && trk.blackout(k))
            consecutive_tracking(ss) = 0;
        end
    end
    epochNIS = [];

    % % ── Real measurement updates (Doppler and Doppler-rate) ─────────────
    % for s = useList
    %     ps = sat.pos(:,1,s,k);
    %     vs = sat.vel(:,1,s,k);
    %     as = sat.acc(:,1,s,k);
    % 
    %     % 1) Doppler
    %     if cfg.meas.useDoppler && ~isnan(meas.dopp.obs(k,s))
    %         z  = meas.dopp.obs(k,s);
    %         Rk = meas.dopp.sigma(k,s)^2;
    % 
    %         if k <= cfg.filter.startupInflateEpochs
    %             Rk = cfg.filter.startupRmult * Rk;
    %         end
    % 
    %         [x, P, accepted, nis, nu, h] = scalar_iekf_update( ...
    %             x, P, z, Rk, ...
    %             @(xi) meas_model_doppler( ...
    %                 xi, r0_guess + R_n2e*xi(iDE:iDU), R_n2e, ps, vs, ...
    %                 k, s, cfg, meas, iCB, nx), ...
    %             cfg, Inx, [iDE iDN iDU]);
    % 
    %         if accepted
    %             epochNIS(end+1,1) = nis; %#ok<AGROW>
    %             nUpdates(k)       = nUpdates(k) + 1;
    %             satUsed(k)        = s;
    %             zPredDoppUsed(k)  = h;
    %             nu_hist(k)        = nu;
    %         end
    %     end
    % 
    %     % 2) Doppler-rate
    %     if cfg.meas.useDopplerRate && ~isnan(meas.rate.obs(k,s))
    %         z  = meas.rate.obs(k,s);
    %         Rk = meas.rate.sigma(k,s)^2;
    % 
    %         [x, P, accepted, nis, nu, h] = scalar_iekf_update( ...
    %             x, P, z, Rk, ...
    %             @(xi) meas_model_doppler_rate( ...
    %                 xi, r0_guess + R_n2e*xi(iDE:iDU), R_n2e, ps, vs, as, ...
    %                 k, s, cfg, meas, iCD, nx), ...
    %             cfg, Inx, [iDE iDN iDU]);
    % 
    %         if accepted
    %             epochNIS(end+1,1) = nis; %#ok<AGROW>
    %             nUpdates(k)       = nUpdates(k) + 1;
    %             satUsed(k)        = s;
    %             zPredRateUsed(k)  = h;
    %             nu_rate_hist(k)   = nu;
    %         end
    %     end
    % end


% % ── Real measurement updates (Doppler and Doppler-rate) ─────────────
    for s = useList
        ps = sat.pos(:,1,s,k);
        vs = sat.vel(:,1,s,k);
        as = sat.acc(:,1,s,k);
        
        hasD  = cfg.meas.useDoppler     && ~isnan(meas.dopp.obs(k,s));
        hasDd = cfg.meas.useDopplerRate && ~isnan(meas.rate.obs(k,s));
        
        % Increment the tracking counter for this specific satellite
        consecutive_tracking(s) = consecutive_tracking(s) + 1;
        
        % ── CASO 1: Doppler + Doppler Rate Simultáneos ──
        if hasD && hasDd
            z  = [meas.dopp.obs(k,s); meas.rate.obs(k,s)];
            
            % Creamos la matriz R diagonal explícitamente
            Rk = diag([meas.dopp.sigma(k,s)^2, meas.rate.sigma(k,s)^2]);
                
            % Startup & Handover Inflation (Lock-in Dampening)
            if consecutive_tracking(s) <= cfg.filter.startupInflateEpochs
                Rk = cfg.filter.startupRmult * Rk;
            end
                
            [x, P, accepted, nis, nu, h] = vector_iekf_update( ...
                x, P, z, Rk, ...
                @(xi) meas_model_dopp_and_rate( ...
                    xi, r0_guess + R_n2e*xi(iDE:iDU), R_n2e, ...
                    ps, vs, as, k, s, cfg, meas, iCB, iCD, nx), ...
                cfg, Inx, [iDE iDN iDU]);
                
            if accepted
                epochNIS(end+1,1) = nis; %#ok<AGROW>
                nUpdates(k)       = nUpdates(k) + 2;
                satUsed(k)        = s;
                zPredDoppUsed(k)  = h(1);
                zPredRateUsed(k)  = h(2);
                nu_hist(k)        = nu(1);
                nu_rate_hist(k)   = nu(2);
            end
            
        % ── CASO 2: Solo Doppler ──
        elseif hasD
            z  = meas.dopp.obs(k,s);
            Rk = meas.dopp.sigma(k,s)^2;
            
            % Startup & Handover Inflation (Lock-in Dampening)
            if consecutive_tracking(s) <= cfg.filter.startupInflateEpochs
                Rk = cfg.filter.startupRmult * Rk;
            end
            
            [x, P, accepted, nis, nu, h] = scalar_iekf_update( ...
                x, P, z, Rk, ...
                @(xi) meas_model_doppler( ...
                    xi, r0_guess + R_n2e*xi(iDE:iDU), R_n2e, ps, vs, ...
                    k, s, cfg, meas, iCB, nx), ...
                cfg, Inx, [iDE iDN iDU]);
                
            if accepted
                epochNIS(end+1,1) = nis; %#ok<AGROW>
                nUpdates(k)       = nUpdates(k) + 1;
                satUsed(k)        = s;
                zPredDoppUsed(k)  = h;
                nu_hist(k)        = nu;
            end
            
        % ── CASO 3: Solo Doppler Rate ──
        elseif hasDd
            z  = meas.rate.obs(k,s);
            Rk = meas.rate.sigma(k,s)^2;
            
            % Startup & Handover Inflation (Lock-in Dampening)
            if consecutive_tracking(s) <= cfg.filter.startupInflateEpochs
                Rk = cfg.filter.startupRmult * Rk;
            end
            
            [x, P, accepted, nis, nu, h] = scalar_iekf_update( ...
                x, P, z, Rk, ...
                @(xi) meas_model_doppler_rate( ...
                    xi, r0_guess + R_n2e*xi(iDE:iDU), R_n2e, ps, vs, as, ...
                    k, s, cfg, meas, iCD, nx), ...
                cfg, Inx, [iDE iDN iDU]);
                
            if accepted
                epochNIS(end+1,1) = nis; %#ok<AGROW>
                nUpdates(k)       = nUpdates(k) + 1;
                satUsed(k)        = s;
                zPredRateUsed(k)  = h;
                nu_rate_hist(k)   = nu;
            end
        end
    end


    % 3) Height pseudo-measurement (once per epoch, not per satellite)
    if cfg.height.enable
        zH  = cfg.height.z0_m;
        RH  = cfg.height.sigma_m^2;

        [x, P] = scalar_iekf_update( ...
            x, P, zH, RH, ...
            @(xi) meas_model_height_constraint(xi, iDU, nx), ...
            cfg, Inx, [iDE iDN iDU]);
        % Height pseudo-measurement is not included in nisHist
    end

    if ~isempty(epochNIS)
        nisHist(k) = mean(epochNIS);
    end

    % Store posterior
    xhist(:,k) = x;
    Pfull{k}   = P;

    rEst = r0_guess + R_n2e * x(iDE:iDU);
    rEstHist(:,k) = rEst;
    errECEF(:,k)  = rEst - recTrue;
    errENU(:,k)   = R_e2n * errECEF(:,k);
    err3d(k)      = norm(errECEF(:,k));

    xhist_prior(:,k+1)    = x;
    Pfull_prior{k+1}      = P;
    rEstHist_prior(:,k+1) = rEst;
    errECEF_prior(:,k+1)  = errECEF(:,k);
    errENU_prior(:,k+1)   = errENU(:,k);
    err3d_prior(k+1)      = err3d(k);
end


%% ── SUMMARY ──────────────────────────────────────────────────────────────

%% ── OUTPUTS ──────────────────────────────────────────────────────────────
out.cfg          = cfg;
out.t            = t;
out.tUTC         = tUTC;
out.sat          = sat;
out.vis          = vis;
out.trk          = trk;

out.xhist        = xhist;
out.Pfull        = Pfull;
out.r0_guess     = r0_guess;
out.rEstHist     = rEstHist;
out.errECEF      = errECEF;
out.errENU       = errENU;
out.err3d        = err3d;

out.nis          = nisHist;
out.nUpdates     = nUpdates;
out.satUsed      = satUsed;

out.zPredDoppUsed = zPredDoppUsed;
out.zPredRateUsed = zPredRateUsed;
out.nu_hist       = nu_hist;
out.nu_rate_hist  = nu_rate_hist;

out.truthAux      = truthAux;
out.meas          = meas;

% Prior-augmented outputs
out.t_prior        = t_prior;
out.xhist_prior    = xhist_prior;
out.Pfull_prior    = Pfull_prior;
out.rEstHist_prior = rEstHist_prior;
out.errECEF_prior  = errECEF_prior;
out.errENU_prior   = errENU_prior;
out.err3d_prior    = err3d_prior;

if cfg.plot.enable
    plot_enu_fastrealistic(t_prior, errENU_prior, Pfull_prior, [iDE iDN iDU]);
    plot_capture_and_measurements(t, trk, meas, satUsed, ...
        zPredDoppUsed, zPredRateUsed, sat, cfg);
end

end

%% ═══════════════════════ VALIDATION / DERIVED CFG ════════════════════════
function cfg = validate_and_derive_cfg(cfg)
require_fields(cfg, {'tleFile','duration_sec','dt','startUTC','stopUTC','startFromTLE', ...
    'maxTLEAgeDays','rx','fc_Hz','c','minElevDeg','coarseStep','minMaxElev', ...
    'maxSats','constellationSelectionWindow_s','useSGP4','allowJ2Fallback', ...
    'init','P0_clkBias','P0_clkDrift','sigma_pos_rw','clk','tropo','meas', ...
    'filter','capture','seedMeas','plot','verbose','height' }, 'cfg.');

require_fields(cfg.rx, {'lat_deg','lon_deg','alt_m','vel_ecef'}, 'cfg.rx.');
require_fields(cfg.init, {'errENU_m','sigENU_m'}, 'cfg.init.');

require_fields(cfg.clk, {'truth','filt','filter_init_bias_mps','filter_init_drift_mps2'}, 'cfg.clk.');
require_fields(cfg.clk.truth, {'initBias_mps','initDrift_mps2', ...
    'sigmaBias_rw_mps_sqrt','sigmaDrift_rw_mps2_sqrt'}, 'cfg.clk.truth.');
require_fields(cfg.clk.filt, {'sigmaBias_rw_mps_sqrt','sigmaDrift_rw_mps2_sqrt'}, 'cfg.clk.filt.');

require_fields(cfg.tropo, {'enable','ztd_truth_m','ztd_model_m','el_floor_deg','map_epsilon_deg'}, 'cfg.tropo.');
% ── Backward-compatible tropo model defaults (added for framework) ──
if ~isfield(cfg.tropo,'model_truth'),  cfg.tropo.model_truth  = 'simple'; end
if ~isfield(cfg.tropo,'model_filter'), cfg.tropo.model_filter = 'simple'; end
if ~isfield(cfg.tropo,'lat_deg'),      cfg.tropo.lat_deg = cfg.rx.lat_deg; end
if ~isfield(cfg.tropo,'alt_m'),        cfg.tropo.alt_m   = cfg.rx.alt_m;  end
if ~isfield(cfg.tropo,'doy'),          cfg.tropo.doy     = 180;           end
if ~isfield(cfg.tropo,'P0_hPa'),       cfg.tropo.P0_hPa  = 1013.25;       end
if ~isfield(cfg.tropo,'T0_K'),         cfg.tropo.T0_K    = 288.15;        end
if ~isfield(cfg.tropo,'e0_hPa'),       cfg.tropo.e0_hPa  = 11.0;          end
require_fields(cfg.meas, {'enable_lighttime','colored_tau_s','sigma_col0_mps','Rinflate', ...
    'useDoppler','useDopplerRate','fdJacStep_m'}, 'cfg.meas.');
require_fields(cfg.filter, {'useNISGate','nisGate','maxIEKFIter','iekfTol', ...
    'maxStep_m','startupInflateEpochs','startupRmult'}, 'cfg.filter.');
require_fields(cfg.capture, {'mode','minDwell_s','switchBlackout_s','elev_margin_deg', ...
    'useServiceAngleFilter','maxServiceOffNadirDeg','switchThreshold_lowElev', ...
    'switchThreshold_highElev','highElevThreshold'}, 'cfg.capture.');
% Backward-compatible defaults for new fields (safe for all existing scripts)
if ~isfield(cfg.capture,'maxSwitches'),   cfg.capture.maxSwitches   = Inf; end
if ~isfield(cfg.capture,'minRaanSepDeg'), cfg.capture.minRaanSepDeg = 0;   end
require_fields(cfg.plot, {'enable'}, 'cfg.plot.');
require_fields(cfg.height, {'enable','sigma_m','z0_m'}, 'cfg.height.');

cfg.lambda = cfg.c / cfg.fc_Hz;

if ~isfield(cfg,'sigma_doppler_Hz')
    error('Missing cfg.sigma_doppler_Hz');
end
if ~isfield(cfg,'sigma_dopplerRate_Hzps')
    error('Missing cfg.sigma_dopplerRate_Hzps');
end

cfg.meas.sigma_rhodot0_mps    = cfg.lambda * cfg.sigma_doppler_Hz;
cfg.meas.sigma_rhoddot0_mps2  = cfg.lambda * cfg.sigma_dopplerRate_Hzps;

if ~isfield(cfg.rx,'acc_ecef') || isempty(cfg.rx.acc_ecef)
    cfg.rx.acc_ecef = [0;0;0];
end
end

function trk = select_capture_stream_online(vis, elev, ~, sat_raan, capture, dt)
% SELECT_CAPTURE_STREAM_ONLINE
%   Sequential satellite tracking with three policies:
%     'stickyBest'     — always follow the highest-elevation visible sat
%     'stickyDiverse'  — always pick the sat with largest RAAN separation
%                        from the current one (no fixed threshold — works
%                        for any plane spacing including OneWeb ~7-20 deg)
%     'stickyPingPong' — lock the globally best-RAAN-separated pair at
%                        first acquisition and alternate strictly A->B->A
%
%   All policies honour minDwell_s, switchBlackout_s and maxSwitches.
%
%   INPUTS
%     vis       [K x Nsat] logical  visibility
%     elev      [K x Nsat] double   elevation [deg]
%     az        [K x Nsat] double   azimuth   [deg, 0-360]  (unused now,
%                                   kept in signature for drop-in compat.)
%     sat_raan  [Nsat x 1] double   RAAN of each satellite [deg]
%     capture   struct  (see validate_and_derive_cfg for required fields)
%     dt        double   time step [s]

[K, ~] = size(vis);
trk.sid      = zeros(K,1);
trk.blackout = false(K,1);
trk.online   = true;

if strcmpi(capture.mode, 'allVisible')
    return;
end

% ── Timing ───────────────────────────────────────────────────────────────
dwellEp    = max(1, round(capture.minDwell_s    / dt));
blackoutEp = max(0, round(capture.switchBlackout_s / dt));

% ── Elevation thresholds ─────────────────────────────────────────────────
elev_margin   = capture.elev_margin_deg;
thresh_low    = capture.switchThreshold_lowElev;
thresh_high   = capture.switchThreshold_highElev;
high_elev_thr = capture.highElevThreshold;

% ── Switch budget (fairness cap) ─────────────────────────────────────────
if isfield(capture,'maxSwitches') && ~isempty(capture.maxSwitches)
    maxSW = capture.maxSwitches;
else
    maxSW = Inf;
end
nSwitches = 0;

% ── State ────────────────────────────────────────────────────────────────
cur           = 0;
dwell         = 0;
blackout      = 0;
switchHistory = [];

% ── PingPong-specific state ───────────────────────────────────────────────
pp_satA   = 0;
pp_satB   = 0;
pp_target = 0;
pp_locked = false;

% ── RAAN angular distance helper (wrapped to [0,180]) ────────────────────
raan_dist_fn = @(a,b) min( mod(abs(sat_raan(a)-sat_raan(b)),360), ...
                            360 - mod(abs(sat_raan(a)-sat_raan(b)),360) );

% ── find_best_pair inline helper (used by stickyPingPong)
%    Given a candidate set, returns sA = highest-elev sat and sB = the
%    candidate with the largest RAAN separation from sA.
%    Call as:  [sA, sB] = p8_best_pair(cands, k, elev, thresh_low, raan_dist_fn)
%    (inlined below at each use site to avoid nested-function syntax)

% ═════════════════════════════════════════════════════════════════════════
%  MAIN LOOP
% ═════════════════════════════════════════════════════════════════════════
for k = 1:K

    % Blackout period — no tracking
    if blackout > 0
        trk.blackout(k) = true;
        blackout = blackout - 1;
        continue;
    end

    visS         = find(vis(k,:));
    currentValid = cur > 0 && ismember(cur,visS) && elev(k,cur) >= elev_margin;

    needSwitch = false;
    newSat     = 0;

    % ── Switch triggers ───────────────────────────────────────────────────
    if cur == 0
        needSwitch = true;

    elseif ~currentValid
        needSwitch = true;                   % dropped below horizon/margin

    elseif dwell >= dwellEp && nSwitches < maxSW

        switch lower(capture.mode)

            % ── stickyBest ───────────────────────────────────────────────
            case 'stickybest'
                if ~isempty(visS)
                    [maxEl, maxIdx] = max(elev(k, visS));
                    bestS     = visS(maxIdx);
                    elevGain  = maxEl - elev(k, cur);
                    thr       = thresh_low + (thresh_high - thresh_low) * ...
                                (elev(k,cur) >= high_elev_thr);
                    recentSw  = ~isempty(switchHistory) && any(switchHistory == bestS);
                    if ~recentSw && elevGain > thr && bestS ~= cur
                        needSwitch = true;
                        newSat     = bestS;
                    end
                end

            % ── stickyDiverse ────────────────────────────────────────────
            %    Rank all visible candidates by RAAN distance from current.
            %    No fixed threshold — the maximum available diversity always
            %    wins, so this works for any plane spacing (7 deg OneWeb or
            %    90 deg Walker).
            case 'stickydiverse'
                if ~isempty(visS)
                    elev_ok = elev(k, visS) >= thresh_low;
                    not_cur = visS ~= cur;
                    not_rec = ~ismember(visS, switchHistory);
                    cands   = visS(elev_ok & not_cur & not_rec);

                    if ~isempty(cands)
                        rd = arrayfun(@(s) raan_dist_fn(cur,s), cands);
                        [best_rd, bi] = max(rd);
                        % Require at least 5 deg RAAN sep to avoid same-pass churn
                        if best_rd >= 5
                            needSwitch = true;
                            newSat     = cands(bi);
                        end
                    end
                end

            % ── stickyPingPong ───────────────────────────────────────────
            %    Lock the globally best-RAAN-separated pair on first
            %    acquisition. Then alternate strictly satA -> satB -> satA.
            %    If the scheduled target is below the horizon, re-lock the
            %    best available replacement.
            case 'stickypingpong'
                if ~isempty(visS)

                    % Initial pair lock (inlined find_best_pair)
                    if ~pp_locked
                        ok_cands = visS(elev(k, visS) >= thresh_low);
                        % sA = highest elevation candidate
                        pp_satA = 0; pp_satB = 0;
                        if ~isempty(ok_cands)
                            [~,ii] = max(elev(k, ok_cands));
                            pp_satA = ok_cands(ii);
                            if numel(ok_cands) >= 2
                                rd_pair = arrayfun(@(s) raan_dist_fn(pp_satA,s), ok_cands);
                                rd_pair(ok_cands == pp_satA) = -1;
                                [~,jj] = max(rd_pair);
                                pp_satB = ok_cands(jj);
                            end
                        end
                        if pp_satA > 0
                            pp_locked = true;
                            pp_target = pp_satA;
                            if pp_satB == 0, pp_satB = pp_satA; end
                            fprintf('  [PingPong] Locked: Sat%d (RAAN=%.1f) <-> Sat%d (RAAN=%.1f)  sep=%.1f deg\n', ...
                                pp_satA, sat_raan(pp_satA), ...
                                pp_satB, sat_raan(pp_satB), ...
                                raan_dist_fn(pp_satA, pp_satB));
                        end
                    end

                    % Attempt scheduled ping-pong switch
                    if pp_locked && pp_target > 0 && pp_target ~= cur
                        if ismember(pp_target, visS) && elev(k, pp_target) >= thresh_low
                            % Target available — execute switch
                            needSwitch = true;
                            newSat     = pp_target;
                        else
                            % Target not visible — re-lock best replacement
                            elev_ok = elev(k, visS) >= thresh_low;
                            not_cur = visS ~= cur;
                            cands   = visS(elev_ok & not_cur);
                            if ~isempty(cands)
                                rd = arrayfun(@(s) raan_dist_fn(cur,s), cands);
                                [~, bi]     = max(rd);
                                replacement = cands(bi);
                                if cur == pp_satA
                                    pp_satB = replacement;
                                else
                                    pp_satA = replacement;
                                end
                                pp_target  = replacement;
                                needSwitch = true;
                                newSat     = replacement;
                                fprintf('  [PingPong] Re-locked: Sat%d (RAAN=%.1f)  sep=%.1f deg\n', ...
                                    replacement, sat_raan(replacement), ...
                                    raan_dist_fn(cur, replacement));
                            end
                        end
                    end
                end

        end % switch
    end % trigger block

    % ── Execute switch ────────────────────────────────────────────────────
    if needSwitch
        if isempty(visS)
            cur   = 0;
            dwell = 0;

        else
            % Fallback: no target chosen yet — pick highest elevation
            if newSat == 0
                [~, ii] = max(elev(k, visS));
                newSat  = visS(ii);

                % PingPong initial acquisition (inlined find_best_pair)
                if strcmpi(capture.mode,'stickyPingPong') && ~pp_locked
                    ok_cands = visS(elev(k, visS) >= thresh_low);
                    pp_satA = 0; pp_satB = 0;
                    if ~isempty(ok_cands)
                        [~,ii2] = max(elev(k, ok_cands));
                        pp_satA = ok_cands(ii2);
                        if numel(ok_cands) >= 2
                            rd_pair = arrayfun(@(s) raan_dist_fn(pp_satA,s), ok_cands);
                            rd_pair(ok_cands == pp_satA) = -1;
                            [~,jj2] = max(rd_pair);
                            pp_satB = ok_cands(jj2);
                        end
                    end
                    if pp_satA > 0
                        pp_locked = true;
                        newSat    = pp_satA;
                        pp_target = pp_satB;
                        if pp_satB == 0, pp_satB = pp_satA; end
                        fprintf('  [PingPong] Initial lock: Sat%d <-> Sat%d  sep=%.1f deg\n', ...
                            pp_satA, pp_satB, raan_dist_fn(pp_satA, pp_satB));
                    end
                end
            end

            % Update PingPong next-target pointer before executing switch
            if strcmpi(capture.mode,'stickyPingPong') && pp_locked && newSat ~= cur
                if newSat == pp_satA
                    pp_target = pp_satB;
                elseif newSat == pp_satB
                    pp_target = pp_satA;
                else
                    if cur == pp_satA
                        pp_target = pp_satB;
                    else
                        pp_target = pp_satA;
                    end
                end
            end

            % Apply blackout if departing a tracked satellite
            if cur > 0 && newSat ~= cur && blackoutEp > 0
                trk.blackout(k) = true;
                blackout        = blackoutEp - 1;
                switchHistory(end+1) = cur; %#ok<AGROW>
                if numel(switchHistory) > 4
                    switchHistory = switchHistory(end-3:end);
                end
                cur       = newSat;
                dwell     = 0;
                nSwitches = nSwitches + 1;
                continue;
            end

            % Instant switch (first acquisition or zero-blackout)
            if cur > 0 && newSat ~= cur
                switchHistory(end+1) = cur; %#ok<AGROW>
                if numel(switchHistory) > 4
                    switchHistory = switchHistory(end-3:end);
                end
                nSwitches = nSwitches + 1;
            end
            cur   = newSat;
            dwell = 0;
        end
    end

    % ── Log tracking epoch ────────────────────────────────────────────────
    if cur > 0 && ismember(cur, visS)
        trk.sid(k) = cur;
        dwell = dwell + 1;
    else
        cur   = 0;
        dwell = 0;
    end

end % k loop
end% (no trailing 'end' — file uses implicit function termination convention)
%% ═══════════════════════ MEASUREMENT SIMULATION ══════════════════════════
