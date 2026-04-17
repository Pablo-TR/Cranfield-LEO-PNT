function m = ekf_extract_metrics(out, cfg)
% EKF_EXTRACT_METRICS   Extract full performance metrics from one solver run.
%
% USAGE
%   m = ekf_extract_metrics(out, cfg)
%
% INPUTS
%   out : output struct from real_strtg_sequential_EKF_fastrealistic
%   cfg : config struct used for that run (needs cfg.filter, cfg.analysis.conv)
%
% OUTPUT
%   m   : struct with all performance metrics (see fields below)
%
% METRIC FIELDS
%   Position/error:     finalEast_m, finalNorth_m, finalHoriz_m,
%                       minHoriz_m, RMSE_E_m, RMSE_N_m, RMSE_h_m
%   Clock:              clkBiasRMSE_mps, clkDriftRMSE_mps2
%   Filter consistency: meanNIS, medianNIS, stdNIS, fractionRejectedByGate,
%                       fracSigACF_dopp, fracSigACF_rate,
%                       excessKurtosis_dopp, excessKurtosis_rate
%   Convergence:        convTime_horiz_s, convTime_east_s, convTime_north_s
%   Run status:         success, failMsg

%% ── Input checks ─────────────────────────────────────────────────────────
if ~isfield(cfg, 'analysis') || ~isfield(cfg.analysis, 'conv')
    cfg.analysis.conv = default_conv_params();
end
conv = cfg.analysis.conv;

t  = out.t(:);
K  = numel(t);

%% ── ENU errors ───────────────────────────────────────────────────────────
eE = out.errENU(1,:).';
eN = out.errENU(2,:).';
eH = sqrt(eE.^2 + eN.^2);

validE = ~isnan(eE);
validN = ~isnan(eN);
validH = ~isnan(eH);

%% ── Position metrics ─────────────────────────────────────────────────────
% Use last valid sample for "final" values
lastH = find(validH, 1, 'last');
lastE = find(validE, 1, 'last');
lastN = find(validN, 1, 'last');

if ~isempty(lastH)
    m.finalEast_m  = eE(lastE);
    m.finalNorth_m = eN(lastN);
    m.finalHoriz_m = eH(lastH);
    m.minHoriz_m   = min(eH(validH));
    m.RMSE_E_m     = sqrt(mean(eE(validE).^2));
    m.RMSE_N_m     = sqrt(mean(eN(validN).^2));
    m.RMSE_h_m     = sqrt(mean(eH(validH).^2));
else
    m.finalEast_m  = NaN;  m.finalNorth_m = NaN;  m.finalHoriz_m = NaN;
    m.minHoriz_m   = NaN;  m.RMSE_E_m = NaN;
    m.RMSE_N_m     = NaN;  m.RMSE_h_m = NaN;
end

%% ── Clock metrics ────────────────────────────────────────────────────────
clkBiasEst   = out.xhist(4,:).';
clkDriftEst  = out.xhist(5,:).';
clkBiasTrue  = out.truthAux.clkBiasTrue(:);
clkDriftTrue = out.truthAux.clkDriftTrue(:);

% Align lengths (should match but guard against off-by-one)
Lb = min(numel(clkBiasEst), numel(clkBiasTrue));
Ld = min(numel(clkDriftEst), numel(clkDriftTrue));

vB = ~isnan(clkBiasEst(1:Lb))  & ~isnan(clkBiasTrue(1:Lb));
vD = ~isnan(clkDriftEst(1:Ld)) & ~isnan(clkDriftTrue(1:Ld));

if any(vB)
    m.clkBiasRMSE_mps  = sqrt(mean((clkBiasEst(1:Lb)  - clkBiasTrue(1:Lb)).^2));
else
    m.clkBiasRMSE_mps  = NaN;
end
if any(vD)
    m.clkDriftRMSE_mps2 = sqrt(mean((clkDriftEst(1:Ld) - clkDriftTrue(1:Ld)).^2));
else
    m.clkDriftRMSE_mps2 = NaN;
end

%% ── NIS metrics ──────────────────────────────────────────────────────────
nis   = out.nis(:);
valid = ~isnan(nis) & nis > 0;
nis_v = nis(valid);

if numel(nis_v) >= 2
    m.meanNIS    = mean(nis_v);
    m.medianNIS  = median(nis_v);
    m.stdNIS     = std(nis_v);
    m.nNIS       = numel(nis_v);
    m.fractionRejectedByGate = mean(nis_v > cfg.filter.nisGate);
else
    m.meanNIS    = NaN;
    m.medianNIS  = NaN;
    m.stdNIS     = NaN;
    m.nNIS       = 0;
    m.fractionRejectedByGate = NaN;
end

%% ── Innovation whiteness & Gaussianity ───────────────────────────────────
[m.fracSigACF_dopp, m.excessKurtosis_dopp] = innovation_stats(out.nu_hist(:));
[m.fracSigACF_rate, m.excessKurtosis_rate] = innovation_stats(out.nu_rate_hist(:));

%% ── Convergence metrics ──────────────────────────────────────────────────
m.convTime_horiz_s = compute_settling_time(t(validH), eH(validH), conv);
m.convTime_east_s  = compute_settling_time(t(validE), abs(eE(validE)), conv);
m.convTime_north_s = compute_settling_time(t(validN), abs(eN(validN)), conv);

%% ── Measurement quality ──────────────────────────────────────────────────
m.nMeasEpochs = sum(out.nUpdates > 0);
m.fracMeasEpochs = m.nMeasEpochs / K;

%% ── Run status ───────────────────────────────────────────────────────────
m.success = isfinite(m.finalHoriz_m) && m.finalHoriz_m < 1e6 && m.nMeasEpochs > 0;
m.failMsg = '';
if ~m.success
    if m.nMeasEpochs == 0
        m.failMsg = 'No measurement epochs processed';
    else
        m.failMsg = sprintf('finalHoriz=%.1f m', m.finalHoriz_m);
    end
end
end

% ═══════════════════════════════════════════════════════════════════════════
function [fracSig, exKurt] = innovation_stats(nu_raw)
% Compute ACF whiteness fraction and excess kurtosis for innovation sequence.
valid = ~isnan(nu_raw) & nu_raw ~= 0;
nu_v  = nu_raw(valid);

if numel(nu_v) < 20
    fracSig = NaN;
    exKurt  = NaN;
    return;
end

sigma = std(nu_v);
if sigma < eps
    fracSig = NaN;
    exKurt  = NaN;
    return;
end

nu_norm = (nu_v - mean(nu_v)) / sigma;

maxLag = min(50, floor(numel(nu_norm)/4));
[acf, lags] = acf_manual(nu_norm, maxLag);
conf_bound = 1.96 / sqrt(numel(nu_norm));
fracSig    = mean(abs(acf(lags ~= 0)) > conf_bound);

% Excess kurtosis: kurt = E[(x-mu)^4] / sigma^4 - 3  (Gaussian → 0)
n    = numel(nu_norm);
m4   = sum(nu_norm.^4) / n;
exKurt = m4 - 3;   % nu_norm already has mean≈0, std≈1
end

% ═══════════════════════════════════════════════════════════════════════════
function tc = compute_settling_time(t, e, conv)
% Robust settling-time metric based on terminal steady-state behaviour.
%
% Algorithm:
%   1. Compute e_ss = median(e over last refWindow_s)
%   2. Tolerance = max(absTol_m, relTol * max(e_ss, errFloor_m))
%   3. For each candidate time t_k, inspect window [t_k, t_k+settleWindow_s]
%   4. Accept if all hold in that window:
%        (a) |mean(e_win) - e_ss| <= tol
%        (b) std(e_win)/max(mean(e_win), errFloor_m) <= maxCoeffVar
%        (c) |linear_slope(e_win)| <= maxSlope_mps
%   5. Return first such t_k; if none, return Inf.

tc = Inf;
t  = t(:);
e  = e(:);
N  = numel(t);

if N < 4
    return;
end

dt = mean(diff(t));
if dt <= 0, return; end

% Step 1: terminal reference
nRef = max(1, round(conv.refWindow_s / dt));
e_ss = median(e(max(1, N-nRef+1):N));

% Step 2: tolerance
tol = max(conv.absTol_m, conv.relTol * max(e_ss, conv.errFloor_m));

% Step 3-5: scan
nSettle = max(2, round(conv.settleWindow_s / dt));

for k = 1:N-1
    i_end = min(N, k + nSettle - 1);
    if i_end - k < 1
        continue;
    end
    ew = e(k:i_end);
    tw = t(k:i_end) - t(k);

    % (a) mean proximity to steady-state
    if abs(mean(ew) - e_ss) > tol
        continue;
    end

    % (b) coefficient of variation
    cv = std(ew) / max(mean(ew), conv.errFloor_m);
    if cv > conv.maxCoeffVar
        continue;
    end

    % (c) slope
    if numel(tw) >= 2
        p = polyfit(tw, ew, 1);
        slope = abs(p(1));
    else
        slope = 0;
    end
    if slope > conv.maxSlope_mps
        continue;
    end

    tc = t(k);
    return;
end
end

% ═══════════════════════════════════════════════════════════════════════════
function conv = default_conv_params()
% Convergence metric default parameters.
% These can be overridden via cfg.analysis.conv.*
conv.refWindow_s    = 15;   % terminal reference window [s]
conv.settleWindow_s = 10;   % settling inspection window [s]
conv.relTol         = 0.15; % 15% relative tolerance on e_ss
conv.absTol_m       = 5.0;  % 5 m absolute tolerance
conv.maxCoeffVar    = 0.25; % 25% max coefficient of variation
conv.maxSlope_mps   = 1.0;  % 1 m/s max linear trend in settling window
conv.errFloor_m     = 5.0;  % floor for normalisation
end
