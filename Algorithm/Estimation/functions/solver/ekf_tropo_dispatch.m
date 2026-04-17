function delay_m = ekf_tropo_dispatch(model, elev_deg, ztd_m, cfg_tropo)
% EKF_TROPO_DISPATCH  Route to the requested troposphere delay model.
%
% INPUTS
%   model      : 'simple' | 'saastamoinen' | 'niell'
%   elev_deg   : elevation angle matrix [K x Nsat] or scalar  [deg]
%   ztd_m      : zenith total delay to use for simple model / override [m]
%                (ignored for saastamoinen/niell which compute ZTD from met params)
%   cfg_tropo  : struct with fields:
%       .el_floor_deg     elevation floor [deg]
%       .map_epsilon_deg  regularisation for simple model [deg]
%       .lat_deg          receiver geodetic latitude [deg]
%       .alt_m            receiver altitude above ellipsoid [m]
%       .doy              day of year (1-365) for Niell seasonal correction
%       .P0_hPa           surface pressure [hPa]
%       .T0_K             surface temperature [K]
%       .e0_hPa           water vapour partial pressure [hPa]
%
% OUTPUT
%   delay_m    : slant tropospheric delay [m], same size as elev_deg
%
% MODELS
%   'simple'       : ZTD / sin(el + eps)  — current baseline
%   'saastamoinen' : Saastamoinen (1972) ZHD + ZWD, 1/sin mapping
%   'niell'        : Niell Mapping Function (Niell 1996) applied to ZHD/ZWD

switch lower(model)
    case 'simple'
        delay_m = tropo_simple(elev_deg, ztd_m, ...
            cfg_tropo.el_floor_deg, cfg_tropo.map_epsilon_deg);

    case 'saastamoinen'
        delay_m = tropo_saastamoinen(elev_deg, ...
            cfg_tropo.P0_hPa, cfg_tropo.T0_K, cfg_tropo.e0_hPa, ...
            cfg_tropo.lat_deg, cfg_tropo.alt_m, cfg_tropo.el_floor_deg);

    case 'niell'
        delay_m = tropo_niell(elev_deg, ...
            cfg_tropo.lat_deg, cfg_tropo.alt_m, cfg_tropo.doy, ...
            cfg_tropo.P0_hPa, cfg_tropo.T0_K, cfg_tropo.e0_hPa, ...
            cfg_tropo.el_floor_deg);

    otherwise
        warning('ekf_tropo_dispatch: unknown model "%s". Falling back to simple.', model);
        delay_m = tropo_simple(elev_deg, ztd_m, ...
            cfg_tropo.el_floor_deg, cfg_tropo.map_epsilon_deg);
end
end

% ═══════════════════════════════════════════════════════════════════════════
% MODEL 1: Simple 1/sin(el) mapping (current baseline)
% ═══════════════════════════════════════════════════════════════════════════
function delay_m = tropo_simple(elev_deg, ztd_m, el_floor, eps_deg)
% Obliquity factor: 1 / sin(el + eps)  with elevation floor
% ztd_m  : zenith total delay [m]
% eps_deg: small regularisation to avoid singularity near horizon
el      = max(elev_deg, el_floor);
delay_m = ztd_m ./ sind(el + eps_deg);
end

% ═══════════════════════════════════════════════════════════════════════════
% MODEL 2: Saastamoinen (1972) — met-param based ZHD/ZWD, 1/sin mapping
% ═══════════════════════════════════════════════════════════════════════════
function delay_m = tropo_saastamoinen(elev_deg, P_hPa, T_K, e_hPa, ...
                                       lat_deg, alt_m, el_floor_deg)
% Saastamoinen zenith delays:
%   ZHD = 0.0022768 * P / f(phi, H)
%   ZWD = 0.002277 * (1255/T + 0.05) * e  (Davis et al. approximation)
%
% Reference:
%   Saastamoinen (1972), IUGG Monograph 15, p.247
%   Davis et al. (1985), Radio Science 20(6):1593-1607
%
% Note: here we use a simple 1/sin(el) mapping to keep the slant formula
%       consistent. Use 'niell' for separate hydrostatic/wet mapping.

% Saastamoinen latitude/height correction factor
phi   = deg2rad(lat_deg);
H_km  = alt_m / 1e3;
fac   = 1 - 0.00266 * cos(2*phi) - 0.00028 * H_km;  % dimensionless

ZHD = 0.0022768 * P_hPa / fac;             % zenith hydrostatic delay [m]
ZWD = 0.002277 * (1255/T_K + 0.05) * e_hPa;  % zenith wet delay [m]
ZTD = ZHD + ZWD;

el      = max(elev_deg, el_floor_deg);
delay_m = ZTD ./ sind(el);
end

% ═══════════════════════════════════════════════════════════════════════════
% MODEL 3: Niell Mapping Function (NMF), Niell (1996)
% ═══════════════════════════════════════════════════════════════════════════
function delay_m = tropo_niell(elev_deg, lat_deg, alt_m, doy, ...
                                P_hPa, T_K, e_hPa, el_floor_deg)
% Applies the Niell (1996) hydrostatic and wet mapping functions.
%
% ZHD and ZWD are computed from surface met parameters using Saastamoinen,
% then each mapped to slant separately with the appropriate NMF function.
%
% Reference:
%   Niell (1996), J. Geophys. Res. 101(B2):3227-3246
%
% Coefficients from Table 1 (hydrostatic) and Table 2 (wet) of Niell (1996).

% ── ZHD and ZWD (Saastamoinen) ───────────────────────────────────────────
phi  = deg2rad(lat_deg);
H_km = alt_m / 1e3;
fac  = 1 - 0.00266 * cos(2*phi) - 0.00028 * H_km;

ZHD = 0.0022768 * P_hPa / fac;
ZWD = 0.002277  * (1255/T_K + 0.05) * e_hPa;

% ── Niell coefficients (interpolated by latitude) ─────────────────────────
lat_abs = abs(lat_deg);

% Hydrostatic mean and seasonal amplitude (Table 1, Niell 1996)
% Reference latitudes: 15, 30, 45, 60, 75 deg
lat_tab = [15, 30, 45, 60, 75];

ah_mean_tab = [1.2769934e-3, 1.2683230e-3, 1.2465397e-3, 1.2196049e-3, 1.2045996e-3];
bh_mean_tab = [2.9153695e-3, 2.9152299e-3, 2.9288445e-3, 2.9022565e-3, 2.9024205e-3];
ch_mean_tab = [62.610505e-3, 62.837393e-3, 63.721774e-3, 63.824265e-3, 64.258455e-3];

ah_amp_tab  = [0.0,          1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5];
bh_amp_tab  = [0.0,          2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 1.1723375e-4];
ch_amp_tab  = [0.0,          9.0128400e-5, 4.3497037e-4, 8.4795348e-4, 1.1207382e-3];

% Wet coefficients (no seasonal variation, Table 2, Niell 1996)
aw_tab = [5.8021897e-4, 5.6794847e-4, 5.8118019e-4, 5.9727542e-4, 6.1641693e-4];
bw_tab = [1.4275268e-3, 1.5138625e-3, 1.4572752e-3, 1.5007428e-3, 1.7599082e-3];
cw_tab = [4.3472961e-2, 4.6729510e-2, 4.3908931e-2, 4.4626982e-2, 5.4736038e-2];

% Interpolate
lat_c = max(15, min(75, lat_abs));

ah_m  = interp1(lat_tab, ah_mean_tab, lat_c);
bh_m  = interp1(lat_tab, bh_mean_tab, lat_c);
ch_m  = interp1(lat_tab, ch_mean_tab, lat_c);
ah_a  = interp1(lat_tab, ah_amp_tab,  lat_c);
bh_a  = interp1(lat_tab, bh_amp_tab,  lat_c);
ch_a  = interp1(lat_tab, ch_amp_tab,  lat_c);

aw    = interp1(lat_tab, aw_tab, lat_c);
bw    = interp1(lat_tab, bw_tab, lat_c);
cw    = interp1(lat_tab, cw_tab, lat_c);

% Seasonal phase: doy=28 is the reference day (winter min in NH)
% Southern hemisphere: shift by half year
if lat_deg < 0
    doy_eff = doy + 182.625;
else
    doy_eff = doy;
end
season_cos = cos(2*pi*(doy_eff - 28)/365.25);

ah = ah_m - ah_a * season_cos;
bh = bh_m - bh_a * season_cos;
ch = ch_m - ch_a * season_cos;

% Height correction for hydrostatic.
aht = 2.53e-5;
bht = 5.49e-3;
cht = 1.14e-3;
mh_corr_top    = (1 + aht/(1 + bht/(1 + cht)));
mh_corr_bottom = @(el) sind(el) + aht./(sind(el) + bht./(sind(el) + cht));
% Height factor: (1/sin(el) - m_ht) * H_km
% where m_ht = continued-fraction form with aht, bht, cht
m_ht = mh_corr_top ./ mh_corr_bottom(max(elev_deg, el_floor_deg));

% ── Apply mapping functions ───────────────────────────────────────────────
el = max(elev_deg, el_floor_deg);

% Hydrostatic mapping
top_h    = 1 + ah/(1 + bh/(1 + ch));
bot_h    = sind(el) + ah./(sind(el) + bh./(sind(el) + ch));
mh       = top_h ./ bot_h + (1./sind(el) - m_ht) * H_km;

% Wet mapping
top_w    = 1 + aw/(1 + bw/(1 + cw));
bot_w    = sind(el) + aw./(sind(el) + bw./(sind(el) + cw));
mw       = top_w ./ bot_w;

delay_m  = ZHD * mh + ZWD * mw;
end
