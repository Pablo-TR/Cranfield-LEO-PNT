function [ZTD, ZHD, ZWD] = ekf_tropo_ztd_from_met(P_hPa, T_K, e_hPa, lat_deg, alt_m)
%EKF_TROPO_ZTD_FROM_MET Compute Saastamoinen zenith delays from met data.

phi  = deg2rad(lat_deg);
H_km = alt_m / 1e3;
fac  = 1 - 0.00266 * cos(2*phi) - 0.00028 * H_km;
ZHD  = 0.0022768 * P_hPa / fac;
ZWD  = 0.002277 * (1255/T_K + 0.05) * e_hPa;
ZTD  = ZHD + ZWD;
end
