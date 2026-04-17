function [p, v] = propagate_j2(tle, tgt)
%PROPAGATE_J2 Lightweight J2 fallback propagation from parsed TLE elements.

sr = tle.satrec_data;
if isempty(sr)
    p = [];
    v = [];
    return;
end

mu = 3.986004418e14;
J2 = 1.08263e-3;
Re = 6378137;

dt_ = seconds(tgt - tle.epochUTC);

a   = sr.a;
e   = sr.ecc;
inc = deg2rad(sr.incl);
O0  = deg2rad(sr.raan);
w0  = deg2rad(sr.argp);
M0  = deg2rad(sr.M0);
n   = sr.n_rad_s;

pp = a*(1-e^2);
ci = cos(inc);

Od = -1.5*n*J2*(Re/pp)^2*ci;
wd =  0.75*n*J2*(Re/pp)^2*(5*ci^2 - 1);
Md =  n + 0.75*n*J2*(Re/pp)^2*sqrt(1-e^2)*(3*ci^2 - 1);

O = O0 + Od*dt_;
w = w0 + wd*dt_;
M = mod(M0 + Md*dt_, 2*pi);

E = M;
for it = 1:20
    dE = (M - E + e*sin(E)) / (1 - e*cos(E));
    E  = E + dE;
    if abs(dE) < 1e-14
        break;
    end
end

nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));
r  = a*(1-e*cos(E));

rp_ = [r*cos(nu); r*sin(nu); 0];
h   = sqrt(mu*pp);
vp_ = (mu/h) * [-sin(nu); e+cos(nu); 0];

cO = cos(O);
sO = sin(O);
cw = cos(w);
sw = sin(w);
ci_ = cos(inc);
si = sin(inc);

Rot = [cO*cw-sO*sw*ci_, -cO*sw-sO*cw*ci_,  sO*si; ...
       sO*cw+cO*sw*ci_, -sO*sw+cO*cw*ci_, -cO*si; ...
       sw*si,            cw*si,            ci_];

p = Rot*rp_;
v = Rot*vp_;

if any(isnan(p)) || norm(p) < 6e6 || norm(p) > 10e6
    p = [];
    v = [];
end
end
