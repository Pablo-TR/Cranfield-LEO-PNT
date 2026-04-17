function [x,P,accepted,nis,nu,h] = vector_iekf_update(x, P, z, Rk, ...
    modelFun, cfg, Inx, posIdx)

xBar = x;
PBar = P;
xi   = xBar;

accepted = true;
nis = NaN;
nu  = nan(size(z));
h   = nan(size(z));

for it = 1:cfg.filter.maxIEKFIter
    [h_i, H_i] = modelFun(xi);

    nu_i = z - h_i;
    S    = H_i * PBar * H_i' + Rk;

    if it == 1
        nis = nu_i' / S * nu_i;
        if cfg.filter.useNISGate && nis > cfg.filter.nisGate
            accepted = false;
            return;
        end
    end

    Kk    = PBar * H_i' / S;
    xiNew = xBar + Kk * (z - h_i - H_i * (xBar - xi));

    dpos = xiNew(posIdx) - xi(posIdx);
    nd   = norm(dpos);
    if nd > cfg.filter.maxStep_m
        xiNew(posIdx) = xi(posIdx) + dpos * (cfg.filter.maxStep_m / nd);
    end

    if norm(xiNew - xi) < cfg.filter.iekfTol
        xi = xiNew;
        break;
    end

    xi = xiNew;
end

[h, H] = modelFun(xi);
nu = z - h;
S  = H * PBar * H' + Rk;
Kk = PBar * H' / S;

x = xi;
P = (Inx - Kk*H) * PBar * (Inx - Kk*H)' + Kk*Rk*Kk';
P = (P + P') / 2;
end

%% ═══════════════════════ MEASUREMENT MODELS ══════════════════════════════
