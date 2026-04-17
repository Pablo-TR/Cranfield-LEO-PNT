function [sat, stats] = build_constellation_fast(tleData, freshIdx, tUTC, recTrue, R_e2n, cfg)
K    = numel(tUTC);

selectionWindow_s = cfg.constellationSelectionWindow_s;
selectionK = max(1, min(K, round(selectionWindow_s / cfg.dt)));
coarseIdx  = 1:cfg.coarseStep:selectionK;

stats.nCoarseReject = 0;
stats.nPropFail     = 0;


cands = struct('name',{},'pos',{},'vel',{},'elev',{},'maxElev',{}, ...
               'raan',{},'propUsed',{},'offNadirDeg',{});

for ii = 1:numel(freshIdx)
    i = freshIdx(ii);

    maxEl = -90;
    okAny = false;
    for ci = coarseIdx
        [p,~,ok,~] = propagate_tle_ecef_fast(tleData(i), tUTC(ci), cfg);
        if ~ok
            continue;
        end
        okAny = true;
        maxEl = max(maxEl, compute_elev(p, recTrue, R_e2n));
    end

    if ~okAny || maxEl < cfg.minMaxElev
        stats.nCoarseReject = stats.nCoarseReject + 1;
        continue;
    end

    pos = nan(3,1,1,K);
    vel = nan(3,1,1,K);
    elev = nan(K,1);
    offNadir = nan(K,1);
    allOk = true;
    used = '';

    for k = 1:K
        [p,v,ok1,used1] = propagate_tle_ecef_fast(tleData(i), tUTC(k), cfg);
        if ~ok1
            allOk = false;
            break;
        end
        pos(:,1,1,k) = p;
        vel(:,1,1,k) = v;
        elev(k)      = compute_elev(p, recTrue, R_e2n);
        offNadir(k)  = off_nadir_deg(p, recTrue);
        used         = used1;
    end

    if ~allOk
        stats.nPropFail = stats.nPropFail + 1;
        continue;
    end

    nc = numel(cands) + 1;
    cands(nc).name        = tleData(i).name;
    cands(nc).pos         = pos;
    cands(nc).vel         = vel;
    cands(nc).elev        = elev;
    cands(nc).maxElev     = maxEl;
    cands(nc).raan        = tleData(i).raan;
    cands(nc).propUsed    = used;
    cands(nc).offNadirDeg = offNadir;
end

if isempty(cands)
    sat = struct('Nsat',0);
    return;
end

[~,ord] = sort([cands.maxElev], 'descend');
cands   = cands(ord(1:min(cfg.maxSats, numel(cands))));
Nsat    = numel(cands);

sat.Nsat        = Nsat;
sat.pos         = zeros(3,1,Nsat,K);
sat.vel         = zeros(3,1,Nsat,K);
sat.elev        = zeros(K,Nsat);
sat.offNadirDeg = zeros(K,Nsat);
sat.name        = strings(Nsat,1);
sat.raan        = zeros(Nsat,1);

for s = 1:Nsat
    sat.pos(:,1,s,:)     = cands(s).pos;
    sat.vel(:,1,s,:)     = cands(s).vel;
    sat.elev(:,s)        = cands(s).elev;
    sat.offNadirDeg(:,s) = cands(s).offNadirDeg;
    sat.name(s)          = string(cands(s).name);
    sat.raan(s)          = cands(s).raan;
end
sat.propagatorUsed = cands(1).propUsed;

end

%% ═══════════════════════ CAPTURE POLICY ══════════════════════════════════
%% ═══════════════════════ CAPTURE POLICY ══════════════════════════════════
