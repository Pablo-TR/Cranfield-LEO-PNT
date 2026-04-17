function [p, v, ok, used] = propagate_tle_ecef_fast(tle, utc, cfg)
%PROPAGATE_TLE_ECEF_FAST Propagate one TLE to ECEF position and velocity.

ok = false;
p = [];
v = [];
used = 'J2 fallback';

if cfg.useSGP4 && ~isempty(tle.satrec)
    try
        tsince = minutes(utc - tle.epochUTC);
        [~, rkm, vkmps] = sgp4(tle.satrec, tsince);
        if numel(rkm) == 3 && all(isfinite(rkm)) && norm(rkm) > 6000 && norm(rkm) < 50000
            [p, v] = eci2ecef(1e3*rkm(:), 1e3*vkmps(:), utc);
            ok = true;
            used = 'SGP4';
            return;
        end
    catch ME
        if ~cfg.allowJ2Fallback
            rethrow(ME);
        end
    end
end

[pi_, vi_] = propagate_j2(tle, utc);
if ~isempty(pi_)
    [p, v] = eci2ecef(pi_, vi_, utc);
    ok = true;
end
end
