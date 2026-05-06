%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    propagateTleEcef.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 14th March 2026
% Last Update: 14th March 2026
% Version: 1.0.0
%
% Description: Standalone TLE propagator extracted from the Estimation
% tooling. It prefers SGP4 when available and falls back to a lightweight
% J2 propagator before converting the state to ECEF.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [positionEcefM, velocityEcefMps, propagationInfo] = ...
    propagateTleEcef(tleEntry, utcTime, varargin)

    parser = inputParser;
    addParameter(parser, 'useSGP4', false);
    addParameter(parser, 'allowJ2Fallback', true);
    parse(parser, varargin{:});

    utcTime = ensureUtcDatetime(utcTime);

    positionEcefM = [];
    velocityEcefMps = [];
    propagationInfo = struct( ...
        'ok', false, ...
        'usedPropagator', "", ...
        'usedSGP4', false);

    if logical(parser.Results.useSGP4) && isfield(tleEntry, 'satrec') && ...
            ~isempty(tleEntry.satrec)
        try
            tsinceMinutes = minutes(utcTime - tleEntry.epochUTC);
            [~, positionKm, velocityKmps] = sgp4(tleEntry.satrec, ...
                tsinceMinutes);

            if numel(positionKm) == 3 && all(isfinite(positionKm)) && ...
                    norm(positionKm) > 6000 && norm(positionKm) < 50000
                [positionEcefM, velocityEcefMps] = eciToEcef( ...
                    1e3 * positionKm(:), 1e3 * velocityKmps(:), utcTime);
                propagationInfo.ok = true;
                propagationInfo.usedPropagator = "SGP4";
                propagationInfo.usedSGP4 = true;
                return
            end
        catch currentError
            if ~logical(parser.Results.allowJ2Fallback)
                rethrow(currentError);
            end
        end
    end

    [positionEciM, velocityEciMps] = propagateJ2(tleEntry, utcTime);
    if isempty(positionEciM)
        return
    end

    [positionEcefM, velocityEcefMps] = eciToEcef(positionEciM, ...
        velocityEciMps, utcTime);
    propagationInfo.ok = true;
    propagationInfo.usedPropagator = "J2";

end

function utcTime = ensureUtcDatetime(utcTime)

    if ~isdatetime(utcTime)
        utcTime = datetime(utcTime, 'TimeZone', 'UTC');
    elseif isempty(utcTime.TimeZone)
        utcTime.TimeZone = 'UTC';
    else
        utcTime.TimeZone = 'UTC';
    end

end

function [positionEciM, velocityEciMps] = propagateJ2(tleEntry, utcTime)

    positionEciM = [];
    velocityEciMps = [];

    if ~isfield(tleEntry, 'satrec_data') || isempty(tleEntry.satrec_data)
        return
    end

    satrecData = tleEntry.satrec_data;
    mu = 3.986004418e14;
    j2 = 1.08263e-3;
    earthRadiusM = 6378137;

    deltaTimeSec = seconds(utcTime - tleEntry.epochUTC);

    semiMajorAxis = satrecData.a;
    eccentricity = satrecData.ecc;
    inclination = deg2rad(satrecData.incl);
    raan0 = deg2rad(satrecData.raan);
    argPerigee0 = deg2rad(satrecData.argp);
    meanAnomaly0 = deg2rad(satrecData.M0);
    meanMotionRadPerSec = satrecData.n_rad_s;

    semiLatusRectum = semiMajorAxis * (1 - eccentricity^2);
    cosInclination = cos(inclination);

    raanRate = -1.5 * meanMotionRadPerSec * j2 * ...
        (earthRadiusM / semiLatusRectum)^2 * cosInclination;
    argPerigeeRate = 0.75 * meanMotionRadPerSec * j2 * ...
        (earthRadiusM / semiLatusRectum)^2 * (5 * cosInclination^2 - 1);
    meanAnomalyRate = meanMotionRadPerSec + 0.75 * meanMotionRadPerSec * ...
        j2 * (earthRadiusM / semiLatusRectum)^2 * ...
        sqrt(1 - eccentricity^2) * (3 * cosInclination^2 - 1);

    raan = raan0 + raanRate * deltaTimeSec;
    argPerigee = argPerigee0 + argPerigeeRate * deltaTimeSec;
    meanAnomaly = mod(meanAnomaly0 + meanAnomalyRate * deltaTimeSec, ...
        2 * pi);

    eccentricAnomaly = meanAnomaly;
    for iteration = 1:20
        deltaE = (meanAnomaly - eccentricAnomaly + ...
            eccentricity * sin(eccentricAnomaly)) / ...
            (1 - eccentricity * cos(eccentricAnomaly));
        eccentricAnomaly = eccentricAnomaly + deltaE;
        if abs(deltaE) < 1e-14
            break;
        end
    end

    trueAnomaly = 2 * atan2( ...
        sqrt(1 + eccentricity) * sin(eccentricAnomaly / 2), ...
        sqrt(1 - eccentricity) * cos(eccentricAnomaly / 2));
    orbitalRadius = semiMajorAxis * ...
        (1 - eccentricity * cos(eccentricAnomaly));

    orbitalPosition = [orbitalRadius * cos(trueAnomaly); ...
        orbitalRadius * sin(trueAnomaly); 0];
    specificAngularMomentum = sqrt(mu * semiLatusRectum);
    orbitalVelocity = (mu / specificAngularMomentum) * ...
        [-sin(trueAnomaly); eccentricity + cos(trueAnomaly); 0];

    cosRaan = cos(raan);
    sinRaan = sin(raan);
    cosArgPerigee = cos(argPerigee);
    sinArgPerigee = sin(argPerigee);
    cosInc = cos(inclination);
    sinInc = sin(inclination);

    rotationMatrix = [ ...
        cosRaan * cosArgPerigee - sinRaan * sinArgPerigee * cosInc, ...
        -cosRaan * sinArgPerigee - sinRaan * cosArgPerigee * cosInc, ...
        sinRaan * sinInc; ...
        sinRaan * cosArgPerigee + cosRaan * sinArgPerigee * cosInc, ...
        -sinRaan * sinArgPerigee + cosRaan * cosArgPerigee * cosInc, ...
        -cosRaan * sinInc; ...
        sinArgPerigee * sinInc, ...
        cosArgPerigee * sinInc, ...
        cosInc];

    positionEciM = rotationMatrix * orbitalPosition;
    velocityEciMps = rotationMatrix * orbitalVelocity;

    if any(~isfinite(positionEciM)) || norm(positionEciM) < 6e6 || ...
            norm(positionEciM) > 1e8
        positionEciM = [];
        velocityEciMps = [];
    end

end

function [positionEcefM, velocityEcefMps] = eciToEcef(positionEciM, ...
    velocityEciMps, utcTime)

    earthRotationRate = 7.2921159e-5;
    julianDate = computeJulianDate(utcTime);
    centuriesSinceJ2000 = (julianDate - 2451545) / 36525;
    gmstDeg = 280.46061837 + ...
        360.98564736629 * (julianDate - 2451545) + ...
        0.000387933 * centuriesSinceJ2000^2 - ...
        centuriesSinceJ2000^3 / 38710000;
    theta = deg2rad(mod(gmstDeg, 360));

    cosTheta = cos(theta);
    sinTheta = sin(theta);

    rotationMatrix = [cosTheta, sinTheta, 0; ...
        -sinTheta, cosTheta, 0; ...
        0, 0, 1];

    rotationDerivative = earthRotationRate * ...
        [-sinTheta, cosTheta, 0; ...
         -cosTheta, -sinTheta, 0; ...
         0, 0, 0];

    positionEcefM = rotationMatrix * positionEciM;
    velocityEcefMps = rotationMatrix * velocityEciMps + ...
        rotationDerivative * positionEciM;

end

function julianDate = computeJulianDate(utcTime)

    julianDate = datenum(utcTime) + 1721058.5;

end

