%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                generateTleDopplerTruth.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 14th March 2026
% Last Update: 14th March 2026
% Version: 1.0.0
%
% Description: Standalone ground-truth generator that produces a Doppler
% curve from TLEs for one selected satellite over a user-defined time
% vector. This extracts the minimal reusable logic from the Estimation
% tooling so Tracking can reuse it with minimal surface-area changes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function truth = generateTleDopplerTruth(varargin)

    parser = inputParser;
    addParameter(parser, 'timeVector', []);
    addParameter(parser, 'tleFile', "");
    addParameter(parser, 'startUTC', []);
    addParameter(parser, 'startFromTLE', 'first');
    addParameter(parser, 'maxTLEAgeDays', 7);
    addParameter(parser, 'rxLatDeg', 52.06869);
    addParameter(parser, 'rxLonDeg', -0.60884);
    addParameter(parser, 'rxAltM', 100);
    addParameter(parser, 'receiverVelocityEcefMps', [0; 0; 0]);
    addParameter(parser, 'carrierFreqHz', 11.325e9);
    addParameter(parser, 'speedOfLight', 299792458);
    addParameter(parser, 'useSGP4', false);
    addParameter(parser, 'allowJ2Fallback', true);
    addParameter(parser, 'selectionMode', 'highestMaxElevation');
    addParameter(parser, 'satelliteNameContains', "");
    addParameter(parser, 'satelliteIndex', []);
    addParameter(parser, 'coarseStep', 10);
    addParameter(parser, 'minElevDeg', 0);
    addParameter(parser, 'minMaxElevDeg', 10);
    addParameter(parser, 'verbose', false);
    parse(parser, varargin{:});

    % Collect all parsed options into a single configuration structure.
    cfg = parser.Results;
    % Store the requested time axis as a column vector.
    timeVector = cfg.timeVector(:);

    if strlength(string(cfg.tleFile)) == 0
        error('generateTleDopplerTruth requires a tleFile.');
    end

    % Read and parse all TLE entries from the requested file.
    [tleData, tleMeta] = readTleFile(cfg.tleFile, ...
        'useSGP4', cfg.useSGP4);
    if isempty(tleData)
        error('No valid TLE entries were found in %s.', ...
            char(string(cfg.tleFile)));
    end

    % Preallocate the array that stores the epoch of each TLE entry.
    epochsUtc = NaT(numel(tleData), 1, 'TimeZone', 'UTC');
    for entryIndex = 1:numel(tleData)
        % Copy the epoch of each parsed TLE into a flat datetime array.
        epochsUtc(entryIndex) = tleData(entryIndex).epochUTC;
    end

    % Mark which TLE entries have a valid epoch.
    validEpochMask = ~isnat(epochsUtc);
    % Mark which TLE entries contain usable propagated orbital data.
    validOrbitalMask = arrayfun(@(entry) ~isempty(entry.satrec_data), ...
        tleData).';
    % Keep only the entries that have both a valid epoch and orbital model.
    validOrbitalIdx = find(validEpochMask & validOrbitalMask);
    if isempty(validOrbitalIdx)
        error('No orbital data could be parsed from the TLE file.');
    end

    % Resolve the absolute UTC start time of the generated truth interval.
    startUtc = resolveStartUtc(cfg.startUTC, cfg.startFromTLE, ...
        epochsUtc(validOrbitalIdx));
    % Build the UTC time vector corresponding to the requested local time grid.
    timeUtc = startUtc + seconds(timeVector - timeVector(1));
    % Force the generated time stamps to remain in UTC.
    timeUtc.TimeZone = 'UTC';

    % Compute the age of each TLE with respect to the start of the scenario.
    tleAgeDays = abs(days(startUtc - epochsUtc));
    % Select the entries that are recent enough according to the requested age limit.
    freshIdx = find(validEpochMask & validOrbitalMask & ...
        (tleAgeDays <= cfg.maxTLEAgeDays));
    if isempty(freshIdx)
        % Fall back to all valid orbital entries when none are recent enough.
        freshIdx = validOrbitalIdx;
    end

    % Build the pool of candidate satellites according to the selection mode.
    candidateIdx = selectCandidatePool(tleData, validOrbitalIdx, ...
        freshIdx, cfg);
    % Evaluate the candidates and keep the best one for this receiver and time span.
    [selectedIdx, selectedMaxElevDeg, selectedPropagator] = ...
        selectBestSatellite(candidateIdx, tleData, timeUtc, cfg);

    % Convert the receiver geodetic position into ECEF coordinates.
    receiverEcefM = geodeticToEcef(cfg.rxLatDeg, cfg.rxLonDeg, ...
        cfg.rxAltM);
    % Store the receiver velocity as a column vector in ECEF coordinates.
    receiverVelocityEcefMps = cfg.receiverVelocityEcefMps(:);
    % Build the constant ECEF-to-ENU rotation for the receiver location.
    ecefToEnu = buildEcefToEnu(cfg.rxLatDeg, cfg.rxLonDeg);

    % Store the total number of simulated epochs.
    numEpochs = numel(timeVector);
    % Preallocate the satellite ECEF position history.
    satPositionEcefM = nan(3, numEpochs);
    % Preallocate the satellite ECEF velocity history.
    satVelocityEcefMps = nan(3, numEpochs);
    % Preallocate the elevation-angle history.
    elevationDeg = nan(numEpochs, 1);
    % Preallocate the line-of-sight range-rate history.
    rangeRateMps = nan(numEpochs, 1);

    for epochIndex = 1:numEpochs
        % Propagate the selected TLE to the current epoch.
        [positionEcefM, velocityEcefMps, propagationInfo] = ...
            propagateTleEcef(tleData(selectedIdx), ...
            timeUtc(epochIndex), ...
            'useSGP4', cfg.useSGP4, ...
            'allowJ2Fallback', cfg.allowJ2Fallback);

        if ~propagationInfo.ok
            error(['TLE propagation failed for selected satellite %s ' ...
                'at epoch %d.'], tleData(selectedIdx).name, epochIndex);
        end

        % Store the propagated satellite position at the current epoch.
        satPositionEcefM(:, epochIndex) = positionEcefM;
        % Store the propagated satellite velocity at the current epoch.
        satVelocityEcefMps(:, epochIndex) = velocityEcefMps;
        % Compute the satellite elevation seen from the receiver.
        elevationDeg(epochIndex) = computeElevation(positionEcefM, ...
            receiverEcefM, ecefToEnu);
        % Compute the radial relative velocity between satellite and receiver.
        rangeRateMps(epochIndex) = computeRangeRate(receiverEcefM, ...
            receiverVelocityEcefMps, positionEcefM, velocityEcefMps);
    end

    % Convert the line-of-sight range rate into Doppler frequency.
    fdHz = -(cfg.carrierFreqHz / cfg.speedOfLight) .* rangeRateMps;
    % Estimate the Doppler rate numerically from the Doppler history.
    fdDotHzPerSec = numericalDerivative(fdHz, timeVector);
    % Estimate the Doppler acceleration numerically from the Doppler-rate history.
    fdDdotHzPerSec2 = numericalDerivative(fdDotHzPerSec, timeVector);
    % Flag the epochs where the satellite is above the requested elevation mask.
    visibilityMask = elevationDeg >= cfg.minElevDeg;

    % Return the Doppler truth, geometric products, and metadata needed by tracking.
    truth = struct( ...
        'time', timeVector, ...
        'timeUTC', timeUtc, ...
        'fdHz', fdHz(:), ...
        'fdDotHzPerSec', fdDotHzPerSec(:), ...
        'fdDdotHzPerSec2', fdDdotHzPerSec2(:), ...
        'rangeRateMps', rangeRateMps(:), ...
        'elevationDeg', elevationDeg(:), ...
        'visibilityMask', visibilityMask(:), ...
        'satelliteName', string(tleData(selectedIdx).name), ...
        'satelliteIndex', selectedIdx, ...
        'selectedMaxElevationDeg', selectedMaxElevDeg, ...
        'propagatorUsed', string(selectedPropagator), ...
        'receiverEcefM', receiverEcefM, ...
        'receiverVelocityEcefMps', receiverVelocityEcefMps, ...
        'satPositionEcefM', satPositionEcefM, ...
        'satVelocityEcefMps', satVelocityEcefMps, ...
        'startUTC', startUtc, ...
        'selectedTleEpochUTC', tleData(selectedIdx).epochUTC, ...
        'tleAgeDaysAtStart', tleAgeDays(selectedIdx), ...
        'tleMeta', tleMeta);

    if logical(cfg.verbose)
        % Print a short summary of the selected satellite and propagation method.
        fprintf(['Selected TLE satellite %s (index %d), max elevation ' ...
            '%.2f deg, propagator %s.\n'], ...
            tleData(selectedIdx).name, selectedIdx, ...
            selectedMaxElevDeg, char(selectedPropagator));
    end

end

function isValid = isValidUtcInput(value)

    isValid = isempty(value) || isdatetime(value) || ischar(value) || ...
        isstring(value);

end

function startUtc = resolveStartUtc(startUtcInput, startFromTle, ...
    validEpochsUtc)

    if ~isempty(startUtcInput)
        if ~isdatetime(startUtcInput)
            % Convert the user-supplied start time into a UTC datetime.
            startUtc = datetime(startUtcInput, 'TimeZone', 'UTC');
        else
            % Reuse the supplied datetime start time.
            startUtc = startUtcInput;
            if isempty(startUtc.TimeZone)
                % Assign UTC when the input datetime has no timezone.
                startUtc.TimeZone = 'UTC';
            else
                % Normalise the input datetime to UTC when a timezone already exists.
                startUtc.TimeZone = 'UTC';
            end
        end
        return
    end

    if isempty(validEpochsUtc)
        error('No valid TLE epochs are available to initialise startUTC.');
    end

    switch lower(string(startFromTle))
        case "median"
            % Start from the median TLE epoch when requested.
            startUtc = median(validEpochsUtc);
        otherwise
            % By default start from the first valid TLE epoch.
            startUtc = validEpochsUtc(1);
    end

end

function candidateIdx = selectCandidatePool(tleData, validOrbitalIdx, ...
    freshIdx, cfg)

    selectionMode = lower(string(cfg.selectionMode));
    % Read the list of all satellite names for name-based selection.
    allNames = string({tleData.name});

    switch selectionMode
        case "highestmaxelevation"
            % Keep all fresh satellites and let the later stage choose the best one.
            candidateIdx = freshIdx;

        case "namecontains"
            % Read the substring used to filter candidate satellite names.
            nameToken = lower(string(cfg.satelliteNameContains));
            if strlength(nameToken) == 0
                error(['selectionMode ''nameContains'' requires ' ...
                    'satelliteNameContains.']);
            end
            % Keep only the valid satellites whose name contains the requested token.
            candidateIdx = validOrbitalIdx(contains( ...
                lower(allNames(validOrbitalIdx)), nameToken));

        case "index"
            if isempty(cfg.satelliteIndex)
                error('selectionMode ''index'' requires satelliteIndex.');
            end
            % Force the candidate set to the requested satellite index.
            candidateIdx = cfg.satelliteIndex;

        otherwise
            error('Unsupported selectionMode: %s', selectionMode);
    end

    % Convert the candidate list to a row vector.
    candidateIdx = candidateIdx(:).';
    % Remove out-of-range satellite indices from the candidate list.
    candidateIdx = candidateIdx(candidateIdx >= 1 & ...
        candidateIdx <= numel(tleData));

    if isempty(candidateIdx)
        error('No TLE candidates matched the requested selection.');
    end

end

function [selectedIdx, selectedMaxElevDeg, selectedPropagator] = ...
    selectBestSatellite(candidateIdx, tleData, timeUtc, cfg)

    % Convert the receiver position to ECEF for geometric calculations.
    receiverEcefM = geodeticToEcef(cfg.rxLatDeg, cfg.rxLonDeg, ...
        cfg.rxAltM);
    % Build the local ENU rotation used for elevation calculations.
    ecefToEnu = buildEcefToEnu(cfg.rxLatDeg, cfg.rxLonDeg);
    % Downsample the time grid to reduce the cost of the coarse candidate search.
    coarseIdx = unique([1:cfg.coarseStep:numel(timeUtc), numel(timeUtc)]);

    % Initialise the best accepted candidate.
    selectedIdx = [];
    selectedMaxElevDeg = -Inf;
    selectedPropagator = "";

    % Initialise a fallback candidate in case none meets the elevation threshold.
    fallbackIdx = [];
    fallbackMaxElevDeg = -Inf;
    fallbackPropagator = "";

    for currentIdx = candidateIdx
        % Reset the maximum elevation seen for the current candidate.
        maxElevDeg = -Inf;
        % Track whether propagation succeeded at least once for this candidate.
        hasPropagation = false;
        % Keep the propagator used on the last successful propagation.
        lastPropagator = "";

        for coarseSampleIdx = coarseIdx
            % Propagate the current candidate to a coarse time sample.
            [positionEcefM, ~, propagationInfo] = propagateTleEcef( ...
                tleData(currentIdx), timeUtc(coarseSampleIdx), ...
                'useSGP4', cfg.useSGP4, ...
                'allowJ2Fallback', cfg.allowJ2Fallback);

            if ~propagationInfo.ok
                continue;
            end

            % Mark this candidate as successfully propagated.
            hasPropagation = true;
            % Store the propagator that produced the current sample.
            lastPropagator = propagationInfo.usedPropagator;
            % Update the best elevation observed for the current candidate.
            maxElevDeg = max(maxElevDeg, computeElevation( ...
                positionEcefM, receiverEcefM, ecefToEnu));
        end

        if ~hasPropagation
            continue;
        end

        if maxElevDeg > fallbackMaxElevDeg
            % Update the fallback candidate with the highest reachable elevation.
            fallbackIdx = currentIdx;
            fallbackMaxElevDeg = maxElevDeg;
            fallbackPropagator = lastPropagator;
        end

        if maxElevDeg >= cfg.minMaxElevDeg && maxElevDeg > selectedMaxElevDeg
            % Update the selected candidate when the elevation threshold is met.
            selectedIdx = currentIdx;
            selectedMaxElevDeg = maxElevDeg;
            selectedPropagator = lastPropagator;
        end
    end

    if isempty(selectedIdx)
        if isempty(fallbackIdx)
            error('No candidate satellite could be propagated from the TLEs.');
        end

        % Fall back to the candidate with the highest elevation even if it
        % does not meet the requested minimum-maximum-elevation threshold.
        selectedIdx = fallbackIdx;
        selectedMaxElevDeg = fallbackMaxElevDeg;
        selectedPropagator = fallbackPropagator;
    end

end

function elevationDeg = computeElevation(satPositionEcefM, ...
    receiverEcefM, ecefToEnu)

    % Express the satellite-receiver geometry in the local ENU frame.
    deltaEnu = ecefToEnu * (satPositionEcefM - receiverEcefM);
    % Convert the ENU geometry into an elevation angle in degrees.
    elevationDeg = rad2deg(atan2(deltaEnu(3), hypot(deltaEnu(1), ...
        deltaEnu(2))));

end

function rangeRateMps = computeRangeRate(receiverEcefM, ...
    receiverVelocityEcefMps, satPositionEcefM, satVelocityEcefMps)

    % Compute the receiver-to-satellite position vector.
    deltaPosition = satPositionEcefM - receiverEcefM;
    % Normalise it to obtain the line-of-sight direction.
    lineOfSight = deltaPosition / norm(deltaPosition);
    % Compute the satellite velocity relative to the receiver.
    deltaVelocity = satVelocityEcefMps - receiverVelocityEcefMps;
    % Project the relative velocity onto the line of sight to get range rate.
    rangeRateMps = deltaVelocity.' * lineOfSight;

end

function receiverEcefM = geodeticToEcef(latDeg, lonDeg, altM)

    % Define the WGS-84 Earth semi-major axis.
    semiMajorAxis = 6378137.0;
    % Define the WGS-84 squared eccentricity.
    eccentricitySquared = 0.00669437999014;
    % Convert latitude to radians.
    latRad = deg2rad(latDeg);
    % Convert longitude to radians.
    lonRad = deg2rad(lonDeg);
    % Compute the prime-vertical radius of curvature.
    primeVerticalRadius = semiMajorAxis / ...
        sqrt(1 - eccentricitySquared * sin(latRad)^2);

    % Convert geodetic latitude, longitude, and altitude into ECEF coordinates.
    receiverEcefM = [(primeVerticalRadius + altM) * cos(latRad) * cos(lonRad); ...
        (primeVerticalRadius + altM) * cos(latRad) * sin(lonRad); ...
        (primeVerticalRadius * (1 - eccentricitySquared) + altM) * ...
        sin(latRad)];

end

function ecefToEnu = buildEcefToEnu(latDeg, lonDeg)

    % Convert latitude to radians.
    latRad = deg2rad(latDeg);
    % Convert longitude to radians.
    lonRad = deg2rad(lonDeg);
    % Precompute the sine of latitude.
    sinLat = sin(latRad);
    % Precompute the cosine of latitude.
    cosLat = cos(latRad);
    % Precompute the sine of longitude.
    sinLon = sin(lonRad);
    % Precompute the cosine of longitude.
    cosLon = cos(lonRad);

    % Build the rotation matrix from ECEF to local ENU coordinates.
    ecefToEnu = [-sinLon, cosLon, 0; ...
        -sinLat * cosLon, -sinLat * sinLon, cosLat; ...
        cosLat * cosLon, cosLat * sinLon, sinLat];

end

function derivative = numericalDerivative(signalValues, timeVector)

    % Estimate the derivative numerically on the provided time grid.
    derivative = gradient(signalValues(:), timeVector(:));

end

