%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                generateDopplerMeasurements.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Generates synthetic Doppler measurements by combining a
% configurable ground-truth model with measurement degradations such as
% Gaussian noise, Doppler-bin quantisation, neighbour-bin jumps, outliers,
% missing measurements, oscillator offset, and slow oscillator drift. The
% ground-truth model can also be driven directly from TLEs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = generateDopplerMeasurements(varargin)

    parser = inputParser;
    addParameter(parser, 'timeVector', []);
    addParameter(parser, 'numEpochs', 100);
    addParameter(parser, 'stepTime', 1e-3);
    addParameter(parser, 'truthModel', 'constant');
    addParameter(parser, 'truthParams', struct());
    addParameter(parser, 'customTruthFcn', []);
    addParameter(parser, 'gaussianNoiseStdHz', 0);
    addParameter(parser, 'dopplerBinWidthHz', 0);
    addParameter(parser, 'dopplerBinOriginHz', 0);
    addParameter(parser, 'neighbourJumpProbability', 0);
    addParameter(parser, 'outlierProbability', 0);
    addParameter(parser, 'outlierMeanHz', 0);
    addParameter(parser, 'outlierStdHz', 0);
    addParameter(parser, 'gapProbability', 0);
    addParameter(parser, 'gapValue', NaN);
    addParameter(parser, 'oscillatorOffsetHz', 0);
    addParameter(parser, 'oscillatorDriftModel', 'none');
    addParameter(parser, 'oscillatorDriftParams', struct());
    addParameter(parser, 'customOscillatorDriftFcn', []);
    addParameter(parser, 'randomSeed', []);
    parse(parser, varargin{:});

    % Collect all parsed options into a single configuration structure.
    cfg = parser.Results;

    if ~isempty(cfg.randomSeed)
        % Make the random draws repeatable when a seed is provided.
        rng(cfg.randomSeed);
    end

    % Build the time base either from the user-supplied vector or from the
    % requested number of epochs and sampling interval.
    timeVector = buildTimeVector(cfg.timeVector, cfg.numEpochs, ...
        cfg.stepTime);
    % Store the number of generated epochs for later array allocation.
    numEpochs = numel(timeVector);

    % Generate the underlying Doppler truth model selected by the user.
    [fdTrue, truthAux] = buildTruthSignal(timeVector, cfg.truthModel, ...
        cfg.truthParams, cfg.customTruthFcn);

    if isstruct(truthAux) && isfield(truthAux, 'fdDotHzPerSec')
        % Reuse the truth-model derivative if it is already available.
        fdDotTrue = truthAux.fdDotHzPerSec(:);
    else
        % Otherwise estimate the Doppler rate numerically from the truth.
        fdDotTrue = numericalDerivative(fdTrue, timeVector);
    end

    if isstruct(truthAux) && isfield(truthAux, 'fdDdotHzPerSec2')
        % Reuse the second derivative if the truth model provides it.
        fdDdotTrue = truthAux.fdDdotHzPerSec2(:);
    else
        % Otherwise estimate the Doppler acceleration numerically.
        fdDdotTrue = numericalDerivative(fdDotTrue, timeVector);
    end

    % Build the slow oscillator drift contribution requested by the user.
    oscillatorDrift = buildOscillatorDrift(timeVector, ...
        cfg.oscillatorDriftModel, cfg.oscillatorDriftParams, ...
        cfg.customOscillatorDriftFcn);
    % Expand the oscillator offset into one value per epoch.
    oscillatorOffset = cfg.oscillatorOffsetHz * ones(numEpochs, 1);

    % Combine the ideal truth, static offset, and slow drift.
    fdIdeal = fdTrue + oscillatorOffset + oscillatorDrift;

    % Draw zero-mean Gaussian measurement noise in Hz.
    gaussianNoise = cfg.gaussianNoiseStdHz * randn(numEpochs, 1);
    % Apply the Gaussian noise to the ideal Doppler track.
    fdAfterNoise = fdIdeal + gaussianNoise;

    if cfg.dopplerBinWidthHz > 0
        % Quantise the noisy Doppler to the closest Doppler-search bin.
        fdQuantised = cfg.dopplerBinOriginHz + ...
            round((fdAfterNoise - cfg.dopplerBinOriginHz) / ...
            cfg.dopplerBinWidthHz) * cfg.dopplerBinWidthHz;
    else
        % Leave the measurements unchanged when no bin quantisation is used.
        fdQuantised = fdAfterNoise;
    end

    % Initialise the neighbour-jump bookkeeping arrays.
    neighbourJumpMask = false(numEpochs, 1);
    neighbourJumpHz = zeros(numEpochs, 1);
    if cfg.dopplerBinWidthHz > 0 && cfg.neighbourJumpProbability > 0
        % Randomly decide which epochs jump to an adjacent Doppler bin.
        neighbourJumpMask = rand(numEpochs, 1) < ...
            cfg.neighbourJumpProbability;
        % Choose whether each selected jump goes to the upper or lower bin.
        jumpDirection = 2 * (rand(numEpochs, 1) > 0.5) - 1;
        % Convert the selected jump direction into a bin-width offset.
        neighbourJumpHz = neighbourJumpMask .* jumpDirection * ...
            cfg.dopplerBinWidthHz;
    end
    % Apply the neighbour-bin jumps to the quantised measurements.
    fdAfterNeighbourJump = fdQuantised + neighbourJumpHz;

    % Randomly mark which epochs will contain outliers.
    outlierMask = rand(numEpochs, 1) < cfg.outlierProbability;
    % Start from zero outlier error at every epoch.
    outlierError = zeros(numEpochs, 1);
    if cfg.outlierProbability > 0 && cfg.outlierStdHz > 0
        % Draw random outlier amplitudes when a stochastic outlier model is used.
        outlierError(outlierMask) = cfg.outlierMeanHz + ...
            cfg.outlierStdHz * randn(nnz(outlierMask), 1);
    elseif cfg.outlierProbability > 0
        % Use a fixed outlier amplitude when only an outlier mean is requested.
        outlierError(outlierMask) = cfg.outlierMeanHz;
    end
    % Add the outlier contribution to the previous measurement stage.
    fdAfterOutliers = fdAfterNeighbourJump + outlierError;

    % Randomly decide which epochs should be treated as missing measurements.
    gapMask = rand(numEpochs, 1) < cfg.gapProbability;
    % Start from the outlier-corrupted measurements.
    fdMeasured = fdAfterOutliers;
    % Replace missing epochs with the configured gap marker, typically NaN.
    fdMeasured(gapMask) = cfg.gapValue;

    % Return the truth, each intermediate perturbation stage, and the final
    % measured Doppler sequence to make downstream analysis easier.
    output = struct( ...
        'time', timeVector, ...
        'fdTrue', fdTrue, ...
        'fdDotTrue', fdDotTrue, ...
        'fdDdotTrue', fdDdotTrue, ...
        'oscillatorOffset', oscillatorOffset, ...
        'oscillatorDrift', oscillatorDrift, ...
        'fdIdeal', fdIdeal, ...
        'gaussianNoise', gaussianNoise, ...
        'fdAfterNoise', fdAfterNoise, ...
        'fdQuantised', fdQuantised, ...
        'neighbourJumpMask', neighbourJumpMask, ...
        'neighbourJumpHz', neighbourJumpHz, ...
        'fdAfterNeighbourJump', fdAfterNeighbourJump, ...
        'outlierMask', outlierMask, ...
        'outlierError', outlierError, ...
        'fdAfterOutliers', fdAfterOutliers, ...
        'gapMask', gapMask, ...
        'fdMeasured', fdMeasured, ...
        'truthAux', truthAux, ...
        'config', cfg);

end

function timeVector = buildTimeVector(inputTimeVector, numEpochs, stepTime)

    if isempty(inputTimeVector)
        % Create a uniformly spaced time axis starting at zero.
        timeVector = (0:numEpochs - 1).' * stepTime;
    else
        % Reformat the supplied time axis as a column vector.
        timeVector = inputTimeVector(:);
    end

end

function [fdTrue, truthAux] = buildTruthSignal(timeVector, truthModel, ...
    truthParams, customTruthFcn)

    truthModel = lower(string(truthModel));
    % Record which truth model has been selected.
    truthAux = struct('model', char(truthModel));

    switch truthModel
        case "constant"
            % Read the constant Doppler level.
            fd0 = getFieldOrDefault(truthParams, 'fd0', 0);
            % Replicate the constant value over the full time axis.
            fdTrue = fd0 * ones(size(timeVector));

        case "linear"
            % Read the initial Doppler value.
            fd0 = getFieldOrDefault(truthParams, 'fd0', 0);
            % Read the constant Doppler rate.
            fdDot = getFieldOrDefault(truthParams, 'fdDot', 0);
            % Shift time so the model starts at zero.
            time0 = timeVector - timeVector(1);
            % Build a first-order Doppler evolution.
            fdTrue = fd0 + fdDot .* time0;

        case "quadratic"
            % Read the initial Doppler value.
            fd0 = getFieldOrDefault(truthParams, 'fd0', 0);
            % Read the initial Doppler rate.
            fdDot = getFieldOrDefault(truthParams, 'fdDot', 0);
            % Read the constant Doppler acceleration.
            fdDdot = getFieldOrDefault(truthParams, 'fdDdot', 0);
            % Shift time so the polynomial starts at zero.
            time0 = timeVector - timeVector(1);
            % Build a second-order Doppler evolution.
            fdTrue = fd0 + fdDot .* time0 + ...
                0.5 * fdDdot .* time0.^2;

        case "sinusoidal"
            % Read the DC offset of the sinusoid.
            bias = getFieldOrDefault(truthParams, 'bias', 0);
            % Read the sinusoid amplitude.
            amplitude = getFieldOrDefault(truthParams, ...
                'amplitude', 80);
            % Read the oscillation frequency.
            frequency = getFieldOrDefault(truthParams, ...
                'frequencyHz', 20);
            % Read the initial phase.
            phase = getFieldOrDefault(truthParams, 'phaseRad', 0);
            % Generate a sinusoidal Doppler variation.
            fdTrue = bias + amplitude * ...
                sin(2 * pi * frequency * timeVector + phase);

        case "satellitepass"
            % Read the carrier frequency used to translate range rate into Doppler.
            fc = getFieldOrDefault(truthParams, 'carrierFreqHz', ...
                2.9e9);
            % Read the speed of light constant used in the Doppler conversion.
            c = getFieldOrDefault(truthParams, 'speedOfLight', 3e8);
            % Read the assumed satellite speed during the pass.
            velocity = getFieldOrDefault(truthParams, ...
                'satelliteVelocityMps', 7500);
            % Read the minimum slant range at closest approach.
            closestApproach = getFieldOrDefault(truthParams, ...
                'closestApproachM', 400e3);
            % Read the time corresponding to the closest approach.
            centreTime = getFieldOrDefault(truthParams, ...
                'centreTimeSec', mean(timeVector));
            % Read an optional constant Doppler bias term.
            bias = getFieldOrDefault(truthParams, 'bias', 0);

            % Shift the time axis so the pass is centred at closest approach.
            timeCentered = timeVector - centreTime;
            % Compute the slant range for a simple straight-line fly-by model.
            range = sqrt(closestApproach.^2 + ...
                (velocity .* timeCentered).^2);
            % Differentiate the range model analytically to get range rate.
            rangeRate = (velocity.^2 .* timeCentered) ./ range;
            % Convert range rate into Doppler frequency.
            fdTrue = -(fc / c) .* rangeRate + bias;

        case {"tle", "tlebased", "tlepass"}
            % Delegate the truth generation to the TLE-based standalone model.
            tleTruth = buildTleTruthSignal(timeVector, truthParams);
            % Extract the Doppler truth from the propagated TLE result.
            fdTrue = tleTruth.fdHz;
            % Keep the extra TLE metadata so downstream code can reuse it.
            truthAux = tleTruth;

        case "custom"
            if isempty(customTruthFcn)
                error(['truthModel is custom but customTruthFcn was ' ...
                    'not provided.']);
            end
            % Evaluate the user-supplied custom truth model.
            fdTrue = evaluateCustomModel(customTruthFcn, ...
                timeVector, truthParams);

        otherwise
            error('Unsupported truthModel: %s', truthModel);
    end

    % Always return the truth signal as a column vector.
    fdTrue = fdTrue(:);
    if numel(fdTrue) ~= numel(timeVector)
        error('The truth model must return one value per epoch.');
    end

end

function tleTruth = buildTleTruthSignal(timeVector, truthParams)

    % Translate alternative field names into the ones expected by the
    % standalone TLE Doppler generator.
    tleParams = normaliseTleTruthParams(truthParams);
    if isfield(tleParams, 'timeVector')
        tleParams = rmfield(tleParams, 'timeVector');
    end

    % Convert the structure into name-value pairs for function forwarding.
    nameValuePairs = structToNameValueCell(tleParams);
    % Generate the TLE-based Doppler truth on the requested time grid.
    tleTruth = generateTleDopplerTruth('timeVector', timeVector, ...
        nameValuePairs{:});

end

function oscillatorDrift = buildOscillatorDrift(timeVector, driftModel, ...
    driftParams, customDriftFcn)

    driftModel = lower(string(driftModel));
    % Shift time so polynomial drift models start at zero.
    time0 = timeVector - timeVector(1);

    switch driftModel
        case "none"
            % Disable oscillator drift completely.
            oscillatorDrift = zeros(size(timeVector));

        case "linear"
            % Read the initial oscillator offset contribution.
            drift0 = getFieldOrDefault(driftParams, 'drift0Hz', 0);
            % Read the constant oscillator drift rate.
            driftRate = getFieldOrDefault(driftParams, ...
                'driftRateHzPerSec', 0);
            % Build a first-order drift model.
            oscillatorDrift = drift0 + driftRate .* time0;

        case "quadratic"
            % Read the initial oscillator offset contribution.
            drift0 = getFieldOrDefault(driftParams, 'drift0Hz', 0);
            % Read the initial oscillator drift rate.
            driftRate = getFieldOrDefault(driftParams, ...
                'driftRateHzPerSec', 0);
            % Read the constant oscillator drift acceleration.
            driftAccel = getFieldOrDefault(driftParams, ...
                'driftAccelHzPerSec2', 0);
            % Build a second-order drift model.
            oscillatorDrift = drift0 + driftRate .* time0 + ...
                0.5 * driftAccel .* time0.^2;

        case "sinusoidal"
            % Read the sinusoidal drift amplitude.
            amplitude = getFieldOrDefault(driftParams, ...
                'amplitudeHz', 0);
            % Read the sinusoidal drift frequency.
            frequency = getFieldOrDefault(driftParams, ...
                'frequencyHz', 0.1);
            % Read the sinusoidal drift phase.
            phase = getFieldOrDefault(driftParams, 'phaseRad', 0);
            % Build a periodic oscillator drift.
            oscillatorDrift = amplitude * ...
                sin(2 * pi * frequency * time0 + phase);

        case "randomwalk"
            % Read the random-walk strength in Hz/sqrt(s).
            stdPerSqrtSecond = getFieldOrDefault(driftParams, ...
                'stdHzPerSqrtSec', 0);
            % Start the random walk from zero.
            oscillatorDrift = zeros(size(timeVector));
            if numel(timeVector) > 1 && stdPerSqrtSecond > 0
                % Compute the time increment at each epoch.
                dt = [0; diff(timeVector)];
                % Draw random-walk increments with variance proportional to dt.
                increments = stdPerSqrtSecond .* sqrt(dt) .* ...
                    randn(numel(timeVector), 1);
                % Accumulate the increments to form the drift sequence.
                oscillatorDrift = cumsum(increments);
            end

        case "custom"
            if isempty(customDriftFcn)
                error(['oscillatorDriftModel is custom but ' ...
                    'customOscillatorDriftFcn was not provided.']);
            end
            % Evaluate the user-supplied custom drift model.
            oscillatorDrift = evaluateCustomModel(customDriftFcn, ...
                timeVector, driftParams);

        otherwise
            error('Unsupported oscillatorDriftModel: %s', driftModel);
    end

    % Always return the oscillator drift as a column vector.
    oscillatorDrift = oscillatorDrift(:);
    if numel(oscillatorDrift) ~= numel(timeVector)
        error('The oscillator drift model must return one value per epoch.');
    end

end

function values = evaluateCustomModel(modelFcn, timeVector, params)

    try
        % First try a custom model that accepts time and parameter inputs.
        values = modelFcn(timeVector, params);
    catch
        % Fall back to a simpler custom model that only needs time.
        values = modelFcn(timeVector);
    end

end

function fieldValue = getFieldOrDefault(inputStruct, fieldName, defaultValue)

    if isfield(inputStruct, fieldName)
        % Use the value stored in the structure when the field exists.
        fieldValue = inputStruct.(fieldName);
    else
        % Fall back to the provided default when the field is absent.
        fieldValue = defaultValue;
    end

end

function nameValuePairs = structToNameValueCell(inputStruct)

    if isempty(inputStruct)
        % Return an empty cell array when there is nothing to forward.
        nameValuePairs = {};
        return
    end

    % Read all structure field names in order.
    fieldNames = fieldnames(inputStruct);
    % Preallocate the name-value output cell array.
    nameValuePairs = cell(1, 2 * numel(fieldNames));
    for fieldIndex = 1:numel(fieldNames)
        % Map each structure field onto a name-value pair position.
        pairIndex = 2 * fieldIndex - 1;
        % Store the field name.
        nameValuePairs{pairIndex} = fieldNames{fieldIndex};
        % Store the corresponding field value.
        nameValuePairs{pairIndex + 1} = inputStruct.(fieldNames{fieldIndex});
    end

end

function tleParams = normaliseTleTruthParams(inputParams)

    % Start from the input structure so only the required fields are adapted.
    tleParams = inputParams;

    if isfield(tleParams, 'fc_Hz') && ~isfield(tleParams, 'carrierFreqHz')
        % Map an alternative carrier-frequency field name to the canonical one.
        tleParams.carrierFreqHz = tleParams.fc_Hz;
    end

    if isfield(tleParams, 'c') && ~isfield(tleParams, 'speedOfLight')
        % Map an alternative speed-of-light field name to the canonical one.
        tleParams.speedOfLight = tleParams.c;
    end

    if isfield(tleParams, 'minMaxElev') && ...
            ~isfield(tleParams, 'minMaxElevDeg')
        % Map the minimum maximum-elevation field onto the expected name.
        tleParams.minMaxElevDeg = tleParams.minMaxElev;
    end

    if isfield(tleParams, 'rx') && isstruct(tleParams.rx)
        if isfield(tleParams.rx, 'lat_deg') && ~isfield(tleParams, 'rxLatDeg')
            % Promote receiver latitude from the nested receiver structure.
            tleParams.rxLatDeg = tleParams.rx.lat_deg;
        end
        if isfield(tleParams.rx, 'lon_deg') && ~isfield(tleParams, 'rxLonDeg')
            % Promote receiver longitude from the nested receiver structure.
            tleParams.rxLonDeg = tleParams.rx.lon_deg;
        end
        if isfield(tleParams.rx, 'alt_m') && ~isfield(tleParams, 'rxAltM')
            % Promote receiver altitude from the nested receiver structure.
            tleParams.rxAltM = tleParams.rx.alt_m;
        end
        if isfield(tleParams.rx, 'vel_ecef') && ...
                ~isfield(tleParams, 'receiverVelocityEcefMps')
            % Promote receiver velocity from the nested receiver structure.
            tleParams.receiverVelocityEcefMps = tleParams.rx.vel_ecef;
        end
    end

    allowedFields = { ...
        'tleFile', 'startUTC', 'startFromTLE', 'maxTLEAgeDays', ...
        'rxLatDeg', 'rxLonDeg', 'rxAltM', 'receiverVelocityEcefMps', ...
        'carrierFreqHz', 'speedOfLight', 'useSGP4', ...
        'allowJ2Fallback', 'selectionMode', 'satelliteNameContains', ...
        'satelliteIndex', 'coarseStep', 'minElevDeg', ...
        'minMaxElevDeg', 'verbose'};

    % List the fields currently present after the renaming stage.
    currentFields = fieldnames(tleParams);
    % Identify any fields that are not accepted by the TLE truth generator.
    discardFields = currentFields(~ismember(currentFields, allowedFields));
    if ~isempty(discardFields)
        % Remove unsupported fields before forwarding the parameter structure.
        tleParams = rmfield(tleParams, discardFields);
    end

end

function derivative = numericalDerivative(signal, timeVector)

    % Work with column vectors to keep the derivative shape consistent.
    signal = signal(:);
    timeVector = timeVector(:);

    if numel(signal) <= 1
        % Return zero when a derivative cannot be estimated from one sample.
        derivative = zeros(size(signal));
        return
    end

    % Estimate the derivative numerically using MATLAB's gradient operator.
    derivative = gradient(signal, timeVector);

end

