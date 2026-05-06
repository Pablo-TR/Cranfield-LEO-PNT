%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%             trackingGapOutlierRobustness.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Computes robustness metrics for tracking performance in the
% presence of missing measurements and outliers.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metric = trackingGapOutlierRobustness(reference, estimate, ...
    gapMask, outlierMask, recoveryWindowSamples)

    % Use five samples as the default recovery window.
    if nargin < 5
        recoveryWindowSamples = 5;
    end

    % Convert the reference sequence into a column vector.
    reference = reference(:);
    % Convert the estimate sequence into a column vector.
    estimate = estimate(:);
    % Convert the gap mask into a logical column vector.
    gapMask = logical(gapMask(:));
    % Convert the outlier mask into a logical column vector.
    outlierMask = logical(outlierMask(:));

    % Mark the samples where both reference and estimate are finite.
    validMask = isfinite(reference) & isfinite(estimate);
    % Build the mask of nominal samples with no impairments.
    nominalMask = validMask & ~gapMask & ~outlierMask;
    % Build the mask of samples directly affected by gaps.
    gapOnlyMask = validMask & gapMask;
    % Build the mask of samples directly affected by outliers.
    outlierOnlyMask = validMask & outlierMask;
    % Build the mask of samples immediately after gaps.
    postGapMask = postEventMask(gapMask, recoveryWindowSamples) & validMask;
    % Build the mask of samples immediately after outliers.
    postOutlierMask = postEventMask(outlierMask, recoveryWindowSamples) & ...
        validMask;

    % Collect the robustness metrics for each impairment region.
    metric = struct( ...
        'numGapSamples', nnz(gapMask), ...
        'numOutlierSamples', nnz(outlierMask), ...
        'recoveryWindowSamples', recoveryWindowSamples, ...
        'rmseNominal', maskedRmse(reference, estimate, nominalMask), ...
        'rmseGapSamples', maskedRmse(reference, estimate, gapOnlyMask), ...
        'rmseOutlierSamples', maskedRmse(reference, estimate, ...
            outlierOnlyMask), ...
        'rmsePostGap', maskedRmse(reference, estimate, postGapMask), ...
        'rmsePostOutlier', maskedRmse(reference, estimate, ...
            postOutlierMask), ...
        'maxErrorNominal', maskedMaxError(reference, estimate, nominalMask), ...
        'maxErrorGapSamples', maskedMaxError(reference, estimate, ...
            gapOnlyMask), ...
        'maxErrorOutlierSamples', maskedMaxError(reference, estimate, ...
            outlierOnlyMask), ...
        'maxErrorPostGap', maskedMaxError(reference, estimate, ...
            postGapMask), ...
        'maxErrorPostOutlier', maskedMaxError(reference, estimate, ...
            postOutlierMask));

    % Compare the post-gap RMSE against the nominal RMSE.
    metric.postGapToNominalRmseRatio = ratioOrNaN(metric.rmsePostGap, ...
        metric.rmseNominal);
    % Compare the post-outlier RMSE against the nominal RMSE.
    metric.postOutlierToNominalRmseRatio = ratioOrNaN( ...
        metric.rmsePostOutlier, metric.rmseNominal);

end

function expandedMask = postEventMask(eventMask, recoveryWindowSamples)

    % Initialise the expanded post-event mask with false values.
    expandedMask = false(size(eventMask));
    % Extract the indices where the event mask is active.
    eventIndices = find(eventMask);
    for idx = 1:numel(eventIndices)
        % Start the post-event window one sample after the event.
        startIdx = eventIndices(idx) + 1;
        % Stop the post-event window after the configured recovery length.
        stopIdx = min(numel(eventMask), ...
            eventIndices(idx) + recoveryWindowSamples);
        if startIdx <= stopIdx
            % Mark the recovery window after the event as active.
            expandedMask(startIdx:stopIdx) = true;
        end
    end

end

function rmseValue = maskedRmse(reference, estimate, mask)

    % Return NaN when the requested mask selects no sample.
    if ~any(mask)
        rmseValue = NaN;
        return
    end

    % Compute the RMSE on the masked subset only.
    rmseValue = trackingRmse(reference(mask), estimate(mask));

end

function maxErrorValue = maskedMaxError(reference, estimate, mask)

    % Return NaN when the requested mask selects no sample.
    if ~any(mask)
        maxErrorValue = NaN;
        return
    end

    % Compute the maximum absolute error on the masked subset.
    maxMetric = trackingMaxError(reference(mask), estimate(mask));
    % Extract the scalar maximum error value from the helper structure.
    maxErrorValue = maxMetric.maxAbsError;

end

function ratioValue = ratioOrNaN(numerator, denominator)

    % Return NaN when the ratio is not numerically meaningful.
    if ~isfinite(numerator) || ~isfinite(denominator) || denominator <= eps
        ratioValue = NaN;
    else
        % Compute the requested ratio when both values are valid.
        ratioValue = numerator / denominator;
    end

end

