%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            trackingInitialErrorSensitivity.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Summarises how the tracking error evolves relative to the
% initial estimation error.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metric = trackingInitialErrorSensitivity(reference, estimate)
    % Convert the reference sequence into a column vector.
    reference = reference(:);
    % Convert the estimate sequence into a column vector.
    estimate = estimate(:);

    % Keep only the samples where both reference and estimate are finite.
    validMask = isfinite(reference) & isfinite(estimate);
    % Return an all-NaN structure when no valid sample pair exists.
    if ~any(validMask)
        metric = struct( ...
            'initialAbsError', NaN, ...
            'peakAbsError', NaN, ...
            'finalAbsError', NaN, ...
            'rmsAbsError', NaN, ...
            'amplificationRatio', NaN, ...
            'finalToInitialRatio', NaN);
        return
    end

    % Compute the signed estimation error on the valid samples.
    signedError = estimate(validMask) - reference(validMask);
    % Compute the absolute error envelope.
    absError = abs(signedError);
    % Read the first absolute error as the initial tracking error.
    initialAbsError = absError(1);
    % Protect the ratios against division by zero.
    denominator = max(initialAbsError, eps);

    % Summarise how the error evolves relative to the initial value.
    metric = struct( ...
        'initialAbsError', initialAbsError, ...
        'peakAbsError', max(absError), ...
        'finalAbsError', absError(end), ...
        'rmsAbsError', sqrt(mean(absError.^2)), ...
        'amplificationRatio', max(absError) / denominator, ...
        'finalToInitialRatio', absError(end) / denominator);

end

