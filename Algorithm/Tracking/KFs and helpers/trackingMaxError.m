%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     trackingMaxError.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Computes the maximum absolute tracking error and the sample
% where it occurs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metric = trackingMaxError(reference, estimate)
    % Convert the reference sequence into a column vector.
    reference = reference(:);
    % Convert the estimate sequence into a column vector.
    estimate = estimate(:);

    % Keep only the samples where both signals are finite.
    validMask = isfinite(reference) & isfinite(estimate);
    % Return an all-NaN metric structure when no valid sample exists.
    if ~any(validMask)
        metric = struct( ...
            'maxAbsError', NaN, ...
            'signedErrorAtMax', NaN, ...
            'index', NaN);
        return
    end

    % Compute the signed estimation error on the valid samples.
    signedError = estimate(validMask) - reference(validMask);
    % Compute the absolute value of the valid errors.
    absError = abs(signedError);
    % Find the largest absolute error and its local index.
    [maxAbsError, localIndex] = max(absError);
    % Map the local valid index back to the original sample index.
    validIndices = find(validMask);

    % Package the maximum-error metrics into one output structure.
    metric = struct( ...
        'maxAbsError', maxAbsError, ...
        'signedErrorAtMax', signedError(localIndex), ...
        'index', validIndices(localIndex));

end

