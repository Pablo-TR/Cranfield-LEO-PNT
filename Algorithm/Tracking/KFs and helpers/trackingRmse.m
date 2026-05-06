%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        trackingRmse.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 13th March 2026
% Last Update: 13th March 2026
% Version: 1.0.0
%
% Description: Computes the root mean square error between a reference and
% an estimate vector.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rmseValue = trackingRmse(reference, estimate)
    % Convert the reference sequence into a column vector.
    reference = reference(:);
    % Convert the estimate sequence into a column vector.
    estimate = estimate(:);

    % Keep only the samples where both reference and estimate are finite.
    validMask = isfinite(reference) & isfinite(estimate);
    % Return NaN when no valid sample pair is available.
    if ~any(validMask)
        rmseValue = NaN;
        return
    end

    % Compute the root mean square error over the valid samples only.
    rmseValue = sqrt(mean((reference(validMask) - estimate(validMask)).^2));

end

