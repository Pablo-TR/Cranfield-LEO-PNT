%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      linearKF.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 5th May 2026
% Last Update: 5th May 2026
% Version: 1.0.0
%
% Description: Runs a 2 or 3-state linear Kalman filter with state(s) [fd; fd_dot; fd_ddot] and a
% scalar Doppler measurement.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = linearKF(fdMeas, stepTime,...
    processNoiseScale, measurementVar, d2fdt2_max, is_3state)
    % Read the total number of Doppler epochs.
    numEpochs = numel(fdMeas);
    % Estimate the initial Doppler rate from the first two samples.
    dfdt1 = (fdMeas(2)-fdMeas(1))/stepTime;
    % Estimate the next Doppler-rate sample from the following pair.
    dfdt2 = (fdMeas(3)-fdMeas(2))/stepTime;
    % Estimate the initial Doppler acceleration from the first two rate estimates.
    d2fdt2 = (dfdt2-dfdt1)/stepTime;

    % Convert the acceleration bound into the expected Doppler increment scale.
    s1 =  0.5*d2fdt2_max*stepTime^2;
    % Convert the acceleration bound into the expected Doppler-rate increment scale.
    s2 = d2fdt2_max*stepTime;

    % Select either the 3-state or the 2-state filter structure.
    if is_3state
        % Use the acceleration bound directly as the third-state scale.
        s3 = d2fdt2_max;
        % Build the initial covariance matrix for the 3-state model.
        P = diag([measurementVar^2, dfdt1^2, d2fdt2^2]);
        % Build a diagonal process-noise matrix from the acceleration scales.
        Q = processNoiseScale*diag([s1^2, s2^2, s3^2]);
        % Build the initial 3-state estimate.
        x = [fdMeas(1);
             dfdt1;
             0];
        % Define the 3-state constant-acceleration transition matrix.
        Fmat = [1, stepTime, 0.5 * stepTime^2;
                0, 1, stepTime;
                0, 0, 1];

        % Define the scalar Doppler measurement model for the 3-state case.
        Hmat = [1, 0, 0];

        % Preallocate the 3-state history array.
        stateHistory = zeros(3, numEpochs);
        % Preallocate the raw innovation history.
        innovation = zeros(1, numEpochs - 1);
        % Preallocate the normalised innovation history.
        normalisedInnovation = zeros(1, numEpochs - 1);
        % Preallocate the covariance history array.
        covarianceHistory = zeros(3, 3, numEpochs);

    else
        % Build the initial covariance matrix for the 2-state model.
        P = diag([measurementVar^2, dfdt1^2]);
        % Build a diagonal process-noise matrix from the acceleration scales.
        Q = processNoiseScale*diag([s1^2, s2^2]);
        % Define the 2-state constant-velocity transition matrix.
        Fmat = [1, stepTime;
                0, 1];
        % Define the scalar Doppler measurement model for the 2-state case.
        Hmat = [1, 0];
        % Build the initial 2-state estimate.
        x = [fdMeas(1);
             dfdt1];

        % Preallocate the 2-state history array.
        stateHistory = zeros(2, numEpochs);
        % Preallocate the raw innovation history.
        innovation = zeros(1, numEpochs - 1);
        % Preallocate the normalised innovation history.
        normalisedInnovation = zeros(1, numEpochs - 1);
        % Preallocate the covariance history array.
        covarianceHistory = zeros(2, 2, numEpochs);
    end 

    % Store the initial state estimate.
    stateHistory(:, 1) = x;
    % Store the initial covariance estimate.
    covarianceHistory(:, :, 1) = P;

    % Run the selected Kalman filter over the full Doppler sequence.
    for epoch = 2:numEpochs
        xPred = Fmat * x;
        PPred = Fmat * P * Fmat' + Q;

        y = fdMeas(epoch);
        S = Hmat * PPred * Hmat' + measurementVar^2;
        Kk = (PPred * Hmat') / S;
        z = y - Hmat * xPred;
        nu = z / sqrt(S);

        x = xPred + Kk * z;
        P = PPred-Kk*Hmat*PPred;

        stateHistory(:, epoch) = x;
        innovation(epoch) = z;
        normalisedInnovation(epoch) = nu;
        covarianceHistory(:, :, epoch) = P;
    end

    % Package the filter estimates and diagnostics into one output structure.
    output = struct( ...
        'fdEst', stateHistory(1, :).', ...
        'fdDotEst', stateHistory(2, :).', ...
        'stateHistory', stateHistory, ...
        'innovation', innovation(:), ...
        'normalisedInnovation', normalisedInnovation(:), ...
        'covarianceHistory', covarianceHistory, ...
        'stepTime', stepTime, ...
        'measurementVar', measurementVar^2, ...
        'processNoiseScale', processNoiseScale);

end
