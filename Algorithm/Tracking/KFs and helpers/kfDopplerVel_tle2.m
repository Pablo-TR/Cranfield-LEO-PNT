%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    kfDopplerVel_tle2.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 15th March 2026
% Last Update: 15th March 2026
% Version: 1.0.0
%
% Description: Runs a 2-state error-state Kalman filter around a Doppler
% reference generated from TLEs. The state is [delta_fd; delta_fd_dot],
% where the estimated absolute outputs are reconstructed as
% fd_est = fd_ref + delta_fd and fd_dot_est = fd_dot_ref + delta_fd_dot.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = kfDopplerVel_tle2(fdMeas, fdRef, stepTime, sigmaFd, sigmaFddot,...
    processNoiseScale, fdDotRef, measurementVar)

    % Read the total number of Doppler epochs.
    numEpochs = numel(fdMeas);
    % Define the 2-state constant-velocity transition matrix for the error state.
    Fmat = [1, stepTime;
            0, 1];
    % Define the scalar measurement model for Doppler error only.
    Hmat = [1, 0];

    % Build the finite-difference Doppler-rate measurement for diagnostics.
    fdDotMeas = [0; diff(fdMeas) / stepTime];
    % Convert the absolute Doppler measurement into a Doppler-error measurement.
    deltaFdMeas = fdMeas - fdRef;
    % Convert the absolute Doppler-rate measurement into a rate-error measurement.
    deltaFdDotMeas = fdDotMeas - fdDotRef;

    % Build the scalar Doppler-error measurement variance.
    R = measurementVar^2;
    % Build the initial 2-state error estimate.
    x = [deltaFdMeas(1);
         deltaFdDotMeas(2)];
    % Build the initial covariance matrix from the user-provided standard deviations.
    P = diag([sigmaFd^2, sigmaFddot^2]);

    % Build the process-noise matrix for the constant-velocity error model.
    Q = processNoiseScale * [stepTime^3 / 3, stepTime^2 / 2;
                             stepTime^2 / 2, stepTime];

    % Preallocate the error-state history array.
    stateHistory = zeros(2, numEpochs);
    % Preallocate the raw innovation history.
    innovation = zeros(1, numEpochs - 1);
    % Preallocate the normalised innovation history.
    normalisedInnovation = zeros(1, numEpochs - 1);
    % Preallocate the covariance history array.
    covarianceHistory = zeros(2, 2, numEpochs);

    % Store the initial state estimate.
    stateHistory(:, 1) = x;
    % Store the initial covariance estimate.
    covarianceHistory(:, :, 1) = P;

    % Run the Kalman filter over the full Doppler sequence.
    for epoch = 1:numEpochs - 1
        xPred = Fmat * x;
        PPred = Fmat * P * Fmat' + Q;

        y = deltaFdMeas(epoch + 1);
        S = Hmat * PPred * Hmat' + R;
        Kk = (PPred * Hmat') / S;
        z = y - Hmat * xPred;
        nu = z / sqrt(S);

        x = xPred + Kk * z;
        josephTerm = eye(2) - Kk * Hmat;
        P = josephTerm * PPred * josephTerm' + Kk * R * Kk';

        stateHistory(:, epoch + 1) = x;
        innovation(epoch) = z;
        normalisedInnovation(epoch) = nu;
        covarianceHistory(:, :, epoch + 1) = P;
    end

    % Extract the estimated Doppler error history.
    deltaFdEst = stateHistory(1, :).';
    % Extract the estimated Doppler-rate error history.
    deltaFdDotEst = stateHistory(2, :).';

    % Package the reconstructed absolute estimates and diagnostics into one output structure.
    output = struct( ...
        'fdEst', fdRef + deltaFdEst, ...
        'fdDotEst', fdDotRef + deltaFdDotEst, ...
        'deltaFdEst', deltaFdEst, ...
        'deltaFdDotEst', deltaFdDotEst, ...
        'fdRef', fdRef, ...
        'fdDotRef', fdDotRef, ...
        'fdDotMeas', fdDotMeas, ...
        'deltaFdMeas', deltaFdMeas, ...
        'deltaFdDotMeas', deltaFdDotMeas, ...
        'stateHistory', stateHistory, ...
        'innovation', innovation(:), ...
        'normalisedInnovation', normalisedInnovation(:), ...
        'covarianceHistory', covarianceHistory, ...
        'measurementCovariance', R, ...
        'stepTime', stepTime, ...
        'sigmaFd', sigmaFd, ...
        'processNoiseScale', processNoiseScale);

end

