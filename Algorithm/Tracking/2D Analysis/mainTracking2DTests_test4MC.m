%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               mainTracking2DTests.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 14th April 2026
% Last Update: 14th April 2026
% Version: 1.0.0
%
% Description: Main file to test the 2D Doppler model and the 2-state and
% 3-state linear Kalman filters in a controlled environment.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

% Fix the random generator for repeatable results.
rng(7);

% Define the baseline model timestep in seconds.
baseStepTime = 1e-3;

% Define the orbital height in metres.
baseH = 1200e3;

% Define the carrier frequency in Hz.
baseF0 = 11.325e9;

% Define the mixing frequency in Hz.
baseFmix = -1;

% Define the initial central angle in rad.
baseLambda0 = -0.6;

% Define the baseline measurement noise standard deviation in Hz.
baseSigmaNoise = 1000;

% Define the baseline process-noise scale.
baseProcessNoiseScale = 1;

% Define the Doppler RMSE threshold in Hz.
dopplerRmseThreshold = 100;

% Define the Doppler-rate RMSE threshold in Hz/s.
dopplerRateRmseThreshold = 1;

% Generate the baseline 2D Doppler model.
[tBase, fdBase, dfdtBase, d2fdt2Base] = generateOWDoppler2DModel(baseStepTime, 0, baseH, baseF0, baseFmix, baseLambda0);

% Convert the time vector to a column vector.
tBase = tBase(:);

% Build a Doppler-rate ground-truth vector with the same length as fdBase.
fdDotBase = [dfdtBase(1); dfdtBase(:)];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1. NOMINAL CASE WITHOUT PERTURBATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the ground truth as the measurement vector.
fdMeasNominal = fdBase;

% Run the 2-state Kalman filter.
output2Nominal = linearKF(fdMeasNominal, baseStepTime, baseProcessNoiseScale, 0, d2fdt2Base, 0);

% Extract the estimated Doppler from the 2-state filter.
fdEst2Nominal = output2Nominal.fdEst;

% Extract the estimated Doppler rate from the 2-state filter.
fdDotEst2Nominal = output2Nominal.fdDotEst;

% Run the 3-state Kalman filter.
output3Nominal = linearKF(fdMeasNominal, baseStepTime, baseProcessNoiseScale, 0, d2fdt2Base, 1);

% Extract the estimated Doppler from the 3-state filter.
fdEst3Nominal = output3Nominal.fdEst;

% Extract the estimated Doppler rate from the 3-state filter.
fdDotEst3Nominal = output3Nominal.fdDotEst;

% Compute the Doppler error of the 2-state filter.
fdError2Nominal = fdEst2Nominal - fdBase;

% Compute the Doppler-rate error of the 2-state filter.
fdDotError2Nominal = fdDotEst2Nominal - fdDotBase;

% Compute the Doppler error of the 3-state filter.
fdError3Nominal = fdEst3Nominal - fdBase;

% Compute the Doppler-rate error of the 3-state filter.
fdDotError3Nominal = fdDotEst3Nominal - fdDotBase;

% Compute the Doppler RMSE of the 2-state filter.
rmseFd2Nominal = sqrt(mean(fdError2Nominal.^2));

% Compute the Doppler-rate RMSE of the 2-state filter.
rmseFdDot2Nominal = sqrt(mean(fdDotError2Nominal.^2));

% Compute the Doppler RMSE of the 3-state filter.
rmseFd3Nominal = sqrt(mean(fdError3Nominal.^2));

% Compute the Doppler-rate RMSE of the 3-state filter.
rmseFdDot3Nominal = sqrt(mean(fdDotError3Nominal.^2));

% Initialise the convergence index of the 2-state filter.
convIndex2Nominal = NaN;

% Search for the first epoch from which the tail RMSE stays below threshold.
for k = 1:length(fdBase)
    fdTailRmse = sqrt(mean((fdEst2Nominal(k:end) - fdBase(k:end)).^2));
    fdDotTailRmse = sqrt(mean((fdDotEst2Nominal(k:end) - fdDotBase(k:end)).^2));
    if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
        convIndex2Nominal = k;
        break
    end
end

% Initialise the convergence index of the 3-state filter.
convIndex3Nominal = NaN;

% Search for the first epoch from which the tail RMSE stays below threshold.
for k = 1:length(fdBase)
    fdTailRmse = sqrt(mean((fdEst3Nominal(k:end) - fdBase(k:end)).^2));
    fdDotTailRmse = sqrt(mean((fdDotEst3Nominal(k:end) - fdDotBase(k:end)).^2));
    if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
        convIndex3Nominal = k;
        break
    end
end

% Convert the 2-state convergence index into time.
if isnan(convIndex2Nominal)
    convTime2Nominal = NaN;
else
    convTime2Nominal = tBase(convIndex2Nominal);
end

% Convert the 3-state convergence index into time.
if isnan(convIndex3Nominal)
    convTime3Nominal = NaN;
else
    convTime3Nominal = tBase(convIndex3Nominal);
end

% Display a header for the first test.
disp(' ')
disp('TEST 1. NOMINAL CASE')

% Display the 2-state metrics.
disp(['2-state KF Doppler RMSE [Hz]: ', num2str(rmseFd2Nominal)])
disp(['2-state KF Doppler-rate RMSE [Hz/s]: ', num2str(rmseFdDot2Nominal)])
disp(['2-state KF convergence time [s]: ', num2str(convTime2Nominal)])

% Display the 3-state metrics.
disp(['3-state KF Doppler RMSE [Hz]: ', num2str(rmseFd3Nominal)])
disp(['3-state KF Doppler-rate RMSE [Hz/s]: ', num2str(rmseFdDot3Nominal)])
disp(['3-state KF convergence time [s]: ', num2str(convTime3Nominal)])

% Create a figure for the nominal Doppler estimation.
figure;

% Plot the ground-truth Doppler.
plot(tBase, fdBase, 'k', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Plot the 2-state estimate.
plot(tBase, fdEst2Nominal, '--', 'LineWidth', 1.2)

% Plot the 3-state estimate.
plot(tBase, fdEst3Nominal, '-.', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Time [s]')

% Label the vertical axis.
ylabel('Doppler [Hz]')

% Add title.
title('Test 1. Nominal case without perturbations')

% Add legend.
legend('Ground truth', '2-state KF', '3-state KF', 'Location', 'best')

% Create a figure for the nominal Doppler-rate estimation.
figure;

% Plot the ground-truth Doppler rate.
plot(tBase, fdDotBase, 'k', 'LineWidth', 1.5)

% Keep the figure active.
hold on

% Plot the 2-state Doppler-rate estimate.
plot(tBase, fdDotEst2Nominal, '--', 'LineWidth', 1.2)

% Plot the 3-state Doppler-rate estimate.
plot(tBase, fdDotEst3Nominal, '-.', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Time [s]')

% Label the vertical axis.
ylabel('Doppler rate [Hz/s]')

% Add title.
title('Test 1. Nominal Doppler-rate estimation')

% Add legend.
legend('Ground truth', '2-state KF', '3-state KF', 'Location', 'best')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4. COMPARISON BETWEEN REGIONS OF THE S-CURVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the representative centre times in seconds.
regionCenterTimeVec = [400 650 850];

% Define the region labels.
regionLabels = {'Early', 'Middle', 'Late'};

% Define the fixed observation-window duration in seconds.
regionWindowDuration = 100;

% Define the number of Monte Carlo runs for this test.
numRunsRegion = 100;

% Define the run that will be used to show example S-curves.
exampleRunIndex = 1;

% Preallocate the 2-state Doppler RMSE by region and run.
rmseFd2RegionRuns = zeros(numRunsRegion, length(regionCenterTimeVec));

% Preallocate the 2-state Doppler-rate RMSE by region and run.
rmseFdDot2RegionRuns = zeros(numRunsRegion, length(regionCenterTimeVec));

% Preallocate the 2-state convergence time by region and run.
convTime2RegionRuns = NaN(numRunsRegion, length(regionCenterTimeVec));

% Preallocate the 3-state Doppler RMSE by region and run.
rmseFd3RegionRuns = zeros(numRunsRegion, length(regionCenterTimeVec));

% Preallocate the 3-state Doppler-rate RMSE by region and run.
rmseFdDot3RegionRuns = zeros(numRunsRegion, length(regionCenterTimeVec));

% Preallocate the 3-state convergence time by region and run.
convTime3RegionRuns = NaN(numRunsRegion, length(regionCenterTimeVec));

% Preallocate the example-run measurement vector.
fdMeasRegionExample = zeros(size(fdBase));

% Create cells to store the example 2-state estimates by region.
fdEst2RegionExample = cell(length(regionCenterTimeVec), 1);

% Create cells to store the example 3-state estimates by region.
fdEst3RegionExample = cell(length(regionCenterTimeVec), 1);

% Create cells to store the example local time vectors by region.
tWindowRegionExample = cell(length(regionCenterTimeVec), 1);

% Create cells to store the example ground-truth Doppler vectors by region.
fdWindowRegionExample = cell(length(regionCenterTimeVec), 1);

% Create a figure to show the analysed regions on the S-curve.
figure;

% Plot the complete Doppler curve.
plot(tBase, fdBase, 'k', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Loop over the different regions to plot them on the S-curve.
for i = 1:length(regionCenterTimeVec)

    % Find the index nearest to the desired centre time.
    [~, centreIndex] = min(abs(tBase - regionCenterTimeVec(i)));

    % Compute half of the window length in samples.
    halfWindowSamples = round(regionWindowDuration / (2 * baseStepTime));

    % Compute the first index of the region.
    startIndex = centreIndex - halfWindowSamples;

    % Compute the last index of the region.
    endIndex = centreIndex + halfWindowSamples;

    % Limit the first index to the available signal.
    startIndex = max(startIndex, 1);

    % Limit the last index to the available signal.
    endIndex = min(endIndex, length(fdBase));

    % Extract the local time vector.
    tWindow = tBase(startIndex:endIndex);

    % Extract the ground-truth Doppler in the region.
    fdWindow = fdBase(startIndex:endIndex);

    % Plot the analysed region on top of the complete S-curve.
    plot(tWindow, fdWindow, 'LineWidth', 3)
end

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Time [s]')

% Label the vertical axis.
ylabel('Doppler [Hz]')

% Add title.
title('Test 4. Regions selected on the S-curve')

% Add legend.
legend('Complete S-curve', 'Early region', 'Middle region', 'Late region', 'Location', 'best')

% Loop over the repeated Monte Carlo runs.
for run = 1:numRunsRegion

    % Build one noisy measurement vector for the current run.
    fdMeasRegion = addNoise(fdBase, baseSigmaNoise);

    % Store the measurement vector of the example run.
    if run == exampleRunIndex
        fdMeasRegionExample = fdMeasRegion;
    end

    % Loop over the different regions.
    for i = 1:length(regionCenterTimeVec)

        % Find the index nearest to the desired centre time.
        [~, centreIndex] = min(abs(tBase - regionCenterTimeVec(i)));

        % Compute half of the window length in samples.
        halfWindowSamples = round(regionWindowDuration / (2 * baseStepTime));

        % Compute the first index of the region.
        startIndex = centreIndex - halfWindowSamples;

        % Compute the last index of the region.
        endIndex = centreIndex + halfWindowSamples;

        % Limit the first index to the available signal.
        startIndex = max(startIndex, 1);

        % Limit the last index to the available signal.
        endIndex = min(endIndex, length(fdBase));

        % Extract the local time vector.
        tWindow = tBase(startIndex:endIndex);

        % Shift the local time to start at zero.
        tWindowLocal = tWindow - tWindow(1);

        % Extract the ground-truth Doppler in the region.
        fdWindow = fdBase(startIndex:endIndex);

        % Extract the noisy Doppler measurements in the region.
        fdMeasWindow = fdMeasRegion(startIndex:endIndex);

        % Extract the ground-truth Doppler rate in the region.
        fdDotWindow = fdDotBase(startIndex:endIndex);

        % Run the 2-state Kalman filter on the region.
        output2 = linearKF(fdMeasWindow, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 0);

        % Extract the 2-state Doppler estimate.
        fdEst2 = output2.fdEst;

        % Extract the 2-state Doppler-rate estimate.
        fdDotEst2 = output2.fdDotEst;

        % Run the 3-state Kalman filter on the region.
        output3 = linearKF(fdMeasWindow, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 1);

        % Extract the 3-state Doppler estimate.
        fdEst3 = output3.fdEst;

        % Extract the 3-state Doppler-rate estimate.
        fdDotEst3 = output3.fdDotEst;

        % Compute the 2-state Doppler RMSE for the region.
        rmseFd2RegionRuns(run, i) = sqrt(mean((fdEst2 - fdWindow).^2));

        % Compute the 2-state Doppler-rate RMSE for the region.
        rmseFdDot2RegionRuns(run, i) = sqrt(mean((fdDotEst2 - fdDotWindow).^2));

        % Compute the 3-state Doppler RMSE for the region.
        rmseFd3RegionRuns(run, i) = sqrt(mean((fdEst3 - fdWindow).^2));

        % Compute the 3-state Doppler-rate RMSE for the region.
        rmseFdDot3RegionRuns(run, i) = sqrt(mean((fdDotEst3 - fdDotWindow).^2));

        % Search for 2-state convergence within the region.
        for k = 1:length(fdWindow)
            fdTailRmse = sqrt(mean((fdEst2(k:end) - fdWindow(k:end)).^2));
            fdDotTailRmse = sqrt(mean((fdDotEst2(k:end) - fdDotWindow(k:end)).^2));
            if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
                convTime2RegionRuns(run, i) = tWindowLocal(k);
                break
            end
        end

        % Search for 3-state convergence within the region.
        for k = 1:length(fdWindow)
            fdTailRmse = sqrt(mean((fdEst3(k:end) - fdWindow(k:end)).^2));
            fdDotTailRmse = sqrt(mean((fdDotEst3(k:end) - fdDotWindow(k:end)).^2));
            if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
                convTime3RegionRuns(run, i) = tWindowLocal(k);
                break
            end
        end

        % Store the example-run data and estimates.
        if run == exampleRunIndex
            tWindowRegionExample{i} = tWindow;
            fdWindowRegionExample{i} = fdWindow;
            fdEst2RegionExample{i} = fdEst2;
            fdEst3RegionExample{i} = fdEst3;
        end
    end
end

% Compute the mean 2-state Doppler RMSE by region.
rmseFd2Region = mean(rmseFd2RegionRuns, 1).';

% Compute the mean 2-state Doppler-rate RMSE by region.
rmseFdDot2Region = mean(rmseFdDot2RegionRuns, 1).';

% Compute the mean 2-state convergence time by region.
convTime2Region = mean(convTime2RegionRuns, 1, 'omitnan').';

% Compute the mean 3-state Doppler RMSE by region.
rmseFd3Region = mean(rmseFd3RegionRuns, 1).';

% Compute the mean 3-state Doppler-rate RMSE by region.
rmseFdDot3Region = mean(rmseFdDot3RegionRuns, 1).';

% Compute the mean 3-state convergence time by region.
convTime3Region = mean(convTime3RegionRuns, 1, 'omitnan').';

% Create a figure for the region comparison metrics.
figure;

% Create the first subplot.
subplot(3,1,1)

% Plot the mean Doppler RMSE by region.
bar(categorical(regionLabels), [rmseFd2Region rmseFd3Region])

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler RMSE [Hz]')

% Add title.
title('Test 4. Comparison between regions of the S-curve')

% Add legend.
legend('2-state KF', '3-state KF', 'Location', 'northwest')

% Create the second subplot.
subplot(3,1,2)

% Plot the mean Doppler-rate RMSE by region.
bar(categorical(regionLabels), [rmseFdDot2Region rmseFdDot3Region])

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler-rate RMSE [Hz/s]')

% Create the third subplot.
subplot(3,1,3)

% Plot the mean convergence time by region.
bar(categorical(regionLabels), [convTime2Region convTime3Region])

% Add grid.
grid on

% Label the vertical axis.
ylabel('Convergence time [s]')

% Create a figure for the example run with the 2-state KF.
figure;

% Plot the complete ground-truth S-curve.
plot(tBase, fdBase, 'k', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the noisy measurements of the example run.
plot(tBase, fdMeasRegionExample, '.', 'MarkerSize', 6)

% Loop over the different regions.
for i = 1:length(regionCenterTimeVec)

    % Plot the 2-state estimate of the current region.
    plot(tWindowRegionExample{i}, fdEst2RegionExample{i}, 'LineWidth', 2)
end

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Time [s]')

% Label the vertical axis.
ylabel('Doppler [Hz]')

% Add title.
title('Test 4. Example run on the S-curve with the 2-state KF')

% Add legend.
legend('Ground truth', 'Noisy measurements', '2-state KF early', '2-state KF middle', '2-state KF late', 'Location', 'best')

% Create a figure for the example run with the 3-state KF.
figure;

% Plot the complete ground-truth S-curve.
plot(tBase, fdBase, 'k', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the noisy measurements of the example run.
plot(tBase, fdMeasRegionExample, '.', 'MarkerSize', 6)

% Loop over the different regions.
for i = 1:length(regionCenterTimeVec)

    % Plot the 3-state estimate of the current region.
    plot(tWindowRegionExample{i}, fdEst3RegionExample{i}, 'LineWidth', 2)
end

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Time [s]')

% Label the vertical axis.
ylabel('Doppler [Hz]')

% Add title.
title('Test 4. Example run on the S-curve with the 3-state KF')

% Add legend.
legend('Ground truth', 'Noisy measurements', '3-state KF early', '3-state KF middle', '3-state KF late', 'Location', 'best')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 5. SENSITIVITY TO THE OBSERVED-WINDOW DURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the centre time of the analysed window in seconds.
windowCentreTime = 650;

% Define the different observation-window durations in seconds.
windowDurationVec = [50 100 150 200 300 400 500 600];

% Build one noisy measurement vector for this test.
fdMeasDuration = addNoise(fdBase, baseSigmaNoise);

% Find the index nearest to the selected centre time.
[~, centreIndex] = min(abs(tBase - windowCentreTime));

% Preallocate the 2-state Doppler RMSE vector.
rmseFd2Duration = zeros(length(windowDurationVec), 1);

% Preallocate the 2-state Doppler-rate RMSE vector.
rmseFdDot2Duration = zeros(length(windowDurationVec), 1);

% Preallocate the 2-state convergence-time vector.
convTime2Duration = zeros(length(windowDurationVec), 1);

% Preallocate the 3-state Doppler RMSE vector.
rmseFd3Duration = zeros(length(windowDurationVec), 1);

% Preallocate the 3-state Doppler-rate RMSE vector.
rmseFdDot3Duration = zeros(length(windowDurationVec), 1);

% Preallocate the 3-state convergence-time vector.
convTime3Duration = zeros(length(windowDurationVec), 1);

% Loop over the different observation-window durations.
for i = 1:length(windowDurationVec)

    % Extract the current window duration.
    currentWindowDuration = windowDurationVec(i);

    % Compute half of the window length in samples.
    halfWindowSamples = round(currentWindowDuration / (2 * baseStepTime));

    % Compute the first index of the window.
    startIndex = centreIndex - halfWindowSamples;

    % Compute the last index of the window.
    endIndex = centreIndex + halfWindowSamples;

    % Limit the first index to the available signal.
    startIndex = max(startIndex, 1);

    % Limit the last index to the available signal.
    endIndex = min(endIndex, length(fdBase));

    % Extract the local time vector.
    tWindow = tBase(startIndex:endIndex);

    % Shift the local time to start at zero.
    tWindowLocal = tWindow - tWindow(1);

    % Extract the ground-truth Doppler in the window.
    fdWindow = fdBase(startIndex:endIndex);

    % Extract the noisy Doppler measurements in the window.
    fdMeasWindow = fdMeasDuration(startIndex:endIndex);

    % Extract the ground-truth Doppler rate in the window.
    fdDotWindow = fdDotBase(startIndex:endIndex);

    % Run the 2-state Kalman filter on the window.
    output2 = linearKF(fdMeasWindow, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 0);

    % Extract the 2-state Doppler estimate.
    fdEst2 = output2.fdEst;

    % Extract the 2-state Doppler-rate estimate.
    fdDotEst2 = output2.fdDotEst;

    % Run the 3-state Kalman filter on the window.
    output3 = linearKF(fdMeasWindow, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 1);

    % Extract the 3-state Doppler estimate.
    fdEst3 = output3.fdEst;

    % Extract the 3-state Doppler-rate estimate.
    fdDotEst3 = output3.fdDotEst;

    % Compute the 2-state Doppler RMSE for the window.
    rmseFd2Duration(i) = sqrt(mean((fdEst2 - fdWindow).^2));

    % Compute the 2-state Doppler-rate RMSE for the window.
    rmseFdDot2Duration(i) = sqrt(mean((fdDotEst2 - fdDotWindow).^2));

    % Compute the 3-state Doppler RMSE for the window.
    rmseFd3Duration(i) = sqrt(mean((fdEst3 - fdWindow).^2));

    % Compute the 3-state Doppler-rate RMSE for the window.
    rmseFdDot3Duration(i) = sqrt(mean((fdDotEst3 - fdDotWindow).^2));

    % Initialise the local convergence time of the 2-state filter.
    convTime2Duration(i) = NaN;

    % Search for 2-state convergence within the current window.
    for k = 1:length(fdWindow)
        fdTailRmse = sqrt(mean((fdEst2(k:end) - fdWindow(k:end)).^2));
        fdDotTailRmse = sqrt(mean((fdDotEst2(k:end) - fdDotWindow(k:end)).^2));
        if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
            convTime2Duration(i) = tWindowLocal(k);
            break
        end
    end

    % Initialise the local convergence time of the 3-state filter.
    convTime3Duration(i) = NaN;

    % Search for 3-state convergence within the current window.
    for k = 1:length(fdWindow)
        fdTailRmse = sqrt(mean((fdEst3(k:end) - fdWindow(k:end)).^2));
        fdDotTailRmse = sqrt(mean((fdDotEst3(k:end) - fdDotWindow(k:end)).^2));
        if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
            convTime3Duration(i) = tWindowLocal(k);
            break
        end
    end
end

% Create a figure for the window-duration sensitivity.
figure;

% Create the first subplot.
subplot(3,1,1)

% Plot the 2-state Doppler RMSE versus window duration.
plot(windowDurationVec, rmseFd2Duration, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state Doppler RMSE versus window duration.
plot(windowDurationVec, rmseFd3Duration, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler RMSE [Hz]')

% Add title.
title('Test 5. Sensitivity to the observation-window duration')

% Add legend.
legend('2-state KF', '3-state KF', 'Location', 'best')

% Create the second subplot.
subplot(3,1,2)

% Plot the 2-state Doppler-rate RMSE versus window duration.
plot(windowDurationVec, rmseFdDot2Duration, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state Doppler-rate RMSE versus window duration.
plot(windowDurationVec, rmseFdDot3Duration, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler-rate RMSE [Hz/s]')

% Create the third subplot.
subplot(3,1,3)

% Plot the 2-state convergence time versus window duration.
plot(windowDurationVec, convTime2Duration, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state convergence time versus window duration.
plot(windowDurationVec, convTime3Duration, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Window duration [s]')

% Label the vertical axis.
ylabel('Convergence time [s]')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 6. SENSITIVITY TO THE SAMPLING PERIOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the list of sampling periods in seconds.
stepTimeVec = [0.25 0.5 1 2.5 5 10];

% Preallocate the 2-state Doppler RMSE vector.
rmseFd2Step = zeros(length(stepTimeVec), 1);

% Preallocate the 2-state Doppler-rate RMSE vector.
rmseFdDot2Step = zeros(length(stepTimeVec), 1);

% Preallocate the 2-state convergence-time vector.
convTime2Step = zeros(length(stepTimeVec), 1);

% Preallocate the 3-state Doppler RMSE vector.
rmseFd3Step = zeros(length(stepTimeVec), 1);

% Preallocate the 3-state Doppler-rate RMSE vector.
rmseFdDot3Step = zeros(length(stepTimeVec), 1);

% Preallocate the 3-state convergence-time vector.
convTime3Step = zeros(length(stepTimeVec), 1);

% Loop over the different sampling periods.
for i = 1:length(stepTimeVec)

    % Extract the current sampling period.
    currentStepTime = stepTimeVec(i);

    % Generate the 2D Doppler model for the current sampling period.
    [tCurrent, fdCurrent, dfdtCurrent, d2fdt2Current] = generateOWDoppler2DModel(currentStepTime, 0, baseH, baseF0, baseFmix, baseLambda0);

    % Convert the current time vector to a column vector.
    tCurrent = tCurrent(:);

    % Build a Doppler-rate ground-truth vector with the same length as fdCurrent.
    fdDotCurrent = [dfdtCurrent(1); dfdtCurrent(:)];

    % Add Gaussian noise to the current measurements.
    fdMeasCurrent = addNoise(fdCurrent, baseSigmaNoise);

    % Run the 2-state Kalman filter for the current sampling period.
    output2 = linearKF(fdMeasCurrent, currentStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Current, 0);

    % Extract the 2-state Doppler estimate.
    fdEst2 = output2.fdEst;

    % Extract the 2-state Doppler-rate estimate.
    fdDotEst2 = output2.fdDotEst;

    % Run the 3-state Kalman filter for the current sampling period.
    output3 = linearKF(fdMeasCurrent, currentStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Current, 1);

    % Extract the 3-state Doppler estimate.
    fdEst3 = output3.fdEst;

    % Extract the 3-state Doppler-rate estimate.
    fdDotEst3 = output3.fdDotEst;

    % Compute the 2-state Doppler RMSE.
    rmseFd2Step(i) = sqrt(mean((fdEst2 - fdCurrent).^2));

    % Compute the 2-state Doppler-rate RMSE.
    rmseFdDot2Step(i) = sqrt(mean((fdDotEst2 - fdDotCurrent).^2));

    % Compute the 3-state Doppler RMSE.
    rmseFd3Step(i) = sqrt(mean((fdEst3 - fdCurrent).^2));

    % Compute the 3-state Doppler-rate RMSE.
    rmseFdDot3Step(i) = sqrt(mean((fdDotEst3 - fdDotCurrent).^2));

    % Initialise the current 2-state convergence time.
    convTime2Step(i) = NaN;

    % Search the 2-state convergence time.
    for k = 1:length(fdCurrent)
        fdTailRmse = sqrt(mean((fdEst2(k:end) - fdCurrent(k:end)).^2));
        fdDotTailRmse = sqrt(mean((fdDotEst2(k:end) - fdDotCurrent(k:end)).^2));
        if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
            convTime2Step(i) = tCurrent(k);
            break
        end
    end

    % Initialise the current 3-state convergence time.
    convTime3Step(i) = NaN;

    % Search the 3-state convergence time.
    for k = 1:length(fdCurrent)
        fdTailRmse = sqrt(mean((fdEst3(k:end) - fdCurrent(k:end)).^2));
        fdDotTailRmse = sqrt(mean((fdDotEst3(k:end) - fdDotCurrent(k:end)).^2));
        if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
            convTime3Step(i) = tCurrent(k);
            break
        end
    end
end

% Create a figure for the sampling-period sensitivity.
figure;

% Create the first subplot.
subplot(3,1,1)

% Plot the 2-state Doppler RMSE versus sampling period.
plot(stepTimeVec, rmseFd2Step, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state Doppler RMSE versus sampling period.
plot(stepTimeVec, rmseFd3Step, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler RMSE [Hz]')

% Add title.
title('Test 6. Sensitivity to the sampling period')

% Add legend.
legend('2-state KF', '3-state KF', 'Location', 'best')

% Create the second subplot.
subplot(3,1,2)

% Plot the 2-state Doppler-rate RMSE versus sampling period.
plot(stepTimeVec, rmseFdDot2Step, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state Doppler-rate RMSE versus sampling period.
plot(stepTimeVec, rmseFdDot3Step, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler-rate RMSE [Hz/s]')

% Create the third subplot.
subplot(3,1,3)

% Plot the 2-state convergence time versus sampling period.
plot(stepTimeVec, convTime2Step, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state convergence time versus sampling period.
plot(stepTimeVec, convTime3Step, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Sampling period [s]')

% Label the vertical axis.
ylabel('Convergence time [s]')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 8. SENSITIVITY TO THE FILTER TUNING USING ONLY THE Q/R FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the list of Q/R fractions to analyse.
qOverRVec = [1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1 3 10 30 100];

% Build one noisy measurement vector for this test.
fdMeasTuning = addNoise(fdBase, baseSigmaNoise);

% Preallocate the 2-state Doppler RMSE vector.
rmseFd2Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 2-state Doppler-rate RMSE vector.
rmseFdDot2Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 2-state convergence-time vector.
convTime2Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 3-state Doppler RMSE vector.
rmseFd3Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 3-state Doppler-rate RMSE vector.
rmseFdDot3Tuning = zeros(length(qOverRVec), 1);

% Preallocate the 3-state convergence-time vector.
convTime3Tuning = zeros(length(qOverRVec), 1);

% Loop over the different Q/R fractions.
for i = 1:length(qOverRVec)

    % Extract the current Q/R ratio.
    currentQoverR = qOverRVec(i);

    % Convert the Q/R ratio into the process-noise scale expected by linearKF.
    currentProcessNoiseScale = currentQoverR * baseSigmaNoise^2;

    % Run the 2-state Kalman filter with the current Q/R ratio.
    output2 = linearKF(fdMeasTuning, baseStepTime, currentProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 0);

    % Extract the 2-state Doppler estimate.
    fdEst2 = output2.fdEst;

    % Extract the 2-state Doppler-rate estimate.
    fdDotEst2 = output2.fdDotEst;

    % Run the 3-state Kalman filter with the current Q/R ratio.
    output3 = linearKF(fdMeasTuning, baseStepTime, currentProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 1);

    % Extract the 3-state Doppler estimate.
    fdEst3 = output3.fdEst;

    % Extract the 3-state Doppler-rate estimate.
    fdDotEst3 = output3.fdDotEst;

    % Compute the 2-state Doppler RMSE.
    rmseFd2Tuning(i) = sqrt(mean((fdEst2 - fdBase).^2));

    % Compute the 2-state Doppler-rate RMSE.
    rmseFdDot2Tuning(i) = sqrt(mean((fdDotEst2 - fdDotBase).^2));

    % Compute the 3-state Doppler RMSE.
    rmseFd3Tuning(i) = sqrt(mean((fdEst3 - fdBase).^2));

    % Compute the 3-state Doppler-rate RMSE.
    rmseFdDot3Tuning(i) = sqrt(mean((fdDotEst3 - fdDotBase).^2));

    % Initialise the 2-state convergence time.
    convTime2Tuning(i) = NaN;

    % Search the 2-state convergence time.
    for k = 1:length(fdBase)
        fdTailRmse = sqrt(mean((fdEst2(k:end) - fdBase(k:end)).^2));
        fdDotTailRmse = sqrt(mean((fdDotEst2(k:end) - fdDotBase(k:end)).^2));
        if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
            convTime2Tuning(i) = tBase(k);
            break
        end
    end

    % Initialise the 3-state convergence time.
    convTime3Tuning(i) = NaN;

    % Search the 3-state convergence time.
    for k = 1:length(fdBase)
        fdTailRmse = sqrt(mean((fdEst3(k:end) - fdBase(k:end)).^2));
        fdDotTailRmse = sqrt(mean((fdDotEst3(k:end) - fdDotBase(k:end)).^2));
        if fdTailRmse <= dopplerRmseThreshold && fdDotTailRmse <= dopplerRateRmseThreshold
            convTime3Tuning(i) = tBase(k);
            break
        end
    end
end

% Create a figure for the tuning sensitivity.
figure;

% Create the first subplot.
subplot(3,1,1)

% Plot the 2-state Doppler RMSE versus Q/R.
semilogx(qOverRVec, rmseFd2Tuning, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state Doppler RMSE versus Q/R.
semilogx(qOverRVec, rmseFd3Tuning, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler RMSE [Hz]')

% Add title.
title('Test 8. Sensitivity to the Q/R fraction')

% Add legend.
legend('2-state KF', '3-state KF', 'Location', 'best')

% Create the second subplot.
subplot(3,1,2)

% Plot the 2-state Doppler-rate RMSE versus Q/R.
semilogx(qOverRVec, rmseFdDot2Tuning, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state Doppler-rate RMSE versus Q/R.
semilogx(qOverRVec, rmseFdDot3Tuning, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the vertical axis.
ylabel('Doppler-rate RMSE [Hz/s]')

% Create the third subplot.
subplot(3,1,3)

% Plot the 2-state convergence time versus Q/R.
semilogx(qOverRVec, convTime2Tuning, '--o', 'LineWidth', 1.2)

% Keep the figure active.
hold on

% Plot the 3-state convergence time versus Q/R.
semilogx(qOverRVec, convTime3Tuning, '-s', 'LineWidth', 1.2)

% Add grid.
grid on

% Label the horizontal axis.
xlabel('Q/R fraction')

% Label the vertical axis.
ylabel('Convergence time [s]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRA TEST. HEATMAP OF REGION VERSUS WINDOW DURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the list of window-centre times in seconds.
heatmapCentreTimeVec = 250:100:1050;

% Define the list of window durations in seconds.
heatmapWindowDurationVec = [50 100 150 200 250 300 350 400];

% Build one noisy measurement vector for this test.
fdMeasHeatmap = addNoise(fdBase, baseSigmaNoise);

% Preallocate the 2-state RMSE heatmap.
rmseFd2Heatmap = zeros(length(heatmapWindowDurationVec), length(heatmapCentreTimeVec));

% Preallocate the 3-state RMSE heatmap.
rmseFd3Heatmap = zeros(length(heatmapWindowDurationVec), length(heatmapCentreTimeVec));

% Loop over the different window durations.
for i = 1:length(heatmapWindowDurationVec)

    % Extract the current window duration.
    currentWindowDuration = heatmapWindowDurationVec(i);

    % Loop over the different centre times.
    for j = 1:length(heatmapCentreTimeVec)

        % Extract the current centre time.
        currentCentreTime = heatmapCentreTimeVec(j);

        % Find the index nearest to the current centre time.
        [~, centreIndex] = min(abs(tBase - currentCentreTime));

        % Compute half of the window length in samples.
        halfWindowSamples = round(currentWindowDuration / (2 * baseStepTime));

        % Compute the first index of the window.
        startIndex = centreIndex - halfWindowSamples;

        % Compute the last index of the window.
        endIndex = centreIndex + halfWindowSamples;

        % Limit the first index to the available signal.
        startIndex = max(startIndex, 1);

        % Limit the last index to the available signal.
        endIndex = min(endIndex, length(fdBase));

        % Extract the ground-truth Doppler in the window.
        fdWindow = fdBase(startIndex:endIndex);

        % Extract the noisy Doppler measurements in the window.
        fdMeasWindow = fdMeasHeatmap(startIndex:endIndex);

        % Run the 2-state Kalman filter on the window.
        output2 = linearKF(fdMeasWindow, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 0);

        % Extract the 2-state Doppler estimate.
        fdEst2 = output2.fdEst;

        % Run the 3-state Kalman filter on the window.
        output3 = linearKF(fdMeasWindow, baseStepTime, baseProcessNoiseScale, baseSigmaNoise, d2fdt2Base, 1);

        % Extract the 3-state Doppler estimate.
        fdEst3 = output3.fdEst;

        % Store the 2-state Doppler RMSE in the heatmap matrix.
        rmseFd2Heatmap(i, j) = sqrt(mean((fdEst2 - fdWindow).^2));

        % Store the 3-state Doppler RMSE in the heatmap matrix.
        rmseFd3Heatmap(i, j) = sqrt(mean((fdEst3 - fdWindow).^2));
    end
end

% Create a figure for the 2-state heatmap.
figure;

% Plot the 2-state RMSE heatmap.
imagesc(heatmapCentreTimeVec, heatmapWindowDurationVec, rmseFd2Heatmap)

% Ensure the vertical axis grows upwards.
set(gca, 'YDir', 'normal')

% Add a colour bar.
colorbar

% Label the horizontal axis.
xlabel('Window centre time [s]')

% Label the vertical axis.
ylabel('Window duration [s]')

% Add title.
title('Extra test. 2-state KF Doppler RMSE heatmap')

% Create a figure for the 3-state heatmap.
figure;

% Plot the 3-state RMSE heatmap.
imagesc(heatmapCentreTimeVec, heatmapWindowDurationVec, rmseFd3Heatmap)

% Ensure the vertical axis grows upwards.
set(gca, 'YDir', 'normal')

% Add a colour bar.
colorbar

% Label the horizontal axis.
xlabel('Window centre time [s]')

% Label the vertical axis.
ylabel('Window duration [s]')

% Add title.
title('Extra test. 3-state KF Doppler RMSE heatmap')
