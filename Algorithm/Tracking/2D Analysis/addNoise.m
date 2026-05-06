%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               addNoise.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 5th May 2026
% Last Update: 5th May 2026
% Version: 1.0.0
%
% Description: Adds Gaussian noise to measurements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fd_noise] = addNoise(fd, sigma_fd)
    % Add zero-mean Gaussian noise with the requested standard deviation.
    fd_noise = fd + sigma_fd * randn(size(fd));
end
