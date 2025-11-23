%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            readSignal.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Algorithm
% Author: Pablo, Marti, Dan & Fredo
% Date: 8th November 2025
% Last Update: 9th November 2025
% Version: 1.1.0
%
% Description: Function that reads the binary data from the raw signal file
% after the SDR, it outputs the detrended signal (I and Q values without DC
% drift, hence zero mean) along with the number of I/Q samples (N) and the
% time array. It takes as arguments the filename of the .bin signal along
% with the sampling frequency Fs [Hz] and t which is the interval of time
% of signal to read [ms].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [signalDetrended, N, timeArray] = readSignal(filename, Fs, t)

    N = floor(t * Fs / 1000); % Total number of I/Q samples to read
    %timeArray = 0:t/1000:(N*t/1000 - t/1000);
    timeArray = 0:t/N:t-t/N;  % Time array in ms
    % File Read Parameters
    dataType = 'int16';

    % Binary File Read
    fprintf('Opening file: %s\n', filename);
    [fid, msg] = fopen(filename, 'r');
    if fid < 0
        error('Could not open file "%s" - Error: %s', filename, msg);
    end
    fprintf('Reading %d samples (%d ms)...\n', N, t);

    % Reading binary data and turning into int16 (-32768 to 32767)
    data = fread(fid, [2, N], dataType);
    fclose(fid);
    fprintf('File read complete.\n');
    
    % Conversion int16 to hexadecimal
    %dataVec = data(:);
    %hexData= dec2hex(typecast(int16(dataVec), 'uint16'), 4);
    
    % Checking if fewer samples than expected were read and adjusting 
    if size(data, 2) < N
        warning('Fewer samples read than expected. End of file?');
        N = size(data, 2);
        if N == 0, error('No data read.'); end
    end

    % Signal Processing (Time Domain)
    dataI = double(data(1, :)) / (2^15 - 1); % Normalising to 2^15 = 32768 (max int16)
    dataQ = double(data(2, :)) / (2^15 - 1);
    signal = dataI + 1i * dataQ;
    signalDetrended = detrend(signal); % Removes DC drift
    fprintf('Signal processed (detrend).\n');

end


