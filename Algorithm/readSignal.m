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
% with the sampling frequency Fs [Hz] and tStart and tEnd which define the 
% interval of the time of signal to read [ms].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [signalDetrended, N, timeArray, binData] = readSignal(filename, Fs, tStart, tEnd)
    
    % Validate inputs
    if tEnd <= tStart
        error('tEnd must be greater than tStart.');
    end
    if Fs <= 0
        error('Fs must be > 0.');
    end

    Nskip = floor(tStart * Fs / 1000); % Nskip: number of I/Q sample pairs to skip from the beginning
    tDur = tEnd - tStart;

    N = floor(tDur * Fs / 1000); % Total number of I/Q samples to read
    timeArray = tStart/1000:1/Fs:tEnd/1000 - 1/Fs; % Time array

    % File Read Parameters
    dataType = 'int16';
    bytesPerValue = 2;
    valuesPerSamplePair = 2; % I and Q
    bytesPerSamplePair = bytesPerValue * valuesPerSamplePair;

    % Binary File Read
    fprintf('Opening file: %s\n', filename);
    [fid, msg] = fopen(filename, 'r');
    if fid < 0
        error('Could not open file "%s" - Error: %s', filename, msg);
    end
    fprintf('Reading %d samples (%d - %d ms)...\n', N, tStart, tEnd);
    
    byteOffset = Nskip * bytesPerSamplePair;
    status = fseek(fid, byteOffset, 'bof');
    if status ~= 0
        fclose(fid);
        error('Failed to seek to byte offset %d. Check file size / tStart.', byteOffset);
    end

    % Reading binary data and turning into int16 (-32768 to 32767)
    data = fread(fid, [2, N], dataType);
    fclose(fid);
    fprintf('File read complete.\n');
    
    % Conversion int16 to hexadecimal
    dataVec = data(:);
    binData= dec2bin(typecast(int16(dataVec), 'uint16'), 4);
    
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


