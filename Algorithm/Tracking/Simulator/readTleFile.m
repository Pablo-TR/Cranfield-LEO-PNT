%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       readTleFile.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Tracking
% Author: Pablo
% Date: 14th March 2026
% Last Update: 14th March 2026
% Version: 1.0.0
%
% Description: Standalone TLE reader extracted from the Estimation
% tooling. It parses named or unnamed TLE entries and returns the orbital
% elements needed by the lightweight Doppler truth generator.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tleData, meta] = readTleFile(tleFile, varargin)

    parser = inputParser;
    addParameter(parser, 'useSGP4', false);
    parse(parser, varargin{:});

    tleFile = char(string(tleFile));
    fid = fopen(tleFile, 'r');
    if fid < 0
        error('Cannot open TLE file: %s', tleFile);
    end

    lines = cell(0, 1);
    while ~feof(fid)
        currentLine = fgetl(fid);
        if ~ischar(currentLine)
            continue;
        end

        currentLine = strtrim(currentLine);
        if ~isempty(currentLine)
            lines{end + 1, 1} = currentLine; %#ok<AGROW>
        end
    end
    fclose(fid);

    hasSGP4 = logical(parser.Results.useSGP4) && ...
        exist('twoline2rv', 'file') == 2 && ...
        exist('sgp4', 'file') == 2;

    tleData = struct('name', {}, 'line1', {}, 'line2', {}, ...
        'epochUTC', {}, 'raan', {}, 'satrec', {}, 'satrec_data', {});

    entryCount = 0;
    lineIndex = 1;
    while lineIndex <= numel(lines)
        [nameLine, line1, line2, nextIndex, isValidEntry] = ...
            getNextTleEntry(lines, lineIndex);

        if ~isValidEntry
            lineIndex = lineIndex + 1;
            continue;
        end

        entryCount = entryCount + 1;
        tleData(entryCount).name = nameLine;
        tleData(entryCount).line1 = line1;
        tleData(entryCount).line2 = line2;
        tleData(entryCount).epochUTC = parseTleEpoch(line1);
        [tleData(entryCount).raan, tleData(entryCount).satrec_data] = ...
            parseTleOrbit(line2);
        tleData(entryCount).satrec = [];

        if hasSGP4
            try
                [~, ~, ~, tleData(entryCount).satrec] = ...
                    twoline2rv(line1, line2, 'c', 'd', 'i', 72);
            catch
                tleData(entryCount).satrec = [];
            end
        end

        lineIndex = nextIndex;
    end

    meta = struct( ...
        'hasSGP4', hasSGP4, ...
        'numEntries', numel(tleData), ...
        'sourceFile', string(tleFile));

end

function [nameLine, line1, line2, nextIndex, isValidEntry] = ...
    getNextTleEntry(lines, lineIndex)

    nameLine = '';
    line1 = '';
    line2 = '';
    nextIndex = lineIndex + 1;
    isValidEntry = false;

    if lineIndex <= numel(lines) - 1 && ...
            isTleLine1(lines{lineIndex}) && ...
            isTleLine2(lines{lineIndex + 1})
        nameLine = sprintf('SAT-%03d', lineIndex);
        line1 = lines{lineIndex};
        line2 = lines{lineIndex + 1};
        nextIndex = lineIndex + 2;
        isValidEntry = true;
        return
    end

    if lineIndex <= numel(lines) - 2 && ...
            isTleLine1(lines{lineIndex + 1}) && ...
            isTleLine2(lines{lineIndex + 2})
        nameLine = lines{lineIndex};
        if startsWith(nameLine, '0 ')
            nameLine = nameLine(3:end);
        end
        line1 = lines{lineIndex + 1};
        line2 = lines{lineIndex + 2};
        nextIndex = lineIndex + 3;
        isValidEntry = true;
    end

end

function isValid = isTleLine1(inputLine)

    isValid = ischar(inputLine) && numel(inputLine) >= 2 && inputLine(1) == '1';

end

function isValid = isTleLine2(inputLine)

    isValid = ischar(inputLine) && numel(inputLine) >= 2 && inputLine(1) == '2';

end

function epochUTC = parseTleEpoch(line1)

    try
        epochYear = str2double(line1(19:20));
        epochDay = str2double(line1(21:32));

        if epochYear < 57
            fullYear = 2000 + epochYear;
        else
            fullYear = 1900 + epochYear;
        end

        epochUTC = datetime(fullYear, 1, 1, 'TimeZone', 'UTC') + ...
            days(epochDay - 1);
    catch
        epochUTC = NaT('TimeZone', 'UTC');
    end

end

function [raanDeg, satrecData] = parseTleOrbit(line2)

    raanDeg = NaN;
    satrecData = [];

    try
        meanMotionRevPerDay = str2double(line2(53:63));
        meanMotionRadPerSec = meanMotionRevPerDay * 2 * pi / 86400;
        mu = 3.986004418e14;

        raanDeg = str2double(line2(18:25));
        satrecData = struct( ...
            'a', (mu / meanMotionRadPerSec^2)^(1 / 3), ...
            'ecc', str2double(['0.' line2(27:33)]), ...
            'incl', str2double(line2(9:16)), ...
            'raan', raanDeg, ...
            'argp', str2double(line2(35:42)), ...
            'M0', str2double(line2(44:51)), ...
            'n_rad_s', meanMotionRadPerSec);
    catch
        satrecData = [];
    end

end

