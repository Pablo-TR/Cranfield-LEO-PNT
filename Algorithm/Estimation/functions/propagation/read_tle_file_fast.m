function [tleData, hasSGP4] = read_tle_file_fast(fn, useSGP4)
%READ_TLE_FILE_FAST Read TLE records and precompute fields used by the report.

fid = fopen(fn, 'r');
if fid < 0
    error('Cannot open: %s', fn);
end

L = {};
while ~feof(fid)
    s = strtrim(fgetl(fid));
    if ~isempty(s) && ischar(s)
        L{end+1} = s; %#ok<AGROW>
    end
end
fclose(fid);

hasSGP4 = useSGP4 && exist('twoline2rv', 'file') == 2 && exist('sgp4', 'file') == 2;

tleData = struct('name', {}, 'line1', {}, 'line2', {}, 'epochUTC', {}, ...
                 'raan', {}, 'satrec', {}, 'satrec_data', {});

i = 1;
c = 0;
n = numel(L);

while i <= n-2
    l1 = L{i+1};
    l2 = L{i+2};

    if length(l1) >= 69 && l1(1) == '1' && length(l2) >= 69 && l2(1) == '2'
        c = c + 1;
        nm = strtrim(L{i});
        if startsWith(nm, '0 ')
            nm = nm(3:end);
        end

        tleData(c).name  = nm;
        tleData(c).line1 = l1;
        tleData(c).line2 = l2;

        try
            ey = str2double(l1(19:20));
            ed = str2double(l1(21:32));
            if ey < 57
                yr = 2000 + ey;
            else
                yr = 1900 + ey;
            end
            tleData(c).epochUTC = datetime(yr, 1, 1, 'TimeZone', 'UTC') + days(ed-1);
        catch
            tleData(c).epochUTC = NaT('TimeZone', 'UTC');
        end

        try
            mm = str2double(l2(53:63));
            nr = mm * 2*pi / 86400;
            mu = 3.986004418e14;

            tleData(c).raan = str2double(l2(18:25));
            tleData(c).satrec_data = struct( ...
                'a',       (mu/nr^2)^(1/3), ...
                'ecc',     str2double(['0.' l2(27:33)]), ...
                'incl',    str2double(l2(9:16)), ...
                'raan',    str2double(l2(18:25)), ...
                'argp',    str2double(l2(35:42)), ...
                'M0',      str2double(l2(44:51)), ...
                'n_rad_s', nr);
        catch
            tleData(c).satrec_data = [];
        end

        tleData(c).satrec = [];
        if hasSGP4
            try
                [~, ~, ~, tleData(c).satrec] = twoline2rv(l1, l2, 'c', 'd', 'i', 72);
            catch
                tleData(c).satrec = [];
            end
        end

        i = i + 3;
    else
        i = i + 1;
    end
end
end
