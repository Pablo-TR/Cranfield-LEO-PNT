%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                            getSignalByID.m
%
% Copyright: Cranfield University
% Project Name: LEO-PNT
% Module: Acquisition
% Author: Alfredo
% Date: 7th March 2026
% Last Update: 7th March 2026
% Version: 1.0.0
%
% Description: Retrieves signal metadata associated with a given signal ID
% from the Excel database file "signalData.xlsx". The Excel file is expected
% to contain a column named "ID" used as a unique identifier for each signal
% entry. The function searches for the row matching the provided ID and
% returns all parameters contained in that row as a MATLAB structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function signalData = getSignalByID(signalID)
%GETROWBYID Read an XLSX file, find a row by ID, and return that row's data.
%
% rowData = GETROWBYID(xlsxFile, signalID)
% rowData = GETROWBYID(xlsxFile, signalID, sheetName)
%
% Inputs:
%   targetID  - ID value to search for
%
% Output:
%   rowData   - Struct containing all columns from the matching row

    T = readtable("signalData.xlsx");

    % Check that ID column exists
    if ~ismember('ID', T.Properties.VariableNames)
        error('The Excel file must contain a column named "ID".');
    end

    % Extract ID column
    idColumn = T.ID;

    % Match target ID robustly for numeric/string/cellstr cases
    if isnumeric(idColumn)
        matchIdx = (idColumn == signalID);
    elseif isstring(idColumn)
        matchIdx = (idColumn == string(signalID));
    elseif iscell(idColumn)
        matchIdx = strcmp(string(idColumn), string(signalID));
    else
        matchIdx = strcmp(string(idColumn), string(signalID));
    end

    % Check results
    if ~any(matchIdx)
        error('No row found with ID = %s.', string(signalID));
    end

    if nnz(matchIdx) > 1
        warning('Multiple rows found with ID = %s. Returning the first match.', string(signalID));
    end

    % Get first matching row
    matchedRow = T(find(matchIdx, 1), :);

    % Convert table row to struct
    signalData = table2struct(matchedRow);

end