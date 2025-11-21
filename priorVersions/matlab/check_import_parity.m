% check_import_parity.m
%
% Simple parity test to verify that the clinical Excel import and the
% processed MATLAB table are consistent on key counts.
%
% Assumes:
%   - Raw Excel file is in ./clinicalData (adjust file name if needed).
%   - Processed table 'tbl' is saved in ./data_processed/procedures_all.mat.

clear; clc;

%% 1. Load raw Excel data (same sheet/range logic as the import)
dataDir  = fullfile(pwd, 'clinicalData');
fileName = 'clinicalData-11-18-2025.xlsx';  % update if your file name changes
filePath = fullfile(dataDir, fileName);
sheetName = 'Sheet1';

if ~isfile(filePath)
    error('Raw Excel file not found at %s', filePath);
end

allCells = readcell(filePath, 'Sheet', sheetName);
headerRow = find(strcmp(allCells(:, 5), 'Case ID'), 1, 'first');
if isempty(headerRow)
    error('Could not find "Case ID" header row in %s.', filePath);
end

range = sprintf('A%d:K1048576', headerRow);
opts = detectImportOptions(filePath, 'Sheet', sheetName, ...
    'DataRange', range, 'PreserveVariableNames', true);
raw = readtable(filePath, opts);
raw = standardizeMissing(raw, {'', 'NA', 'N/A'});
raw = raw(~ismissing(raw.("Case ID")), :);

%% 2. Load processed MATLAB table
procPath = fullfile(pwd, 'data_processed', 'procedures_all.mat');
if ~isfile(procPath)
    error('Processed file not found at %s. Run the import to create it.', procPath);
end

S = load(procPath);
tbl = S.tbl;

%% 3. Basic row-count parity
n_raw_rows = height(raw);
n_tbl_rows = height(tbl);

fprintf('Raw rows: %d, Processed rows: %d\n', n_raw_rows, n_tbl_rows);

if n_raw_rows ~= n_tbl_rows
    warning('Row count mismatch between raw and processed data.');
end

%% 4. Unique ID parity (caseID vs. procedure_id)
if ismember('Case ID', raw.Properties.VariableNames) && ismember('procedure_id', tbl.Properties.VariableNames)
    n_raw_ids = numel(unique(raw.("Case ID")));
    n_tbl_ids = numel(unique(tbl.procedure_id));
    fprintf('Unique caseID (raw): %d, Unique procedure_id (processed): %d\n', n_raw_ids, n_tbl_ids);
    if n_raw_ids ~= n_tbl_ids
        warning('Unique ID count mismatch between raw and processed data.');
    end
else
    warning('"Case ID" or procedure_id column not found; skipping ID parity check.');
end

%% 5. Parity by procedure type (Slices by Codes vs. procedure_type)
codeVar = 'Slices by Codes';
if ismember(codeVar, raw.Properties.VariableNames) && ismember('procedure_type', tbl.Properties.VariableNames)
    fprintf('\nParity check: procedure type (PVI vs PVI+)\n');

    raw_code_counts = groupsummary(raw, codeVar);
    tbl_type_counts = groupsummary(tbl, 'procedure_type');

    disp('Raw counts by Code:');
    disp(raw_code_counts(:, {codeVar, 'GroupCount'}));

    disp('Processed counts by procedure_type:');
    disp(tbl_type_counts(:, {'procedure_type', 'GroupCount'}));
else
    warning('Procedure type column not found in raw or processed data; skipping.');
end

%% 6. Parity for Affera indicator (Supplies vs. is_affera)
suppliesVar = 'Slices by Supplies and Implants';
if ismember(suppliesVar, raw.Properties.VariableNames) && ismember('is_affera', tbl.Properties.VariableNames)
    fprintf('\nParity check: Affera vs non-Affera (Sphere-9 Catheter).\n');

    is_affera_raw = contains(string(raw.(suppliesVar)), 'Sphere-9 Catheter', 'IgnoreCase', true);
    n_raw_affera  = sum(is_affera_raw);
    n_tbl_affera  = sum(tbl.is_affera);

    fprintf('Raw Affera (Supplies contains Sphere-9 Catheter): %d\n', n_raw_affera);
    fprintf('Processed Affera (is_affera == true): %d\n', n_tbl_affera);

    if n_raw_affera ~= n_tbl_affera
        warning('Affera count mismatch between raw and processed data.');
    end
else
    warning('Supplies or is_affera column not found; skipping Affera parity check.');
end

fprintf('\nParity checks complete.\n');
