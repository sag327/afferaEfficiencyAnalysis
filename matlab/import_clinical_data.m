function tbl = import_clinical_data(rawFile, sheetName)
% IMPORT_CLINICAL_DATA  Load clinical Excel export and build analysis-ready table.
%
%   tbl = import_clinical_data() opens a file-selection dialog so you can
%   choose the clinical Excel file, then extracts the procedure-level rows
%   and returns a tidy MATLAB table with derived Affera indicators.
%
%   tbl = import_clinical_data(rawFile, sheetName) lets you specify the
%   source file path and sheet name directly (defaults to Sheet1).
%
%   Baseline-learning metrics are computed in the analysis script, not here.
%
% OUTPUT VARIABLES (per procedure row)
%   - procedure_id: string (from 'Case ID')
%   - procedure_date: datetime
%   - operator_name: string (from 'Primary Surgeon')
%   - operator_id: categorical version of operator_name
%   - procedure_type: categorical ('PVI', 'PVI+', etc. from 'Slices by Codes')
%   - supplies_raw: string (from 'Slices by Supplies and Implants')
%   - is_affera: logical (true if supplies mention 'Sphere-9 Catheter')
%   - tool_type: categorical ('affera' vs 'other')
%   - procedure_duration_min: double (from 'Procedure Start to Procedure Complete (Minutes)')
%   - log_duration: log transform of duration (useful for mixed models)
%   - procedure_primary, service, case_location: optional metadata columns
%
% The resulting table is also saved to ./data_processed/procedures_all.mat for reuse.

if nargin < 1 || isempty(rawFile)
    rawFile = "";
else
    rawFile = string(rawFile);
end

if nargin < 2 || isempty(sheetName)
    sheetName = "Sheet1";
else
    sheetName = string(sheetName);
end

if rawFile == ""
    defaultDir = fullfile(pwd, 'clinicalData');
    if usejava('desktop')
        [fname, fpath] = uigetfile({'*.xlsx;*.xls;*.csv', 'Clinical data files (*.xlsx, *.xls, *.csv)'}, ...
            'Select clinical data file', defaultDir);
        if isequal(fname, 0)
            error('No clinical data file selected.');
        end
        rawFile = fullfile(fpath, fname);
    else
        error('No clinical data file specified and GUI selection is unavailable. Provide a file path.');
    end
end

if ~isfile(rawFile)
    % Allow passing just a file name while assuming ./clinicalData/ as base.
    candidate = fullfile('clinicalData', rawFile);
    if isfile(candidate)
        rawFile = candidate;
    else
        error('Raw Excel file not found: %s', rawFile);
    end
end

rawFile = char(rawFile);
sheetName = char(sheetName);

fprintf('Importing clinical data from %s (sheet %s)...\n', rawFile, sheetName);

%% 1. Detect the header row that contains "Case ID"
% The exported worksheet has metadata in the first several rows; the actual
% header row contains "Case ID" in column 5. We detect it automatically so
% the script stays robust to future exports.
allCells = readcell(rawFile, 'Sheet', sheetName);
caseIdCol = 5;
headerRow = find(strcmp(allCells(:, caseIdCol), 'Case ID'), 1, 'first');
if isempty(headerRow)
    error('Could not find "Case ID" header in column %d of %s.', caseIdCol, rawFile);
end

% Build a range that starts at the header row and spans all 11 columns.
excelRange = sprintf('A%d:K1048576', headerRow);
opts = detectImportOptions(rawFile, 'Sheet', sheetName, ...
    'DataRange', excelRange, 'PreserveVariableNames', true);

% Read the rectangular range; the first row becomes variable names.
rawTbl = readtable(rawFile, opts);
rawTbl = standardizeMissing(rawTbl, {'', 'NA', 'N/A'});

%% 2. Remove empty rows and ensure Case ID is present
caseIdRaw = rawTbl.("Case ID");
maskMissing = ismissing(caseIdRaw);
rawTbl = rawTbl(~maskMissing, :);
caseIdRaw = caseIdRaw(~maskMissing);
caseIdStr = strtrim(string(caseIdRaw));
maskEmpty = caseIdStr == "";
rawTbl = rawTbl(~maskEmpty, :);
caseIdStr = caseIdStr(~maskEmpty);

%% 3. Extract and transform relevant columns
% Procedure date (use the 'Date' column, not 'Start Date' / 'End Date').
dateRaw = rawTbl.("Date");
if isdatetime(dateRaw)
    procedureDate = dateRaw;
elseif isnumeric(dateRaw)
    procedureDate = datetime(dateRaw, 'ConvertFrom', 'excel');
else
    procedureDate = datetime(string(dateRaw));
end
procedureDate.Format = 'yyyy-MM-dd';

% Operator names
operatorName = strtrim(string(rawTbl.("Primary Surgeon")));
operatorName(operatorName == "") = "UnknownOperator";

% Procedure type (PVI vs PVI+)
procedureType = strtrim(string(rawTbl.("Slices by Codes")));
procedureType(procedureType == "") = "UnknownCode";

% Supplies (used to flag Affera/PFA/RF and detect undefined catheters).
suppliesRaw = strtrim(string(rawTbl.("Slices by Supplies and Implants")));
suppliesRaw(suppliesRaw == "") = "UnknownSupply";

% Map supplies to catheter labels and exposure flags (handles new and legacy formats).
[catheter_primary, is_affera, is_pfa_catheter, is_rf_catheter, catheter_defined] = ...
    parse_supplies_column(suppliesRaw);

% Procedure duration in minutes
durationRaw = rawTbl.("Procedure Start to Procedure Complete (Minutes)");
procedureDuration = double(durationRaw);

% Additional metadata (optional)
procedurePrimary = strtrim(string(rawTbl.("Procedure (Primary)")));
procedurePrimary(procedurePrimary == "") = "UnknownProcedure";
service          = strtrim(string(rawTbl.Service));
service(service == "") = "UnknownService";
caseLocation     = strtrim(string(rawTbl.("Case Location")));
caseLocation(caseLocation == "") = "UnknownLocation";

% Drop undefined catheter rows (e.g., "None of the above") to avoid length mismatches.
if any(~catheter_defined)
    maskKeep = catheter_defined;
    fprintf('Excluding %d procedures with undefined catheter (e.g., "None of the above").\n', ...
        sum(~maskKeep));
    caseIdStr         = caseIdStr(maskKeep);
    procedureDate     = procedureDate(maskKeep);
    operatorName      = operatorName(maskKeep);
    procedureType     = procedureType(maskKeep);
    suppliesRaw       = suppliesRaw(maskKeep);
    catheter_primary  = catheter_primary(maskKeep);
    is_affera         = is_affera(maskKeep);
    is_pfa_catheter   = is_pfa_catheter(maskKeep);
    is_rf_catheter    = is_rf_catheter(maskKeep);
    catheter_defined  = catheter_defined(maskKeep);
    durationRaw       = durationRaw(maskKeep);
    procedureDuration = procedureDuration(maskKeep);
    procedurePrimary  = procedurePrimary(maskKeep);
    service           = service(maskKeep);
    caseLocation      = caseLocation(maskKeep);
end

%% 4. Assemble the analysis-ready table
tbl = table;
tbl.procedure_id          = caseIdStr;
tbl.procedure_date        = procedureDate;
tbl.operator_name         = operatorName;
tbl.operator_id           = categorical(operatorName);
tbl.procedure_type        = categorical(procedureType);
tbl.supplies_raw          = suppliesRaw;
tbl.catheter_primary      = categorical(catheter_primary);
tbl.is_affera             = is_affera;
tbl.is_pfa_catheter       = is_pfa_catheter;
tbl.is_rf_catheter        = is_rf_catheter;
tbl.catheter_defined      = catheter_defined;
tbl.exclude_from_analysis = ~catheter_defined;

toolType = repmat("other", height(tbl), 1);
toolType(tbl.is_affera) = "affera";
tbl.tool_type = categorical(toolType);

% Affera launch timing (global earliest Affera case).
affera_launch_date_global = NaT;
if any(tbl.is_affera)
    affera_launch_date_global = min(tbl.procedure_date(tbl.is_affera));
    tbl.time_since_affera_launch = days(tbl.procedure_date - affera_launch_date_global);
    tbl.is_pre_affera_global = tbl.procedure_date < affera_launch_date_global;
    tbl.is_post_affera_global = ~tbl.is_pre_affera_global;
else
    warning('No Affera cases detected; affera_launch_date_global is undefined.');
    tbl.time_since_affera_launch = nan(height(tbl), 1);
    tbl.is_pre_affera_global = false(height(tbl), 1);
    tbl.is_post_affera_global = false(height(tbl), 1);
end

tbl.procedure_duration_min = procedureDuration;
tbl.log_duration           = log(tbl.procedure_duration_min);

tbl.procedure_primary = procedurePrimary;
tbl.service           = service;
tbl.case_location     = caseLocation;

% Operator-level Affera/non-Affera counts for reporting.
ops = categories(tbl.operator_id);
numOps = numel(ops);
operator_case_counts = table;
operator_case_counts.operator_id = categorical(ops);
operator_case_counts.n_affera_cases = zeros(numOps, 1);
operator_case_counts.n_non_affera_cases = zeros(numOps, 1);
operator_case_counts.n_rf_cases = zeros(numOps, 1);
for k = 1:numOps
    op = ops{k};
    idx = (tbl.operator_id == op);
    operator_case_counts.n_affera_cases(k) = sum(idx & tbl.is_affera);
    operator_case_counts.n_non_affera_cases(k) = sum(idx & ~tbl.is_affera);
    operator_case_counts.n_rf_cases(k) = sum(idx & tbl.is_rf_catheter);
end

fprintf('\nOperator case counts (Affera vs non-Affera):\n');
disp(operator_case_counts);

% Sort rows for downstream analyses (per-operator chronological order).
tbl = sortrows(tbl, {'operator_id', 'procedure_date', 'procedure_id'});

% Metadata for downstream scripts.
meta = struct();
meta.affera_launch_date_global = affera_launch_date_global;
meta.operator_case_counts = operator_case_counts;

%% 5. Persist the table for reuse
outDir = fullfile('data_processed');
if ~isfolder(outDir)
    mkdir(outDir);
end
outFile = fullfile(outDir, 'procedures_all.mat');
save(outFile, 'tbl', 'meta');

fprintf('Saved %d procedures (%d operators, %d Affera cases) to %s\n', ...
    height(tbl), numel(categories(tbl.operator_id)), sum(tbl.is_affera), outFile);

end

% -------------------------------------------------------------------------
function [catheter_primary, is_affera, is_pfa_catheter, is_rf_catheter, catheter_defined] = parse_supplies_column(suppliesRaw)
% PARSE_SUPPLIES_COLUMN  Map supplies strings to catheter labels and flags.
%
% Handles both legacy and new coding:
%   - New: 'Sphere-9', 'RF', 'PulseSelect', 'Farawave', 'Cryo', 'None of the above'
%   - Legacy: 'Sphere-9 Catheter', 'Other catheters (RF)', 'PulseSelect', 'Farawave'

delims = {',', ';', '|', newline};
n = numel(suppliesRaw);
catheter_primary   = repmat("UNKNOWN", n, 1);
is_affera          = false(n, 1);
is_pfa_catheter    = false(n, 1);
is_rf_catheter     = false(n, 1);
catheter_defined   = false(n, 1);

for ii = 1:n
    entry = string(suppliesRaw(ii));
    lowEntry = lower(entry);

    % Exclude undefined catheter rows.
    if contains(lowEntry, 'none of the above')
        catheter_primary(ii) = "UNKNOWN";
        continue;
    end

    tokens = split(entry, delims);
    tokens = strtrim(tokens);
    tokens(tokens == "") = [];
    tokensLow = lower(tokens);

    recognized = false;

    if any(contains(tokensLow, 'sphere-9'))
        catheter_primary(ii) = "SPHERE-9";
        is_affera(ii) = true;
        is_pfa_catheter(ii) = true;
        recognized = true;
    elseif any(contains(tokensLow, 'pulseselect'))
        catheter_primary(ii) = "PULSESELECT";
        is_pfa_catheter(ii) = true;
        recognized = true;
    elseif any(contains(tokensLow, 'farawave'))
        catheter_primary(ii) = "FARAWAVE";
        is_pfa_catheter(ii) = true;
        recognized = true;
    elseif any(contains(tokensLow, 'other catheters (rf)') | strcmpi(tokensLow, 'rf'))
        catheter_primary(ii) = "RF";
        is_rf_catheter(ii) = true;
        recognized = true;
    elseif any(contains(tokensLow, 'cryo'))
        catheter_primary(ii) = "CRYO";
        recognized = true;
    end

    if ~recognized && ~isempty(tokens)
        % Legacy fallback: first token as label.
        catheter_primary(ii) = upper(tokens(1));
        recognized = true;
    end

    catheter_defined(ii) = recognized;
end

end
