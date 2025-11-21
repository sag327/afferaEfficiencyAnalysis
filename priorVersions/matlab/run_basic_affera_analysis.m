function results = run_basic_affera_analysis(tbl, varargin)
% RUN_BASIC_AFFERA_ANALYSIS  Fit bare-bones Affera learning-curve models.
%
%   RUN_BASIC_AFFERA_ANALYSIS(TBL) accepts a procedure-level table produced
%   by import_clinical_data (one row per procedure) and computes the
%   minimal derived variables plus four mixed-effects models addressing:
%       1) Overall Affera vs. non-Affera duration difference.
%       2) Affera learning curve vs. within-operator case index.
%       3) Interaction between baseline operator speed and the learning curve.
%       4) Interaction between baseline operator speed and Affera vs. non-Affera difference.
%
%   results = RUN_BASIC_AFFERA_ANALYSIS(TBL) returns a struct containing the
%   fitted LinearMixedModel objects (fields lme_q1 ... lme_q4; models 3/4 may
%   be empty if baseline data are unavailable).
%
%   RUN_BASIC_AFFERA_ANALYSIS(..., 'ProcedureType', TYPE) filters the input
%   data to the specified procedure type before modeling. TYPE can be:
%       - 'all'  : PVI and PVI+ together (default)
%       - 'PVI'  : PVI cases only
%       - 'PVI+' : PVI+ cases only
%       - 'sequential' : run the analysis three times in sequence:
%                        (all), then (PVI only), then (PVI+ only),
%                        returning a struct with fields .all, .PVI, .PVIplus.
%
%   RUN_BASIC_AFFERA_ANALYSIS() with no input attempts to load the processed
%   table saved at ./data_processed/procedures_all.mat (variable 'tbl').
%
%   Key arguments:
%     'MinCases' – single per-operator minimum used for question-specific
%                  eligibility thresholds:
%                    Q1: at least MinCases pre-Affera non-Affera cases AND
%                        at least MinCases post-Affera Affera cases.
%                    Q2: at least MinCases Affera cases.
%                    Q3: a valid baseline (from baselineCaseNumberRange)
%                        AND at least MinCases Affera cases.
%                    Q4: a valid baseline AND at least MinCases Affera
%                        AND at least MinCases non-Affera cases.
%                  Operators failing a question’s threshold are excluded
%                  from that question only (not globally).

parser = inputParser;
parser.addParameter('ProcedureType', "all");
parser.addParameter('MinAfferaCases', 10, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('NonPFABaselineOnly', false);
parser.addParameter('BaselinePVIOnly', false);
parser.addParameter('baselineCaseNumberRange', Inf, @(x) isnumeric(x) && isvector(x));
parser.addParameter('Verbose', false, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});
procTypeFilter = parser.Results.ProcedureType;
minAfferaCases = parser.Results.MinAfferaCases;
nonPFABaselineOnly = logical(parser.Results.NonPFABaselineOnly);
baselinePVIOnly = logical(parser.Results.BaselinePVIOnly);
baselineRangeInput = double(parser.Results.baselineCaseNumberRange(:)');
verboseOutput = parser.Results.Verbose;
if isempty(baselineRangeInput)
    baselineRangeInput = [0 Inf];
elseif numel(baselineRangeInput) == 1
    baselineRangeInput = [0 baselineRangeInput];
else
    baselineRangeInput = baselineRangeInput(1:2);
end
baselineMinCasesBaseline = max(0, floor(baselineRangeInput(1)));
baselineMaxCasesBaseline = baselineRangeInput(2);
if ~(isfinite(baselineMaxCasesBaseline) && baselineMaxCasesBaseline > 0)
    baselineMaxCasesBaseline = Inf;
end
if isfinite(baselineMaxCasesBaseline) && baselineMaxCasesBaseline < baselineMinCasesBaseline
    baselineMaxCasesBaseline = baselineMinCasesBaseline;
end
baselineRangeRequested = [baselineMinCasesBaseline baselineMaxCasesBaseline];
if ischar(procTypeFilter) || isstring(procTypeFilter)
    procTypeFilter = string(procTypeFilter);
end
if iscell(procTypeFilter)
    procTypeFilter = string(procTypeFilter);
end
procTypeFilter = strtrim(procTypeFilter);
procTypeFilter(procTypeFilter=="") = "all";

% If ProcedureType == 'sequential', run three analyses in sequence and return.
if ischar(procTypeFilter) || isstring(procTypeFilter)
    if any(strcmpi(string(procTypeFilter), "sequential"))
        if verboseOutput
            fprintf('\nRunning analyses sequentially for: ALL (PVI+PVI+), PVI only, PVI+ only.\n');
        end
        results = struct();
        results.all = run_basic_affera_analysis(tbl, ...
            'ProcedureType', 'all', ...
            'MinAfferaCases', minAfferaCases, ...
            'NonPFABaselineOnly', nonPFABaselineOnly, ...
            'BaselinePVIOnly', baselinePVIOnly, ...
            'baselineCaseNumberRange', baselineRangeRequested, ...
            'Verbose', verboseOutput);
        results.PVI = run_basic_affera_analysis(tbl, ...
            'ProcedureType', 'PVI', ...
            'MinAfferaCases', minAfferaCases, ...
            'NonPFABaselineOnly', nonPFABaselineOnly, ...
            'BaselinePVIOnly', baselinePVIOnly, ...
            'baselineCaseNumberRange', baselineRangeRequested, ...
            'Verbose', verboseOutput);
        results.PVIplus = run_basic_affera_analysis(tbl, ...
            'ProcedureType', 'PVI+', ...
            'MinAfferaCases', minAfferaCases, ...
            'NonPFABaselineOnly', nonPFABaselineOnly, ...
            'BaselinePVIOnly', baselinePVIOnly, ...
            'baselineCaseNumberRange', baselineRangeRequested, ...
            'Verbose', verboseOutput);
        return;
    end
end

meta = struct();
metaAvailable = false;

% Prepare naming info for optional text report.
resultsDir = 'results';
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end
timeStamp = datestr(now, 'yyyymmdd_HHMMSS');

if ~verboseOutput
    warning('Verbose output disabled; only summary and warnings will be shown.');
end

if nargin < 1 || isempty(tbl)
    defaultMat = fullfile('data_processed', 'procedures_all.mat');
    if ~isfile(defaultMat)
        error('No input table provided and default file %s not found. Run import_clinical_data or pass tbl explicitly.', defaultMat);
    end
    S = load(defaultMat);
    if ~isfield(S, 'tbl')
        error('File %s does not contain a variable named ''tbl''.', defaultMat);
    end
    tbl = S.tbl;
    if isfield(S, 'meta')
        meta = S.meta;
        metaAvailable = true;
    end
    if verboseOutput
        fprintf('Loaded %d procedures from %s\n', height(tbl), defaultMat);
    end
else
    if verboseOutput
        fprintf('Running analysis on provided table with %d procedures.\n', height(tbl));
    end
end

% -------------------------------------------------------------------------
% Analysis settings header
% -------------------------------------------------------------------------
fprintf('\n================ Affera Analysis Settings ================\n');
fprintf('Procedure type filter       : %s\n', strjoin(string(procTypeFilter), ', '));
fprintf('MinAfferaCases (per-operator Affera threshold): %g\n', minAfferaCases);
fprintf('  Q1: >=MinAfferaCases post-Affera Affera AND baseline requirement (>=%g pre-Affera non-Affera)\n', baselineMinCasesBaseline);
fprintf('  Q2: >=MinAfferaCases Affera cases\n');
fprintf('  Q3: baseline available AND >=MinAfferaCases Affera cases\n');
fprintf('  Q4: baseline available AND >=MinAfferaCases Affera cases AND >=%g non-Affera cases\n', baselineMinCasesBaseline);
fprintf('Baseline restricted to PVI only          : %s\n', ...
    ternary(baselinePVIOnly, 'YES', 'NO'));
if metaAvailable
    if isfield(meta, 'affera_launch_date_global') && ~isempty(meta.affera_launch_date_global) && ~isnat(meta.affera_launch_date_global)
        fprintf('Global Affera launch date                     : %s\n', datestr(meta.affera_launch_date_global, 'yyyy-mm-dd'));
    end
end
fprintf('Non-PFA baseline only                        : %s\n', ...
    ternary(nonPFABaselineOnly, 'YES', 'NO'));
fprintf('Baseline min cases/operator                  : %g\n', baselineMinCasesBaseline);
if isfinite(baselineMaxCasesBaseline)
    fprintf('Baseline max cases/operator                  : %g\n', baselineMaxCasesBaseline);
else
    fprintf('Baseline max cases/operator                  : no cap (Inf)\n');
end
fprintf('=========================================================\n\n');

    function out = ternary(cond, a, b)
        if cond
            out = a;
        else
            out = b;
        end
end

%% 1. Basic validation and housekeeping
requiredVars = {'procedure_duration_min', 'tool_type', 'operator_id', 'procedure_date', ...
    'procedure_id', 'is_affera', 'is_pre_affera_global', 'is_post_affera_global'};
missingVars = setdiff(requiredVars, tbl.Properties.VariableNames);
if ~isempty(missingVars)
    error('Missing required variables in table: %s', strjoin(missingVars, ', '));
end

if ~iscategorical(tbl.operator_id)
    tbl.operator_id = categorical(tbl.operator_id);
end

if ~iscategorical(tbl.tool_type)
    tbl.tool_type = categorical(tbl.tool_type);
end

if ~isdatetime(tbl.procedure_date)
    tbl.procedure_date = datetime(tbl.procedure_date);
end

tbl.procedure_duration_min = double(tbl.procedure_duration_min);
invalidDuration = isnan(tbl.procedure_duration_min) | ~(tbl.procedure_duration_min > 0);
if any(invalidDuration)
    warning('Removing %d procedures with missing or non-positive duration.', sum(invalidDuration));
    tbl(invalidDuration, :) = [];
end
if isempty(tbl)
    error('No valid procedures remain after filtering invalid durations.');
end
tbl.log_duration = log(tbl.procedure_duration_min);

if ~ismember('is_affera', tbl.Properties.VariableNames)
    tbl.is_affera = (tbl.tool_type == 'affera');
end

% Keep an unfiltered, cleaned copy of the data for baseline calculations.
tbl_all = tbl;

tbl_pre_mask = tbl_all.is_pre_affera_global & ~tbl_all.is_affera;
tbl_post_affera_mask = tbl_all.is_post_affera_global & tbl_all.is_affera;

ops_all = categories(tbl_all.operator_id);
nOps = numel(ops_all);
operator_case_summary = table('Size', [nOps 5], ...
    'VariableTypes', {'categorical', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'operator_id', 'n_pre_non_affera_cases', 'n_post_affera_cases', ...
    'n_affera_cases_total', 'n_non_affera_cases_total'});
operator_case_summary.operator_id = categorical(ops_all);

for k = 1:nOps
    op = ops_all{k};
    idx = (tbl_all.operator_id == op);
    operator_case_summary.n_pre_non_affera_cases(k) = sum(idx & tbl_pre_mask);
    operator_case_summary.n_post_affera_cases(k) = sum(idx & tbl_post_affera_mask);
    operator_case_summary.n_affera_cases_total(k) = sum(idx & tbl_all.is_affera);
    operator_case_summary.n_non_affera_cases_total(k) = sum(idx & ~tbl_all.is_affera);
end

% Baseline operator statistics from equal-sized pre-Affera samples (computed here).
% These are computed from the full (unfiltered) dataset tbl_all so that,
% for example, PVI-only baseline cases can still be used when modeling is
% restricted to PVI+ procedures.
tbl_baseline = tbl_all;
if ~iscategorical(tbl_baseline.operator_id)
    tbl_baseline.operator_id = categorical(tbl_baseline.operator_id);
end
if ~ismember('is_pfa_catheter', tbl_baseline.Properties.VariableNames)
    tbl_baseline.is_pfa_catheter = false(height(tbl_baseline), 1);
end

baseline_eligible = tbl_baseline.is_pre_affera_global;
if baselinePVIOnly && ismember('procedure_type', tbl_baseline.Properties.VariableNames)
    baseline_eligible = baseline_eligible & (tbl_baseline.procedure_type == categorical("PVI"));
end
if nonPFABaselineOnly
    baseline_eligible = baseline_eligible & ~tbl_baseline.is_pfa_catheter;
else
    baseline_eligible = baseline_eligible & ~tbl_baseline.is_affera;
end

baseline_counts = zeros(nOps, 1);
for k = 1:nOps
    op = ops_all{k};
    baseline_counts(k) = sum(baseline_eligible & (tbl_baseline.operator_id == op));
end

eligible_for_threshold = baseline_counts >= baselineMinCasesBaseline;
if any(eligible_for_threshold)
    minPerOperator = min(baseline_counts(eligible_for_threshold));
    baseline_n_cases_global = min(minPerOperator, baselineMaxCasesBaseline);
else
    baseline_n_cases_global = 0;
end

baseline_mean = nan(nOps, 1);

if baseline_n_cases_global > 0
    for k = 1:nOps
        op = ops_all{k};
        idxElig = find(baseline_eligible & (tbl_baseline.operator_id == op));
        if numel(idxElig) >= baseline_n_cases_global
            sub = tbl_baseline(idxElig, {'procedure_date', 'procedure_id'});
            [~, order] = sortrows(sub, {'procedure_date', 'procedure_id'}, {'descend', 'descend'});
            idxSorted = idxElig(order);
            idxTake = idxSorted(1:baseline_n_cases_global);
            baseline_mean(k) = mean(tbl_baseline.procedure_duration_min(idxTake), 'omitnan');
        end
    end
end

baseline_speed = nan(nOps, 1);
if any(~isnan(baseline_mean))
    mu = mean(baseline_mean, 'omitnan');
    sigma = std(baseline_mean, 'omitnan');
    if sigma == 0
        baseline_speed(~isnan(baseline_mean)) = 0;
    else
        baseline_speed = (baseline_mean - mu) ./ sigma;
    end
end

baseline_summary = table;
baseline_summary.operator_id = categorical(ops_all);
baseline_summary.pre_baseline_case_count = baseline_counts;
baseline_summary.baseline_n_cases_used = repmat(baseline_n_cases_global, nOps, 1);
baseline_summary.baseline_duration_mean_operator = baseline_mean;
baseline_summary.baseline_speed_operator = baseline_speed;

% Augment operator summary with baseline availability (global, before any filtering).
baseline_lookup = table;
baseline_lookup.operator_id = baseline_summary.operator_id;
baseline_lookup.baseline_speed_value = baseline_summary.baseline_speed_operator;
operator_case_summary = outerjoin(operator_case_summary, baseline_lookup, ...
    'Keys', 'operator_id', 'MergeKeys', true);
operator_case_summary.has_baseline_speed = ~isnan(operator_case_summary.baseline_speed_value);

eligible_ops_global_mask = operator_case_summary.has_baseline_speed & ...
    (operator_case_summary.n_pre_non_affera_cases >= baselineMinCasesBaseline) & ...
    (operator_case_summary.n_post_affera_cases >= minAfferaCases) & ...
    (operator_case_summary.n_affera_cases_total >= minAfferaCases) & ...
    (operator_case_summary.n_non_affera_cases_total >= baselineMinCasesBaseline);
eligible_ops_global = operator_case_summary.operator_id(eligible_ops_global_mask);

if isempty(eligible_ops_global)
    error('No operators meet the global cohort requirements: baseline cases >= %g, non-Affera cases >= %g, Affera cases >= %g. Adjust thresholds or baseline settings.', ...
        baselineMinCasesBaseline, baselineMinCasesBaseline, minAfferaCases);
end

excluded_global = operator_case_summary(~eligible_ops_global_mask, :);

if verboseOutput
    fprintf('\nOperator-level case summary:\n');
    disp(operator_case_summary);
    fprintf('Global operator cohort size (all questions): %d operators.\n', numel(eligible_ops_global));
end

% Apply procedure-type restriction to the modeling dataset after defining the cohort.
tbl_model = tbl_all;
if ~ismember('procedure_type', tbl_model.Properties.VariableNames)
    warning('procedure_type column not found; cannot filter by PVI/PVI+. Using all procedures.');
    procTypeFilter = "all";
end

if ~strcmpi(procTypeFilter, "all")
    keepTypes = categorical(procTypeFilter);
    tbl_model = tbl_model(ismember(tbl_model.procedure_type, keepTypes), :);
    if verboseOutput
        fprintf('Filtered data to %s only (%d procedures remain before operator restriction).\n', ...
            strjoin(string(procTypeFilter), ', '), height(tbl_model));
    end
end

% Restrict to the globally eligible operator cohort.
tbl = tbl_model(ismember(tbl_model.operator_id, eligible_ops_global), :);
if isempty(tbl)
    warning('No procedures remain after applying the operator cohort and procedure-type filter (%s).', ...
        strjoin(string(procTypeFilter), ', '));
end

% Attach baseline metrics to each row for downstream modeling.
tbl = outerjoin(tbl, baseline_summary, 'Keys', 'operator_id', 'MergeKeys', true);

if ~isempty(tbl)
    minDate = min(tbl.procedure_date);
    tbl.time_days = days(tbl.procedure_date - minDate);
else
    tbl.time_days = [];
end

results = struct('lme_q1', [], 'lme_q2', [], 'lme_q3', [], 'lme_q4', [], ...
    'operator_case_summary', operator_case_summary, ...
    'eligible_operator_ids', eligible_ops_global, ...
    'excluded_q1', table(), 'excluded_q2', table(), ...
    'excluded_q3', table(), 'excluded_q4', table(), ...
    'excluded_global', excluded_global, ...
    'baseline_summary', baseline_summary, ...
    'baseline_n_cases_global', baseline_n_cases_global, ...
    'baseline_non_pfa_only', nonPFABaselineOnly, ...
    'baseline_range_requested', baselineRangeRequested);

if verboseOutput
    fprintf('Prepared %d procedures across %d eligible operators for modeling.\n', ...
        height(tbl), numel(unique(tbl.operator_id)));
end

%% 2. Question 1 – Overall effect of Affera
% For Q1, Affera cases are compared specifically against pre-Affera baseline
% cases that meet the user-specified baseline criteria (pre-Affera era,
% optional PVI-only baseline, optional non-PFA-only baseline). This uses the
% same global operator cohort, but the comparison group is restricted to
% those baseline-eligible rows rather than all non-Affera cases.

results.excluded_q1 = operator_case_summary(~eligible_ops_global_mask, ...
    {'operator_id', 'n_pre_non_affera_cases', 'n_post_affera_cases'});

% Baseline (control) group: pre-Affera, non-Affera baseline-eligible cases
% from globally eligible operators, using the baseline eligibility computed
% earlier (independent of ProcedureType filter).
baseline_mask_q1 = baseline_eligible & ismember(tbl_baseline.operator_id, eligible_ops_global);
tbl_q1_baseline = tbl_baseline(baseline_mask_q1, :);

% Affera (exposed) group: post-Affera Affera cases from globally eligible
% operators, optionally restricted by procedure type (all vs PVI vs PVI+).
affera_mask_q1 = tbl_all.is_affera & tbl_all.is_post_affera_global & ...
    ismember(tbl_all.operator_id, eligible_ops_global);
if ~strcmpi(procTypeFilter, "all") && ismember('procedure_type', tbl_all.Properties.VariableNames)
    keepTypes_q1 = categorical(procTypeFilter);
    affera_mask_q1 = affera_mask_q1 & ismember(tbl_all.procedure_type, keepTypes_q1);
end
tbl_q1_affera = tbl_all(affera_mask_q1, :);

tbl_q1 = [tbl_q1_baseline; tbl_q1_affera];

if isempty(tbl_q1)
    warning('Q1: No baseline and/or Affera procedures available after applying baseline criteria, operator cohort, and procedure-type filter. Skipping model.');
else
    % Recompute time_days for the Q1 dataset only.
    minDate_q1 = min(tbl_q1.procedure_date);
    tbl_q1.time_days = days(tbl_q1.procedure_date - minDate_q1);

    if verboseOutput
        fprintf('\nQ1: fitting overall Affera effect with %d operators (%d procedures).\n', ...
            numel(unique(tbl_q1.operator_id)), height(tbl_q1));
        if ~isempty(results.excluded_q1)
            fprintf('Operators excluded from the global cohort (applies to all questions):\n');
            disp(results.excluded_q1);
        end
        fprintf('Q1 baseline group: %d procedures (pre-Affera, baseline-eligible).\n', ...
            height(tbl_q1_baseline));
        fprintf('Q1 Affera group:   %d procedures (post-Affera Affera, ProcedureType=%s).\n', ...
            height(tbl_q1_affera), strjoin(string(procTypeFilter), ', '));
        fprintf('Fitting Q1 model (Affera vs. pre-Affera baseline)...\n');
    end
    formula_q1 = 'log_duration ~ is_affera + time_days + (1|operator_id)';
    results.lme_q1 = fitlme(tbl_q1, formula_q1);
    if verboseOutput
        disp(results.lme_q1);
    end
end

%% 3. Question 2 – Affera learning curve
if verboseOutput
    fprintf('\nComputing within-operator Affera case index...\n');
end
tbl = sortrows(tbl, {'operator_id', 'procedure_date', 'procedure_id'});
tbl.affera_index = zeros(height(tbl), 1);

ops = categories(tbl.operator_id);
for k = 1:numel(ops)
    op = ops{k};
    idx = (tbl.operator_id == op) & tbl.is_affera;
    tbl.affera_index(idx) = (1:sum(idx)).';
end

results.excluded_q2 = operator_case_summary(~eligible_ops_global_mask, ...
    {'operator_id', 'n_affera_cases_total'});

tbl_affera = tbl(tbl.is_affera & tbl.affera_index > 0, :);

if isempty(tbl_affera)
    warning('Q2: No eligible Affera procedures available after filtering. Skipping model.');
else
    if verboseOutput
        fprintf('Q2: fitting learning curve with %d operators (%d Affera procedures).\n', ...
            numel(unique(tbl_affera.operator_id)), height(tbl_affera));
        if ~isempty(results.excluded_q2)
            fprintf('Operators excluded from the global cohort (applies to all questions):\n');
            disp(results.excluded_q2);
        end
        fprintf('Fitting Q2 model (Affera learning curve)...\n');
    end
    formula_q2 = 'log_duration ~ affera_index + time_days + (1 + affera_index|operator_id)';
    results.lme_q2 = fitlme(tbl_affera, formula_q2);
    if verboseOutput
        disp(results.lme_q2);
    end
end

%% 4. Baseline availability summary
baseline_valid_ops = operator_case_summary.operator_id(operator_case_summary.has_baseline_speed);
if verboseOutput
    fprintf('\nBaseline speed available for %d of %d operators.\n', ...
        numel(baseline_valid_ops), height(operator_case_summary));
end

%% 5. Question 3 – Baseline speed vs. learning curve
results.excluded_q3 = operator_case_summary(~eligible_ops_global_mask, ...
    {'operator_id', 'has_baseline_speed', 'n_affera_cases_total'});

tbl_affera_bs = tbl(tbl.is_affera & tbl.affera_index > 0 & ~isnan(tbl.baseline_speed_operator), :);

if isempty(tbl_affera_bs)
    warning('Q3: No eligible Affera procedures with baseline speed available after filtering. Skipping model.');
else
    if verboseOutput
        fprintf('\nQ3: fitting baseline speed × learning model with %d operators (%d Affera procedures).\n', ...
            numel(unique(tbl_affera_bs.operator_id)), height(tbl_affera_bs));
        if ~isempty(results.excluded_q3)
            fprintf('Operators excluded from the global cohort (applies to all questions):\n');
            disp(results.excluded_q3);
        end
        fprintf('Fitting Q3 model (baseline speed modifies learning)...\n');
    end
    formula_q3 = 'log_duration ~ affera_index*baseline_speed_operator + time_days + (1 + affera_index|operator_id)';
    results.lme_q3 = fitlme(tbl_affera_bs, formula_q3);
    if verboseOutput
        disp(results.lme_q3);
    end
end

%% 6. Question 4 – Baseline speed vs. Affera vs. non-Affera difference
results.excluded_q4 = operator_case_summary(~eligible_ops_global_mask, ...
    {'operator_id', 'has_baseline_speed', 'n_affera_cases_total', 'n_non_affera_cases_total'});

tbl_bs = tbl(~isnan(tbl.baseline_speed_operator), :);

if isempty(tbl_bs)
    warning('Q4: No operators meet the baseline/exposure requirements. Skipping model.');
else
    if verboseOutput
        fprintf('\nQ4: fitting baseline speed × Affera effect with %d operators (%d procedures).\n', ...
            numel(unique(tbl_bs.operator_id)), height(tbl_bs));
        if ~isempty(results.excluded_q4)
            fprintf('Operators excluded from the global cohort (applies to all questions):\n');
            disp(results.excluded_q4);
        end
        fprintf('Fitting Q4 model (baseline speed modifies Affera effect)...\n');
    end
    formula_q4 = 'log_duration ~ is_affera*baseline_speed_operator + time_days + (1|operator_id)';
    results.lme_q4 = fitlme(tbl_bs, formula_q4);
    if verboseOutput
        disp(results.lme_q4);
    end
end

if verboseOutput
    fprintf('\nBasic analysis complete.\n');
end

% -------------------------------------------------------------------------
% Write concise text report to results/ directory
% -------------------------------------------------------------------------
ptypeLabel = strjoin(string(procTypeFilter), '_');
ptypeLabel = strtrim(ptypeLabel);
if any(strcmpi(ptypeLabel, "PVI+"))
    ptypeLabel = "PVIplus";
end
ptypeLabel = regexprep(ptypeLabel, '[^A-Za-z0-9]', '');
if isempty(ptypeLabel)
    ptypeLabel = 'all';
end
logFile = fullfile(resultsDir, sprintf('affera_analysis_%s_%s.txt', ptypeLabel, timeStamp));

fid = fopen(logFile, 'w');
if fid == -1
    warning('Could not open log file for writing: %s', logFile);
else
    fprintf(fid, 'Affera Analysis Report\n');
    fprintf(fid, '======================\n\n');
    fprintf(fid, 'Procedure type filter       : %s\n', strjoin(string(procTypeFilter), ', '));
    fprintf(fid, 'MinAfferaCases (per-operator Affera threshold): %g\n', minAfferaCases);
    fprintf(fid, '  Q1: >=MinAfferaCases post-Affera Affera cases AND >=%g pre-Affera non-Affera cases\n', baselineMinCasesBaseline);
    fprintf(fid, '  Q2: >=MinAfferaCases Affera cases\n');
    fprintf(fid, '  Q3: baseline available AND >=MinAfferaCases Affera cases\n');
    fprintf(fid, '  Q4: baseline available AND >=MinAfferaCases Affera cases AND >=%g non-Affera cases\n', baselineMinCasesBaseline);
    fprintf(fid, 'Baseline restricted to PVI  : %s\n', ternary(baselinePVIOnly, 'YES', 'NO'));
    fprintf(fid, 'Non-PFA baseline only       : %s\n', ternary(nonPFABaselineOnly, 'YES', 'NO'));
    fprintf(fid, 'Baseline min cases/operator  : %g\n', baselineMinCasesBaseline);
    if isfinite(baselineMaxCasesBaseline)
        fprintf(fid, 'Baseline max cases/operator  : %g\n', baselineMaxCasesBaseline);
    else
        fprintf(fid, 'Baseline max cases/operator  : no cap (Inf)\n');
    end
    if metaAvailable && isfield(meta, 'affera_launch_date_global') && ~isempty(meta.affera_launch_date_global) && ~isnat(meta.affera_launch_date_global)
        fprintf(fid, 'Global Affera launch date   : %s\n', datestr(meta.affera_launch_date_global, 'yyyy-mm-dd'));
    end
    fprintf(fid, '\nBaseline summary:\n');
    fprintf(fid, '  Operators with baseline   : %d of %d\n', ...
        sum(operator_case_summary.has_baseline_speed), height(operator_case_summary));
    fprintf(fid, '  baseline_n_cases_global   : %g\n', baseline_n_cases_global);
    % Narrative description of baseline cohort
    if baseline_n_cases_global > 0
        % Describe case type used for baseline
        if baselinePVIOnly && nonPFABaselineOnly
            caseDesc = 'PVI, non-PFA (non-Affera) cases';
        elseif baselinePVIOnly && ~nonPFABaselineOnly
            caseDesc = 'PVI cases (all energy types)';
        elseif ~baselinePVIOnly && nonPFABaselineOnly
            caseDesc = 'non-PFA (non-Affera) cases';
        else
            caseDesc = 'non-Affera cases (all procedures)';
        end
        fprintf(fid, '  Baseline is computed from the %g most recent pre-Affera %s per operator,\n', ...
            baseline_n_cases_global, caseDesc);
        fprintf(fid, '  among operators with at least %g eligible baseline cases.\n', baselineMinCasesBaseline);
    else
        fprintf(fid, '  Baseline could not be computed (no operators met the baseline case requirements).\n');
    end
    fprintf(fid, '\nQuestion 1 – Overall Affera effect\n');
    if isempty(results.lme_q1)
        fprintf(fid, '  Model not fitted (insufficient eligible operators).\n');
    else
        coef = results.lme_q1.Coefficients;
        row = strcmp(coef.Name, 'is_affera_1');
        if any(row)
            est = coef.Estimate(row);
            se  = coef.SE(row);
            t   = coef.tStat(row);
            p   = coef.pValue(row);
            ci  = coefCI(results.lme_q1);
            ci_lo = ci(row,1);
            ci_hi = ci(row,2);
            fprintf(fid, '  log-duration diff (Affera vs non-Affera): %.4f (95%% CI %.4f, %.4f), p = %.3g\n', ...
                est, ci_lo, ci_hi, p);
        else
            fprintf(fid, '  Affera coefficient not found in model.\n');
        end
    end
    fprintf(fid, '\nQuestion 2 – Affera learning curve\n');
    if isempty(results.lme_q2)
        fprintf(fid, '  Model not fitted (insufficient Affera data).\n');
    else
        coef = results.lme_q2.Coefficients;
        row = strcmp(coef.Name, 'affera_index');
        if any(row)
            est = coef.Estimate(row);
            p   = coef.pValue(row);
            ci  = coefCI(results.lme_q2);
            ci_lo = ci(row,1);
            ci_hi = ci(row,2);
            fprintf(fid, '  Per-case change in log-duration (index): %.4f (95%% CI %.4f, %.4f), p = %.3g\n', ...
                est, ci_lo, ci_hi, p);
        else
            fprintf(fid, '  affera_index term not found in model.\n');
        end
    end
    fprintf(fid, '\nQuestion 3 – Baseline speed × learning\n');
    if isempty(results.lme_q3)
        fprintf(fid, '  Model not fitted (baseline and/or Affera case thresholds not met).\n');
    else
        coef = results.lme_q3.Coefficients;
        row = strcmp(coef.Name, 'affera_index:baseline_speed_operator');
        if any(row)
            est = coef.Estimate(row);
            p   = coef.pValue(row);
            ci  = coefCI(results.lme_q3);
            ci_lo = ci(row,1);
            ci_hi = ci(row,2);
            fprintf(fid, '  Interaction (index × baseline speed): %.4f (95%% CI %.4f, %.4f), p = %.3g\n', ...
                est, ci_lo, ci_hi, p);
        else
            fprintf(fid, '  index × baseline_speed term not found in model.\n');
        end
    end
    fprintf(fid, '\nQuestion 4 – Baseline speed × Affera effect\n');
    if isempty(results.lme_q4)
        fprintf(fid, '  Model not fitted (baseline and/or exposure thresholds not met).\n');
    else
        coef = results.lme_q4.Coefficients;
        row = strcmp(coef.Name, 'is_affera_1:baseline_speed_operator');
        if any(row)
            est = coef.Estimate(row);
            p   = coef.pValue(row);
            ci  = coefCI(results.lme_q4);
            ci_lo = ci(row,1);
            ci_hi = ci(row,2);
            fprintf(fid, '  Interaction (Affera × baseline speed): %.4f (95%% CI %.4f, %.4f), p = %.3g\n', ...
                est, ci_lo, ci_hi, p);
        else
            fprintf(fid, '  Affera × baseline_speed term not found in model.\n');
        end
    end

    % Detailed section: include tables and full model output.
    fprintf(fid, '\n\nDetailed results\n');
    fprintf(fid, '----------------\n\n');
    fprintf(fid, 'Operator-level case summary:\n');
    fprintf(fid, '%-25s %8s %8s %8s %8s %8s %12s\n', ...
        'operator_id', 'pre_non', 'post_aff', 'affera', 'non_aff', 'baseline', 'baseline_z');
    for i = 1:height(operator_case_summary)
        fprintf(fid, '%-25s %8d %8d %8d %8d %8d %12.4f\n', ...
            string(operator_case_summary.operator_id(i)), ...
            operator_case_summary.n_pre_non_affera_cases(i), ...
            operator_case_summary.n_post_affera_cases(i), ...
            operator_case_summary.n_affera_cases_total(i), ...
            operator_case_summary.n_non_affera_cases_total(i), ...
            operator_case_summary.has_baseline_speed(i), ...
            operator_case_summary.baseline_speed_value(i));
    end

	    fprintf(fid, '\nBaseline summary table:\n');
	    fprintf(fid, '%-25s %8s %8s %12s %12s\n', ...
	        'operator_id', 'pre_base', 'N_used', 'mean_dur', 'baseline_z');
    for i = 1:height(baseline_summary)
	        fprintf(fid, '%-25s %8d %8d %12.4f %12.4f\n', ...
	            string(baseline_summary.operator_id(i)), ...
	            baseline_summary.pre_baseline_case_count(i), ...
	            baseline_summary.baseline_n_cases_used(i), ...
	            baseline_summary.baseline_duration_mean_operator(i), ...
	            baseline_summary.baseline_speed_operator(i));
	    end
	
	    fprintf(fid, '\nGlobal cohort exclusions (applies to Q1–Q4):\n');
    if isempty(excluded_global)
        fprintf(fid, '    (none)\n');
    else
        fprintf(fid, '    %-25s %8s %8s %8s %8s %8s\n', ...
            'operator_id', 'pre_non', 'post_aff', 'affera', 'non_aff', 'baseline');
	        for i = 1:height(excluded_global)
	            fprintf(fid, '    %-25s %8d %8d %8d %8d %8d\n', ...
	                string(excluded_global.operator_id(i)), ...
	                excluded_global.n_pre_non_affera_cases(i), ...
	                excluded_global.n_post_affera_cases(i), ...
	                excluded_global.n_affera_cases_total(i), ...
	                excluded_global.n_non_affera_cases_total(i), ...
	                excluded_global.has_baseline_speed(i));
	        end
	    end

	    % -----------------------------------------------------------------
	    % Case-level tables for each question (by operator)
	    % -----------------------------------------------------------------
	    fprintf(fid, '\n\nIncluded cases by question and operator:\n');
	    fprintf(fid, '----------------------------------------\n\n');

	    % Q1 cases: tbl_q1 (baseline + Affera) if available.
	    if exist('tbl_q1', 'var') && ~isempty(tbl_q1)
	        fprintf(fid, 'Q1 included cases (baseline + Affera):\n');
	        ops_print = categories(tbl_q1.operator_id);
	        for k = 1:numel(ops_print)
	            op = ops_print{k};
	            opRows = tbl_q1(tbl_q1.operator_id == op, :);
	            fprintf(fid, '\n  Operator: %s (n = %d)\n', string(op), height(opRows));
	            fprintf(fid, '    %-15s %-12s %-10s %-8s %-8s %-8s %10s\n', ...
	                'procedure_id', 'date', 'type', 'Affera', 'preEra', 'postEra', 'duration');
	            for i = 1:height(opRows)
	                fprintf(fid, '    %-15s %-12s %-10s %-8d %-8d %-8d %10.2f\n', ...
	                    string(opRows.procedure_id(i)), ...
	                    datestr(opRows.procedure_date(i), 'yyyy-mm-dd'), ...
	                    char(opRows.procedure_type(i)), ...
	                    opRows.is_affera(i), ...
	                    opRows.is_pre_affera_global(i), ...
	                    opRows.is_post_affera_global(i), ...
	                    opRows.procedure_duration_min(i));
	            end
	        end
	        fprintf(fid, '\n');
	    end

	    % Q2 cases: tbl_affera (Affera-only learning-curve dataset).
	    if exist('tbl_affera', 'var') && ~isempty(tbl_affera)
	        fprintf(fid, 'Q2 included Affera cases (learning curve):\n');
	        ops_print = categories(tbl_affera.operator_id);
	        for k = 1:numel(ops_print)
	            op = ops_print{k};
	            opRows = tbl_affera(tbl_affera.operator_id == op, :);
	            fprintf(fid, '\n  Operator: %s (n = %d)\n', string(op), height(opRows));
	            fprintf(fid, '    %-15s %-12s %-10s %-10s %-8s %10s\n', ...
	                'procedure_id', 'date', 'type', 'affera_idx', 'postEra', 'duration');
	            for i = 1:height(opRows)
	                fprintf(fid, '    %-15s %-12s %-10s %10.0f %-8d %10.2f\n', ...
	                    string(opRows.procedure_id(i)), ...
	                    datestr(opRows.procedure_date(i), 'yyyy-mm-dd'), ...
	                    char(opRows.procedure_type(i)), ...
	                    opRows.affera_index(i), ...
	                    opRows.is_post_affera_global(i), ...
	                    opRows.procedure_duration_min(i));
	            end
	        end
	        fprintf(fid, '\n');
	    end

	    % Q3 cases: tbl_affera_bs (Affera with baseline_speed for learning).
	    if exist('tbl_affera_bs', 'var') && ~isempty(tbl_affera_bs)
	        fprintf(fid, 'Q3 included Affera cases (learning × baseline speed):\n');
	        ops_print = categories(tbl_affera_bs.operator_id);
	        for k = 1:numel(ops_print)
	            op = ops_print{k};
	            opRows = tbl_affera_bs(tbl_affera_bs.operator_id == op, :);
	            fprintf(fid, '\n  Operator: %s (n = %d)\n', string(op), height(opRows));
	            fprintf(fid, '    %-15s %-12s %-10s %-10s %12s %10s\n', ...
	                'procedure_id', 'date', 'type', 'affera_idx', 'baseline_z', 'duration');
	            for i = 1:height(opRows)
	                fprintf(fid, '    %-15s %-12s %-10s %10.0f %12.4f %10.2f\n', ...
	                    string(opRows.procedure_id(i)), ...
	                    datestr(opRows.procedure_date(i), 'yyyy-mm-dd'), ...
	                    char(opRows.procedure_type(i)), ...
	                    opRows.affera_index(i), ...
	                    opRows.baseline_speed_operator(i), ...
	                    opRows.procedure_duration_min(i));
	            end
	        end
	        fprintf(fid, '\n');
	    end

	    % Q4 cases: tbl_bs (all cases with baseline_speed used in Affera vs non-Affera × baseline).
	    if exist('tbl_bs', 'var') && ~isempty(tbl_bs)
	        fprintf(fid, 'Q4 included cases (Affera vs non-Affera × baseline speed):\n');
	        ops_print = categories(tbl_bs.operator_id);
	        for k = 1:numel(ops_print)
	            op = ops_print{k};
	            opRows = tbl_bs(tbl_bs.operator_id == op, :);
	            fprintf(fid, '\n  Operator: %s (n = %d)\n', string(op), height(opRows));
	            fprintf(fid, '    %-15s %-12s %-10s %-8s %12s %10s\n', ...
	                'procedure_id', 'date', 'type', 'Affera', 'baseline_z', 'duration');
	            for i = 1:height(opRows)
	                fprintf(fid, '    %-15s %-12s %-10s %-8d %12.4f %10.2f\n', ...
	                    string(opRows.procedure_id(i)), ...
	                    datestr(opRows.procedure_date(i), 'yyyy-mm-dd'), ...
	                    char(opRows.procedure_type(i)), ...
	                    opRows.is_affera(i), ...
	                    opRows.baseline_speed_operator(i), ...
	                    opRows.procedure_duration_min(i));
	            end
	        end
	        fprintf(fid, '\n');
	    end

	    if ~isempty(results.lme_q1)
	        fprintf(fid, '\nFull Q1 model output:\n');
        q1ModelText = evalc('disp(results.lme_q1);');
        q1ModelText = regexprep(q1ModelText, '<[^>]+>', '');
        fprintf(fid, '%s\n', q1ModelText);
    end
    if ~isempty(results.lme_q2)
        fprintf(fid, '\nFull Q2 model output:\n');
        q2ModelText = evalc('disp(results.lme_q2);');
        q2ModelText = regexprep(q2ModelText, '<[^>]+>', '');
        fprintf(fid, '%s\n', q2ModelText);
    end
    if ~isempty(results.lme_q3)
        fprintf(fid, '\nFull Q3 model output:\n');
        q3ModelText = evalc('disp(results.lme_q3);');
        q3ModelText = regexprep(q3ModelText, '<[^>]+>', '');
        fprintf(fid, '%s\n', q3ModelText);
    end
    if ~isempty(results.lme_q4)
        fprintf(fid, '\nFull Q4 model output:\n');
        q4ModelText = evalc('disp(results.lme_q4);');
        q4ModelText = regexprep(q4ModelText, '<[^>]+>', '');
        fprintf(fid, '%s\n', q4ModelText);
    end

    fclose(fid);
    fprintf('Analysis report written to %s\n', logFile);
end
end
