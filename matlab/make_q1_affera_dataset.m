function [tbl_q1, operator_counts_all] = make_q1_affera_dataset(tbl, meta, baseline_days, min_cases_per_group)
% MAKE_Q1_AFFERA_DATASET  Build analysis-ready table for Q1 Affera model.
%
%   tbl_q1 = make_q1_affera_dataset(tbl, meta) uses the imported procedure-
%   level table and metadata from import_clinical_data to construct the
%   subset and variables required by the Q1 mixed-effects model
%   specification (Affera vs non-Affera efficiency analysis).
%
%   tbl_q1 = make_q1_affera_dataset(tbl, meta, baseline_days) allows
%   overriding the default baseline window length (in days) before the
%   first Affera case (default: 180).

if nargin < 3 || isempty(baseline_days)
    baseline_days = 180;
end

if nargin < 4 || isempty(min_cases_per_group)
    min_cases_per_group = 15;
end

if ~isfield(meta, 'affera_launch_date_global') || isempty(meta.affera_launch_date_global)
    error('make_q1_affera_dataset:MissingAfferaLaunchDate', ...
        'meta.affera_launch_date_global is missing or empty.');
end

Affera_start = meta.affera_launch_date_global;

baseline_start = Affera_start - days(baseline_days);
baseline_end   = Affera_start;  % exclusive upper bound

% Logical flags from imported table.
if ~ismember('is_affera', tbl.Properties.VariableNames)
    error('make_q1_affera_dataset:MissingVariable', ...
        'Table is missing required variable ''is_affera''.');
end
if ~ismember('is_pfa_catheter', tbl.Properties.VariableNames)
    error('make_q1_affera_dataset:MissingVariable', ...
        'Table is missing required variable ''is_pfa_catheter''.');
end

isAffera_raw = tbl.is_affera;

isBaseline = ~tbl.is_pfa_catheter & ...
             ~isAffera_raw & ...
             tbl.procedure_date >= baseline_start & ...
             tbl.procedure_date <  baseline_end;

% v2: non-PFA post-Affera cases (used to improve identifiability of time trend).
isPostNonPFA = ~tbl.is_pfa_catheter & ...
               ~isAffera_raw & ...
               tbl.procedure_date >= Affera_start;

include = isAffera_raw | isBaseline | isPostNonPFA;

% Work first on the full baseline+post-Affera-non-PFA+Affera set for
% operator summaries and mask construction.
tbl_inc           = tbl(include, :);
baselineMaskInc   = isBaseline(include);
postNonPFAMaskInc = isPostNonPFA(include);
afferaMaskInc     = isAffera_raw(include);

% Duration in minutes (spec name).
if ~ismember('procedure_duration_min', tbl_inc.Properties.VariableNames)
    error('make_q1_affera_dataset:MissingVariable', ...
        'Table is missing required variable ''procedure_duration_min''.');
end
tbl_inc.duration_minutes = tbl_inc.procedure_duration_min;

% Affera indicator (analysis variable name).
tbl_inc.isAffera = double(tbl_inc.is_affera);

% PVI+ indicator derived from procedure_type (PVI vs PVI+).
if ~ismember('procedure_type', tbl_inc.Properties.VariableNames)
    error('make_q1_affera_dataset:MissingVariable', ...
        'Table is missing required variable ''procedure_type''.');
end
tbl_inc.isPVIplus = tbl_inc.procedure_type == 'PVI+';

% Time since earliest procedure date in the analysis cohort (using all
% baseline+Affera rows before operator filtering).
if ~ismember('procedure_date', tbl_inc.Properties.VariableNames)
    error('make_q1_affera_dataset:MissingVariable', ...
        'Table is missing required variable ''procedure_date''.');
end
tbl_inc.time_days = days(tbl_inc.procedure_date - min(tbl_inc.procedure_date));

% Ensure log_duration is present; reuse import value.
if ~ismember('log_duration', tbl_inc.Properties.VariableNames)
    if ismember('duration_minutes', tbl_inc.Properties.VariableNames)
        tbl_inc.log_duration = log(tbl_inc.duration_minutes);
    else
        error('make_q1_affera_dataset:MissingVariable', ...
            'Table is missing variables needed to construct log_duration.');
    end
end

% Optional baseline indicators for reporting (on full included set).
tbl_inc.isBaselineEra   = baselineMaskInc;
tbl_inc.isPostNonPFAEra = postNonPFAMaskInc;

% Operator-level inclusion based on minimum baseline and Affera cases, and
% construction of an operator summary table that includes all operators in
% the included set (baseline window + post-Affera non-PFA + Affera).
ops = categories(tbl_inc.operator_id);
nOps = numel(ops);
keepOp = false(nOps, 1);

operator_id_all = categorical(ops);
nBaseline_all   = zeros(nOps, 1);
nAffera_all     = zeros(nOps, 1);

for k = 1:nOps
    op = ops{k};
    idx = (tbl_inc.operator_id == op);
    nBaseline_all(k) = sum(idx & tbl_inc.isBaselineEra);
    nAffera_all(k)   = sum(idx & (tbl_inc.isAffera ~= 0));
    keepOp(k) = (nBaseline_all(k) >= min_cases_per_group) && ...
                (nAffera_all(k)   >= min_cases_per_group);
end

includedFlag = keepOp;

operator_counts_all = table(operator_id_all, nBaseline_all, nAffera_all, includedFlag, ...
    'VariableNames', {'operator_id', 'n_baseline', 'n_affera', 'included_in_model'});

% Apply operator filtering to build the final analysis dataset.
keep_rows = ismember(tbl_inc.operator_id, ops(keepOp));
tbl_q1 = tbl_inc(keep_rows, :);

% Drop unused operator_id categories so downstream summaries only reflect
% operators actually included in the analysis dataset.
tbl_q1.operator_id = removecats(tbl_q1.operator_id);

end
