function [tbl_cohort, operator_summary] = build_affera_operator_cohort(tbl, varargin)
% BUILD_AFFERA_OPERATOR_COHORT  Define an operator cohort based on Affera/baseline counts.
%
%   [TBL_COHORT, OP_SUMMARY] = BUILD_AFFERA_OPERATOR_COHORT() loads
%   data_processed/procedures_all.mat (variable 'tbl'), applies user-specified
%   thresholds on Affera and baseline cases, and returns:
%       - TBL_COHORT: subset of tbl containing only operators in the cohort.
%       - OP_SUMMARY: per-operator table with counts and a logical 'in_cohort'.
%
%   [TBL_COHORT, OP_SUMMARY] = BUILD_AFFERA_OPERATOR_COHORT(TBL, ...) uses the
%   provided procedure-level table instead of loading from disk.
%
%   Criteria are defined separately for PVI and PVI+:
%
%     Affera exposure counts (per operator):
%       - MinAfferaPVI     : minimum number of Affera PVI cases.
%       - MinAfferaPVIplus : minimum number of Affera PVI+ cases.
%       (Affera cases are defined as is_affera == true; by construction these
%        are post-Affera-era procedures.)
%
%     Baseline counts (per operator):
%       - MinBaselinePVI      : minimum number of baseline PVI cases.
%       - MinBaselinePVIplus  : minimum number of baseline PVI+ cases.
%       Baseline is defined as:
%           is_pre_affera_global == true
%           & ~is_pfa_catheter           (non-PFA)
%           & ~is_affera                 (non-Affera)
%
%   Name/value arguments:
%       'MinAfferaPVI'        (default 0)
%       'MinAfferaPVIplus'    (default 0)
%       'MinBaselinePVI'      (default 0)
%       'MinBaselinePVIplus'  (default 0)
%       'Verbose'             (default true)

parser = inputParser;
parser.addParameter('MinAfferaPVI', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('MinAfferaPVIplus', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('MinBaselinePVI', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('MinBaselinePVIplus', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parser.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
parser.parse(varargin{:});

minAfferaPVI = parser.Results.MinAfferaPVI;
minAfferaPVIplus = parser.Results.MinAfferaPVIplus;
minBaselinePVI = parser.Results.MinBaselinePVI;
minBaselinePVIplus = parser.Results.MinBaselinePVIplus;
verboseOutput = parser.Results.Verbose;

% Load default table if none provided.
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
    if verboseOutput
        fprintf('Loaded %d procedures from %s\n', height(tbl), defaultMat);
    end
else
    if verboseOutput
        fprintf('Building operator cohort from provided table with %d procedures.\n', height(tbl));
    end
end

%% Basic validation and housekeeping
requiredVars = {'procedure_type', 'is_affera', 'operator_id', ...
    'procedure_date', 'is_pre_affera_global', 'is_pfa_catheter'};
missingVars = setdiff(requiredVars, tbl.Properties.VariableNames);
if ~isempty(missingVars)
    error('Missing required variables in table: %s', strjoin(missingVars, ', '));
end

if ~iscategorical(tbl.operator_id)
    tbl.operator_id = categorical(tbl.operator_id);
end
if ~iscategorical(tbl.procedure_type)
    tbl.procedure_type = categorical(tbl.procedure_type);
end

% Ensure boolean/logical masks.
tbl.is_affera = logical(tbl.is_affera);
tbl.is_pre_affera_global = logical(tbl.is_pre_affera_global);
tbl.is_pfa_catheter = logical(tbl.is_pfa_catheter);

%% Define masks for procedure types
isPVI = tbl.procedure_type == categorical("PVI");
isPVIplus = tbl.procedure_type == categorical("PVI+");

%% Affera counts (PVI and PVI+)
afferaPVI_mask = tbl.is_affera & isPVI;
afferaPVIplus_mask = tbl.is_affera & isPVIplus;

%% Baseline counts: pre-Affera, non-PFA, non-Affera (PVI and PVI+)
baseline_mask = tbl.is_pre_affera_global & ~tbl.is_pfa_catheter & ~tbl.is_affera;
baselinePVI_mask = baseline_mask & isPVI;
baselinePVIplus_mask = baseline_mask & isPVIplus;

%% Per-operator counts
ops = categories(tbl.operator_id);
nOps = numel(ops);

operator_summary = table('Size', [nOps 6], ...
    'VariableTypes', {'categorical', 'double', 'double', 'double', 'double', 'logical'}, ...
    'VariableNames', {'operator_id', 'n_affera_PVI', 'n_affera_PVIplus', ...
    'n_baseline_PVI', 'n_baseline_PVIplus', 'in_cohort'});
operator_summary.operator_id = categorical(ops);

for k = 1:nOps
    op = ops{k};
    idxOp = (tbl.operator_id == op);

    operator_summary.n_affera_PVI(k) = sum(idxOp & afferaPVI_mask);
    operator_summary.n_affera_PVIplus(k) = sum(idxOp & afferaPVIplus_mask);
    operator_summary.n_baseline_PVI(k) = sum(idxOp & baselinePVI_mask);
    operator_summary.n_baseline_PVIplus(k) = sum(idxOp & baselinePVIplus_mask);
end

% Apply cohort criteria.
operator_summary.in_cohort = ...
    (operator_summary.n_affera_PVI >= minAfferaPVI) & ...
    (operator_summary.n_affera_PVIplus >= minAfferaPVIplus) & ...
    (operator_summary.n_baseline_PVI >= minBaselinePVI) & ...
    (operator_summary.n_baseline_PVIplus >= minBaselinePVIplus);

cohort_ops = operator_summary.operator_id(operator_summary.in_cohort);

%% Build cohort table
tbl_cohort = tbl(ismember(tbl.operator_id, cohort_ops), :);

if verboseOutput
    fprintf('\nCohort criteria:\n');
    fprintf('  MinAfferaPVI        : %g\n', minAfferaPVI);
    fprintf('  MinAfferaPVIplus    : %g\n', minAfferaPVIplus);
    fprintf('  MinBaselinePVI      : %g\n', minBaselinePVI);
    fprintf('  MinBaselinePVIplus  : %g\n', minBaselinePVIplus);
    fprintf('Operators meeting criteria: %d of %d\n', sum(operator_summary.in_cohort), nOps);
    fprintf('Cohort procedures retained: %d of %d\n', height(tbl_cohort), height(tbl));

    % Print per-operator table of counts.
    fprintf('\nOperator-level Affera and baseline counts (PVI / PVI+):\n');
    fprintf('%-20s %10s %12s %12s %14s %10s\n', ...
        'operator_id', 'base_PVI', 'base_PVIplus', 'affera_PVI', 'affera_PVIplus', 'in_cohort');
    for k = 1:nOps
        fprintf('%-20s %10d %12d %12d %14d %10d\n', ...
            string(operator_summary.operator_id(k)), ...
            operator_summary.n_baseline_PVI(k), ...
            operator_summary.n_baseline_PVIplus(k), ...
            operator_summary.n_affera_PVI(k), ...
            operator_summary.n_affera_PVIplus(k), ...
            operator_summary.in_cohort(k));
    end
end

end
