function [lme, results, tbl_q1] = run_q1_affera_model(clinicalData, varargin)
% RUN_Q1_AFFERA_MODEL  Q1 Affera vs non-Affera efficiency analysis.
%
%   [lme, results, tbl_q1] = run_q1_affera_model(clinicalData) takes the
%   procedure-level table produced by import_clinical_data (or loaded from
%   clinicalData/importedClinicalData.mat), constructs the Q1 analysis
%   dataset, fits the mixed-effects model, prints a summary, and returns
%   the fitted model and results.
%
%   [...] = run_q1_affera_model(clinicalData, 'BaselineDays', X, 'MinCasesPerGroup', Y, 'EnableLearningCurve', TF)
%   allows overriding:
%     - BaselineDays: length of the pre-Affera baseline window (default 180)
%     - MinCasesPerGroup: minimum number of baseline and Affera cases per operator (default 15)
%     - EnableLearningCurve: logical flag to include the Affera case index term (default false)
%
%   For backward compatibility, the legacy positional arguments
%   (baseline_days, min_cases_per_group, opts struct) are still accepted.

% Defaults.
baseline_days = 180;
min_cases_per_group = 15;
opts = struct('enableLearningCurve', false);

% Backward-compatible positional usage: (clinicalData, baseline_days, min_cases_per_group, opts)
if ~isempty(varargin) && isnumeric(varargin{1})
    baseline_days = varargin{1};
    varargin(1) = [];
    if ~isempty(varargin) && isnumeric(varargin{1})
        min_cases_per_group = varargin{1};
        varargin(1) = [];
    end
    if ~isempty(varargin) && isstruct(varargin{1})
        tmpOpts = varargin{1};
        varargin(1) = [];
        if isfield(tmpOpts, 'enableLearningCurve')
            opts.enableLearningCurve = logical(tmpOpts.enableLearningCurve);
        end
    end
end

% Name-value parsing.
p = inputParser;
p.addParameter('BaselineDays', baseline_days, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinCasesPerGroup', min_cases_per_group, @(x) isnumeric(x) && isscalar(x));
p.addParameter('EnableLearningCurve', opts.enableLearningCurve, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

baseline_days = p.Results.BaselineDays;
min_cases_per_group = p.Results.MinCasesPerGroup;
opts.enableLearningCurve = logical(p.Results.EnableLearningCurve);
opts.baseline_days = baseline_days;
opts.min_cases_per_group = min_cases_per_group;

tbl = clinicalData;

if ~ismember('procedure_date', tbl.Properties.VariableNames) || ...
   ~ismember('is_affera', tbl.Properties.VariableNames)
    error('run_q1_affera_model:MissingVariables', ...
        'Input table must contain procedure_date and is_affera variables.');
end

meta = struct();
afferaIdx = tbl.is_affera ~= 0;
if ~any(afferaIdx)
    error('run_q1_affera_model:NoAfferaCases', ...
        'Input table contains no Affera cases (is_affera == true).');
end
meta.affera_launch_date_global = min(tbl.procedure_date(afferaIdx));

[tbl_q1, operator_counts_all] = make_q1_affera_dataset(tbl, meta, baseline_days, min_cases_per_group);
diagnostics = [];

if opts.enableLearningCurve
    tbl_q1 = add_affera_learning_curve_variables(tbl_q1);
    diagnostics = compute_learning_curve_diagnostics(tbl_q1);
end

[lme, results] = fit_q1_affera_model(tbl_q1, opts);

% Attach operator summary (all operators with included flag) and post-Affera
% non-PFA case count.
results.operator_counts = operator_counts_all;
if ismember('isPostNonPFAEra', tbl_q1.Properties.VariableNames)
    results.nPostNonPFA = sum(tbl_q1.isPostNonPFAEra);
end
if ~isempty(diagnostics)
    results.diagnostics = diagnostics;
end

print_q1_affera_summary(lme, results, tbl_q1, opts, diagnostics);

% Save outputs for downstream use.
outDir = 'results';
if ~isfolder(outDir)
    mkdir(outDir);
end

save(fullfile(outDir, 'q1_affera_model_results.mat'), 'tbl_q1', 'lme', 'results');

end
