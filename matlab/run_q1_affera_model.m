function [lme, results, tbl_q1] = run_q1_affera_model(clinicalData, baseline_days, min_cases_per_group)
% RUN_Q1_AFFERA_MODEL  Q1 Affera vs non-Affera efficiency analysis.
%
%   [lme, results, tbl_q1] = run_q1_affera_model(clinicalData) takes the
%   procedure-level table produced by import_clinical_data (or loaded from
%   clinicalData/importedClinicalData.mat), constructs the Q1 analysis
%   dataset, fits the mixed-effects model, prints a summary, and returns
%   the fitted model and results.
%
%   [...] = run_q1_affera_model(clinicalData, baseline_days, min_cases_per_group)
%   allows overriding:
%     - baseline_days: length of the pre-Affera baseline window (default 180)
%     - min_cases_per_group: minimum number of baseline and Affera cases
%       required per operator (default 15 for each group).

if nargin < 2 || isempty(baseline_days)
    baseline_days = 180;
end

if nargin < 3 || isempty(min_cases_per_group)
    min_cases_per_group = 15;
end

tbl = clinicalData;

if ~ismember('procedure_date', tbl.Properties.VariableNames) || ...
   ~ismember('is_affera', tbl.Properties.VariableNames)
    error('run_q1_affera_model:MissingVariables', ...
        'Input table must contain procedure_date and is_affera variables.');
end

meta = struct();
afferaIdx = tbl.is_affera;
if ~any(afferaIdx)
    error('run_q1_affera_model:NoAfferaCases', ...
        'Input table contains no Affera cases (is_affera == true).');
end
meta.affera_launch_date_global = min(tbl.procedure_date(afferaIdx));

[tbl_q1, operator_counts_all] = make_q1_affera_dataset(tbl, meta, baseline_days, min_cases_per_group);
[lme, results] = fit_q1_affera_model(tbl_q1);

% Attach operator summary (all operators with included flag) and post-Affera
% non-PFA case count.
results.operator_counts = operator_counts_all;
if ismember('isPostNonPFAEra', tbl_q1.Properties.VariableNames)
    results.nPostNonPFA = sum(tbl_q1.isPostNonPFAEra);
end

print_q1_affera_summary(lme, results, tbl_q1);

% Save outputs for downstream use.
outDir = 'results';
if ~isfolder(outDir)
    mkdir(outDir);
end

save(fullfile(outDir, 'q1_affera_model_results.mat'), 'tbl_q1', 'lme', 'results');

end
