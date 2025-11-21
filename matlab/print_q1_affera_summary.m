function print_q1_affera_summary(lme, results, tbl_q1, opts, diagnostics)
% PRINT_Q1_AFFERA_SUMMARY  Display a human-readable summary of Q1 results.
%
%   print_q1_affera_summary(lme, results, tbl_q1) prints key sample counts
%   and interpretable fixed-effect estimates, including Affera effects in
%   PVI-only and PVI+ cases.
%
%   [...] = print_q1_affera_summary(lme, results, tbl_q1, opts, diagnostics)
%   also reports learning-curve and collinearity diagnostics when provided.

if nargin < 4 || isempty(opts)
    opts = struct();
end
if nargin < 5
    diagnostics = [];
end

fprintf('\n=== Q1 Affera vs Non-Affera Mixed-Effects Model ===\n');

% Report settings if available.
if isfield(opts, 'baseline_days') || isfield(opts, 'min_cases_per_group') || isfield(opts, 'enableLearningCurve')
    fprintf('Settings:\n');
    if isfield(opts, 'baseline_days')
        fprintf('  Baseline window (days): %d\n', opts.baseline_days);
    end
    if isfield(opts, 'min_cases_per_group')
        fprintf('  Minimum cases per group (baseline and Affera): %d\n', opts.min_cases_per_group);
    end
    if isfield(opts, 'enableLearningCurve')
        fprintf('  Learning curve enabled: %s\n', string(opts.enableLearningCurve));
    end
end

% Sample and operator counts.
nTotal        = height(tbl_q1);
nBaseline     = sum(tbl_q1.isBaselineEra);
nAffera       = sum(tbl_q1.isAffera ~= 0);
hasPostNonPFA = ismember('isPostNonPFAEra', tbl_q1.Properties.VariableNames);
if hasPostNonPFA
    nPostNonPFA = sum(tbl_q1.isPostNonPFAEra);
else
    nPostNonPFA = NaN;
end

isPVIplus = logical(tbl_q1.isPVIplus);
nPVI      = sum(~isPVIplus);
nPVIplus  = sum(isPVIplus);

nOps = numel(categories(tbl_q1.operator_id));

fprintf('Total procedures included: %d\n', nTotal);
fprintf('  Baseline (non-PFA, pre-Affera): %d\n', nBaseline);
fprintf('  Affera cases: %d\n', nAffera);
if hasPostNonPFA
    fprintf('  Post-Affera non-PFA cases: %d\n', nPostNonPFA);
end
fprintf('PVI-only cases: %d\n', nPVI);
fprintf('PVI+ cases: %d\n', nPVIplus);
fprintf('Number of operators: %d\n', nOps);

% Operator-level case counts (all operators with inclusion flag).
if isfield(results, 'operator_counts') && ~isempty(results.operator_counts)
    fprintf('\nOperator case counts (baseline window + Affera), with inclusion flag:\n');
    opTbl = results.operator_counts;
    hasIncludedFlag = ismember('included_in_model', opTbl.Properties.VariableNames);
    if hasIncludedFlag
        fprintf('%-20s %12s %12s %10s\n', 'Operator', 'Baseline n', 'Affera n', 'Included');
    else
        fprintf('%-20s %12s %12s\n', 'Operator', 'Baseline n', 'Affera n');
    end
    for i = 1:height(opTbl)
        if hasIncludedFlag
            incStr = 'no';
            if opTbl.included_in_model(i)
                incStr = 'yes';
            end
            fprintf('%-20s %12d %12d %10s\n', ...
                char(opTbl.operator_id(i)), opTbl.n_baseline(i), opTbl.n_affera(i), incStr);
        else
            fprintf('%-20s %12d %12d\n', ...
                char(opTbl.operator_id(i)), opTbl.n_baseline(i), opTbl.n_affera(i));
        end
    end
end

% Descriptive duration summaries.
bs  = results.baseline_stats;
as  = results.affera_stats;
hasPostStats = isfield(results, 'post_non_pfa_stats');
if hasPostStats
    pnp = results.post_non_pfa_stats;
else
    pnp = [];
end

% Case-mix counts by era and complexity.
fprintf('\nCase mix by era and complexity (procedure counts):\n');
fprintf('  Baseline (non-PFA, pre-Affera):\n');
fprintf('    Total: %d  |  PVI-only: %d  |  PVI+: %d\n', ...
    bs.overall.n, bs.pvi.n, bs.pvi_plus.n);

fprintf('  Affera era (Affera procedures only):\n');
fprintf('    Total: %d  |  PVI-only: %d  |  PVI+: %d\n', ...
    as.overall.n, as.pvi.n, as.pvi_plus.n);

if hasPostStats && pnp.overall.n > 0
    fprintf('  Post-Affera non-PFA:\n');
    fprintf('    Total: %d  |  PVI-only: %d  |  PVI+: %d\n', ...
        pnp.overall.n, pnp.pvi.n, pnp.pvi_plus.n);
end

% Operator share of baseline and post-Affera non-PFA procedures.
fprintf('\nOperator share of baseline and post-Affera non-PFA procedures:\n');
opsInModel = categories(tbl_q1.operator_id);
for k = 1:numel(opsInModel)
    op = opsInModel{k};
    maskOp = (tbl_q1.operator_id == op);
    nBaseOp = sum(maskOp & tbl_q1.isBaselineEra);
    if hasPostNonPFA
        nPostNonPFAOp = sum(maskOp & tbl_q1.isPostNonPFAEra);
    else
        nPostNonPFAOp = 0;
    end
    pctBase = 100 * nBaseOp / max(nBaseline, 1);
    pctPost = 100 * nPostNonPFAOp / max(nPostNonPFA, 1 + (nPostNonPFA == 0));
    fprintf('  %-20s Baseline: %3d (%.1f%% of baseline)', op, nBaseOp, pctBase);
    if hasPostNonPFA && nPostNonPFA > 0
        fprintf('  |  Post-Affera non-PFA: %3d (%.1f%% of post-Affera non-PFA)\n', ...
            nPostNonPFAOp, pctPost);
    else
        fprintf('  |  Post-Affera non-PFA:    0 (n/a)\n');
    end
end

fprintf('\nBaseline duration (non-PFA, pre-Affera):\n');
fprintf('  All baseline procedures: n = %d, mean = %.1f min, median = %.1f min\n', ...
    bs.overall.n, bs.overall.mean_duration, bs.overall.median_duration);
fprintf('  Baseline PVI-only:       n = %d, mean = %.1f min, median = %.1f min\n', ...
    bs.pvi.n, bs.pvi.mean_duration, bs.pvi.median_duration);
fprintf('  Baseline PVI+:           n = %d, mean = %.1f min, median = %.1f min\n', ...
    bs.pvi_plus.n, bs.pvi_plus.mean_duration, bs.pvi_plus.median_duration);

if ~isempty(bs.by_catheter) && height(bs.by_catheter) > 0
    fprintf('\nBaseline duration by catheter_primary (all PVI statuses):\n');
    fprintf('%-30s %6s %12s %12s\n', 'Catheter', 'n', 'Mean (min)', 'Median (min)');
    for i = 1:height(bs.by_catheter)
        row = bs.by_catheter(i, :);
        fprintf('%-30s %6d %12.1f %12.1f\n', ...
            char(row.catheter_primary), row.n, row.mean_duration, row.median_duration);
    end
end

if hasPostStats && pnp.overall.n > 0
    fprintf('\nPost-Affera non-PFA duration:\n');
    fprintf('  All post-Affera non-PFA procedures: n = %d, mean = %.1f min, median = %.1f min\n', ...
        pnp.overall.n, pnp.overall.mean_duration, pnp.overall.median_duration);
    fprintf('  Post-Affera non-PFA PVI-only:       n = %d, mean = %.1f min, median = %.1f min\n', ...
        pnp.pvi.n, pnp.pvi.mean_duration, pnp.pvi.median_duration);
    fprintf('  Post-Affera non-PFA PVI+:           n = %d, mean = %.1f min, median = %.1f min\n', ...
        pnp.pvi_plus.n, pnp.pvi_plus.mean_duration, pnp.pvi_plus.median_duration);

    if ~isempty(pnp.by_catheter) && height(pnp.by_catheter) > 0
        fprintf('\nPost-Affera non-PFA duration by catheter_primary:\n');
        fprintf('%-30s %6s %12s %12s\n', 'Catheter', 'n', 'Mean (min)', 'Median (min)');
        for i = 1:height(pnp.by_catheter)
            row = pnp.by_catheter(i, :);
            fprintf('%-30s %6d %12.1f %12.1f\n', ...
                char(row.catheter_primary), row.n, row.mean_duration, row.median_duration);
        end
    end
end

fprintf('\nAffera duration:\n');
fprintf('  All Affera procedures: n = %d, mean = %.1f min, median = %.1f min\n', ...
    as.overall.n, as.overall.mean_duration, as.overall.median_duration);
fprintf('  Affera PVI-only:       n = %d, mean = %.1f min, median = %.1f min\n', ...
    as.pvi.n, as.pvi.mean_duration, as.pvi.median_duration);
fprintf('  Affera PVI+:           n = %d, mean = %.1f min, median = %.1f min\n', ...
    as.pvi_plus.n, as.pvi_plus.mean_duration, as.pvi_plus.median_duration);

% Fixed-effect summary table (percent scale).
fprintf('\nFixed effects (percent change in duration):\n');
fprintf('%-25s %26s %26s %12s\n', 'Term', 'Percent change', '95% CI (percent)', 'p-value');

for i = 1:numel(results.names)
    name = results.names{i};
    pe   = results.pct_est(i);
    peLo = results.pct_lo(i);
    peHi = results.pct_hi(i);
    pVal = NaN;
    if isfield(results, 'pValue') && numel(results.pValue) >= i
        pVal = results.pValue(i);
    end

    fprintf('%-25s %10.1f%% %15s %6.1f%%, %6.1f%%] %12.3g\n', ...
        name, pe, '[', peLo, peHi, pVal);
end

% Highlight key Affera effects.
if ~isempty(results.idxAffera)
    i = results.idxAffera;
    fprintf('\nAffera effect in PVI-only cases (isAffera):\n');
    fprintf('  Percent change in duration = %.1f%% [%.1f, %.1f]%%\n', ...
        results.pct_est(i), results.pct_lo(i), results.pct_hi(i));
    if isfield(results, 'pValue') && numel(results.pValue) >= i
        fprintf('  p-value = %.3g\n', results.pValue(i));
    end
end

if ~isempty(results.idxAfferaPVIplus)
    i = results.idxAfferaPVIplus;
    fprintf('\nAdditional Affera effect in PVI+ vs PVI (isAffera:isPVIplus):\n');
    fprintf('  Percent change (additional) = %.1f%% [%.1f, %.1f]%%\n', ...
        results.pct_est(i), results.pct_lo(i), results.pct_hi(i));
    if isfield(results, 'pValue') && numel(results.pValue) >= i
        fprintf('  p-value = %.3g\n', results.pValue(i));
    end
end

% Overall Affera effect across PVI/PVI+ mix (linear contrast).
if isfield(results, 'overallAffera') && ~isempty(results.overallAffera)
    oa = results.overallAffera;
    fprintf('\nOverall Affera effect across PVI/PVI+ case mix:\n');
    if isfield(oa, 'pPVIplus')
        fprintf('  Weighted by mean isPVIplus = %.3f\n', oa.pPVIplus);
    end
    fprintf('  Percent change in duration = %.1f%% [%.1f, %.1f]%%\n', ...
        oa.pct_est, oa.pct_lo, oa.pct_hi);
    if isfield(oa, 'pValue') && ~isnan(oa.pValue)
        fprintf('  p-value = %.3g\n', oa.pValue);
    end
end

% Optional learning-curve term.
if isfield(results, 'idxAfferaIndex') && ~isempty(results.idxAfferaIndex)
    i = results.idxAfferaIndex;
    fprintf('\nLearning curve (Affera case index, centered):\n');
    fprintf('  Percent change per additional Affera case = %.2f%% [%.2f, %.2f]%%\n', ...
        results.pct_est(i), results.pct_lo(i), results.pct_hi(i));
    if isfield(results, 'pValue') && numel(results.pValue) >= i
        fprintf('  p-value = %.3g\n', results.pValue(i));
    end

    % Simple interpretation: case #1 vs case #20 (difference of 19 cases).
    deltaCases = 19;
    pct_delta = (exp(results.beta(i) * deltaCases) - 1) * 100;
    pct_delta_lo = (exp(results.ci(i, 1) * deltaCases) - 1) * 100;
    pct_delta_hi = (exp(results.ci(i, 2) * deltaCases) - 1) * 100;
    fprintf('  Approximate change from Affera case #1 to #20: %.1f%% [%.1f, %.1f]%%\n', ...
        pct_delta, pct_delta_lo, pct_delta_hi);
end

fprintf('\nModel fit summary (AIC/BIC): AIC = %.1f, BIC = %.1f\n', ...
    lme.ModelCriterion.AIC, lme.ModelCriterion.BIC);

% Collinearity diagnostics (printed only when provided).
if ~isempty(diagnostics)
    fprintf('\nCollinearity diagnostics (learning curve enabled):\n');
    preds = diagnostics.predictors;
    if ~isempty(diagnostics.corrMatrix)
        fprintf('Correlation matrix (predictors):\n');
        fprintf('%15s', '');
        for j = 1:numel(preds)
            fprintf('%12s', preds{j});
        end
        fprintf('\n');
        for r = 1:numel(preds)
            fprintf('%15s', preds{r});
            for c = 1:numel(preds)
                fprintf('%12.2f', diagnostics.corrMatrix(r, c));
            end
            fprintf('\n');
        end
    end
    if ~isempty(diagnostics.vif)
        fprintf('\nVariance inflation factors:\n');
        for j = 1:numel(preds)
            vifVal = diagnostics.vif(j);
            if isnan(vifVal)
                fprintf('  %-20s VIF: NaN\n', preds{j});
            else
                fprintf('  %-20s VIF: %.2f\n', preds{j}, vifVal);
            end
        end
    end
    if isfield(diagnostics, 'messages') && ~isempty(diagnostics.messages)
        fprintf('\nDiagnostics messages:\n');
        for i = 1:numel(diagnostics.messages)
            fprintf('  - %s\n', diagnostics.messages{i});
        end
    end
end

end
