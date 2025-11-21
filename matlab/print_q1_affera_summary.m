function print_q1_affera_summary(lme, results, tbl_q1)
% PRINT_Q1_AFFERA_SUMMARY  Display a human-readable summary of Q1 results.
%
%   print_q1_affera_summary(lme, results, tbl_q1) prints key sample counts
%   and interpretable fixed-effect estimates, including Affera effects in
%   PVI-only and PVI+ cases.

fprintf('\n=== Q1 Affera vs Non-Affera Mixed-Effects Model ===\n');

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

% Fixed-effect summary table.
fprintf('\nFixed effects (log-duration scale and percent effect):\n');
fprintf('%-25s %10s %22s %26s\n', 'Term', 'Beta', '95% CI (log scale)', 'Percent effect [95% CI]');

for i = 1:numel(results.names)
    name = results.names{i};
    b    = results.beta(i);
    ci   = results.ci(i, :);
    pe   = results.pct_est(i);
    peLo = results.pct_lo(i);
    peHi = results.pct_hi(i);

    fprintf('%-25s %10.3f [%7.3f, %7.3f] %10.1f%% [%6.1f, %6.1f]%%\n', ...
        name, b, ci(1), ci(2), pe, peLo, peHi);
end

% Highlight key Affera effects.
if ~isempty(results.idxAffera)
    i = results.idxAffera;
    fprintf('\nAffera effect in PVI-only cases (isAffera):\n');
    fprintf('  Beta = %.3f, 95%% CI = [%.3f, %.3f]\n', ...
        results.beta(i), results.ci(i, 1), results.ci(i, 2));
    fprintf('  Percent change in duration = %.1f%% [%.1f, %.1f]%%\n', ...
        results.pct_est(i), results.pct_lo(i), results.pct_hi(i));
end

if ~isempty(results.idxAfferaPVIplus)
    i = results.idxAfferaPVIplus;
    fprintf('\nAdditional Affera effect in PVI+ vs PVI (isAffera:isPVIplus):\n');
    fprintf('  Beta = %.3f, 95%% CI = [%.3f, %.3f]\n', ...
        results.beta(i), results.ci(i, 1), results.ci(i, 2));
    fprintf('  Percent change (additional) = %.1f%% [%.1f, %.1f]%%\n', ...
        results.pct_est(i), results.pct_lo(i), results.pct_hi(i));
end

fprintf('\nModel fit summary (AIC/BIC): AIC = %.1f, BIC = %.1f\n', ...
    lme.ModelCriterion.AIC, lme.ModelCriterion.BIC);

end
