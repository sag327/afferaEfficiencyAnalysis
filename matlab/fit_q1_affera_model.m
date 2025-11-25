function [lme, results] = fit_q1_affera_model(tbl_q1, opts)
% FIT_Q1_AFFERA_MODEL  Fit the specified mixed-effects model for Q1.
%
%   [lme, results] = fit_q1_affera_model(tbl_q1) fits the linear mixed-
%   effects model:
%
%       log_duration ~ isAffera * isPVIplus + time_days + (1|operator_id)
%
%   and returns the LinearMixedModel object (lme) along with a results
%   struct containing fixed-effect estimates, confidence intervals, and
%   percent effects, as well as descriptive duration statistics for
%   baseline and Affera cohorts.
%
%   [...] = fit_q1_affera_model(tbl_q1, opts) enables an optional learning-
%   curve term (Affera case index) when opts.enableLearningCurve is true.

if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'enableLearningCurve')
    opts.enableLearningCurve = false;
end
if ~isfield(opts, 'EnableOperatorEffect')
    opts.EnableOperatorEffect = false;
end

requiredVars = {'log_duration', 'isAffera', 'isPVIplus', 'time_days', 'operator_id'};
missingVars = setdiff(requiredVars, tbl_q1.Properties.VariableNames);
if ~isempty(missingVars)
    error('fit_q1_affera_model:MissingVariables', ...
        'Analysis table is missing required variables: %s', strjoin(missingVars, ', '));
end

fixedFormula = 'log_duration ~ isAffera * isPVIplus + time_days';

if opts.enableLearningCurve
    if ~ismember('affera_index_ctr', tbl_q1.Properties.VariableNames)
        error('fit_q1_affera_model:MissingAfferaIndex', ...
            'Learning curve enabled but affera_index_ctr is missing.');
    end
    fixedFormula = [fixedFormula, ' + affera_index_ctr'];
end

if opts.EnableOperatorEffect
    if ~ismember('baseline_speed_ctr', tbl_q1.Properties.VariableNames)
        error('fit_q1_affera_model:MissingBaselineSpeed', ...
            'Operator effect enabled but baseline_speed_ctr is missing.');
    end
    fixedFormula = [fixedFormula, ' + baseline_speed_ctr'];
    fixedFormula = [fixedFormula, ' + isAffera:baseline_speed_ctr'];
end

fprintf('final fixed formula for lme: %s', fixedFormula);

lme = fitlme(tbl_q1, sprintf('%s + (1|operator_id)', fixedFormula));

beta  = fixedEffects(lme);
ci    = coefCI(lme);
names = lme.CoefficientNames(:);
pVals = lme.Coefficients.pValue;

pct_est = (exp(beta)     - 1) * 100;
pct_lo  = (exp(ci(:, 1)) - 1) * 100;
pct_hi  = (exp(ci(:, 2)) - 1) * 100;

idxAffera = find(strcmp(names, 'isAffera'));
idxAfferaPVIplus = find(strcmp(names, 'isAffera:isPVIplus'));
if isempty(idxAfferaPVIplus)
    base = 'isAffera:isPVIplus';
    mask = strncmp(names, base, numel(base));
    idxAfferaPVIplus = find(mask, 1);
end
idxAfferaIndex = find(strcmp(names, 'affera_index_ctr'));
idxBaselineSpeed = find(strcmp(names, 'baseline_speed_ctr'));
idxAfferaBaseline = find(strcmp(names, 'isAffera:baseline_speed_ctr'));
if isempty(idxAfferaBaseline)
    base = 'isAffera:baseline_speed_ctr';
    mask = strncmp(names, base, numel(base));
    idxAfferaBaseline = find(mask, 1);
end

% Overall Affera effect across PVI/PVI+ case mix (linear contrast).
overallAffera = struct();
if ~isempty(idxAffera)
    % Case-mix weight for PVI+ (using analysis dataset).
    pPVIplus = mean(logical(tbl_q1.isPVIplus));

    c = zeros(numel(beta), 1);
    c(idxAffera) = 1;
    if ~isempty(idxAfferaPVIplus)
        c(idxAfferaPVIplus) = pPVIplus;
    end

    beta_overall = c' * beta;
    covB = lme.CoefficientCovariance;
    var_overall = c' * covB * c;
    se_overall = sqrt(max(var_overall, 0));

    df = lme.DFE;
    if isfinite(df) && df > 0 && se_overall > 0
        tcrit = tinv(0.975, df);
        ci_overall = [beta_overall - tcrit * se_overall, ...
                      beta_overall + tcrit * se_overall];
    else
        ci_overall = [NaN NaN];
    end

    p_overall = NaN;
    try
        p_overall = coefTest(lme, c');
    catch
        % Leave p_overall as NaN if contrast test fails.
    end

    pct_est_overall = (exp(beta_overall)      - 1) * 100;
    pct_lo_overall  = (exp(ci_overall(1))     - 1) * 100;
    pct_hi_overall  = (exp(ci_overall(2))     - 1) * 100;

    overallAffera.beta     = beta_overall;
    overallAffera.ci       = ci_overall;
    overallAffera.pct_est  = pct_est_overall;
    overallAffera.pct_lo   = pct_lo_overall;
    overallAffera.pct_hi   = pct_hi_overall;
    overallAffera.pValue   = p_overall;
    overallAffera.pPVIplus = pPVIplus;
    overallAffera.hasInteraction = ~isempty(idxAfferaPVIplus);
end

results = struct();
results.names        = names;
results.beta         = beta;
results.ci           = ci;
results.pct_est      = pct_est;
results.pct_lo       = pct_lo;
results.pct_hi       = pct_hi;
results.pValue       = pVals;
results.idxAffera    = idxAffera;
results.idxAfferaPVIplus = idxAfferaPVIplus;
results.idxAfferaIndex   = idxAfferaIndex;
results.idxBaselineSpeed = idxBaselineSpeed;
results.idxAfferaBaseline = idxAfferaBaseline;
results.overallAffera    = overallAffera;

% -------------------------------------------------------------------------
% Descriptive duration summaries (baseline vs Affera, PVI vs PVI+,
% and baseline by catheter_primary).
% -------------------------------------------------------------------------

% Determine duration variable to use.
if ismember('duration_minutes', tbl_q1.Properties.VariableNames)
    duration = tbl_q1.duration_minutes;
elseif ismember('procedure_duration_min', tbl_q1.Properties.VariableNames)
    duration = tbl_q1.procedure_duration_min;
else
    error('fit_q1_affera_model:MissingDuration', ...
        'Table must contain duration_minutes or procedure_duration_min.');
end

baselineMask   = logical(tbl_q1.isBaselineEra);
afferaMask     = tbl_q1.isAffera ~= 0;
isPVIplus      = logical(tbl_q1.isPVIplus);
hasComparisonMask = ismember('isComparisonEra', tbl_q1.Properties.VariableNames);
hasPostNonPFAMask = ismember('isPostNonPFAEra', tbl_q1.Properties.VariableNames);
if hasComparisonMask
    comparisonMask = logical(tbl_q1.isComparisonEra);
elseif hasPostNonPFAMask
    comparisonMask = logical(tbl_q1.isPostNonPFAEra);
else
    comparisonMask = false(size(baselineMask));
end

statsFn = @(m) struct( ...
    'n', sum(m), ...
    'mean_duration', mean(duration(m), 'omitnan'), ...
    'median_duration', median(duration(m), 'omitnan'));

baseline_stats = struct();
baseline_stats.overall   = statsFn(baselineMask);
baseline_stats.pvi       = statsFn(baselineMask & ~isPVIplus);
baseline_stats.pvi_plus  = statsFn(baselineMask &  isPVIplus);

% Baseline by catheter_primary (all PVI statuses combined).
if ismember('catheter_primary', tbl_q1.Properties.VariableNames)
    cats = categories(tbl_q1.catheter_primary);
    nCats = numel(cats);
    catNames = strings(nCats, 1);
    nVec     = zeros(nCats, 1);
    meanVec  = nan(nCats, 1);
    medVec   = nan(nCats, 1);
    for k = 1:nCats
        cat = cats{k};
        mask = baselineMask & (tbl_q1.catheter_primary == cat);
        catNames(k) = string(cat);
        nVec(k)     = sum(mask);
        meanVec(k)  = mean(duration(mask), 'omitnan');
        medVec(k)   = median(duration(mask), 'omitnan');
    end
    baseline_stats.by_catheter = table(catNames, nVec, meanVec, medVec, ...
        'VariableNames', {'catheter_primary', 'n', 'mean_duration', 'median_duration'});
else
    baseline_stats.by_catheter = table();
end

affera_stats = struct();
affera_stats.overall  = statsFn(afferaMask);
affera_stats.pvi      = statsFn(afferaMask & ~isPVIplus);
affera_stats.pvi_plus = statsFn(afferaMask &  isPVIplus);

comparison_stats = struct();
if any(comparisonMask)
    comparison_stats.overall  = statsFn(comparisonMask);
    comparison_stats.pvi      = statsFn(comparisonMask & ~isPVIplus);
    comparison_stats.pvi_plus = statsFn(comparisonMask &  isPVIplus);
    % Catheter breakdown for post-Affera comparison-group cases.
    if ismember('catheter_primary', tbl_q1.Properties.VariableNames)
        cats = categories(tbl_q1.catheter_primary);
        nCats = numel(cats);
        catNames = strings(nCats, 1);
        nVec     = zeros(nCats, 1);
        meanVec  = nan(nCats, 1);
        medVec   = nan(nCats, 1);
        for k = 1:nCats
            cat = cats{k};
            mask = comparisonMask & (tbl_q1.catheter_primary == cat);
            catNames(k) = string(cat);
            nVec(k)     = sum(mask);
            meanVec(k)  = mean(duration(mask), 'omitnan');
            medVec(k)   = median(duration(mask), 'omitnan');
        end
        comparison_stats.by_catheter = table(catNames, nVec, meanVec, medVec, ...
            'VariableNames', {'catheter_primary', 'n', 'mean_duration', 'median_duration'});
    else
        comparison_stats.by_catheter = table();
    end
else
    comparison_stats.overall  = statsFn(false(size(baselineMask)));
    comparison_stats.pvi      = comparison_stats.overall;
    comparison_stats.pvi_plus = comparison_stats.overall;
    comparison_stats.by_catheter = table();
end

% -------------------------------------------------------------------------
% Operator-level baseline duration and fast/slow group comparisons
% (only when baseline_speed_ctr and Affera×baseline are present).
% -------------------------------------------------------------------------

fastSlow = struct();
if ~isempty(idxBaselineSpeed) && ~isempty(idxAfferaBaseline) && ...
        ismember('baseline_speed_ctr', tbl_q1.Properties.VariableNames)

    % Operator-level baseline duration (minutes) and baseline_speed_ctr.
    ops = categories(tbl_q1.operator_id);
    nOps = numel(ops);
    opBaselineMeanMin   = nan(nOps, 1);
    opBaselineMedianMin = nan(nOps, 1);
    opBaselineSpeedCtr  = nan(nOps, 1);

    for k = 1:nOps
        op = ops{k};
        maskOpBase = (tbl_q1.operator_id == op) & baselineMask & ~afferaMask;
        if any(maskOpBase)
            opBaselineMeanMin(k)   = mean(duration(maskOpBase), 'omitnan');
            opBaselineMedianMin(k) = median(duration(maskOpBase), 'omitnan');
            opBaselineSpeedCtr(k)  = mean(tbl_q1.baseline_speed_ctr(maskOpBase), 'omitnan');
        end
    end

    validOps = isfinite(opBaselineMeanMin) & isfinite(opBaselineSpeedCtr);
    if any(validOps)
        fastSlow.baseline_mean_per_op   = opBaselineMeanMin;
        fastSlow.baseline_median_per_op = opBaselineMedianMin;
        fastSlow.baseline_speed_ctr_per_op = opBaselineSpeedCtr;

        % Quartile-based groups on operator-level baseline mean duration.
        q = quantile(opBaselineMeanMin(validOps), [0.25 0.50 0.75]);
        q1 = q(1);
        q2 = q(2);
        q3 = q(3);
        maskQ1 = validOps & (opBaselineMeanMin <= q1);
        maskQ2 = validOps & (opBaselineMeanMin > q1  & opBaselineMeanMin <= q2);
        maskQ3 = validOps & (opBaselineMeanMin > q2  & opBaselineMeanMin <= q3);
        maskQ4 = validOps & (opBaselineMeanMin > q3);

        fastMaskOp = maskQ1;
        slowMaskOp = maskQ4;

        if any(fastMaskOp) && any(slowMaskOp)
            fastSlow.fast_operator_ids = categorical(ops(fastMaskOp));
            fastSlow.slow_operator_ids = categorical(ops(slowMaskOp));
            fastSlow.n_fast_ops = sum(fastMaskOp);
            fastSlow.n_slow_ops = sum(slowMaskOp);

            fastSlow.baseline_mean_duration_fast   = mean(opBaselineMeanMin(fastMaskOp), 'omitnan');
            fastSlow.baseline_median_duration_fast = median(opBaselineMeanMin(fastMaskOp), 'omitnan');
            fastSlow.baseline_mean_duration_slow   = mean(opBaselineMeanMin(slowMaskOp), 'omitnan');
            fastSlow.baseline_median_duration_slow = median(opBaselineMeanMin(slowMaskOp), 'omitnan');

            % Build contrasts for each quartile group (PVI-only Affera effect).
            covB = lme.CoefficientCovariance;
            df = lme.DFE;
            tcrit = tinv(0.975, df);

            quartMasks = {maskQ1, maskQ2, maskQ3, maskQ4};
            quartiles(1:4) = struct('n_ops', 0);
            refContrast = [];

            for qIdx = 1:4
                m = quartMasks{qIdx};
                if ~any(m)
                    quartiles(qIdx).n_ops = 0;
                    continue;
                end
                b_q = mean(opBaselineSpeedCtr(m), 'omitnan');
                mean_q = mean(opBaselineMeanMin(m), 'omitnan');
                med_q  = median(opBaselineMeanMin(m), 'omitnan');

                c_q = zeros(numel(beta), 1);
                c_q(idxAffera) = 1;
                c_q(idxAfferaBaseline) = b_q;

                beta_q = c_q' * beta;
                var_q  = c_q' * covB * c_q;
                se_q   = sqrt(max(var_q, 0));
                ci_q_log = [beta_q - tcrit * se_q, ...
                            beta_q + tcrit * se_q];
                pct_q     = (exp(beta_q)      - 1) * 100;
                pct_q_lo  = (exp(ci_q_log(1)) - 1) * 100;
                pct_q_hi  = (exp(ci_q_log(2)) - 1) * 100;
                p_q = NaN;
                if se_q > 0 && isfinite(df) && df > 0
                    t_q = beta_q / se_q;
                    p_q = 2 * tcdf(-abs(t_q), df);
                end

                quartiles(qIdx).n_ops = sum(m);
                quartiles(qIdx).baseline_mean_duration = mean_q;
                quartiles(qIdx).baseline_median_duration = med_q;
                quartiles(qIdx).b = b_q;
                quartiles(qIdx).beta = beta_q;
                quartiles(qIdx).ci = ci_q_log;
                quartiles(qIdx).pct_est = pct_q;
                quartiles(qIdx).pct_lo = pct_q_lo;
                quartiles(qIdx).pct_hi = pct_q_hi;
                quartiles(qIdx).pValue = p_q;

                if qIdx == 1
                    refContrast = c_q;
                end
            end

            fastSlow.quartiles = quartiles;

            % Also keep explicit fast (Q1) and slow (Q4) entries and their difference.
            if quartiles(1).n_ops > 0 && quartiles(4).n_ops > 0
                fastSlow.fast = quartiles(1);
                fastSlow.slow = quartiles(4);

                c_diff = zeros(numel(beta), 1);
                c_diff = (quartiles(4).b - quartiles(1).b) * 0;
                % Directly use difference of contrasts for slow vs fast.
                c_fast = refContrast;
                c_slow = zeros(numel(beta), 1);
                c_slow(idxAffera) = 1;
                c_slow(idxAfferaBaseline) = quartiles(4).b;
                c_diff = c_slow - c_fast;

                beta_diff = c_diff' * beta;
                var_diff  = c_diff' * covB * c_diff;
                se_diff   = sqrt(max(var_diff, 0));
                ci_diff_log = [beta_diff - tcrit * se_diff, ...
                               beta_diff + tcrit * se_diff];
                pct_diff     = (exp(beta_diff)      - 1) * 100;
                pct_diff_lo  = (exp(ci_diff_log(1)) - 1) * 100;
                pct_diff_hi  = (exp(ci_diff_log(2)) - 1) * 100;
                p_diff = NaN;
                if se_diff > 0 && isfinite(df) && df > 0
                    t_diff = beta_diff / se_diff;
                    p_diff = 2 * tcdf(-abs(t_diff), df);
                end

                fastSlow.delta_slow_vs_fast = struct('beta', beta_diff, 'ci', ci_diff_log, ...
                    'pct_est', pct_diff, 'pct_lo', pct_diff_lo, 'pct_hi', pct_diff_hi, ...
                    'pValue', p_diff);
            end
        end
    end
end

results.baseline_stats = baseline_stats;
results.affera_stats   = affera_stats;
results.comparison_stats = comparison_stats;
% Backward-compatible field name (legacy non-PFA naming).
results.post_non_pfa_stats = comparison_stats;
if ~isempty(fieldnames(fastSlow))
    results.fastSlow = fastSlow;
end

end
