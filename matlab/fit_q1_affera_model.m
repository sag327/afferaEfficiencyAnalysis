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
hasPostNonPFAMask = ismember('isPostNonPFAEra', tbl_q1.Properties.VariableNames);
if hasPostNonPFAMask
    postNonPFAMask = logical(tbl_q1.isPostNonPFAEra);
else
    postNonPFAMask = false(size(baselineMask));
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

post_non_pfa_stats = struct();
if hasPostNonPFAMask
    post_non_pfa_stats.overall  = statsFn(postNonPFAMask);
    post_non_pfa_stats.pvi      = statsFn(postNonPFAMask & ~isPVIplus);
    post_non_pfa_stats.pvi_plus = statsFn(postNonPFAMask &  isPVIplus);
    % Catheter breakdown for post-Affera non-PFA cases.
    if ismember('catheter_primary', tbl_q1.Properties.VariableNames)
        cats = categories(tbl_q1.catheter_primary);
        nCats = numel(cats);
        catNames = strings(nCats, 1);
        nVec     = zeros(nCats, 1);
        meanVec  = nan(nCats, 1);
        medVec   = nan(nCats, 1);
        for k = 1:nCats
            cat = cats{k};
            mask = postNonPFAMask & (tbl_q1.catheter_primary == cat);
            catNames(k) = string(cat);
            nVec(k)     = sum(mask);
            meanVec(k)  = mean(duration(mask), 'omitnan');
            medVec(k)   = median(duration(mask), 'omitnan');
        end
        post_non_pfa_stats.by_catheter = table(catNames, nVec, meanVec, medVec, ...
            'VariableNames', {'catheter_primary', 'n', 'mean_duration', 'median_duration'});
    else
        post_non_pfa_stats.by_catheter = table();
    end
else
    post_non_pfa_stats.overall  = statsFn(false(size(baselineMask)));
    post_non_pfa_stats.pvi      = post_non_pfa_stats.overall;
    post_non_pfa_stats.pvi_plus = post_non_pfa_stats.overall;
    post_non_pfa_stats.by_catheter = table();
end

results.baseline_stats = baseline_stats;
results.affera_stats   = affera_stats;
results.post_non_pfa_stats = post_non_pfa_stats;

end
