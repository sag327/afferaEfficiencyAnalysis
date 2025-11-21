function [lme, results] = fit_q1_affera_model(tbl_q1)
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

requiredVars = {'log_duration', 'isAffera', 'isPVIplus', 'time_days', 'operator_id'};
missingVars = setdiff(requiredVars, tbl_q1.Properties.VariableNames);
if ~isempty(missingVars)
    error('fit_q1_affera_model:MissingVariables', ...
        'Analysis table is missing required variables: %s', strjoin(missingVars, ', '));
end

lme = fitlme(tbl_q1, ...
    'log_duration ~ isAffera * isPVIplus + time_days + (1|operator_id)');

beta  = fixedEffects(lme);
ci    = coefCI(lme);
names = lme.CoefficientNames(:);

pct_est = (exp(beta)     - 1) * 100;
pct_lo  = (exp(ci(:, 1)) - 1) * 100;
pct_hi  = (exp(ci(:, 2)) - 1) * 100;

idxAffera = find(strcmp(names, 'isAffera'));
idxAfferaPVIplus = find(strcmp(names, 'isAffera:isPVIplus'));

results = struct();
results.names        = names;
results.beta         = beta;
results.ci           = ci;
results.pct_est      = pct_est;
results.pct_lo       = pct_lo;
results.pct_hi       = pct_hi;
results.idxAffera    = idxAffera;
results.idxAfferaPVIplus = idxAfferaPVIplus;

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
