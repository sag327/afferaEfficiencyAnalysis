function diagnostics = compute_learning_curve_diagnostics(tbl_q1)
% COMPUTE_LEARNING_CURVE_DIAGNOSTICS  Collinearity checks for learning term.
%
%   diagnostics = compute_learning_curve_diagnostics(tbl_q1) computes
%   correlations and variance inflation factors (VIF) for the predictors
%   used in the learning-curve-augmented model.

predictors = {'isAffera', 'isPVIplus', 'time_days', 'affera_index_ctr'};
missing = setdiff(predictors, tbl_q1.Properties.VariableNames);
if ~isempty(missing)
    error('compute_learning_curve_diagnostics:MissingVariables', ...
        'Table is missing required variables: %s', strjoin(missing, ', '));
end

X = tbl_q1{:, predictors};
rows_ok = all(isfinite(X), 2);
X_clean = X(rows_ok, :);
diagnostics = struct();
diagnostics.predictors = predictors;

if isempty(X_clean)
    warning('compute_learning_curve_diagnostics:NoCompleteRows', ...
        'No complete rows available for diagnostics.');
    diagnostics.corrMatrix = [];
    diagnostics.vif = [];
    diagnostics.highCorrMask = [];
    diagnostics.vifWarnings = struct('gt5', [], 'gt10', []);
    diagnostics.messages = "No complete rows; diagnostics skipped.";
    return;
end

corrMatrix = corr(X_clean, 'Rows', 'complete');
diagnostics.corrMatrix = corrMatrix;

absCorr = abs(corrMatrix);
highCorrMask = absCorr > 0.7 & ~eye(size(absCorr));
diagnostics.highCorrMask = highCorrMask;

p = size(X_clean, 2);
vifValues = nan(p, 1);
messages = strings(0, 1);

if size(X_clean, 1) > p
    for j = 1:p
        y = X_clean(:, j);
        Z = X_clean(:, setdiff(1:p, j));
        try
            mdl = fitlm(Z, y);
            R2 = mdl.Rsquared.Ordinary;
            if R2 < 1
                vifValues(j) = 1 / (1 - R2);
            else
                vifValues(j) = Inf;
            end
        catch e %#ok<NASGU>
            vifValues(j) = NaN;
            messages(end+1, 1) = sprintf('VIF computation failed for %s.', predictors{j}); %#ok<AGROW>
        end
    end
else
    messages(end+1, 1) = "Insufficient rows for VIF computation (n <= p).";
end

diagnostics.vif = vifValues;
diagnostics.vifWarnings = struct('gt5', vifValues > 5, 'gt10', vifValues > 10);

if any(highCorrMask(:))
    messages(end+1, 1) = "High correlation detected (|r| > 0.7).";
end
if any(diagnostics.vifWarnings.gt10, 'all')
    messages(end+1, 1) = "VIF > 10 detected.";
elseif any(diagnostics.vifWarnings.gt5, 'all')
    messages(end+1, 1) = "VIF > 5 detected.";
end
diagnostics.messages = messages;

% Emit brief warnings to the command window.
for i = 1:numel(messages)
    fprintf('[Diagnostics] %s\n', messages{i});
end

end
