function tbl_q1 = add_operator_baseline_speed(tbl_q1)
% ADD_OPERATOR_BASELINE_SPEED  Add operator-level baseline speed covariates.
%
%   tbl_q1 = add_operator_baseline_speed(tbl_q1) computes, for each
%   operator, the mean log-duration across baseline non-PFA cases
%   (isBaselineEra == 1) and attaches:
%     - baseline_speed_raw: operator's mean baseline log-duration
%     - baseline_speed_ctr: centered baseline_speed_raw (mean zero across ops)

requiredVars = {'operator_id', 'log_duration', 'isBaselineEra'};
missingVars = setdiff(requiredVars, tbl_q1.Properties.VariableNames);
if ~isempty(missingVars)
    error('add_operator_baseline_speed:MissingVariables', ...
        'Table is missing required variables: %s', strjoin(missingVars, ', '));
end

ops = categories(tbl_q1.operator_id);
nOps = numel(ops);
baseline_speed_raw_op = nan(nOps, 1);

for k = 1:nOps
    op = ops{k};
    maskBaselineOp = (tbl_q1.operator_id == op) & logical(tbl_q1.isBaselineEra);
    if any(maskBaselineOp)
        baseline_speed_raw_op(k) = mean(tbl_q1.log_duration(maskBaselineOp), 'omitnan');
    end
end

mu_baseline = mean(baseline_speed_raw_op, 'omitnan');
baseline_speed_ctr_op = baseline_speed_raw_op - mu_baseline;

n = height(tbl_q1);
tbl_q1.baseline_speed_raw = nan(n, 1);
tbl_q1.baseline_speed_ctr = nan(n, 1);

for k = 1:nOps
    maskOp = (tbl_q1.operator_id == ops{k});
    tbl_q1.baseline_speed_raw(maskOp) = baseline_speed_raw_op(k);
    tbl_q1.baseline_speed_ctr(maskOp) = baseline_speed_ctr_op(k);
end

end
