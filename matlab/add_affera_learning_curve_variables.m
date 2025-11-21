function tbl_q1 = add_affera_learning_curve_variables(tbl_q1)
% ADD_AFFERA_LEARNING_CURVE_VARIABLES  Add Affera case index variables.
%
%   tbl_q1 = add_affera_learning_curve_variables(tbl_q1) takes the Q1
%   analysis table and adds:
%     - affera_index: sequential Affera case number per operator (1,2,3,...)
%     - affera_index_ctr: centered version of affera_index

requiredVars = {'operator_id', 'procedure_date', 'isAffera'};
missingVars = setdiff(requiredVars, tbl_q1.Properties.VariableNames);
if ~isempty(missingVars)
    error('add_affera_learning_curve_variables:MissingVariables', ...
        'Table is missing required variables: %s', strjoin(missingVars, ', '));
end

n = height(tbl_q1);
tbl_q1.affera_index     = zeros(n, 1);
tbl_q1.affera_index_ctr = zeros(n, 1);

isA = tbl_q1.isAffera == 1;
if ~any(isA)
    warning('add_affera_learning_curve_variables:NoAfferaCases', ...
        'No Affera cases found; affera_index columns remain zero.');
    return;
end

ops = categories(tbl_q1.operator_id);
hasProcId = ismember('procedure_id', tbl_q1.Properties.VariableNames);
for k = 1:numel(ops)
    op = ops{k};
    maskA = isA & tbl_q1.operator_id == op;
    if ~any(maskA)
        continue;
    end
    rows = find(maskA);
    if hasProcId
        [~, order] = sortrows(tbl_q1(rows, :), {'procedure_date', 'procedure_id'});
    else
        [~, order] = sort(tbl_q1.procedure_date(rows));
    end
tbl_q1.affera_index(rows(order)) = (1:numel(rows)).';
end

mean_idx = mean(tbl_q1.affera_index(isA), 'omitnan');
tbl_q1.affera_index_ctr(isA) = tbl_q1.affera_index(isA) - mean_idx;
tbl_q1.affera_index_ctr(~isA) = 0;

end
