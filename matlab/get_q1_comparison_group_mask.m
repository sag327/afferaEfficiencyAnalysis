function mask = get_q1_comparison_group_mask(tbl, comparisonGroup)
% GET_Q1_COMPARISON_GROUP_MASK  Build mask for the chosen comparison group.
%
%   mask = get_q1_comparison_group_mask(tbl, comparisonGroup) returns a
%   logical mask identifying comparison-group cases (non-Affera) according
%   to the requested comparisonGroup string:
%     - 'nonPFA' (default): non-PFA, non-Affera cases.
%     - 'RF': RF cases (is_rf_catheter), non-Affera.
%
%   This centralizes the definition of the comparison group so that all
%   downstream code (dataset construction, summaries) uses the same logic.

if nargin < 2 || isempty(comparisonGroup)
    comparisonGroup = "nonPFA";
end

if ~ismember('is_affera', tbl.Properties.VariableNames)
    error('get_q1_comparison_group_mask:MissingVariable', ...
        'Table is missing required variable ''is_affera''.');
end

compOpt = lower(string(comparisonGroup));

switch compOpt
    case "nonpfa"
        if ~ismember('is_pfa_catheter', tbl.Properties.VariableNames)
            error('get_q1_comparison_group_mask:MissingVariable', ...
                'Table is missing required variable ''is_pfa_catheter'' for nonPFA comparison.');
        end
        mask = ~tbl.is_pfa_catheter & ~tbl.is_affera;

    case "rf"
        if ~ismember('is_rf_catheter', tbl.Properties.VariableNames)
            % Compute RF indicator on the fly if possible.
            hasSupplies = ismember('supplies_raw', tbl.Properties.VariableNames);
            hasCatheter = ismember('catheter_primary', tbl.Properties.VariableNames);
            if ~(hasSupplies || hasCatheter)
                error('get_q1_comparison_group_mask:MissingVariable', ...
                    'Table is missing required variable ''is_rf_catheter'' for RF comparison.');
            end
            is_rf = false(height(tbl), 1);
            if hasCatheter
                is_rf = is_rf | strcmpi(tbl.catheter_primary, "OTHER CATHETERS (RF)");
            end
            if hasSupplies
                lowSupLocal = lower(string(tbl.supplies_raw));
                is_rf = is_rf | contains(lowSupLocal, 'other catheters (rf)');
            end
            mask = is_rf & ~tbl.is_affera;
        else
            mask = tbl.is_rf_catheter & ~tbl.is_affera;
        end

    otherwise
        error('get_q1_comparison_group_mask:InvalidOption', ...
            'Unknown comparisonGroup option: %s (use ''nonPFA'' or ''RF'')', comparisonGroup);
end

if ~any(mask)
    error('get_q1_comparison_group_mask:EmptyMask', ...
        'No comparison-group cases found for option %s.', comparisonGroup);
end

end
