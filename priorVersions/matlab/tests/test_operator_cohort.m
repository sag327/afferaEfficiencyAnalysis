% test_operator_cohort.m
%
% Command-line regression test that verifies run_basic_affera_analysis uses
% a shared operator cohort across ALL / PVI / PVI+ analyses and that the
% learning-curve models (including Q3/Q4) are fit successfully on synthetic
% data.

testDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(testDir);
addpath(matlabDir); %#ok<*MCAP>

tbl = build_synthetic_table();

seqResults = run_basic_affera_analysis(tbl, ...
    'ProcedureType', 'sequential', ...
    'MinAfferaCases', 1, ...
    'baselineCaseNumberRange', [1 2], ...
    'Verbose', true);

ops_all = sort(string(seqResults.all.eligible_operator_ids));
ops_pvi = sort(string(seqResults.PVI.eligible_operator_ids));
ops_pvip = sort(string(seqResults.PVIplus.eligible_operator_ids));

assert(isequal(ops_all, ops_pvi) && isequal(ops_all, ops_pvip), ...
    'Operator cohorts differ between ALL/PVI/PVI+ analyses.');

assert(~isempty(seqResults.all.lme_q3) && ~isempty(seqResults.all.lme_q4), ...
    'ALL dataset did not produce Q3/Q4 models.');
assert(~isempty(seqResults.PVI.lme_q3) && ~isempty(seqResults.PVI.lme_q4), ...
    'PVI dataset did not produce Q3/Q4 models.');
assert(~isempty(seqResults.PVIplus.lme_q3) && ~isempty(seqResults.PVIplus.lme_q4), ...
    'PVI+ dataset did not produce Q3/Q4 models.');

disp('test_operator_cohort: all checks passed.');

function tbl = build_synthetic_table()
operators = ["OpA", "OpB", "OpC"];
nOps = numel(operators);

rows = [];
idx = 0;

for oo = 1:nOps
    op = operators(oo);
    baseDate = datetime(2024, 1, 1) + days(oo);

    entries = [
        struct('dt', baseDate - days(30), 'type', "PVI",  'is_affera', false, 'is_pre', true),  ...
        struct('dt', baseDate - days(28), 'type', "PVI+", 'is_affera', false, 'is_pre', true),  ...
        struct('dt', baseDate - days(26), 'type', "PVI",  'is_affera', false, 'is_pre', true),  ...
        struct('dt', baseDate + days(2),  'type', "PVI",  'is_affera', true,  'is_pre', false), ...
        struct('dt', baseDate + days(4),  'type', "PVI+", 'is_affera', true,  'is_pre', false), ...
        struct('dt', baseDate + days(7),  'type', "PVI",  'is_affera', true,  'is_pre', false), ...
        struct('dt', baseDate + days(12), 'type', "PVI+", 'is_affera', true,  'is_pre', false), ...
        struct('dt', baseDate + days(16), 'type', "PVI",  'is_affera', true,  'is_pre', false), ...
        struct('dt', baseDate + days(21), 'type', "PVI+", 'is_affera', true,  'is_pre', false), ...
        struct('dt', baseDate + days(8),  'type', "PVI+", 'is_affera', false, 'is_pre', false), ...
        struct('dt', baseDate + days(10), 'type', "PVI",  'is_affera', false, 'is_pre', false)
    ];

    for ee = 1:numel(entries)
        idx = idx + 1;
        rows(idx).procedure_id = sprintf('%s_case_%02d', op, ee); %#ok<AGROW>
        rows(idx).operator_id = op;
        rows(idx).operator_name = op;
        rows(idx).procedure_date = entries(ee).dt;
        rows(idx).procedure_type = entries(ee).type;
        rows(idx).tool_type = ternary(entries(ee).is_affera, "affera", "other");
        rows(idx).is_affera = entries(ee).is_affera;
        rows(idx).is_pre_affera_global = entries(ee).is_pre;
        rows(idx).is_post_affera_global = ~entries(ee).is_pre;
        rows(idx).procedure_duration_min = 90 - 5*ee + 2*oo;
        rows(idx).is_pfa_catheter = false;
    end
end

tbl = struct2table(rows);
tbl.tool_type = categorical(tbl.tool_type);
tbl.operator_id = categorical(tbl.operator_id);
tbl.operator_name = categorical(tbl.operator_name);
tbl.procedure_type = categorical(tbl.procedure_type);
tbl.procedure_id = string(tbl.procedure_id);
tbl.log_duration = log(tbl.procedure_duration_min);
tbl.time_days = days(tbl.procedure_date - min(tbl.procedure_date));
tbl.is_post_affera_global = logical(tbl.is_post_affera_global);
tbl.is_pre_affera_global = logical(tbl.is_pre_affera_global);
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
