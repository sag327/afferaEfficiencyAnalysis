# Affera Comparison Group Refactor – Implementation Plan

This document describes how to refactor the Q1 Affera efficiency analysis so that the **comparison group** for Affera can be selected by the user (e.g., **non‑PFA** vs **RF only**), with a **single point of definition** for that group. It also outlines command‑line tests to verify correct behavior and parity.

The goals are:
- Maintain current behavior by default (**non‑PFA comparison**).
- Allow analysts to set the comparison group to **RF only** via a simple option.
- Centralize comparison‑group logic in one helper function.
- Ensure all downstream components (dataset construction, baseline speed, summaries, diagnostics) respect this choice.

---

## 1. Current Behavior (Baseline)

### 1.1 Imported variables

`matlab/import_clinical_data.m` currently constructs:
- `tbl.is_affera` – true when the supplies mention “Sphere-9 Catheter” (Affera).
- `tbl.is_pfa_catheter` – true when `supplies_raw` contains any of:
  - `'sphere-9'`, `'farawave'`, `'pulseselect'`.
- `tbl.catheter_primary` – first token from `supplies_raw`, upper‑cased (e.g., “SPHERE-9 CATHETER”, “OTHER CATHETERS (RF)”, etc.).

RF ablations in this dataset are identifiable via:
- `supplies_raw` containing `'Other catheters (RF)'`, and/or
- `catheter_primary == 'OTHER CATHETERS (RF)'` (exact form to be confirmed).

### 1.2 Comparison group = “non‑PFA”

`matlab/make_q1_affera_dataset.m` currently defines:
- Raw Affera indicator:
  - `isAffera_raw = tbl.is_affera;`
- Baseline “non‑PFA, non‑Affera” cohort:
  - `isBaseline = ~tbl.is_pfa_catheter & ~isAffera_raw & baseline_date_window;`
- Post‑Affera “non‑PFA, non‑Affera” cohort:
  - `isPostNonPFA = ~tbl.is_pfa_catheter & ~isAffera_raw & date >= Affera_start;`
- Inclusion mask:
  - `include = isAffera_raw | isBaseline | isPostNonPFA;`

Downstream:
- `isBaselineEra` and `isPostNonPFAEra` are derived from these masks and used by:
  - operator inclusion logic in `make_q1_affera_dataset.m`,
  - `add_operator_baseline_speed.m` (baseline operator speed),
  - descriptive summaries in `fit_q1_affera_model.m`,
  - printing in `print_q1_affera_summary.m`.

Thus, “non‑PFA” is currently “all non‑Affera cases that are not flagged as PFA by `is_pfa_catheter`”.

---

## 2. High‑Level Refactor Design

### 2.1 New comparison‑group options

We introduce a user‑facing option with two initial modes:
- `ComparisonGroup = 'nonPFA'` (default, matches existing behavior).
- `ComparisonGroup = 'RF'` (restrict comparison group to RF ablations only).

Future extensions (if needed) could add other groups (e.g., “all non‑Affera”, “non‑PFA PVI‑only”), but the refactor should keep the logic centralized so these are easy to add.

### 2.2 Single source of truth

All logic defining the comparison group will be moved into a dedicated helper, e.g.:

```matlab
% matlab/get_q1_comparison_group_mask.m
function mask = get_q1_comparison_group_mask(tbl, comparisonGroup)
```

Responsibilities:
- Validate required variables (e.g., `is_affera`, `is_pfa_catheter`, `is_rf_catheter`).
- Implement mapping from `comparisonGroup` string → logical mask on `tbl`:
  - `'nonPFA'`: non‑PFA, non‑Affera.
  - `'RF'`: RF only, non‑Affera.
- Provide clear errors/warnings when an option is incompatible with the available data (e.g., no RF cases, or `is_rf_catheter` not present).

All components that need the comparison group will use this helper instead of re‑encoding the logic.

---

## 3. Step‑by‑Step Code Changes

### 3.1 Add RF indicator in import script

**File:** `matlab/import_clinical_data.m`

**Goal:** Explicitly flag RF cases so they can be used as a comparison group.

Planned changes:
- After constructing `catheter_primary` and `supplies_raw`, add:
  - `tbl.is_rf_catheter`: logical, true for RF ablations.
- Base definition on the current data:
  - Use `catheter_primary == "OTHER CATHETERS (RF)"` as the primary rule.
  - Optionally also check `contains(lower(supplies_raw), 'other catheters (rf)')` for robustness.
- Save `is_rf_catheter` in `procedures_all.mat` alongside existing fields.

Sanity checks:
- Print a short RF summary after the existing operator counts, e.g.:
  - Total RF procedures and RF count by operator.
- Confirm that the total RF count matches expectations from the raw Excel file (via a one‑off exploratory script if desired).

### 3.2 New helper for comparison group

**File (new):** `matlab/get_q1_comparison_group_mask.m`

**Signature:**

```matlab
function mask = get_q1_comparison_group_mask(tbl, comparisonGroup)
```

**Responsibilities:**
- Input validation:
  - Ensure `comparisonGroup` is a char/string scalar.
  - Ensure `is_affera` exists.
  - Ensure `is_pfa_catheter` exists for `'nonPFA'`.
  - Ensure `is_rf_catheter` exists for `'RF'`.
- Implement logic:
  - For `'nonPFA'`:
    - `mask = ~tbl.is_pfa_catheter & ~tbl.is_affera;`
  - For `'RF'`:
    - `mask = tbl.is_rf_catheter & ~tbl.is_affera;`
- Edge cases:
  - If `sum(mask) == 0`, throw a descriptive error (no comparison cases found).
  - Optionally, warn if the mask size is small (e.g., `< 20` cases).

This function becomes the only place where the definition of “comparison group” lives.

### 3.3 Extend Q1 dataset builder to accept comparison group

**File:** `matlab/make_q1_affera_dataset.m`

**Goal:** Replace hard‑coded “non‑PFA” logic with a call to `get_q1_comparison_group_mask`, and expose the comparison‑group choice via an argument.

Planned changes:
- Extend function signature, e.g.:

```matlab
function [tbl_q1, operator_counts_all] = make_q1_affera_dataset(tbl, meta, baseline_days, min_cases_per_group, comparisonGroup)
```

with:
- `comparisonGroup` as a required or optional argument (if omitted, default to `'nonPFA'` to preserve current behavior).

Inside the function:
- Use the helper to build the base comparison mask:

```matlab
isAffera_raw = tbl.is_affera;
compMask = get_q1_comparison_group_mask(tbl, comparisonGroup);  % non‑PFA or RF
```

- Redefine cohort flags:
  - `isBaseline = compMask & date-in-baseline-window;`
  - `isPostComp = compMask & date >= Affera_start;`
- Replace:
  - `isPostNonPFA` → `isPostComparison` internally.
  - `include = isAffera_raw | isBaseline | isPostComp;`

For reporting columns on the included table `tbl_inc`:
- Keep semantically clear names:
  - `tbl_inc.isBaselineEra = baselineMaskInc;`
  - `tbl_inc.isComparisonEra = postCompMaskInc;`
- For backward compatibility with existing code:
  - Option 1 (preferred for clarity): adjust downstream code to use `isComparisonEra` instead of `isPostNonPFAEra`.
  - Option 2: set `tbl_inc.isPostNonPFAEra = tbl_inc.isComparisonEra;` with a comment noting that, depending on `comparisonGroup`, this may be RF‑only rather than strictly “non‑PFA”.

Operator inclusion logic:
- No structural change needed; it already counts:
  - `n_baseline`: using `isBaselineEra`.
  - `n_affera`: using `isAffera ~= 0`.
- These counts will now reflect the chosen comparison group (non‑PFA or RF‑only).

### 3.4 Wire comparison group from top‑level API

**File:** `matlab/run_q1_affera_model.m`

**Goal:** Expose a clean user‑facing option and pass it down to the dataset builder and summary.

Planned changes:
- Extend `inputParser` with:

```matlab
p.addParameter('ComparisonGroup', 'nonPFA', @(s) ischar(s) || isstring(s));
```

- After parsing:

```matlab
opts.comparisonGroup = string(p.Results.ComparisonGroup);
```

- When calling `make_q1_affera_dataset`:

```matlab
[tbl_q1, operator_counts_all] = make_q1_affera_dataset(tbl, meta, baseline_days, min_cases_per_group, opts.comparisonGroup);
```

- Attach `opts.comparisonGroup` to the options struct passed into:
  - `fit_q1_affera_model`.
  - `print_q1_affera_summary`.

Backward compatibility:
- Positional arguments remain unchanged.
- If `ComparisonGroup` is not specified, behavior defaults to the existing non‑PFA based definition.

### 3.5 Update summary printing to reflect comparison group

**File:** `matlab/print_q1_affera_summary.m`

**Goals:**
- Report which comparison group was used.
- Rename textual labels so they are comparison‑group aware.

Planned changes:
- In the settings section:

```matlab
if isfield(opts, 'comparisonGroup')
    fprintf('  Comparison group: %s\n', opts.comparisonGroup);
end
```

- Where the code currently prints:
  - “Post‑Affera non‑PFA” counts/durations, switch to a more general label:
    - “Post‑Affera comparison group (non‑PFA)” when `comparisonGroup == 'nonPFA'`.
    - “Post‑Affera comparison group (RF only)” when `comparisonGroup == 'RF'`.

Implementation approach:
- Use `opts.comparisonGroup` to build a human‑readable label:

```matlab
if strcmp(opts.comparisonGroup, 'RF')
    compLabel = 'RF comparison group';
else
    compLabel = 'non-PFA comparison group';
end
```

- Substitute `compLabel` in all summary text where the comparison group is referenced.

No change is required to the fixed‑effects model structure; the refactor only changes which procedures are in the comparison cohort.

### 3.6 Adjust descriptive stats to use generic comparison‑group naming

**File:** `matlab/fit_q1_affera_model.m`

Current behavior:
- Uses `isBaselineEra` and `isPostNonPFAEra` to compute:
  - `baseline_stats` (baseline non‑PFA era).
  - `post_non_pfa_stats` (post‑Affera non‑PFA).

Planned changes:
- Internal logic:
  - Continue to compute baseline stats using `isBaselineEra`.
  - Replace dependence on `isPostNonPFAEra` with the more general `isComparisonEra` (or keep `isPostNonPFAEra` as alias if Option 2 in §3.3 is chosen).
- Naming:
  - Optionally rename `post_non_pfa_stats` → `post_comparison_stats` in the `results` struct for clarity.
  - If keeping the old field name for backward compatibility, clearly document that it now represents the active comparison group (non‑PFA or RF).

The key requirement is that the **logic** follows the centralized comparison‑group mask; labels will be harmonized in `print_q1_affera_summary`.

---

## 4. Command‑Line Testing and Parity Checks

This section summarizes tests to verify the refactor and ensure parity with the pre‑refactor behavior when `ComparisonGroup = 'nonPFA'`.

### 4.1 Minimal synthetic test (sanity check)

Create a small synthetic dataset in MATLAB with:
- A handful of procedures spanning:
  - Baseline non‑PFA RF cases.
  - Baseline non‑PFA non‑RF cases.
  - Post‑Affera RF and non‑RF non‑PFA cases.
  - Affera cases.
- Explicit settings of:
  - `is_affera`, `is_pfa_catheter`, `is_rf_catheter`.
  - `procedure_date`, `procedure_duration_min`, `procedure_type`, `operator_id`.

Run from the project root:

```bash
matlab -batch "addpath('matlab'); \
tbl = your_synthetic_table_constructor(); \
[lme1, results1, tbl_q1_1] = run_q1_affera_model(tbl, ...
    'BaselineDays', 180, ...
    'MinCasesPerGroup', 0, ...
    'ComparisonGroup', 'nonPFA'); \
[lme2, results2, tbl_q1_2] = run_q1_affera_model(tbl, ...
    'BaselineDays', 180, ...
    'MinCasesPerGroup', 0, ...
    'ComparisonGroup', 'RF'); \
disp(results1.operator_counts); \
disp(results2.operator_counts);"
```

Checks:
- `tbl_q1_1` includes all non‑PFA non‑Affera cases, both RF and non‑RF.
- `tbl_q1_2` includes only RF non‑Affera cases.
- Affera counts and operator IDs included are consistent between the two runs, except for operators whose baseline RF counts drop below any `MinCasesPerGroup` threshold when using RF‑only.

### 4.2 Parity test vs. pre‑refactor behavior (non‑PFA)

Before refactoring (or using a saved pre‑refactor version), run:

```bash
matlab -batch "addpath('matlab'); \
load('data_processed/procedures_all.mat', 'tbl'); \
[lme_old, results_old, tbl_q1_old] = run_q1_affera_model(tbl);"
```

After refactoring, run with explicit non‑PFA comparison:

```bash
matlab -batch "addpath('matlab'); \
load('data_processed/procedures_all.mat', 'tbl'); \
[lme_new, results_new, tbl_q1_new] = run_q1_affera_model(tbl, ...
    'ComparisonGroup', 'nonPFA');"
```

Parity checks:
- **Row counts**:
  - `height(tbl_q1_new)` should exactly match `height(tbl_q1_old)`.
  - Counts by `isBaselineEra`, `isAffera`, and (if present) `isPostNonPFAEra` should match.
- **Operator inclusion**:
  - `results_old.operator_counts` vs `results_new.operator_counts` should have identical `n_baseline`, `n_affera`, and `included_in_model`.
- **Fixed effects**:
  - `lme_old.Coefficients` and `lme_new.Coefficients` should match numerically (up to floating‑point tolerance).
  - `results_old.pct_est`/`pct_lo`/`pct_hi` vs `results_new` should be identical.

Any differences here indicate a bug in the refactor’s non‑PFA logic and should be investigated before proceeding.

### 4.3 Real‑data RF comparison run

Once parity is confirmed for non‑PFA:

```bash
matlab -batch "addpath('matlab'); \
load('data_processed/procedures_all.mat', 'tbl'); \
[lme_rf, results_rf, tbl_q1_rf] = run_q1_affera_model(tbl, ...
    'ComparisonGroup', 'RF');"
```

Checks:
- Summary output should clearly indicate:
  - `Comparison group: RF` in the settings section.
  - Labels such as “Post‑Affera comparison group (RF only)” where appropriate.
- Baseline and post‑comparison group counts should:
  - Be ≤ the non‑PFA counts.
  - Match manual RF counts from `tbl.is_rf_catheter` within the defined time windows.
- Operator inclusion:
  - Operators with few RF baseline cases may drop out due to `MinCasesPerGroup`, as expected.

### 4.4 Optional: quick internal consistency checks

Optionally, add a small assertion script (not required for core functionality) that:
- Confirms `is_rf_catheter => ~is_pfa_catheter` (RF should not be flagged as PFA).
- Checks that:
  - With `ComparisonGroup = 'RF'`, every non‑Affera baseline/comparison row has `is_rf_catheter == true`.
  - With `ComparisonGroup = 'nonPFA'`, every non‑Affera baseline/comparison row has `is_pfa_catheter == false`.

These checks can be part of a simple diagnostic script or run interactively during analysis.

---

## 5. Summary

This refactor:
- Adds a clear `ComparisonGroup` API to `run_q1_affera_model`.
- Centralizes comparison‑group definition in `get_q1_comparison_group_mask`.
- Introduces an explicit `is_rf_catheter` flag in the imported data.
- Leaves the mixed‑effects model structure unchanged while allowing analysts to choose between:
  - A broad non‑PFA comparison (current default).
  - A more focused RF‑only comparison.
- Includes concrete command‑line tests to verify:
  - Parity with existing non‑PFA behavior.
  - Correct behavior of the RF‑only option.

Once implemented and tested, this will make it much easier to align the analysis with manuscript language (e.g., “Affera vs RF” vs “Affera vs non‑PFA”) without touching the core modeling code.

