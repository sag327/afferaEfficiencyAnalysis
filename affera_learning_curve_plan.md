
# Affera Learning Curve Analysis Plan

## Overview

This document outlines the rationale, objectives, and and methodological plan for evaluating learning curves associated with adoption of the Affera pulsed-field ablation (PFA) system. The goal is to determine whether procedural efficiency improves with increasing Affera experience and whether such learning effects differ between pulmonary vein isolation (PVI) and PVI-plus (PVI+) procedures.

This specification is intended for coding assistants (e.g., Codex) to guide implementation in MATLAB. The goal is to generate a reproducible workflow that quantifies Affera-specific learning while controlling for operator variation, procedural complexity, and secular trends in laboratory performance.

---

## Rationale

The currently implemented mixed-effects model estimates the **average** effect of Affera compared with non–PFA ablation. However, it does **not** distinguish:

- the **instantaneous** effect of adopting Affera (e.g., first case vs. RF), from  
- the **experience-related** improvement as operators gain proficiency.

Without modeling Affera-specific learning curves, early Affera cases and late Affera cases are treated as equivalent. This may obscure meaningful changes in efficiency and potentially inflate or shrink the estimated Affera effect.

To address this, the analysis should incorporate an **Affera case index**, representing the sequence number of Affera cases performed by each operator (1, 2, 3, …). This index enables explicit measurement of whether (and how quickly) performance improves with Affera use.

---

## Methodological Challenges: Collinearity

A key risk in adding an Affera case index is **collinearity with calendar time**. Because operators typically begin using Affera around the same date and continue over time, Affera case number often increases in parallel with calendar time. This may lead to:

- inflated standard errors  
- unstable coefficient estimates  
- difficulty separating **technology learning** from **overall laboratory efficiency trends**

This problem was mitigated in the Q1 model by including non–PFA cases after Affera introduction, creating temporal overlap. However, further protection is needed.

---

## Plan to Detect and Mitigate Collinearity

The script should include explicit tests for collinearity between:

- `affera_case_index`
- `time_days`
- `isAffera`
- `isPVIplus`
- interaction terms

### 1. Correlation Diagnostics

Compute Pearson and Spearman correlations:

```matlab
corrMatrix = corr(table{:, predictors}, 'Rows', 'complete');
```

Any absolute correlation above **0.7** should trigger a warning.

### 2. Variance Inflation Factors (VIF)

Compute VIF for fixed-effect predictors:

```matlab
vifValues = computeVIF(designMatrix);
```

Warnings should be triggered at:

- VIF > 5: moderate concern  
- VIF > 10: serious concern  

### 3. Centered Affera Case Index

Compute an operator-specific centered index:

```matlab
tbl.affera_index_centered = tbl.affera_index - mean(tbl.affera_index(tbl.isAffera==1));
```

Centering reduces correlation with both the Affera main effect and calendar time.

### 4. Optional Orthogonalization

Residualize Affera index against time:

```matlab
mdl_temp = fitlm(tbl.time_days, tbl.affera_index);
tbl.affera_index_orth = mdl_temp.Residuals.Raw;
```

This yields an Affera index that is mathematically independent of calendar time, ensuring interpretability of both effects.

---

## Proposed Learning Curve Model

\[
\log(\text{duration}) =
\beta_0 +
\beta_1(\text{Affera}) +
\beta_2(\text{PVI+}) +
\beta_3(\text{time}) +
\beta_4(\text{Affera} \times \text{PVI+}) +
\beta_5(\text{Affera index (centered/orthogonalized)}) +
(1|\text{operator})
\]

Optional random slopes may be included:

```matlab
'(1 + affera_index_centered | operator_id)'
```

---

## Outputs and Interpretation

The script should report:

- Effect of Affera use (holding experience at mean level)  
- Learning curve slope (change in log-duration per Affera case)  
- Predicted duration improvement between early vs. late Affera experience  
- Whether PVI+ cases exhibit different learning dynamics  
- Diagnostic outputs: correlation matrix, VIF table, and warnings  

Graphical outputs may include:

- Operator-specific or pooled learning curves  
- Marginal effects plots  
- Forest plots of fixed effects  

---

## Implementation Plan (MATLAB, DRY/KISS)

The learning-curve implementation should be a **minimal extension** of the existing Q1 pipeline. After the Affera case index is added, the model simply gains **one extra fixed-effect term**; the core workflow, outputs, and entry point remain the same.

### 1. High-level structure

1. Keep the existing Q1 code path intact:
   - `matlab/import_clinical_data.m`
   - `matlab/make_q1_affera_dataset.m`
   - `matlab/fit_q1_affera_model.m`
   - `matlab/print_q1_affera_summary.m`
   - `matlab/run_q1_affera_model.m`
2. Add **one new helper** to compute Affera case indices:
   - `matlab/add_affera_learning_curve_variables.m`
3. Optionally add **one small diagnostics helper**:
   - `matlab/compute_learning_curve_diagnostics.m`
4. Lightly extend `run_q1_affera_model.m` and `fit_q1_affera_model.m` to:
   - Accept an optional `opts` struct.
   - When `opts.enableLearningCurve == true`, call the helper to add the index and include it as an additional fixed effect.

No separate “learning-curve model” function or parallel result object is needed; we still fit **one** LME, just with an extra predictor.

### 2. Extend `run_q1_affera_model.m` (optional learning-curve flag)

Update the function signature to accept an optional `opts` struct while keeping the same outputs:

```matlab
function [lme, results, tbl_q1] = run_q1_affera_model( ...
    clinicalData, baseline_days, min_cases_per_group, opts)
```

Inputs (existing + new):
- `clinicalData`: table from `import_clinical_data` (unchanged).
- `baseline_days` (optional): default `180` (unchanged).
- `min_cases_per_group` (optional): default `15` (unchanged).
- `opts` (optional struct), with key field:
  - `opts.enableLearningCurve` (logical, default `false`).

Defaults:
- If `opts` is omitted or empty, internally set `opts.enableLearningCurve = false`.

Internal steps (pseudocode):
1. Construct `meta` as in the Q1 spec:
   ```matlab
   meta.affera_launch_date_global = min(clinicalData.procedure_date(clinicalData.is_affera));
   ```
2. Build the core Q1 dataset (unchanged):
   ```matlab
   tbl_q1 = make_q1_affera_dataset(clinicalData, meta, baseline_days, min_cases_per_group);
   ```
3. If `opts.enableLearningCurve` is `true`, add the Affera index variables and run diagnostics:
   ```matlab
   if opts.enableLearningCurve
       tbl_q1 = add_affera_learning_curve_variables(tbl_q1);
       diagnostics = compute_learning_curve_diagnostics(tbl_q1);
   end
   ```
4. Fit the mixed-effects model, passing `opts` through so the formula can include the index term when enabled:
   ```matlab
   if nargin < 4
       opts = struct();
   end
   [lme, results] = fit_q1_affera_model(tbl_q1, opts);
   print_q1_affera_summary(lme, results, tbl_q1, opts);
   ```

Behavior:
- When `opts.enableLearningCurve == false` (default), this reproduces the existing Q1 model exactly.
- When `opts.enableLearningCurve == true`, the same function builds the same dataset but augments it with an Affera index and calls a slightly richer model.

### 3. Helper: `add_affera_learning_curve_variables.m`

Create:

```matlab
function tbl_q1 = add_affera_learning_curve_variables(tbl_q1)
```

Responsibilities:
1. Ensure required variables exist: `operator_id`, `procedure_date`, `isAffera`, `time_days`, `log_duration`, `isPVIplus`.
2. Initialize new columns:
   ```matlab
   n = height(tbl_q1);
   tbl_q1.affera_index     = nan(n, 1);
   tbl_q1.affera_index_ctr = nan(n, 1);
   ```
3. Compute Affera case index per operator:
   - Restrict to Affera rows:
     ```matlab
     isA = tbl_q1.isAffera == 1;
     ```
   - For each `operator_id` with any Affera cases:
     - Subset rows for that operator where `isA` is true.
     - Sort those rows by `procedure_date` (and `procedure_id` as a tiebreaker if needed).
     - Assign indices `1, 2, 3, ...` in chronological order to `tbl_q1.affera_index` for those rows.
   - Non-Affera rows remain `0` in `affera_index`.

   Implementation hint (simple and clear):
   - Use `findgroups` on `operator_id(isA)` and loop over groups, or a straightforward loop over `categories(tbl_q1.operator_id)`.

4. Center the Affera index:
   - Compute mean index among Affera cases:
     ```matlab
     mean_idx = mean(tbl_q1.affera_index(isA), 'omitnan');
     ```
   - Center:
     ```matlab
     tbl_q1.affera_index_ctr(isA) = tbl_q1.affera_index(isA) - mean_idx;
     ```
   - Non-Affera rows are set to `0` so non-Affera rows remain in the model.

This helper is the **only** place where the learning-curve index is constructed; it is small, testable, and independent of the modeling step.

### 4. Optional helper: `compute_learning_curve_diagnostics.m`

Create:

```matlab
function diagnostics = compute_learning_curve_diagnostics(tbl_q1)
```

Responsibilities (kept simple):
1. Assemble predictors to check for collinearity:
   ```matlab
   predictors = {'isAffera', 'isPVIplus', 'time_days', 'affera_index_ctr'};
   X = tbl_q1{:, predictors};
   rows_ok = all(isfinite(X), 2);
   X_clean = X(rows_ok, :);
   ```
2. Correlation matrix:
   ```matlab
   corrMatrix = corr(X_clean, 'Rows', 'complete');
   ```
   - Compute absolute correlations and flag any `> 0.7`.
3. Simple VIF calculation:
   - For each predictor column `j`:
     - Regress that column on the others and compute `R^2`.
     - `VIF_j = 1 / (1 - R^2)`.
4. Package into a `diagnostics` struct:
   ```matlab
   diagnostics.predictors = predictors;
   diagnostics.corrMatrix = corrMatrix;
   diagnostics.vif        = vifValues;  % if implemented
   ```
5. Print brief warnings to the command window if:
   - Any pairwise `|corr| > 0.7`.
   - Any `VIF > 5` or `VIF > 10` (if VIF is computed).

This helper is called automatically whenever the learning-curve term is enabled; it does not change the model, only reports potential collinearity.

### 5. Extend `fit_q1_affera_model.m` (add one term)

Update the function signature to optionally accept `opts`:

```matlab
function [lme, results] = fit_q1_affera_model(tbl_q1, opts)
```

Defaults:
- If `opts` is not provided, treat `opts.enableLearningCurve` as `false`.

Responsibilities:
1. Build the base fixed-effects formula (current Q1 spec):
   ```matlab
   baseFormula = 'log_duration ~ isAffera * isPVIplus + time_days';
   ```
2. If `opts.enableLearningCurve == true`, append the centered Affera index:
   ```matlab
   if isfield(opts, 'enableLearningCurve') && opts.enableLearningCurve
       baseFormula = [baseFormula, ' + affera_index_ctr'];
   end
   ```
3. Random-effects term remains unchanged:
   ```matlab
   reTerm = '(1|operator_id)';
   formula = sprintf('%s + %s', baseFormula, reTerm);
   ```
4. Fit the model:
   ```matlab
   lme = fitlme(tbl_q1, formula);
   ```
5. Extract fixed effects and CIs as already done in the Q1 spec:
   ```matlab
   beta    = fixedEffects(lme);
   ci      = coefCI(lme);
   names   = lme.CoefficientNames;

   pct_est = (exp(beta)     - 1) * 100;
   pct_lo  = (exp(ci(:, 1)) - 1) * 100;
   pct_hi  = (exp(ci(:, 2)) - 1) * 100;
   ```
6. Package into the existing `results` struct/table format (unchanged), with one additional row where `names == 'affera_index_ctr'` when learning is enabled.

Thus the “learning-curve model” is literally the original Q1 model with `+ affera_index_ctr` added to the fixed effects.

### 6. Lightly extend `print_q1_affera_summary.m`

Without changing its overall behavior, extend the summary function so that, **if** a coefficient named `affera_index_ctr` is present, it:
- Prints the estimate, 95% CI, and percent change per one-unit increase in Affera index.
- Optionally translates this to a more intuitive comparison, e.g.:
  - “Expected percent change in duration between Affera case #1 vs #20” (by multiplying the slope by 19 and exponentiating).

If the coefficient is absent (learning curve disabled), the summary remains exactly as currently specified.

### 7. Minimal changes to existing code

To respect DRY and KISS:
- Do **not** modify `import_clinical_data.m` or `make_q1_affera_dataset.m`.
- Keep the **same entry point** and outputs:
  - `run_q1_affera_model` still returns `[lme, results, tbl_q1]`.
- Add only:
  - `add_affera_learning_curve_variables.m`
  - (Optionally) `compute_learning_curve_diagnostics.m`
- Make small, backwards-compatible extensions:
  - `run_q1_affera_model.m` accepts `opts` and, when `opts.enableLearningCurve` is true, calls the new helper before fitting.
  - `fit_q1_affera_model.m` conditionally adds `+ affera_index_ctr` to the fixed-effects formula.
  - `print_q1_affera_summary.m` detects and reports the `affera_index_ctr` coefficient if present.

With this structure, the learning-curve analysis:
- Reuses the existing Q1 pipeline end-to-end.
- Adds only one new predictor to the model when requested.
- Keeps the public interface simple (`opts.enableLearningCurve` on/off).
- Avoids duplicate models, extra result types, or unnecessary complexity.

---

## Summary

This analysis plan allows us to:

- Separate instant Affera benefit from operator learning  
- Control for secular efficiency gains  
- Quantify how experience shapes performance  
- Detect and address collinearity risks  

This document should be provided to Codex to guide implementation of a robust Affera learning curve model.
