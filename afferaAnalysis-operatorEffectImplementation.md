# Affera Operator Baseline Efficiency Effect — Implementation Plan

This document describes a concrete, modular implementation plan for extending the existing Affera Q1 mixed-effects model to evaluate whether **baseline operator speed** modifies the **efficiency benefit of Affera**. It assumes familiarity with the core model and learning-curve extension as specified in `afferaAnalysis-modelSpecification.md` (especially Sections 4, 10, and 11).

The goal is to add a small, optional extension to the existing MATLAB pipeline that:

- Derives a baseline speed metric for each operator from **non-PFA baseline cases**.
- Adds this metric and its interaction with Affera use as fixed effects in the model.
- Keeps all modeling on the **log-duration** scale internally, while reporting effects as **percent change in duration**.
- Is enabled via a simple flag (e.g. `'EnableOperatorEffect', true`) without disrupting existing workflows.

---

## 1. High-Level Behavior and User API

The operator-baseline effect should be an **optional extension**, similar to the existing learning-curve flag.

### 1.1 Entry point: `run_q1_affera_model`

`run_q1_affera_model` already supports name-value options (e.g. `'EnableLearningCurve', true`). We will:

- Add a new boolean option:

  - `'EnableOperatorEffect'` (default: `false`)

- Typical calls:

  ```matlab
  % Core Q1 model only
  [lme, results, tbl_q1] = run_q1_affera_model(clinicalData);

  % Q1 + learning curve
  [lme, results, tbl_q1] = run_q1_affera_model( ...
      clinicalData, 'EnableLearningCurve', true);

  % Q1 + operator baseline effect
  [lme, results, tbl_q1] = run_q1_affera_model( ...
      clinicalData, 'EnableOperatorEffect', true);

  % Q1 + learning curve + operator baseline effect
  [lme, results, tbl_q1] = run_q1_affera_model( ...
      clinicalData, 'EnableLearningCurve', true, ...
                   'EnableOperatorEffect', true);
  ```

### 1.2 Behavioral guarantees

- When `'EnableOperatorEffect'` is `false` (default), behavior of the Q1 and/or learning-curve model is unchanged.
- When `'EnableOperatorEffect'` is `true`:
  - The same dataset `tbl_q1` is used, but augmented with per-operator baseline speed variables.
  - The fixed-effects formula in `fit_q1_affera_model` gains the operator-speed main effect and its interaction with Affera.
  - The summary (`print_q1_affera_summary`) reports these effects in **percent change** terms, with log-scale values still stored internally in `results`.

---

## 2. Baseline Operator Speed Metric

### 2.1 Conceptual definition

For each operator, we define a **baseline speed** metric based on that operator’s **baseline non-PFA** cases (pre-Affera window):

- Let `B_op` be the set of baseline procedures for operator `op` (`isBaselineEra == 1` and `operator_id == op` in `tbl_q1`).
- Using `log_duration` from those cases, define:

  ```matlab
  baseline_speed_raw(op) = mean(tbl_q1.log_duration(B_op), 'omitnan');
  ```

  This is the operator’s mean log-duration in the baseline non-PFA era.

- To aid interpretation and reduce collinearity, we then **center** this metric across operators:

  ```matlab
  mu_baseline = mean(baseline_speed_raw, 'omitnan');
  baseline_speed_ctr(op) = baseline_speed_raw(op) - mu_baseline;
  ```

### 2.2 Attaching to the analysis table

For all rows of `tbl_q1` belonging to operator `op` (baseline, post-Affera non-PFA, and Affera cases), we assign:

- `tbl_q1.baseline_speed_raw(row) = baseline_speed_raw(op)`
- `tbl_q1.baseline_speed_ctr(row) = baseline_speed_ctr(op)`

Operators with no baseline cases (which should not occur if the Q1 inclusion rules are enforced) can be guarded against by:

- Checking for empty `B_op`; if empty, either:
  - Skipping operator-baseline effect for that operator (not recommended), or
  - Falling back to the global mean (preferred if needed).
  - In practice, the existing `min_cases_per_group` logic should ensure every included operator has some baseline cases.

---

## 3. Helper Function: `add_operator_baseline_speed.m`

To keep the implementation DRY and modular, introduce a dedicated helper:

```matlab
function tbl_q1 = add_operator_baseline_speed(tbl_q1)
```

### 3.1 Responsibilities

1. **Preconditions**
   - Verify that `tbl_q1` contains:
     - `operator_id` (categorical)
     - `log_duration`
     - `isBaselineEra` (logical or numeric 0/1)
   - If any are missing, throw a descriptive error.

2. **Compute baseline metrics per operator**

   - Get operator categories:

     ```matlab
     ops = categories(tbl_q1.operator_id);
     ```

   - Initialize arrays `baseline_speed_raw` and `baseline_speed_ctr` indexed by operator.
   - For each operator `op`:
     - Identify baseline rows:

       ```matlab
       maskBaselineOp = (tbl_q1.operator_id == op) & tbl_q1.isBaselineEra;
       ```

     - Compute `baseline_speed_raw(op)` as mean `log_duration` over `maskBaselineOp`.

   - Compute global mean across operators and create centered values:

     ```matlab
     mu_baseline = mean(baseline_speed_raw, 'omitnan');
     baseline_speed_ctr = baseline_speed_raw - mu_baseline;
     ```

3. **Attach to `tbl_q1`**

   - Add new columns (one row per procedure):

     ```matlab
     n = height(tbl_q1);
     tbl_q1.baseline_speed_raw = nan(n, 1);
     tbl_q1.baseline_speed_ctr = nan(n, 1);
     ```

   - For each operator `op`, assign the appropriate value to all its rows.

4. **Return `tbl_q1`**

   - The table returned is identical to the input except for two new operator-level covariates:
     - `baseline_speed_raw`
     - `baseline_speed_ctr`

### 3.2 Call site

In `run_q1_affera_model`, after building `tbl_q1` (and after optionally adding `affera_index_ctr` when learning is enabled):

```matlab
if opts.EnableOperatorEffect
    tbl_q1 = add_operator_baseline_speed(tbl_q1);
end
```

This keeps the new behavior isolated and simple.

---

## 4. Model Formula Changes in `fit_q1_affera_model`

`fit_q1_affera_model` currently:

- Builds a base formula (Q1 model):

  ```matlab
  baseFormula = 'log_duration ~ isAffera * isPVIplus + time_days';
  ```

- Optionally adds the learning-curve term:

  ```matlab
  if opts.EnableLearningCurve
      baseFormula = [baseFormula, ' + affera_index_ctr'];
  end
  ```

### 4.1 Add operator baseline effect (main + interaction)

Extend the formula construction to recognize `opts.EnableOperatorEffect`:

1. **Main effect of baseline speed**

   - When `EnableOperatorEffect` is true, ensure that `baseline_speed_ctr` exists in `tbl_q1`.
   - Append to the fixed-effects formula:

     ```matlab
     if opts.EnableOperatorEffect
         baseFormula = [baseFormula, ' + baseline_speed_ctr'];
     end
     ```

2. **Interaction with Affera use**

   - Add the interaction term `isAffera:baseline_speed_ctr`:

     ```matlab
     if opts.EnableOperatorEffect
         baseFormula = [baseFormula, ' + isAffera:baseline_speed_ctr'];
     end
     ```

3. **Random effects**

   - Keep the random-effects structure unchanged:

     ```matlab
     formula = sprintf('%s + (1|operator_id)', baseFormula);
     ```

4. **Result extraction**

   - After fitting the model and extracting `beta`, `ci`, and `names`, identify:
     - `idxBaselineSpeed = find(strcmp(names, 'baseline_speed_ctr'));`
     - `idxAfferaBaseline =` index where `names` starts with `'isAffera:baseline_speed_ctr'`.
   - Add these indices to the `results` struct:

     ```matlab
     results.idxBaselineSpeed   = idxBaselineSpeed;
     results.idxAfferaBaseline  = idxAfferaBaseline;
     ```

   - Percent-scale values and p-values will already be populated via the existing transformation and coefficient extraction logic:
     - `results.pct_est`, `pct_lo`, `pct_hi`, and `pValue`.

Internally, the operator-speed terms are still estimated on the log-duration scale. Only the console reporting uses percent change.

---

## 5. Summary Output Changes in `print_q1_affera_summary`

`print_q1_affera_summary` already:

- Prints a fixed-effects table in percent scale.
- Highlights:
  - Affera effect in PVI-only (`isAffera`).
  - Additional Affera effect in PVI+ (`isAffera:isPVIplus`).
  - Overall Affera effect (linear contrast).
  - Learning-curve term (if present).

### 5.1 Add operator baseline effect reporting

When `idxBaselineSpeed` and/or `idxAfferaBaseline` exist in `results`, add a new section:

1. **Baseline operator speed main effect**

   - Identify index `i = results.idxBaselineSpeed`.
   - Print:

     ```text
     Baseline operator efficiency (non-PFA era):
       Percent change in duration per unit increase in baseline speed = XX.X% [L, U]%, p = ...
     ```

   - “Per unit” can be clarified in an inline note:
     - e.g., “per 1-unit increase in centered mean log-duration.”
   - Optionally, compute approximate percent change between “slow” and “fast” operators (e.g., 25th vs 75th percentile of `baseline_speed_ctr`), but this can be deferred or added later.

2. **Interaction: Affera × baseline speed**

   - Identify index `j = results.idxAfferaBaseline`.
   - Print:

     ```text
     Effect modification by baseline operator efficiency (Affera × baseline_speed):
       Additional percent change in duration per unit higher baseline speed = YY.Y% [L, U]%, p = ...
     ```

   - Interpretation example (inline sentence):
     - If the interaction is negative (on percent scale), note that operators with slower baseline performance (higher baseline duration) experience a **larger percent reduction** in procedure duration when using Affera, and vice versa.

All values reported here are taken from `results.pct_est`, `pct_lo`, `pct_hi`, and `pValue`, consistent with the rest of the summary.

---

## 6. Diagnostics and Saved Outputs

### 6.1 Collinearity diagnostics

The existing `compute_learning_curve_diagnostics` function currently focuses on:

- `{'isAffera', 'isPVIplus', 'time_days', 'affera_index_ctr'}`.

For the operator baseline effect, we can:

- Either extend that helper, or introduce a small new helper, to also include `baseline_speed_ctr` when `EnableOperatorEffect` is true.
- Minimal extension (DRY approach):
  - Reuse the same diagnostics machinery, simply adding `baseline_speed_ctr` to the predictors list when present.
  - Emit additional warnings if high correlations or elevated VIFs involve the operator-speed terms.

### 6.2 Saved results

No changes to saving behavior are required beyond the expanded `results` struct:

- `run_q1_affera_model` already saves `lme`, `results`, and `tbl_q1` to `results/q1_affera_model_results.mat`.
- With this extension, `results` will additionally contain:
  - `idxBaselineSpeed`, `idxAfferaBaseline`.
  - Percent-scale effects and p-values for these coefficients.
  - Log-scale `beta` and `ci` values for possible use in technical appendices.

---

## 7. Summary of Changes (DRY & KISS)

To implement the operator baseline efficiency effect:

1. **Add one helper**:
   - `add_operator_baseline_speed.m` to compute per-operator baseline speed from non-PFA baseline cases and attach `baseline_speed_raw` and `baseline_speed_ctr` to `tbl_q1`.

2. **Extend the entry point**:
   - In `run_q1_affera_model`, add a new name-value option `'EnableOperatorEffect'` and, when true, call `add_operator_baseline_speed` after dataset construction (and after the Affera learning index, if enabled).

3. **Extend model formula construction**:
   - In `fit_q1_affera_model`, conditionally add `+ baseline_speed_ctr + isAffera:baseline_speed_ctr` to the fixed-effects part of the formula when `EnableOperatorEffect` is true.
   - Track indices of these new coefficients in `results`.

4. **Extend summary reporting**:
   - In `print_q1_affera_summary`, detect `idxBaselineSpeed` and `idxAfferaBaseline` and print percent-scale interpretations and p-values for these terms.

5. **Diagnostics**:
   - Optionally extend the existing diagnostics helper to include `baseline_speed_ctr` when present.

This approach:

- Keeps all modeling on the log-duration scale internally.
- Reports **all user-facing effects** (including operator baseline and its interaction with Affera) as **percent changes in duration** with 95% CIs and p-values.
- Adds minimal, well-localized additions to the existing Q1 + learning-curve pipeline, preserving DRY and KISS principles.

---

## 8. Command-Line MATLAB Testing

After implementing the operator baseline efficiency extension, basic verification should be performed using **command-line MATLAB**. The goal is to confirm:

- The new options parse correctly (`EnableOperatorEffect`, with and without `EnableLearningCurve`).
- `tbl_q1` contains the expected new variables (`baseline_speed_raw`, `baseline_speed_ctr`).
- The fitted model includes the new fixed effects and interaction.
- The `results` struct exposes the indices and percent-scale summaries.
- Console output reports operator baseline effects in percent, with no log-scale values shown.

### 8.1 Minimal synthetic test

Create a small synthetic dataset directly in a `matlab -batch` call to verify wiring without depending on external files:

```bash
matlab -batch "addpath('matlab'); \
dates = (datetime(2025,1,1) + days([0 10 20 30 40 50])).'; \
is_affera = [0;0;0;1;1;1]; \
is_pfa_catheter = [0;0;0;1;1;1]; \
procedure_duration_min = [130;120;110;100;95;90]; \
procedure_type = categorical({'PVI';'PVI+';'PVI';'PVI';'PVI+';'PVI'}); \
operator_id = categorical({'A';'A';'B';'A';'B';'B'}); \
log_duration = log(procedure_duration_min); \
tbl = table((1:6).', dates, operator_id, procedure_type, is_affera, is_pfa_catheter, ...
    procedure_duration_min, log_duration, ...
    'VariableNames', {'procedure_id','procedure_date','operator_id','procedure_type', ...
                      'is_affera','is_pfa_catheter','procedure_duration_min','log_duration'}); \
[lme, results, tbl_q1] = run_q1_affera_model(tbl, ...
    'BaselineDays', 180, 'MinCasesPerGroup', 0, ...
    'EnableOperatorEffect', true); \
disp(results.idxBaselineSpeed); disp(results.idxAfferaBaseline); \
disp(tbl_q1(:, {'operator_id','isBaselineEra','baseline_speed_raw','baseline_speed_ctr'}));"
```

Sanity checks from this run:

- `tbl_q1` shows non-NaN `baseline_speed_raw` and `baseline_speed_ctr` values per operator.
- `results.idxBaselineSpeed` and `results.idxAfferaBaseline` are non-empty and correspond to coefficients named `baseline_speed_ctr` and `isAffera:baseline_speed_ctr` (or similar).
- The console summary includes a “Baseline operator efficiency” section and an “Effect modification by baseline operator efficiency” section, with percent-scale effects and p-values.

### 8.2 Test on real imported data

Assuming `clinicalData/importedClinicalData.mat` contains a `clinicalData` table produced by `import_clinical_data.m`, run:

```bash
matlab -batch "addpath('matlab'); \
load('clinicalData/importedClinicalData.mat','clinicalData'); \
[lme_base, results_base, tbl_q1_base] = run_q1_affera_model(clinicalData); \
[lme_op, results_op, tbl_q1_op] = run_q1_affera_model(clinicalData, ...
    'EnableOperatorEffect', true); \
[lme_full, results_full, tbl_q1_full] = run_q1_affera_model(clinicalData, ...
    'EnableLearningCurve', true, 'EnableOperatorEffect', true); \
disp(results_op.idxBaselineSpeed); disp(results_op.idxAfferaBaseline);"
```

Check that:

- `tbl_q1_op` and `tbl_q1_full` include `baseline_speed_raw` and `baseline_speed_ctr`.
- The fixed-effect names in `lme_op` / `lme_full` include the new terms.
- The console summaries show the operator-baseline sections in percent scale only, with no log-scale coefficients printed.
- The saved MAT file `results/q1_affera_model_results.mat` (from the last run) contains the expanded `results` struct with the new indices and percent-scale fields.

These command-line tests provide a quick, reproducible way to confirm that the operator baseline efficiency extension is wired correctly, behaves as an optional feature, and respects the log-scale modeling / percent-scale reporting conventions described in the main specification.

