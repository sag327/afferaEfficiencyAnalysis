# Q1 Model Specification — Affera vs Non-Affera Efficiency Analysis
This document defines the exact data requirements, preprocessing, and mixed-effects model structure for the new Q1 analysis. It is written as a technical specification for an automated coding assistant (e.g., Codex) to generate MATLAB code.

---

## 1. Objective
Estimate the effect of **Affera** catheter use on **procedure duration** for AF ablation, using all applicable operators and all relevant AF ablation procedures.

The model should:
- Include **all operators**
- Include **all non-PFA AF ablation procedures** from the **6 months prior** to the introduction of Affera
- Include **all Affera cases**
- Allow separate estimation of Affera’s effect in **PVI-only** vs **PVI+** (complex PVI) cases
- Adjust for **calendar time**, to control for secular improvements in lab efficiency
- Account for operator-to-operator differences via **random effects**

---

## 2. Dataset Requirements

### Included procedure types:
- Any AF ablation procedure that is either:
  - **Affera** (PFA), performed at/after the Affera start date
  - **Non-PFA** (e.g., RF) AF ablation within the **6 months prior to Affera introduction**
  - **Non-PFA** (e.g., RF) AF ablation performed **on or after** Affera introduction (post-Affera non-PFA cases)

### Variables required in the table:
The MATLAB table `tbl` must include:

| column name | description |
|-------------|-------------|
| `duration_minutes` | Total procedure duration in minutes |
| `isAffera` | 1 = Affera case, 0 = non-PFA case |
| `isPVIplus` | 1 = PVI+ (complex PVI), 0 = PVI-only |
| `operator_id` | Categorical operator identifier |
| `procedure_date` | Datetime of procedure |
| `time_days` | Numeric days since reference date (computed from procedure_date) |

### Optional but useful:
- `case_id`
- `procedure_type` (string or categorical)
- `isBaselineEra` (0/1)

Codex should compute `time_days` as:

```matlab
tbl.time_days = days(tbl.procedure_date - min(tbl.procedure_date));
```

---

## 3. Data Inclusion Logic

1. **Affera cohort**
   - Include all procedures where `isAffera == 1`, performed at or after `Affera_start`.

2. **Non-Affera baseline cohort (pre-Affera, non-PFA)**
   - Include procedures where:
     - `isAffera == 0`
     - Procedure is a non-PFA AF ablation
     - `procedure_date` is **between (Affera_start - 180 days)** and **Affera_start**

3. **Non-Affera post-Affera cohort (post-Affera, non-PFA) — v2 update**
   - Include procedures where:
     - `isAffera == 0`
     - Procedure is a non-PFA AF ablation
     - `procedure_date` is **on or after `Affera_start`** (through the end of the observation period)

4. **Operator inclusion**
   - Conceptually, include *all operators* who appear in any of the above cohorts (Affera, baseline non-PFA, or post-Affera non-PFA).
   - No minimum case count threshold is required at the model-specification level (implementation may optionally enforce thresholds as described in Section 7).

5. **PVI vs PVI+**
   - Coded via `isPVIplus`, used as both a main effect and an interaction term.

---

## 4. Model Specification

### Outcome:
Log-transformed duration:

```matlab
tbl.log_duration = log(tbl.duration_minutes);
```

### Mixed-effects model formula:

log(duration) =
beta0 +
beta1 * isAffera +
beta2 * isPVIplus +
beta3 * time_days +
beta4 * (isAffera * isPVIplus) +
(1 | operator_id)

### MATLAB syntax:

```matlab
lme = fitlme(tbl, ...
 'log_duration ~ isAffera * isPVIplus + time_days + (1|operator_id)');
```

### Interpretation:
- beta1 (isAffera): Effect of Affera vs non-Affera for **PVI-only** (reference group).
- beta4 (isAffera:isPVIplus): Additional effect of Affera in **PVI+** compared to **PVI**.
- beta2 (isPVIplus): Baseline difference between PVI and PVI+ in the non-Affera era.
- beta3 (time_days): Secular efficiency trend independent of Affera.
- Random intercept (operator_id): Each operator has their own baseline duration; Affera effects are pooled across operators.

Internally, all fixed effects are estimated on the **log-duration scale**, but all **reported** effects in the console summary are expressed as **percent change in duration** (with 95% CIs and p-values) via the transformation:

```matlab
pct = (exp(beta) - 1) * 100;
```

The saved `results` struct should contain both log-scale (`beta`, `ci`) and percent-scale (`pct_est`, `pct_lo`, `pct_hi`) quantities so that log-scale values remain available for technical appendices or downstream analyses.

---

## 5. Outputs Required

Codex should generate MATLAB code that:

1. Loads the dataset.
2. Filters baseline and Affera cases correctly.
3. Constructs `time_days`, `log_duration`, `isAffera`, and `isPVIplus`.
4. Fits the mixed-effects model exactly as specified.
5. Extracts:
   - beta1 (Affera effect for PVI)
   - beta4 (additional Affera effect for PVI+)
   - 95% CIs for all coefficients
   - Percentage effects:

     ```matlab
     pct = (exp(beta) - 1) * 100;
     ```

6. Prints a clean, human-readable summary.

---

## 6. Additional Notes

- Do **not** subsample or equalize baseline cases.
- Do **not** exclude operators based on case mix.
- All AF ablations in the baseline window are valid.
- Affera cases may include PVI and/or PVI+.
- The model must be fit using **all available cases**, not randomly sampled subsets.

---

## 7. Implementation Plan

This section documents the concrete, modular MATLAB implementation that consumes the output of `import_clinical_data.m` and runs the Q1 mixed-effects model.

### 7.1 Data source and assumptions

- The import script `matlab/import_clinical_data.m` produces a table `tbl` (or `clinicalData`) with at least the following variables:
  - `procedure_id`
  - `procedure_date`
  - `operator_name`
  - `operator_id`
  - `procedure_type` (categorical, currently with categories `PVI` and `PVI+`)
  - `supplies_raw`
  - `catheter_primary`
  - `is_affera` (logical, true for Affera/Sphere-9 cases)
  - `is_pfa_catheter` (logical, true for any PFA catheter including Affera, Farawave, PulseSelect)
  - `tool_type`
  - `procedure_duration_min`
  - `log_duration` (already computed as `log(procedure_duration_min)`)
  - `procedure_primary`
  - `service`
  - `case_location`
- A companion struct `meta` contains:
  - `meta.affera_launch_date_global` — the first Affera case date (Affera_start).
- All rows in the imported data are AF ablation procedures; no extra arrhythmia filtering is required.

### 7.2 High-level script structure (DRY and modular)

New MATLAB entry point and helper functions:

1. `matlab/run_q1_affera_model.m` (function)
   - Signature:
     ```matlab
     [lme, results, tbl_q1] = run_q1_affera_model(clinicalData, baseline_days, min_cases_per_group)
     ```
   - Inputs:
     - `clinicalData`: table produced by `import_clinical_data` (or loaded from `clinicalData/importedClinicalData.mat`), with the variables listed in 7.1.
     - `baseline_days` (optional): length of the pre-Affera baseline window in days (default: 180).
     - `min_cases_per_group` (optional): minimum number of baseline cases **and** Affera cases required per operator to be included in the analysis (default: 15 for each group).
   - Outputs:
     - `lme`: fitted `LinearMixedModel` object.
     - `results`: struct with fixed-effect coefficients, CIs, and percent effects.
     - `tbl_q1`: final analysis dataset used to fit the model.
   - Internally:
     - Derives `meta.affera_launch_date_global` from `clinicalData` (earliest Affera case).
     - Calls a helper to construct the Q1 analysis dataset.
     - Calls a helper to fit the mixed-effects model.
     - Calls a helper to print (and optionally save) a clean summary.

2. `matlab/make_q1_affera_dataset.m` (function)
   - Input: `tbl`, `meta`, `baseline_days`, and `min_cases_per_group`.
   - Output: `tbl_q1`, a table that matches the variable requirements in Sections 2–3.
   - This function does *not* re-implement any import logic; it only derives analysis-specific variables from the existing `tbl` and applies optional operator-level inclusion rules.

3. `matlab/fit_q1_affera_model.m` (function)
   - Input: `tbl_q1`.
   - Output: `lme` (the `LinearMixedModel` object) and a compact `results` struct or table with betas, 95% CIs, and percent effects.

4. `matlab/print_q1_affera_summary.m` (function)
   - Input: `lme`, `results`, `tbl_q1`.
   - Prints a human-readable textual summary and, optionally, writes outputs into `results/` (e.g., `.mat` or `.csv` files).

### 7.3 Construction of the Q1 analysis dataset

This logic lives in `make_q1_affera_dataset.m` and operates only on the already-imported `tbl` and `meta`.

1. Affera launch date:
   - `Affera_start = meta.affera_launch_date_global;`

2. Baseline window definition:
   - Baseline window is the 180 days prior to the first Affera case:

     ```matlab
     baseline_start = Affera_start - days(180);
     baseline_end   = Affera_start;  % exclusive upper bound
     ```

3. Cohort flags:
   - Affera cohort:

     ```matlab
     isAffera_raw = tbl.is_affera;
     ```

   - Non-Affera, non-PFA baseline cohort (pre-Affera):

     ```matlab
     isBaseline = ~tbl.is_pfa_catheter & ...
                  ~isAffera_raw & ...
                  tbl.procedure_date >= baseline_start & ...
                  tbl.procedure_date <  baseline_end;
     ```

   - Non-Affera, non-PFA post-Affera cohort (v2 update):

     ```matlab
     isPostNonPFA = ~tbl.is_pfa_catheter & ...
                    ~isAffera_raw & ...
                    tbl.procedure_date >= Affera_start;
     ```

   - All procedures are AF ablations by design; no extra AF filter is required.

4. Inclusion mask and subset:

   ```matlab
   include = isAffera_raw | isBaseline | isPostNonPFA;
   tbl_q1  = tbl(include, :);
   ```

5. Operator-level inclusion (optional min cases per group):
   - For each operator:
     - Count baseline cases: `nBaseline_op = sum(tbl_q1.isBaselineEra & tbl_q1.operator_id == op)`.
     - Count Affera cases: `nAffera_op = sum(tbl_q1.isAffera ~= 0 & tbl_q1.operator_id == op)`.
   - Keep only operators with:

     ```matlab
     nBaseline_op >= min_cases_per_group && nAffera_op >= min_cases_per_group
     ```

   - Filter `tbl_q1` to include only rows belonging to retained operators.

6. Variable construction/renaming to match the spec:

   - `duration_minutes`:

     ```matlab
     tbl_q1.duration_minutes = tbl_q1.procedure_duration_min;
     ```

   - `isAffera` (analysis-ready):

     ```matlab
     tbl_q1.isAffera = double(tbl_q1.is_affera);  % or logical, as long as it is usable in fitlme
     ```

   - `isPVIplus` (PVI vs PVI+):
     - Derived from `procedure_type`, which currently has categories `PVI` and `PVI+`:

       ```matlab
       tbl_q1.isPVIplus = tbl_q1.procedure_type == 'PVI+';
       ```

     - This implies:
       - PVI-only reference group: `procedure_type == 'PVI'` → `isPVIplus == 0`.
       - PVI+ (complex PVI): `procedure_type == 'PVI+'` → `isPVIplus == 1`.

   - `operator_id` and `procedure_date` are already present and reused directly.

   - `time_days`:

     ```matlab
     tbl_q1.time_days = days(tbl_q1.procedure_date - min(tbl_q1.procedure_date));
     ```

   - `log_duration`:
     - Reuse the already computed column:

       ```matlab
       tbl_q1.log_duration = tbl_q1.log_duration;
       ```

   - Optional but useful:
     - `tbl_q1.isBaselineEra = isBaseline(include);` for reporting.

### 7.4 Mixed-effects model fitting

Implemented in `fit_q1_affera_model.m`:

1. Model specification:

   ```matlab
   lme = fitlme(tbl_q1, ...
       'log_duration ~ isAffera * isPVIplus + time_days + (1|operator_id)');
   ```

2. Coefficient extraction:
   - Retrieve fixed effects and their 95% confidence intervals:

     ```matlab
     beta    = fixedEffects(lme);
     ci      = coefCI(lme);  % same order as beta
     names   = lme.CoefficientNames;
     ```

   - Identify key coefficients:
     - `beta1` (Affera effect in PVI-only): coefficient with name `isAffera`.
     - `beta4` (additional Affera effect in PVI+): coefficient with name `isAffera:isPVIplus`.

3. Percent effect transformation:

   ```matlab
   pct_est = (exp(beta)      - 1) * 100;
   pct_lo  = (exp(ci(:, 1))  - 1) * 100;
   pct_hi  = (exp(ci(:, 2))  - 1) * 100;
   ```

4. Results packaging:
   - Construct a `results` struct or table containing:
     - `names`
     - `beta`, `ci_lo`, `ci_hi`
     - `pct_est`, `pct_lo`, `pct_hi`
     - indices or fields specifically marking `beta1` and `beta4`.

### 7.5 Summary and output

Implemented in `print_q1_affera_summary.m`:

1. Sample and operator counts:
   - Total number of included procedures.
   - Counts by baseline vs Affera (`isBaselineEra` vs `isAffera`).
   - Counts by PVI vs PVI+ (`isPVIplus`).
   - Number of operators (`numel(categories(tbl_q1.operator_id))`).

2. Coefficient summary:
   - For each fixed effect:
     - Name
     - Estimate (log scale), 95% CI.
     - Percent effect and 95% CI.

3. Highlighted interpretation:
   - Affera effect for PVI-only (beta1): print estimate and percent effect with clear label.
   - Additional Affera effect for PVI+ vs PVI (beta4): print estimate and percent effect with clear label.

4. Descriptive duration summaries:
   - Baseline (non-PFA, pre-Affera) procedure duration statistics:
     - Mean and median `duration_minutes` for:
       - All baseline procedures combined.
       - Baseline PVI-only (`isPVIplus == 0`).
       - Baseline PVI+ (`isPVIplus == 1`).
     - Mean and median `duration_minutes` by `catheter_primary` among baseline procedures (all PVI statuses combined).
   - Affera procedure duration statistics:
     - Mean and median `duration_minutes` for:
       - All Affera procedures combined.
       - Affera PVI-only (`isPVIplus == 0`).
       - Affera PVI+ (`isPVIplus == 1`).

5. Operator summary table:
   - For each operator appearing in the combined baseline-window + Affera cohort:
     - `operator_id`
     - `n_baseline`: number of baseline procedures for that operator.
     - `n_affera`: number of Affera procedures for that operator.
     - `included_in_model`: logical flag indicating whether the operator meets the minimum case thresholds and is retained in the final analysis dataset.
   - This table should be printed to the command line and saved as part of the `results` struct.

6. Optional file outputs:
   - Write a summary table or struct (including fixed effects, descriptive duration summaries, and operator counts) to `results/q1_affera_model_summary.mat` and/or `.csv`. The `results` struct saved by `run_q1_affera_model` should contain:
     - Fixed-effect coefficients and percent effects.
     - Baseline duration statistics (overall, by PVI status, and by catheter_primary).
     - Affera duration statistics (overall and by PVI status).
     - Operator-level baseline and Affera case counts.

### 7.6 Entry-point orchestration

The function `run_q1_affera_model.m` will:

1. Accept a `clinicalData` table (from `import_clinical_data` or `clinicalData/importedClinicalData.mat`).
2. Internally construct a `meta` struct containing:

   ```matlab
   meta.affera_launch_date_global = min(clinicalData.procedure_date(clinicalData.is_affera));
   ```

3. Call:

   ```matlab
   tbl_q1 = make_q1_affera_dataset(clinicalData, meta, baseline_days, min_cases_per_group);
   [lme, results] = fit_q1_affera_model(tbl_q1);
   print_q1_affera_summary(lme, results, tbl_q1);
   ```

4. Optionally save `lme`, `results`, and `tbl_q1` into `results/` for downstream exploration.

This design keeps the import step (`import_clinical_data.m`) separate and unchanged, adheres to DRY principles, and localizes all Q1-specific logic in small, testable functions. It also allows the user to adjust the baseline window length and the minimum number of baseline and Affera cases per operator (with defaults of 180 days and 15 cases, respectively).

---

## 8. Future extensions (optional)

Codex may optionally add:
- Random slopes for `isAffera` or `time_days`
- Nonlinear time components (splines)
- Operator-level clustering diagnostics

By default, all operators can be included; however, for robustness, analysis code may enforce user-configurable minimum case-count thresholds for baseline and Affera cases per operator (e.g., 15 each) as described in Section 7.3. The primary model structure (fixed effects and random effects) must still be implemented exactly as specified.

---

## 9. v2 Update: Include Post-Affera Non-PFA Cases

### Rationale

In the original Q1 model, all **baseline** (non-PFA) cases were drawn from the 6 months **before** Affera introduction, and **all cases after** Affera introduction were Affera. This created a structural confounding problem:

- Before Affera: only non-PFA cases (`isAffera == 0`).
- After Affera: only Affera cases (`isAffera == 1`).

As a result, **Affera use and calendar time were nearly perfectly collinear**, and the mixed-effects model could not reliably distinguish:

- the true **Affera effect**, vs.
- **secular efficiency improvements** over time.

To address this identifiability issue, v2 of the Q1 specification explicitly includes **non-PFA AF ablation cases that occur after Affera introduction** (Section 3.3). This creates **temporal overlap** between Affera and non-Affera cases, allowing the model to estimate:

- a **time trend** (secular change) based on both RF/cryo and Affera cases, and
- an **Affera-specific effect** that is not just a proxy for “later in time.”

Practically, this means the analysis dataset now consists of:

- All Affera cases (PFA) from `Affera_start` onward.
- All non-PFA AF ablations in the 180 days prior to `Affera_start` (baseline).
- All non-PFA AF ablations on or after `Affera_start` (post-Affera non-PFA).

The model formula and interpretation in Sections 4–5 remain unchanged, but the expanded non-PFA cohort improves identifiability of both the **time_days** term and the **Affera effect**.

---

## 10. Affera Learning-Curve Extension (Implemented)

### 10.1 Rationale

The primary Q1 mixed-effects model estimates the **average** effect of Affera compared with non-PFA ablation, but it does not distinguish:

- the **instantaneous** effect of adopting Affera (first few cases vs RF/cryo), from
- the **experience-related** improvement as operators gain proficiency with Affera.

Without explicitly modeling learning, early and late Affera cases are treated as exchangeable. This can obscure meaningful within-operator improvements and potentially bias the average Affera effect (up or down).

To address this, the analysis incorporates an **Affera case index**, `affera_index`, representing the sequence number of Affera cases performed by each operator (1, 2, 3, …). This allows:

- estimation of a **learning curve slope** (change in duration with Affera experience), and
- separation of “technology effect” from “learning effect.”

### 10.2 Model extension

Starting from the core Q1 model, the learning-curve extension adds a single fixed-effect term for a centered Affera case index:

- Construct per-operator Affera index:

  - For each operator_id, sort that operator’s Affera cases by `procedure_date` (and `procedure_id` as tiebreaker if needed).
  - Assign `affera_index = 1, 2, 3, …` for those Affera rows.
  - Non-Affera rows are given `affera_index = 0`.

- Center the Affera index across all Affera cases:

  ```matlab
  isA = tbl_q1.isAffera == 1;
  mean_idx = mean(tbl_q1.affera_index(isA), 'omitnan');
  tbl_q1.affera_index_ctr(isA) = tbl_q1.affera_index(isA) - mean_idx;
  tbl_q1.affera_index_ctr(~isA) = 0;
  ```

- Extended model (implemented when `EnableLearningCurve` is true):

  \[
  \log(\text{duration}) =
  \beta_0 +
  \beta_1(\text{Affera}) +
  \beta_2(\text{PVI+}) +
  \beta_3(\text{time}) +
  \beta_4(\text{Affera} \times \text{PVI+}) +
  \beta_5(\text{Affera index, centered}) +
  (1|\text{operator}).
  \]

  In MATLAB syntax this corresponds to

  ```matlab
  'log_duration ~ isAffera * isPVIplus + time_days + affera_index_ctr + (1|operator_id)'
  ```

  (Random intercept only; any random slopes are considered future extensions.)

### 10.3 Collinearity diagnostics

Because Affera experience tends to increase with calendar time, there is a risk that `affera_index_ctr` is collinear with `time_days` and/or `isAffera`. To monitor this, the implementation includes diagnostics (via `compute_learning_curve_diagnostics.m`) whenever the learning term is enabled:

- Predictors checked:

  ```matlab
  predictors = {'isAffera', 'isPVIplus', 'time_days', 'affera_index_ctr'};
  ```

- Diagnostics:
  - Correlation matrix (Pearson), with warnings if any |r| > 0.7.
  - Simple variance inflation factors (VIF) for each predictor, with warnings for VIF > 5 and VIF > 10.

These diagnostics are printed to the console but do not change the fitted model.

### 10.4 Outputs and interpretation

When `EnableLearningCurve` is enabled, the summary output and `results` struct additionally report:

- **Learning curve slope**:
  - Percent change in duration per additional Affera case:

    ```matlab
    pct_est_learning = results.pct_est(idxAfferaIndex);
    ```

  - 95% CI and p-value on the percent scale.

- **Practical contrast**:
  - Approximate percent change in duration from Affera case #1 to case #20, using:

    ```matlab
    deltaCases = 19;
    pct_delta = (exp(beta_learning * deltaCases) - 1) * 100;
    ```

- All Affera, PVI vs PVI+, and overall Affera effects remain defined and reported as in the core model, but now conditional on the mean centered Affera index.

Internally, the learning term is estimated on the log-duration scale; the console summary reports only percent changes and their CIs/p-values, while the log-scale coefficient and covariance are retained in `lme` and `results`.

---

## 11. Operator Baseline Efficiency Effect (Planned Extension)

### 11.1 Motivation

The current models estimate:

- an average Affera effect vs non-PFA,
- how that effect differs between PVI-only and PVI+, and
- how Affera experience (learning curve) influences procedure duration.

However, they do not yet address whether an **operator’s baseline efficiency** (on non-PFA cases) predicts the **magnitude of the Affera efficiency benefit**. Clinically, one might expect:

- Slower baseline operators (longer non-PFA durations) could see larger relative improvements with Affera.
- Faster baseline operators might see smaller incremental benefit.

To evaluate this, we plan to introduce an interaction term between Affera use and a continuous **baseline operator speed** metric.

### 11.2 Baseline operator speed definition

For each operator, we will derive a baseline efficiency measure using that operator’s **non-PFA baseline cases** (pre-Affera, within the defined baseline window):

- Candidate definition:
  - Mean log-duration of baseline non-PFA cases:

    ```matlab
    % Conceptual only; implementation TBD
    baseline_log_dur_op = mean(log(duration_minutes_baseline_op));
    ```

  - Or, equivalently, mean baseline duration (possibly log-transformed and standardized).

- This baseline speed will then be:
  - Attached to each procedure for that operator as an operator-level covariate, e.g. `baseline_speed_op`.
  - Centered across operators (and potentially scaled) to reduce collinearity and aid interpretation.

### 11.3 Planned model extension

Starting from the Q1 (and optionally learning-curve) model, we plan to add:

- A main effect for baseline operator speed.
- An interaction between Affera use and baseline speed.

Conceptual model:

\[
\log(\text{duration}) =
\beta_0 +
\beta_1(\text{Affera}) +
\beta_2(\text{PVI+}) +
\beta_3(\text{time}) +
\beta_4(\text{Affera} \times \text{PVI+}) +
\beta_5(\text{baseline speed}) +
\beta_6(\text{Affera} \times \text{baseline speed}) +
(1|\text{operator}) \quad [+ \; \beta_7(\text{Affera index}) \text{ if learning curve enabled}].
\]

In MATLAB formula notation this would look like:

```matlab
'log_duration ~ isAffera * isPVIplus + time_days + baseline_speed + isAffera:baseline_speed + ... + (1|operator_id)'
```

where `baseline_speed` is an operator-level covariate derived from non-PFA baseline cases, and the `...` indicates optional inclusion of `affera_index_ctr` when the learning-curve term is enabled.

### 11.4 Interpretation of the operator-speed interaction

- The **baseline_speed** main effect captures cross-operator differences in duration in the non-PFA era (faster vs slower operators).
- The key interaction term **Affera × baseline_speed** tests whether the Affera effect depends on an operator’s baseline efficiency:
  - A negative interaction coefficient (on the log scale) implies that **slower baseline operators** (higher baseline duration) experience **larger relative reductions in duration** when using Affera.
  - A positive interaction implies the opposite.

Reported quantities (once implemented) will follow the same percent-scale convention as the core model:

- Percent change in duration per 1-unit (or 1–SD) increase in baseline speed.
- Difference in Affera percent effect comparing “slow” vs “fast” operators (e.g., between the 25th and 75th percentiles of baseline speed).

### 11.5 Collinearity and reporting

Because baseline speed is derived solely from **pre-Affera non-PFA** cases, it should be less correlated with `isAffera` and `time_days` than naive measures that mix eras. Nevertheless, once implemented:

- `baseline_speed` will be included in the same correlation/VIF diagnostics as the existing predictors.
- Console output will report:
  - Percent-scale effects and CIs for baseline speed and its interaction with Affera.
  - No log-scale coefficients in the user-facing summary.

As with all other fixed effects, the underlying log-scale estimates and covariance for the operator-speed terms will be preserved in `lme` and `results` for potential technical appendix use.

---

# END OF SPECIFICATION
