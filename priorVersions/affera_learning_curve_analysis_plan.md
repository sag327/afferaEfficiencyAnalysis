# Affera Learning Curve Analysis – Codebase Plan

## 1. Study context and goals

- **Clinical context**: Atrial fibrillation ablations, with the recent introduction of the Affera tool.
- **Sample**: ~1300 Affera cases over ~10 months, performed by 16 operators, plus historical cases using other tools.
- **Primary outcome**: Procedure duration (minutes).
- **High-level goal**: Quantify how procedure duration changes with Affera use over time, at both the global and operator levels, accounting for baseline operator speed and case mix.

## 2. Key questions to address

1. **Effect of Affera on procedure duration (overall effect)**  
   - How does procedure duration differ between Affera and non-Affera procedures, after adjusting for operator, temporal trends, and case mix?

2. **Learning curve with Affera (within-operator over time)**  
   - For each operator, how does procedure duration change as they gain experience with Affera (e.g., case number 1 vs. 10 vs. 50)?

3. **Impact of baseline operator speed on learning curve**  
   - Do faster baseline operators (before Affera) learn Affera faster or slower than slower baseline operators?

4. **Impact of baseline operator speed on Affera vs. non-Affera speed difference**  
   - Is the Affera vs. non-Affera duration difference larger for slow operators than for fast operators, after learning has occurred?

These questions will guide both the **statistical models** and the **codebase structure**.

---

## 3. Data structure assumptions

We will treat the data as **procedure-level observations**, with one row per procedure. We expect at least the following variables (names can be modified later to match the raw data):

### 3.1. Core identifiers
- `patient_id` (optional): Unique patient identifier (for clustering/control if repeat procedures).
- `operator_id`: Identifier for the operator (1–16).
- `procedure_id`: Unique procedure identifier.

### 3.2. Timing
- `procedure_date`: Calendar date of the procedure.
- `procedure_start_time` / `procedure_end_time` (optional): To compute duration if not already available.
- `affera_launch_date`: Date Affera became available at the center (may be global or operator-specific if adoption staggered).

### 3.3. Exposure / tool indicators
- `tool_type`: Categorical (e.g., `affera`, `other_tool_A`, `other_tool_B`, etc.).
- `is_affera`: Binary indicator (1 = Affera, 0 = non-Affera).
- `is_post_affera_era`: Binary indicator whether procedure occurred after Affera introduction (useful for diff-in-diff style models).

### 3.4. Outcome
- `procedure_duration_min`: Continuous outcome, total duration in minutes.

### 3.5. Case mix / confounders (to refine later)
- `procedure_type` (e.g., PAF vs. persistent AF, redo vs. first-time).
- `patient_age`, `sex`, key comorbidities.
- `left_atrial_size`, other anatomical or imaging measures if available.
- `anesthesia_type`, `mapping_system`, and other operational factors that may affect duration.

### 3.6. Derived learning-curve variables (to be created in code)
For Affera procedures:
- `affera_case_index_operator`: Within-operator sequence number of Affera cases (1, 2, 3, …).
- `affera_case_index_center`: Overall Affera case index at the center (if relevant).
- `time_since_first_affera_operator`: Days or months since each operator’s first Affera case.

Affera launch timing (global and operator-level):
- `affera_launch_date_global`: The earliest `procedure_date` among all Affera cases in the dataset. This defines the global Affera launch date for primary analyses.
- `time_since_affera_launch`: Continuous variable (e.g., days) = `procedure_date - affera_launch_date_global`. Can be negative for pre-launch cases if desired, or truncated at 0.
- `is_pre_affera_global`: Indicator (1 = `procedure_date` < `affera_launch_date_global`).
- `is_post_affera_global`: Indicator (1 = `procedure_date` ≥ `affera_launch_date_global`).
- `is_pfa_catheter`: Indicator for PFA tools (e.g., Sphere-9, Farawave, PulseSelect) based on `supplies_raw` or catheter labels, used to optionally exclude PFA from baseline calculations.

Optional (for later refinement, not required immediately):
- `first_affera_date_operator`: For each operator, date of their first Affera case.
- `is_pre_affera_operator`: Indicator (1 = `procedure_date` < `first_affera_date_operator` for that operator).
- `is_post_affera_operator`: Indicator (1 = `procedure_date` ≥ `first_affera_date_operator`).

For baseline speed:
- `baseline_n_cases_global`: The largest integer \( N \) such that **all operators** used in baseline-based analyses have at least \( N \) eligible pre-Affera cases available, subject to an analyst-specified case-range constraint (e.g., `baselineCaseNumberRange = [N_{\min}, N_{\max}]`). Concretely: for each operator \( j \), compute \( C_j \) = number of pre-Affera cases eligible for baseline; set `baseline_n_cases_global = min( \min_j C_j, N_{\max} )` over operators with \( C_j \ge N_{\min} \).
- `baseline_sample_case_ids_operator`: For each operator, the IDs (or row indices) of the \( N \) pre-Affera cases selected to define baseline (e.g., the **most recent** \( N \) pre-Affera cases, sorted by procedure date within operator).
- `baseline_duration_mean_operator`: Average procedure duration for each operator using exactly `baseline_n_cases_global` pre-Affera (global) cases from the baseline sample.
- `baseline_duration_median_operator` (optional): Median duration over the same baseline sample.
- `baseline_speed_operator`: A transformed measure (e.g., centered/scaled `baseline_duration_mean_operator`) to use in interactions.

### 3.7. Operator eligibility (global cohort)

To keep all downstream comparisons fair, define a **single global operator cohort** before any procedure-type filtering or modeling:

- Compute operator-level counts and baseline availability using the full dataset (pre- and post-Affera, all procedure types).
- Require every included operator to satisfy the strictest criteria needed for any question (e.g., valid baseline speed, at least `baselineCaseNumberRange(1)` pre-Affera non-Affera cases, at least `MinAfferaCases` Affera cases, and sufficient non-Affera totals for Q4).
- Use this same cohort for **all** analyses (Questions 1–4) and for every procedure-type subset (`all`, `PVI`, `PVI+`). Operators that are globally eligible remain in the models even when baseline data rely on non-PVI procedures, ensuring comparable cohorts across subsets.
---

## 4. Statistical analysis blueprint

### 4.1. Notation
- Let \( i \) index procedures, \( j \) index operators.
- Outcome: \( Y_{ij} \) = procedure duration (e.g., \( \log(\text{minutes}) \) if needed).
- Main exposure: `is_affera` indicator and learning-curve index for Affera cases.

We will likely use **mixed-effects models** (or hierarchical Bayesian models) to capture:
- Operator-specific baselines.
- Operator-specific learning curves.
- Correlation within operators and across time.

#### Transformations
- Consider modeling \( \log(Y_{ij}) \) if the distribution is skewed.
- Center and scale continuous predictors (e.g., baseline speed, time, case index).

---

### 4.2. Question 1 – Overall effect of Affera

**Goal**: Estimate the average difference in procedure duration between Affera and non-Affera procedures, adjusting for covariates.

**Core model idea** (mixed-effects regression):

- Model (conceptual form):
  \[
  Y_{ij} = \beta_0 + \beta_1 \text{Affera}_{ij} + \beta_2 \text{Time}_{ij} + \gamma^\top X_{ij}
    + b_{0j} + \epsilon_{ij}
  \]
  - `Affera_ij`: 1 if Affera, else 0.
  - `Time_ij`: calendar time or time since Affera introduction (to capture secular trends).
  - `X_ij`: vector of case-mix covariates.
  - `b_{0j}`: random intercept for operator \( j \).

**Key outputs**:
- \( \beta_1 \): adjusted mean difference (or ratio, if log-scale) in duration between Affera and non-Affera.
- Sensitivity analyses: different time adjustments, case-mix specification, alternative outcomes (e.g., log vs. raw minutes).

---

### 4.3. Question 2 – Learning curve with Affera

**Goal**: Characterize how procedure duration changes with Affera experience for each operator.

**Key predictor**:
- `affera_case_index_operator` (1, 2, 3, …), possibly transformed as `log(1 + index)` or modeled with splines.

**Model sketch**:

  \[
  Y_{ij}^{\text{(Affera)}} = \beta_0 + f(\text{Index}_{ij}) + \gamma^\top X_{ij}
    + b_{0j} + b_{1j} \text{Index}_{ij} + \epsilon_{ij}
  \]

- \( f(\text{Index}_{ij}) \): global learning curve function (e.g., spline or log-linear).
- \( b_{1j} \): random slope on index to allow operator-specific learning rates.

**Outputs**:
- Global learning curve (average operator): predicted duration vs. Affera case number.
- Operator-specific curves: shrinkage estimates of each operator’s curve.
- Time to “plateau” or to reach a predefined target reduction.

---

### 4.4. Question 3 – Baseline operator speed vs. learning curve

**Goal**: Assess whether baseline operator speed modifies the learning curve with Affera.

**Key variables**:
- `baseline_speed_operator` (e.g., centered/scaled baseline average duration from **pre-Affera** cases only).
- Interaction between `baseline_speed_operator` and Affera learning index.

**Model extension**:

  \[
  Y_{ij}^{\text{(Affera)}} = \beta_0 + f(\text{Index}_{ij}) + \beta_2 \text{BaselineSpeed}_j
    + \beta_3 \left( \text{BaselineSpeed}_j \times \text{Index}_{ij} \right)
    + \gamma^\top X_{ij} + b_{0j} + b_{1j} \text{Index}_{ij} + \epsilon_{ij}
  \]

**Interpretation**:
- \( \beta_3 \neq 0 \): evidence that baseline speed modifies the learning curve slope.
- Visualizations: learning curves stratified by baseline speed tertiles (fast / medium / slow).

---

### 4.5. Question 4 – Baseline speed vs. Affera vs. non-Affera difference

**Goal**: Evaluate whether the Affera vs. non-Affera duration difference depends on baseline operator speed.

**Approach**:
- Use all procedures (Affera and non-Affera), but define `baseline_speed_operator` strictly from **pre-Affera** cases.
- Include interaction between `is_affera` and `baseline_speed_operator`.
- Optionally allow different learning curves for Affera vs. non-Affera over time.

**Model sketch**:

  \[
  Y_{ij} = \beta_0 + \beta_1 \text{Affera}_{ij}
    + \beta_2 \text{BaselineSpeed}_j
    + \beta_3 \left( \text{Affera}_{ij} \times \text{BaselineSpeed}_j \right)
    + \beta_4 \text{Time}_{ij}
    + \gamma^\top X_{ij}
    + b_{0j} + \epsilon_{ij}
  \]

**Outputs**:
- Predicted Affera vs. non-Affera difference at different levels of baseline speed.
- Operator-level summaries: e.g., predicted duration difference for each operator at their observed baseline speed.

---

### 4.6. Comparison groups and selection bias (pre- vs. post-Affera)

#### Key issue: post-launch non-Affera cases are selected

Once Affera is available, operators tend to:
- Use Affera for “routine” cases.
- Reserve non-Affera for special situations (redo/complex anatomy, device conflicts, equipment downtime, specific workflow constraints, early adoption phase, or days Affera is unavailable).

Consequences:
- Post-launch non-Affera cases are not representative of earlier non-Affera cases.
- Comparing post-launch Affera (often routine) to post-launch non-Affera (often atypical) introduces confounding by indication and selection bias that modeling alone cannot fully solve.

This is why relying solely on within-era non-Affera controls is risky.

#### Options for comparison groups

**Option A (preferred for “Affera vs traditional”) – Pre-Affera non-Affera vs. post-Affera Affera**
- Use pre-Affera, non-Affera cases as the primary comparison group.
- Rationale:
  - Pre-Affera cases are not subject to Affera-related selection.
  - They provide an unbiased “traditional” baseline with larger sample size and full operator variation.
  - Time trends can be modeled to transport the baseline forward.
- Conceptual model:
  \[
  \log(\text{duration}) \sim \text{Affera} + \text{time\_since\_affera\_launch}
    + \text{Affera} \times \text{time\_since\_affera\_launch (optional)}
    + (1|\text{operator})
  \]
- Advantages:
  - Cleaner estimate of the Affera effect.
  - Time accounts for secular improvements.
  - Pre-Affera baseline is fair and robust; typically more compelling to reviewers.

**Option B (secondary) – Post-launch Affera vs. post-launch non-Affera**
- Restrict to post-Affera era and compare contemporary Affera vs. non-Affera.
- Caveats:
  - Remaining non-Affera cases are selectively chosen.
  - Must adjust for as many covariates as possible (procedure type, redo status, anatomy, etc.).
  - Treat as supportive/sensitivity analysis rather than the primary estimate.

**Option C (not recommended as primary) – Post-launch Affera vs. non-Affera only**
- Direct comparison within the post-launch window only.
- This is the most biased approach if non-Affera cases are rare and selected.
- Still useful as an additional robustness check, but not for the main claim.

#### Recommended hybrid strategy

To align with the four core questions while minimizing bias:

- **Primary analysis for Q1 (overall Affera effect)**  
  - Compare **pre-Affera non-Affera** to **post-Affera Affera**, adjusting for time trends and operator effects.  
  - Exclude post-Affera non-Affera cases from this primary contrast.

- **Secondary analysis for Q1 (contemporary comparison)**  
  - Restrict to **post-Affera era only** and compare Affera vs. non-Affera, fully adjusted and clearly labeled as sensitivity analysis.

- **Learning-curve analyses for Q2–Q3**  
  - Use **post-Affera Affera-era data only**, since learning is defined with respect to Affera exposure.

- **Baseline operator speed for Q3–Q4**  
  - Compute **baseline operator speed using pre-Affera non-Affera data only** to avoid post-launch selection bias.
  - Use this baseline measure in:
    - Q3: interaction with Affera learning index.
    - Q4: interaction with Affera vs. non-Affera effect.

Implementation details for baseline speed (as implemented in the MATLAB analysis script):
- Define `affera_launch_date_global` as the earliest Affera `procedure_date` across all operators.
- For each operator \( j \), identify their **pre-Affera** cases as those with `procedure_date` < `affera_launch_date_global`, then define a baseline-eligible subset according to the chosen rule:
  - **All pre-Affera non-Affera**: `is_pre_affera_global & ~is_affera`.
  - **Non-PFA-only baseline** (if requested): `is_pre_affera_global & ~is_pfa_catheter`, where `is_pfa_catheter` flags Sphere-9, Farawave, PulseSelect, etc.
  - **PVI-only baseline** (if requested): additionally restrict to `procedure_type == 'PVI'` so that baseline is computed only from PVI cases.
- Let \( C_j \) = number of such baseline-eligible pre-Affera cases for operator \( j \).
- Compute `baseline_n_cases_global = min( \min_j C_j , N_{\max} )`, where `baselineCaseNumberRange = [N_{\min}, N_{\max}]` is a user-specified range (e.g., `[10 50]` or `[0 Inf]` for no cap). Only operators with \( C_j \ge N_{\min} \) are considered when taking the minimum. This is the maximum number of baseline cases that **every included operator** can contribute.
  - If some operators have very few pre-Affera cases (e.g., \( C_j = 0 \) or 1), the analyst may:
    - Lower the requirement (e.g., include only operators with \( C_j \geq N_{\text{min}} \)), or
    - Restrict baseline-based analyses (Q3–Q4) to operators meeting a minimum threshold.
- For each operator, sort their **baseline-eligible** pre-Affera cases in **descending** order by date (most recent first) and select the first `baseline_n_cases_global` cases as the **baseline sample**.
- Compute operator-level baseline statistics from this equal-sized sample:
  - `baseline_duration_mean_operator` (primary).
  - Optionally: `baseline_duration_median_operator`, `baseline_duration_sd_operator`.
- Transform `baseline_duration_mean_operator` (e.g., center/scale) to create `baseline_speed_operator` for use in Q3–Q4 models.

#### Operator inclusion thresholds should be **question-specific**

Do not globally exclude operators before modeling. Instead, determine eligibility per question based on the data each model requires:
- In the MATLAB implementation this is controlled by two complementary arguments:
  - `MinAfferaCases` – per-operator minimum number of Affera cases required for inclusion. Used in:
    - **Q1**: operator must have ≥`MinAfferaCases` post-Affera Affera cases (in addition to the baseline requirement below).
    - **Q2**: operator must have ≥`MinAfferaCases` Affera cases.
    - **Q3**: operator must have ≥`MinAfferaCases` Affera cases *and* a valid baseline (per `baselineCaseNumberRange`).
    - **Q4**: operator must have ≥`MinAfferaCases` Affera cases and a valid baseline.
  - `baselineCaseNumberRange = [N_{\min}, N_{\max}]` – governs baseline eligibility:
    - Operators must have ≥`N_{\min}` eligible pre-Affera baseline cases to contribute baseline metrics or be included in Q1/Q3/Q4.
    - Baseline uses the `baseline_n_cases_global = min( \min_j C_j , N_{\max} )` most recent eligible cases per operator.
    - **Q4** additionally requires ≥`N_{\min}` non-Affera cases (mirroring the baseline minimum) to ensure both exposures are represented post-baseline.
  Operators failing a question’s `MinAfferaCases`/baseline-range rule are excluded from that question only; they remain eligible for other questions if their data support it.

When reporting results, print the operator-level counts and explicitly list which operators were excluded for each question (and why). This preserves as much data as possible while keeping each model internally valid.

Rationale:
- Pre-Affera data provide the cleanest non-Affera comparison.
- Post-Affera non-Affera data are biased but still informative for sensitivity checks.
- Modeling time trends captures lab-wide process improvements.
- Learning curves and Affera exposure inherently relate to the post-launch period.
- Baseline operator speed must be defined from an unbiased pre-Affera period to be meaningful.

This hybrid design maximizes internal validity while still exploiting contemporary data for robustness.

---

## 5. Proposed codebase structure

This section outlines how to organize the analysis code specifically for a **MATLAB** workflow using scripts and functions (no GUI).

### 5.1. High-level project layout

- `data_raw/`  
  - Raw input files (Affera-era data, historical data, operator metadata).
  - Read-only; do not modify.

- `data_processed/`  
  - Cleaned and merged datasets ready for analysis.
  - Example: `procedures_all.parquet`, `operator_baseline_summary.csv`.

- `matlab/`  
  - Core MATLAB `.m` files organized by responsibility (see below).

- `notebooks/` or `scratch/`  
  - Optional exploratory scripts for EDA, sanity checks, and figures (can also live directly under `matlab/`).

- `results/`  
  - Model outputs, tables, and figures used in the manuscript.

- `config/`  
  - Configuration files for paths, variable name mappings, and modeling options.

- `tests/` (optional but recommended)  
  - Unit tests for key data processing and modeling utilities (can be simple MATLAB scripts checking key functions).

### 5.2. Core MATLAB modules (within `matlab/`)

- `load_procedure_data.m`  
  - Functions to read raw datasets from `data_raw/` (e.g., CSV, Excel, or MAT files).
  - Centralizes paths and basic integrity checks (e.g., uniqueness of IDs).

- `clean_procedure_data.m`  
  - Cleans variable names, handles missing values, harmonizes date formats.
  - Constructs `procedure_duration_min` if needed from start/end times.

- `engineer_affera_features.m`  
  - Creates derived variables:
    - `is_affera`, `is_post_affera_era`.
    - `affera_case_index_operator`, `time_since_first_affera_operator`.
    - `is_pfa_catheter` and other indicators used by later modeling.

- `fit_q1_overall_effect.m`  
  - Implements mixed-effects models for Question 1 using `fitlme` (Statistics and Machine Learning Toolbox).
  - Example formula (on a log scale):  
    - `lme = fitlme(tbl, 'log_duration ~ is_affera + time + covariates + (1|operator_id)');`
  - Saves fitted model objects and key summaries to `results/`.

- `fit_q2_q3_learning_curve.m`  
  - Implements learning-curve models with operator-specific slopes.  
  - Uses formulas like:  
    - `lme = fitlme(tbl_affera, 'log_duration ~ index + baseline_speed + index:baseline_speed + covariates + (1 + index|operator_id)');`
  - Produces outputs for Questions 2 and 3.

- `fit_q4_baseline_interaction.m`  
  - Models the Affera vs. non-Affera difference as a function of baseline speed.  
  - Uses an interaction term `is_affera:baseline_speed` in `fitlme`.

- `plot_learning_curves.m`  
  - Generates learning-curve plots (overall and by operator).

- `plot_operator_effects.m`  
  - Visualizes operator-level random effects and Affera vs. non-Affera differences (e.g., forest plots).

- `affera_utils.m` (or a `+affera` package folder)  
  - Shared helper functions (e.g., path handling, transformation utilities, default plotting styles).

### 5.3. High-level MATLAB analysis scripts

Create a small set of top-level `.m` scripts in `matlab/` to run the pipeline end-to-end:

- `run_prepare_data.m`  
  - Adds necessary paths.  
  - Calls `load_procedure_data`, `clean_procedure_data`, and `engineer_affera_features`.  
  - Saves processed datasets to `data_processed/` (e.g., `procedures_all.mat` and `operator_baseline_summary.mat`).

- `run_models_q1.m`  
  - Loads processed data from `data_processed/`.  
  - Calls `fit_q1_overall_effect`.  
  - Saves model objects and summary tables to `results/`.

- `run_models_q2_q3.m`  
  - Loads processed data.  
  - Calls `fit_q2_q3_learning_curve`.  
  - Saves learning-curve estimates and diagnostic outputs.

- `run_models_q4.m`  
  - Loads processed data.  
  - Calls `fit_q4_baseline_interaction`.  
  - Saves predicted Affera vs. non-Affera differences vs. baseline speed.

- `run_generate_figures_and_tables.m`  
  - Loads saved model outputs.  
  - Calls plotting and table-generation functions to produce manuscript-ready outputs in `results/`.

Sequential procedure-type analyses:
- In the MATLAB implementation, allow an option to **run all procedure-type strata sequentially** from a single call (e.g., `ProcedureType = 'sequential'`), which:
  - Runs the full analysis for:
    - Pooled PVI + PVI+ (`ProcedureType = 'all'`),
    - PVI-only (`ProcedureType = 'PVI'`),
    - PVI+-only (`ProcedureType = 'PVI+'`).
  - Returns a structured result with separate components for each stratum (e.g., `.all`, `.PVI`, `.PVIplus`), each containing its own models, baseline summaries, and operator inclusion/exclusion tables.
 - Each analysis call also writes a plain-text report (`affera_analysis_<ptype>_YYYYMMDD_HHMMSS.txt`) under `results/`, with:
   - A high-level summary of settings and key Q1–Q4 estimates.
   - A detailed section containing operator-level counts, baseline summaries, per-question exclusions, and full mixed-model output (`fitlme` summaries).

---

## 6. Implementation checklist

Use this section as a living checklist while building the codebase.

### 6.0. Phase 1 – Bare-bones MATLAB script
- [ ] Create `data_raw/` and place a single procedure-level data file (e.g., `procedures.csv`) there.
- [ ] Implement `matlab/run_basic_affera_analysis.m` that:
  - [ ] Loads the raw table with one row per procedure using `readtable`.
  - [ ] Creates minimal derived variables: `is_affera`, `log_duration`, `time_days`, simple `affera_index`, and per-operator baseline speed from non-Affera cases.
  - [ ] Defines a single global operator cohort that satisfies the strictest baseline/affera requirements and reuses it across Questions 1–4 and across `all` vs. `PVI` vs. `PVI+` subsets.
  - [ ] Fits four `fitlme` models corresponding to Questions 1–4 (overall Affera effect, simple linear learning curve, baseline speed vs. learning, baseline speed vs. Affera vs. non-Affera difference).
  - [ ] Prints model summaries to the command window.
- [ ] Implement and run a simple **parity test** script (`matlab/check_import_parity.m`) that:
  - [ ] Re-reads the original Excel file from `clinicalData/`.
  - [ ] Loads the processed `tbl` from `data_processed/procedures_all.mat`.
  - [ ] Confirms that the number of rows and unique `caseID`/`procedure_id` in the raw data match the processed table.
  - [ ] Confirms that counts by `Code` (PVI vs. PVI+) in the raw data match counts by `procedure_type` in the processed table.
  - [ ] Confirms that the number of rows with `Sphere-9 Catheter` in `Supplies` in the raw data matches the number of `is_affera == true` rows in the processed table.

### 6.1. Data preparation
- [ ] Inventory all raw files and document their structure.
- [ ] Define final variable names and types in a data dictionary.
- [ ] Implement `data_loading` module and basic integrity checks.
- [ ] Implement `data_cleaning` to produce a clean procedure-level dataset.
- [ ] Implement `feature_engineering` to compute:
  - [ ] `affera_launch_date_global` based on the earliest Affera `procedure_date`.
  - [ ] Global pre/post flags: `is_pre_affera_global`, `is_post_affera_global`, `time_since_affera_launch`.
  - [ ] Affera learning indices (`affera_case_index_operator`, `time_since_first_affera_operator`).
  - [ ] PFA and catheter indicators (`is_pfa_catheter`, `catheter_primary`) and any other fields needed to define baseline-eligible cohorts.

### 6.2. Modeling
- **Question 1 (overall effect)**  
  - [ ] Implement mixed-effects model(s) for Affera vs. non-Affera duration.  
  - [ ] Decide on outcome scale (raw vs. log-transformed).  
  - [ ] Perform sensitivity analyses for time trends and case mix.

- **Question 2 (learning curve)**  
  - [ ] Define learning curve functional form (e.g., spline, log index).  
  - [ ] Implement models with random intercepts and slopes by operator.  
  - [ ] Extract and visualize operator-specific and overall curves.

- **Question 3 (baseline speed vs. learning)**  
  - [ ] Compute baseline-speed metrics per operator from pre-Affera data.  
  - [ ] Add baseline-speed × learning-index interaction to learning-curve models.  
  - [ ] Visualize learning curves stratified by baseline speed.

- **Question 4 (baseline speed vs. Affera vs. non-Affera difference)**  
  - [ ] Implement models with Affera × baseline-speed interaction.  
  - [ ] Generate predicted effects at different baseline-speed levels.  
  - [ ] Summarize operator-level differences.

### 6.3. Reporting and visualization
- [ ] Define key figures (e.g., learning curves, operator forest plots, time trends).  
- [ ] Implement plotting functions with consistent aesthetics.  
- [ ] Create table-generating code for summary statistics and model results.

### 6.4. Reproducibility and documentation
- [ ] Add a `README.md` describing how to run the analysis.  
- [ ] Document versions of key packages and software.  
- [ ] (Optional) Set up tests for critical data-processing and modeling functions.

---

## 7. Open design decisions / to discuss

Use this section to capture choices that still need to be finalized:

- Confirm MATLAB toolboxes (e.g., Statistics and Machine Learning Toolbox for `fitlme`, and any spline/curve-fitting tools) are available and documented.
- Exact functional form of the learning curve (spline vs. log vs. piecewise linear).
- Handling of outliers in procedure duration (e.g., extreme complications, unusual cases).
- Strategy for case-mix adjustment and which covariates to include in the primary vs. sensitivity models.
- Whether to allow operator-specific time trends separate from Affera learning (e.g., center-wide process improvements).

As the codebase develops, this document should be updated to reflect final decisions and any deviations from the initial plan.
