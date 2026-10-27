# Ki-67, 1p36 loss and 9p21/CDKN2A loss as prognostic factors in skull base chordoma

Data, code and analysis report for the systematic review and meta-analysis *"Ki-67, 1p36 loss and 9p21/CDKN2A loss as prognostic factors in skull base chordoma: a systematic review and meta-analysis"* (PROSPERO CRD420261390460).

**Authors:** Cássio da Silva Minuzzo, Pedro Malfussi Horst, Nicolas de Cesaro Schpchacki, Francisco Martins Lamas, Francisco Ratto Fabris, Erion Junior de Andrade.

**Contact:** Cássio da Silva Minuzzo (cassio.minuzzo@ufcspa.edu.br), Universidade Federal de Ciências da Saúde de Porto Alegre (UFCSPA), Porto Alegre, Brazil.

## Analysis report

`index.html` is a self-contained report with every analysis of the review:
- study selection;
- characteristics of the included studies;
- all extracted and derived hazard ratios;
- pooled and display-only syntheses with forest plots;
- subgroup and sensitivity analyses;
- small-study effects;
- risk of bias (QUIPS), with the justification of each judgement;
- certainty of evidence (GRADE).

View it online at https://cassiominuzzo.github.io/skull-base-chordoma-biomarkers/, or download the file and open it in any web browser. All numbers are recomputed from `data/estimates.csv` when the report is rendered.

## Contents

| Path | Content |
|---|---|
| `index.html` | Rendered analysis report |
| `index.qmd` | Source of the report (Quarto) |
| `R/meta_analysis.R` | Meta-analysis code: pooled estimates, subgroup analysis, sensitivity analyses and Egger's test |
| `results/` | Tables written by `R/meta_analysis.R` |
| `data/estimates.csv` | All 94 hazard ratios extracted or derived from the 22 included studies |
| `data/cohort_attributes.csv` | Cohort-level variables used in subgroup and sensitivity analyses |
| `data/risk_of_bias_quips.csv` | QUIPS ratings and their justifications (40 assessments) |
| `data/study_characteristics.csv` | Characteristics of the included studies |
| `data/biomarker_assessment.csv` | Assay, specimen and cut-off of each biomarker in each study |
| `data/summary_of_findings.csv`, `data/grade_evidence_profile.csv` | Certainty of evidence (GRADE) |
| `data/excluded_reports.csv`, `data/prisma_counts.csv` | Study selection |
| `figures/prisma_flow_diagram.png` | PRISMA 2020 flow diagram |

## Main variables in `data/estimates.csv`

| Variable | Description |
|---|---|
| `estimate_id`, `study_id`, `study`, `year` | Identifiers of the estimate and the study |
| `biomarker`, `outcome`, `comparison`, `adjustment` | What the hazard ratio compares and whether it is adjusted |
| `n`, `events` | Patients and events in the analysis (blank when not reported) |
| `hr`, `ci_lower`, `ci_upper`, `se_log_hr` | Hazard ratio, 95% confidence interval and standard error of the log hazard ratio. All estimates are oriented so that values above 1 indicate worse prognosis with higher Ki-67 or with loss of 1p36 or 9p21 |
| `source_code`, `source` | How the estimate was obtained: reported; confidence interval derived from the p value; Cox model fitted to individual patient data; reconstructed Kaplan–Meier curve; or log-rank statistics (Tierney method) |
| `synthesis_id`, `role`, `use` | The synthesis the estimate belongs to and its role (main analysis, display, sensitivity analysis or not used) |
| `sensitivity_label` | Label of alternative estimates used in sensitivity analyses |
| `quips_overall` | Overall QUIPS risk of bias of the estimate for its outcome |

## Methods in brief

Adjusted and unadjusted hazard ratios were analysed separately, and binary and continuous scales were never combined. The primary analysis was the synthesis of adjusted binary Ki-67 estimates for recurrence or progression. Pooling required at least three independent compatible cohorts.

Pooled estimates come from inverse-variance random-effects models. τ² was estimated by restricted maximum likelihood, and confidence intervals used the Hartung–Knapp–Sidik–Jonkman method with the ad hoc variance correction (metafor, `test = "adhoc"`). Heterogeneity is reported as I², τ² and 95% prediction intervals.

Small-study effects were assessed only with at least 10 cohorts. Full methods are described in the article and its supplementary material.

## How to reproduce

Requirements:
- R 4.3 or later, with the packages `metafor`, `DT`, `ggplot2`, `knitr` and `rmarkdown`;
- [Quarto](https://quarto.org).

```bash
Rscript R/meta_analysis.R      # analyses only; writes results/
quarto render index.qmd        # full report; writes index.html
```

## Citation

The citation of the published article will be added here after publication.

## Licence

- Code (`R/`, `index.qmd`): MIT licence (`LICENSE`).
- Data, results, figure and report (`data/`, `results/`, `figures/`, `index.html`): CC BY 4.0 (`LICENSE-DATA.md`).
