# Supplementary data and code

Conflicts caused by cut-in behaviors: A five-type framework of drivers' reactions and conflict escalation

This file contains the anonymized case-level data, the code used for the analyses and the outputs behind the tables, figures and Appendices A-C of the manuscript.

## Data

The 1,003 cases were coded from judicial documents published on China Judgments Online. `case_id` is a sequential study number. Court case numbers, names and other identifying details are not included. Variable definitions are given in `codebook.csv` and Table 1 of the manuscript.

- `reproducibility/case_data.csv`: coded variables for the 1,003 cases (input to the analyses).
- `reproducibility/final_analysis_data.csv`: the same records with the derived variables and sample indicators used in the models.

## Where the reported results are

| Manuscript | File |
|---|---|
| Table 4, Table 5 | `reproducibility/tables/descriptives.csv`, `reproducibility/tables/injury.csv` |
| Table 6, Table B1 (Appendix B) | `reproducibility/tables/regression_categorical_adjusted.csv` |
| Tables 7-9, chi-square tests and residuals | `reproducibility/association_results.json`; time and road residual matrices in `figure_data/Figure_5_residuals.csv` and `figure_data/Figure_6_residuals.csv` |
| Table 10 | `reproducibility/tables/lca_fit_primary.csv` |
| Table 11 | `reproducibility/tables/lca_profiles.csv`, `reproducibility/tables/lca_classes.csv` |
| Figure 7 | `reproducibility/tables/lca_type_weighted.csv`; case-level posterior probabilities in `reproducibility/tables/lca_posterior.csv` |
| Figures 1-3, 5-7 | `figure_data/` |
| LCA without emotional expression (Appendix C) | `reproducibility/tables/lca_fit_no_emotion.csv` |
| Comparison of three and four classes (Appendix C) | `lca_comparison/summary.json`; full profiles of both models and bootstrap ranges in `lca_comparison/profiles_and_bootstrap_ranges.csv` |

## Running the analyses

R 4.5.3 with brglm2 1.1.1, poLCA 1.6.0.2 and jsonlite 2.0.0 (full session information in `reproducibility/R_session_info.txt`). From the `reproducibility` folder, run in order:

1. `Rscript analyze.R`: derived variables, descriptive summaries, chi-square tests and residuals, the adjusted multinomial model and the one- to five-class LCA.
2. `Rscript finalize_statistics.R`: pairwise LCA fit diagnostics and an independent recalculation of the observed-data likelihood.
3. `Rscript export_figure_data.R`: data underlying the figures.

The comparison of three and four classes is in `lca_comparison`. It requires the outputs of step 1 and Python 3 with numpy, scipy and numba, in addition to R and poLCA. From that folder, run `Rscript refit.R`, `python engine.py`, `python resample.py`, `python complete_cases.py`, `python check_simulations.py`, `Rscript check_simulations.R` and `python summarize.py`. The bootstrap and cross-validation steps are computationally intensive.

Random seeds are set in the scripts. Different software versions or random number generators may change the last digits of simulation-based results.
