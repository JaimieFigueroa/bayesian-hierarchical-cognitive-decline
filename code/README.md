# Analysis Code

This folder contains the R and JAGS code used for the Bayesian cognitive aging analysis.

## Analyses

- `q1_baseline_analysis.R` — Baseline descriptive analysis and comparison between centenarian offspring and controls.
- `q2_baseline_bayesian.R` — Bayesian regression model evaluating baseline TICS scores.
- `q3_hierarchical_model.R` — Bayesian hierarchical model evaluating longitudinal change in TICS scores.
- `q4_latent_trajectory_model.R` — Bayesian latent trajectory mixture model identifying distinct cognitive trajectories.
- `q5_missing_data.R` — Missing-data assessment and sensitivity analysis.

The models were implemented in R using JAGS through the `rjags` package, with posterior analysis conducted using `coda`.
