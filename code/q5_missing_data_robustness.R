# Q5: Missing Data and Robustness Analysis
#
# Purpose:
# 1. Describe missingness in longitudinal TICS measurements.
# 2. Examine whether missingness is associated with observed
#    participant characteristics.
# 3. Refit the Q3 hierarchical model after mean imputation.
# 4. Compare the imputed results with the complete-case analysis.
#
# Software: R + JAGS


library(rjags)
library(coda)


# ---------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------

TICS <- read.csv("tics.data.2025.csv")


# ---------------------------------------------------------
# 2. Summarize missingness
# ---------------------------------------------------------

tics_vars <- c(
  "TICS01",
  "TICS02",
  "TICS03",
  "TICS04",
  "TICS05"
)

miss_summary <- colSums(
  is.na(
    TICS[, tics_vars]
  )
)

miss_percent <- round(
  miss_summary / nrow(TICS) * 100,
  1
)

missing_table <- data.frame(
  Variable = names(miss_summary),
  Missing = as.numeric(miss_summary),
  Percent_Missing = as.numeric(miss_percent)
)

print(missing_table)


# ---------------------------------------------------------
# 3. Create missingness indicator for TICS05
# ---------------------------------------------------------

TICS$miss_TICS05 <- as.numeric(
  is.na(TICS$TICS05)
)


# ---------------------------------------------------------
# 4. Logistic regression for missingness
# ---------------------------------------------------------
#
# Outcome:
#   miss_TICS05 = 1 if TICS05 is missing
#
# Predictors:
#   group
#   baseline TICS
#   baseline age
#   sex
#
# This examines whether missingness is associated
# with observed participant characteristics.


model_mar <- glm(
  miss_TICS05 ~ group + TICS01 + Age01 + sex,
  data = TICS,
  family = binomial
)

summary(model_mar)


# Odds ratios and confidence intervals
odds_ratios <- exp(
  coef(model_mar)
)

conf_intervals <- exp(
  confint(model_mar)
)

mar_results <- data.frame(
  Odds_Ratio = odds_ratios,
  CI_Lower = conf_intervals[, 1],
  CI_Upper = conf_intervals[, 2]
)

print(mar_results)


# ---------------------------------------------------------
# 5. Mean imputation
# ---------------------------------------------------------
#
# Mean imputation is used here as a sensitivity analysis
# to examine how conclusions from the complete-case model
# change when missing values are replaced by observed means.
#
# This is NOT treated as a preferred missing-data method.


TICS_imp <- TICS


# ---------------------------------------------------------
# 6. Impute age variables
# ---------------------------------------------------------

age_vars <- c(
  "Age02",
  "Age03",
  "Age04",
  "Age05"
)

for (col in age_vars) {

  TICS_imp[[col]][
    is.na(TICS_imp[[col]])
  ] <- mean(
    TICS_imp[[col]],
    na.rm = TRUE
  )
}


# ---------------------------------------------------------
# 7. Impute follow-up TICS scores
# ---------------------------------------------------------

tics_followup <- c(
  "TICS02",
  "TICS03",
  "TICS04",
  "TICS05"
)

for (col in tics_followup) {

  TICS_imp[[col]][
    is.na(TICS_imp[[col]])
  ] <- mean(
    TICS_imp[[col]],
    na.rm = TRUE
  )
}


# ---------------------------------------------------------
# 8. Impute BMI
# ---------------------------------------------------------

TICS_imp$BMI[
  is.na(TICS_imp$BMI)
] <- mean(
  TICS_imp$BMI,
  na.rm = TRUE
)


# ---------------------------------------------------------
# 9. Impute baseline TICS if necessary
# ---------------------------------------------------------

TICS_imp$TICS01[
  is.na(TICS_imp$TICS01)
] <- mean(
  TICS_imp$TICS01,
  na.rm = TRUE
)


# ---------------------------------------------------------
# 10. Calculate follow-up time using imputed ages
# ---------------------------------------------------------

Time_imp <- numeric(4)

Time_imp[1] <- mean(
  TICS_imp$Age02 - TICS_imp$Age01,
  na.rm = TRUE
)

Time_imp[2] <- Time_imp[1] +
  mean(
    TICS_imp$Age03 - TICS_imp$Age02,
    na.rm = TRUE
  )

Time_imp[3] <- Time_imp[2] +
  mean(
    TICS_imp$Age04 - TICS_imp$Age03,
    na.rm = TRUE
  )

Time_imp[4] <- Time_imp[3] +
  mean(
    TICS_imp$Age05 - TICS_imp$Age04,
    na.rm = TRUE
  )


# ---------------------------------------------------------
# 11. Create longitudinal outcome matrix
# ---------------------------------------------------------

y_imp <- as.matrix(
  TICS_imp[, tics_followup]
)


# ---------------------------------------------------------
# 12. Prepare imputed data for JAGS
# ---------------------------------------------------------

data_imp <- list(
  N = nrow(TICS_imp),

  y = y_imp,

  Time = Time_imp,

  group = TICS_imp$group,

  TICS01 = TICS_imp$TICS01,

  Age01 = TICS_imp$Age01,

  BMI = TICS_imp$BMI,

  Sex = TICS_imp$sex,

  Smoking = TICS_imp$SH.Ever.Smoked.,

  HTN = TICS_imp$MC.HTN,

  Diabetes = TICS_imp$MC.Diabetes.Mellitus,

  Stroke = TICS_imp$MC.Stroke,

  CAD = TICS_imp$MC.Coronary.Artery.Disease
)


# ---------------------------------------------------------
# 13. Hierarchical model using imputed data
# ---------------------------------------------------------

model_q5 <- "

model {

  for (i in 1:N) {

    for (j in 1:4) {

      y[i,j] ~ dnorm(mu[i,j], tau)

      mu[i,j] <- alpha[i]
                   + (beta + delta * group[i]) * Time[j]
                   + gamma * TICS01[i]
                   + b_age * (Age01[i] - mean(Age01[]))
                   + b_bmi * (BMI[i] - mean(BMI[]))
                   + b_sex * Sex[i]
                   + b_smoke * Smoking[i]
                   + b_htn * HTN[i]
                   + b_diab * Diabetes[i]
                   + b_stroke * Stroke[i]
                   + b_cad * CAD[i]
    }

    alpha[i] ~ dnorm(mu_alpha, tau_alpha)
  }


  mu_alpha ~ dnorm(0, 0.001)

  beta ~ dnorm(0, 0.001)

  delta ~ dnorm(0, 0.001)

  gamma ~ dnorm(0, 0.001)


  b_age ~ dnorm(0, 0.001)

  b_bmi ~ dnorm(0, 0.001)

  b_sex ~ dnorm(0, 0.001)

  b_smoke ~ dnorm(0, 0.001)

  b_htn ~ dnorm(0, 0.001)

  b_diab ~ dnorm(0, 0.001)

  b_stroke ~ dnorm(0, 0.001)

  b_cad ~ dnorm(0, 0.001)


  tau ~ dgamma(1, 1)

  tau_alpha ~ dgamma(1, 1)


  sigma <- 1 / sqrt(tau)

  sigma_alpha <- 1 / sqrt(tau_alpha)
}

"


# ---------------------------------------------------------
# 14. Fit imputed model
# ---------------------------------------------------------

jags_imp <- jags.model(
  textConnection(model_q5),
  data = data_imp,
  n.adapt = 3000,
  n.chains = 4
)


# Burn-in
update(
  jags_imp,
  3000
)


# ---------------------------------------------------------
# 15. Posterior sampling
# ---------------------------------------------------------

samples_imp <- coda.samples(
  jags_imp,

  c(
    "mu_alpha",
    "beta",
    "delta",
    "gamma",
    "b_age",
    "b_bmi",
    "b_sex",
    "b_smoke",
    "b_htn",
    "b_diab",
    "b_stroke",
    "b_cad",
    "sigma",
    "sigma_alpha"
  ),

  n.iter = 30000
)


# ---------------------------------------------------------
# 16. Summarize imputed model
# ---------------------------------------------------------

summary(
  samples_imp
)


# ---------------------------------------------------------
# 17. Gelman-Rubin diagnostics
# ---------------------------------------------------------

gelman.diag(
  samples_imp
)


# ---------------------------------------------------------
# 18. Compare key parameters with complete-case model
# ---------------------------------------------------------

summary(
  samples_imp[, c(
    "beta",
    "delta",
    "gamma",
    "b_age",
    "b_bmi",
    "b_sex"
  )]
)
