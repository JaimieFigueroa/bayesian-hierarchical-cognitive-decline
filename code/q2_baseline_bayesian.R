# Q2: Bayesian Baseline Cognitive Model
# Bayesian Cognitive Aging Project
#
# Purpose:
# Estimate the association between centenarian offspring status
# and baseline cognitive function (TICS01) after adjustment for
# potential confounders.
#
# Model:
# TICS01 ~ Normal(mu, tau)
#
# The primary parameter of interest is beta1, representing the
# adjusted difference in baseline TICS score between centenarian
# offspring and controls.

# ---------------------------------------------------------
# 1. Load packages and data
# ---------------------------------------------------------

library(rjags)
library(coda)

data <- read.csv("tics.data.2025.csv")

TICS <- data

# Convert study group to numeric indicator
# 0 = Control
# 1 = Centenarian Offspring

TICS$ptype.num <- as.numeric(TICS$ptype) - 1


# ---------------------------------------------------------
# 2. Prepare covariates
# ---------------------------------------------------------

# Mean imputation was used for missing covariates in the
# original analysis.

impute_mean <- function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

Age01_imp <- impute_mean(TICS$Age01)
BMI_imp <- impute_mean(TICS$BMI)
ptype_b <- impute_mean(TICS$ptype.num)
sex_b <- impute_mean(TICS$sex)
smoked_b <- impute_mean(TICS$SH.Ever.Smoked.)
stroke_b <- impute_mean(TICS$MC.Stroke)
dm_b <- impute_mean(TICS$MC.Diabetes.Mellitus)
htn_b <- impute_mean(TICS$MC.HTN)
cad_b <- impute_mean(TICS$MC.Coronary.Artery.Disease)


# ---------------------------------------------------------
# 3. Bayesian model specification
# ---------------------------------------------------------

model.base <- "

model {

  # Likelihood
  for(i in 1:N) {

    TICS01[i] ~ dnorm(mu[i], tau)

    mu[i] <- beta0
             + beta1 * ptype[i]
             + beta2 * age[i]
             + beta3 * sex[i]
             + beta4 * bmi[i]
             + beta5 * smoked[i]
             + beta6 * stroke[i]
             + beta7 * dm[i]
             + beta8 * htn[i]
             + beta9 * cad[i]
  }

  # Weakly informative priors
  beta0 ~ dnorm(0, 0.001)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  beta3 ~ dnorm(0, 0.001)
  beta4 ~ dnorm(0, 0.001)
  beta5 ~ dnorm(0, 0.001)
  beta6 ~ dnorm(0, 0.001)
  beta7 ~ dnorm(0, 0.001)
  beta8 ~ dnorm(0, 0.001)
  beta9 ~ dnorm(0, 0.001)

  # Residual precision
  tau ~ dgamma(0.01, 0.01)

  # Residual standard deviation
  sigma <- 1 / sqrt(tau)
}
"


# ---------------------------------------------------------
# 4. Prepare data for JAGS
# ---------------------------------------------------------

jags_data <- list(
  N = nrow(TICS),
  TICS01 = TICS$TICS01,
  ptype = ptype_b,
  age = Age01_imp,
  sex = sex_b,
  bmi = BMI_imp,
  smoked = smoked_b,
  stroke = stroke_b,
  dm = dm_b,
  htn = htn_b,
  cad = cad_b
)


# ---------------------------------------------------------
# 5. Fit Bayesian model
# ---------------------------------------------------------

set.seed(123)

jags_model <- jags.model(
  textConnection(model.base),
  data = jags_data,
  n.chains = 3,
  n.adapt = 3000
)

# Burn-in
update(
  jags_model,
  5000
)


# ---------------------------------------------------------
# 6. Posterior sampling
# ---------------------------------------------------------

samples_q2 <- coda.samples(
  jags_model,
  variable.names = c(
    "beta0",
    "beta1",
    "beta2",
    "beta3",
    "beta4",
    "beta5",
    "beta6",
    "beta7",
    "beta8",
    "beta9",
    "sigma"
  ),
  n.iter = 30000
)


# ---------------------------------------------------------
# 7. Posterior summaries
# ---------------------------------------------------------

summary(samples_q2)


# ---------------------------------------------------------
# 8. Convergence diagnostics
# ---------------------------------------------------------

# Gelman-Rubin diagnostic
gelman.diag(
  samples_q2,
  multivariate = TRUE
)

# Geweke diagnostic
geweke.diag(samples_q2)


# ---------------------------------------------------------
# 9. Diagnostic plots
# ---------------------------------------------------------

png(
  "q2_trace_density.png",
  width = 1200,
  height = 800,
  res = 150
)

plot(
  samples_q2[, c(
    "beta0",
    "beta1",
    "beta2",
    "beta3"
  )]
)

dev.off()


# Gelman-Rubin diagnostic plot
png(
  "q2_gelman_rubin.png",
  width = 1000,
  height = 800,
  res = 150
)

gelman.plot(samples_q2)

dev.off()
