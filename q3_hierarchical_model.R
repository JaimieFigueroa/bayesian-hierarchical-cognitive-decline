# Q3: Bayesian Hierarchical Model
# Rate of Change in TICS Among Centenarian Offspring and Controls
#
# Outcome: TICS02-TICS05
# Exposure: Centenarian offspring vs. control
# Primary parameter: group-by-time interaction (delta)
#
# Model estimates whether the rate of cognitive change differs
# between centenarian offspring and controls after adjustment
# for baseline cognition and baseline risk factors.

library(rjags)
library(coda)

# ---------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------

TICS <- read.csv("tics.data.2025.csv")

# ---------------------------------------------------------
# 2. Calculate average follow-up time
# ---------------------------------------------------------

Time <- numeric(4)

Time[1] <- mean(TICS$Age02 - TICS$Age01, na.rm = TRUE)

Time[2] <- Time[1] +
  mean(TICS$Age03 - TICS$Age02, na.rm = TRUE)

Time[3] <- Time[2] +
  mean(TICS$Age04 - TICS$Age03, na.rm = TRUE)

Time[4] <- Time[3] +
  mean(TICS$Age05 - TICS$Age04, na.rm = TRUE)

Time


# ---------------------------------------------------------
# 3. Create longitudinal outcome matrix
# ---------------------------------------------------------

y <- as.matrix(
  TICS[, c("TICS02", "TICS03", "TICS04", "TICS05")]
)


# ---------------------------------------------------------
# 4. Select complete cases for model covariates
# ---------------------------------------------------------

cc_idx_q3 <- which(
  complete.cases(
    TICS[, c(
      "TICS01",
      "Age01",
      "BMI",
      "sex",
      "SH.Ever.Smoked.",
      "MC.HTN",
      "MC.Diabetes.Mellitus",
      "MC.Stroke",
      "MC.Coronary.Artery.Disease"
    )]
  )
)


# ---------------------------------------------------------
# 5. Prepare JAGS data
# ---------------------------------------------------------

data_q3 <- list(
  N = length(cc_idx_q3),
  y = y[cc_idx_q3, ],
  Time = Time,

  group = TICS$group[cc_idx_q3],

  TICS01 = TICS$TICS01[cc_idx_q3],
  Age01 = TICS$Age01[cc_idx_q3],
  BMI = TICS$BMI[cc_idx_q3],

  Sex = TICS$sex[cc_idx_q3],
  Smoking = TICS$SH.Ever.Smoked.[cc_idx_q3],
  HTN = TICS$MC.HTN[cc_idx_q3],
  Diabetes = TICS$MC.Diabetes.Mellitus[cc_idx_q3],
  Stroke = TICS$MC.Stroke[cc_idx_q3],
  CAD = TICS$MC.Coronary.Artery.Disease[cc_idx_q3]
)


# ---------------------------------------------------------
# 6. Bayesian hierarchical model
# ---------------------------------------------------------

model_q3 <- "

model {

  # Longitudinal outcome
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

    # Subject-specific random intercept
    alpha[i] ~ dnorm(mu_alpha, tau_alpha)
  }


  # Population-level parameters

  mu_alpha ~ dnorm(0, 0.001)

  # Average slope among controls
  beta ~ dnorm(0, 0.001)

  # Difference in slope between offspring and controls
  delta ~ dnorm(0, 0.001)

  # Effect of baseline TICS
  gamma ~ dnorm(0, 0.001)


  # Covariate effects

  b_age ~ dnorm(0, 0.001)
  b_bmi ~ dnorm(0, 0.001)
  b_sex ~ dnorm(0, 0.001)
  b_smoke ~ dnorm(0, 0.001)
  b_htn ~ dnorm(0, 0.001)
  b_diab ~ dnorm(0, 0.001)
  b_stroke ~ dnorm(0, 0.001)
  b_cad ~ dnorm(0, 0.001)


  # Precision parameters

  tau ~ dgamma(1, 1)
  tau_alpha ~ dgamma(1, 1)

  # Convert precision to standard deviation

  sigma <- 1 / sqrt(tau)
  sigma_alpha <- 1 / sqrt(tau_alpha)
}

"


# ---------------------------------------------------------
# 7. Fit the JAGS model
# ---------------------------------------------------------

jags_q3 <- jags.model(
  textConnection(model_q3),
  data = data_q3,
  n.adapt = 3000,
  n.chains = 4
)

# Burn-in
update(
  jags_q3,
  3000
)


# ---------------------------------------------------------
# 8. Obtain posterior samples
# ---------------------------------------------------------

samples_q3 <- coda.samples(
  jags_q3,

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
# 9. Model summary
# ---------------------------------------------------------

summary(samples_q3)


# ---------------------------------------------------------
# 10. Gelman-Rubin convergence diagnostics
# ---------------------------------------------------------

gelman.diag(samples_q3)


# ---------------------------------------------------------
# 11. Combine posterior samples
# ---------------------------------------------------------

samples_q3_combined <- as.mcmc(
  do.call(rbind, samples_q3)
)


# ---------------------------------------------------------
# 12. Trace and density plots
# ---------------------------------------------------------

png(
  "q3_trace_density.png",
  width = 1200,
  height = 800,
  res = 150
)

par(
  mar = c(2, 2, 2, 1)
)

plot(
  samples_q3[, c("beta", "delta", "gamma")]
)

dev.off()


# ---------------------------------------------------------
# 13. Autocorrelation plots
# ---------------------------------------------------------

png(
  "q3_acf.png",
  width = 800,
  height = 900,
  res = 150
)

par(
  mfrow = c(3, 1),
  mar = c(4, 4, 2, 1)
)

acf(
  as.numeric(samples_q3_combined[, "beta"]),
  main = "",
  ylab = "beta"
)

acf(
  as.numeric(samples_q3_combined[, "delta"]),
  main = "",
  ylab = "delta"
)

acf(
  as.numeric(samples_q3_combined[, "gamma"]),
  main = "",
  ylab = "gamma"
)

dev.off()
