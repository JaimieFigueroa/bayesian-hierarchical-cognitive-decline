# Q4: Bayesian Latent Trajectory Mixture Model
# Identifying subgroups with different rates of cognitive change
#
# Outcome: TICS02-TICS05
# Baseline cognition: TICS01
# Model: Two-class Bayesian latent trajectory mixture model
# Software: R + JAGS

library(rjags)
library(coda)


# ---------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------

dat <- read.csv("tics.data.2025.csv")


# ---------------------------------------------------------
# 2. Calculate average follow-up time
# ---------------------------------------------------------

Time <- numeric(4)

Time[1] <- mean(
  dat$Age02 - dat$Age01,
  na.rm = TRUE
)

Time[2] <- Time[1] +
  mean(
    dat$Age03 - dat$Age02,
    na.rm = TRUE
  )

Time[3] <- Time[2] +
  mean(
    dat$Age04 - dat$Age03,
    na.rm = TRUE
  )

Time[4] <- Time[3] +
  mean(
    dat$Age05 - dat$Age04,
    na.rm = TRUE
  )

Time


# ---------------------------------------------------------
# 3. Complete-case analysis
# ---------------------------------------------------------

use <- complete.cases(
  dat[, c(
    "TICS01",
    "TICS02",
    "TICS03",
    "TICS04",
    "TICS05"
  )]
)


# Longitudinal outcome matrix
y <- as.matrix(
  dat[use, c(
    "TICS02",
    "TICS03",
    "TICS04",
    "TICS05"
  )]
)


# Baseline TICS
x0 <- dat$TICS01[use]


# Number of participants
n.subj <- nrow(y)


# ---------------------------------------------------------
# 4. Visualize individual trajectories
# ---------------------------------------------------------

par(
  mfrow = c(1, 1),
  mar = c(3, 3, 1.5, 0.5),
  mgp = c(1.6, 0.5, 0),
  cex.lab = 0.9,
  cex.axis = 0.8,
  tck = -0.02
)

plot(
  c(min(Time), max(Time)),
  range(y),
  type = "n",
  xlab = "Time",
  ylab = "TICS"
)

for (i in 1:n.subj) {

  lines(
    Time,
    y[i, ],
    col = rgb(0, 0, 0, 0.2)
  )
}


# ---------------------------------------------------------
# 5. Bayesian latent trajectory mixture model
# ---------------------------------------------------------

model_string <- "

model {

  # Assign each participant to one of two latent classes
  for (i in 1:n.subj) {

    w[i] ~ dcat(theta[])

    # Repeated TICS measurements
    for (j in 1:4) {

      y[i,j] ~ dnorm(mu[i,j], tau)

      mu[i,j] <-
          b0[w[i]]
        + g0[w[i]] * x0[i]
        + b1[w[i]] * Time[j]
    }
  }


  # Class probabilities
  theta[1:2] ~ ddirch(alpha[])


  # Class-specific parameters
  for (k in 1:2) {

    # Class-specific slope
    b1[k] ~ dnorm(0, 0.01)

    # Effect of baseline TICS
    g0[k] ~ dnorm(0, 0.001)

    # Dirichlet prior parameters
    alpha[k] <- 1
  }


  # Class-specific intercepts
  #
  # Identifiability constraint:
  # Class 2 intercept is defined relative to Class 1

  b0[1] ~ dnorm(0, 0.001)

  b0[2] <- b0[1] + delta

  delta ~ dgamma(1, 1)


  # Residual precision
  tau ~ dgamma(1, 1)

}

"


# ---------------------------------------------------------
# 6. Prepare data for JAGS
# ---------------------------------------------------------

jags_data <- list(
  y = y,
  x0 = x0,
  Time = Time,
  n.subj = n.subj
)


# ---------------------------------------------------------
# 7. Fit model
# ---------------------------------------------------------

set.seed(123)

mod <- jags.model(
  textConnection(model_string),
  data = jags_data,
  n.chains = 3,
  n.adapt = 2000
)


# ---------------------------------------------------------
# 8. Burn-in
# ---------------------------------------------------------

update(
  mod,
  5000
)


# ---------------------------------------------------------
# 9. Posterior sampling
# ---------------------------------------------------------

samp <- coda.samples(
  mod,

  variable.names = c(
    "b0",
    "b1",
    "g0",
    "theta",
    "tau",
    "w"
  ),

  n.iter = 10000
)


# ---------------------------------------------------------
# 10. Posterior summary
# ---------------------------------------------------------

summary(
  samp[, c(
    "b0[1]",
    "b0[2]",
    "b1[1]",
    "b1[2]",
    "g0[1]",
    "g0[2]",
    "theta[1]",
    "theta[2]",
    "tau"
  )]
)


# ---------------------------------------------------------
# 11. Posterior diagnostic plots
# ---------------------------------------------------------

plot(
  samp[, c(
    "b0[1]",
    "b0[2]",
    "b1[1]",
    "b1[2]",
    "theta[1]",
    "theta[2]"
  )]
)


# ---------------------------------------------------------
# 12. Estimate latent class membership
# ---------------------------------------------------------

mat <- as.matrix(samp)

w_cols <- grep(
  "^w\\[",
  colnames(mat)
)

w_post <- mat[, w_cols]


# Posterior median class assignment
group <- apply(
  w_post,
  2,
  median
)


# Class counts
table(group)


# Class proportions
prop.table(
  table(group)
)


# ---------------------------------------------------------
# 13. Visualize trajectories by latent class
# ---------------------------------------------------------

par(
  mfrow = c(1, 1),
  mar = c(3, 3, 1.5, 0.5),
  mgp = c(1.6, 0.5, 0),
  cex.lab = 0.9,
  cex.axis = 0.8,
  tck = -0.02
)

plot(
  c(min(Time), max(Time)),
  range(y),
  type = "n",
  xlab = "Time",
  ylab = "TICS"
)


# Class-specific colors
cols <- ifelse(
  group == 1,
  rgb(0, 0, 1, 0.2),
  rgb(1, 0, 0, 0.2)
)


# Individual trajectories
for (i in 1:n.subj) {

  lines(
    Time,
    y[i, ],
    col = cols[i]
  )
}


# Mean trajectory for Class 1
mean1 <- colMeans(
  y[group == 1, , drop = FALSE]
)


# Mean trajectory for Class 2
mean2 <- colMeans(
  y[group == 2, , drop = FALSE]
)


# Add class-specific mean trajectories
lines(
  Time,
  mean1,
  col = "blue",
  lwd = 3
)

lines(
  Time,
  mean2,
  col = "red",
  lwd = 3
)


legend(
  "topright",
  legend = c(
    "Class 1",
    "Class 2"
  ),
  col = c(
    "blue",
    "red"
  ),
  lty = 1,
  lwd = 3,
  bty = "n"
)


# ---------------------------------------------------------
# 14. Posterior membership probabilities
# ---------------------------------------------------------

prob_class1 <- colMeans(
  w_post == 1
)

prob_class2 <- colMeans(
  w_post == 2
)


membership <- data.frame(
  group = group,
  prob_class1 = prob_class1,
  prob_class2 = prob_class2
)

head(membership)


# ---------------------------------------------------------
# 15. Convergence diagnostics
# ---------------------------------------------------------

params_conv1 <- c(
  "b0[1]",
  "b0[2]",
  "b1[1]",
  "b1[2]"
)

params_conv2 <- c(
  "g0[1]",
  "g0[2]",
  "theta[1]",
  "tau"
)


# Geweke diagnostics
geweke.diag(
  samp[, params_conv1],
  frac1 = 0.1,
  frac2 = 0.5
)

geweke.diag(
  samp[, params_conv2],
  frac1 = 0.1,
  frac2 = 0.5
)


# Gelman-Rubin diagnostics
gelman.diag(
  samp[, params_conv1],
  multivariate = FALSE
)

gelman.diag(
  samp[, params_conv2],
  multivariate = FALSE
)


# Gelman-Rubin diagnostic plots
par(
  mar = c(2, 2, 1, 1)
)

gelman.plot(
  samp[, params_conv1]
)

gelman.plot(
  samp[, params_conv2]
)
