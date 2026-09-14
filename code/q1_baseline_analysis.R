# Q1: Baseline Characteristics
# Bayesian Cognitive Aging Project
#
# Purpose:
# Compare baseline characteristics between centenarian offspring
# and control participants and identify potential confounders
# for subsequent analyses.

# ---------------------------------------------------------
# 1. Load data and packages
# ---------------------------------------------------------

library(gtsummary)
library(dplyr)

data <- read.csv("tics.data.2025.csv")

TICS <- data

# Convert group indicator to factor
TICS$ptype <- as.factor(TICS$ptype)
levels(TICS$ptype) <- c("Control", "Centenarian Offspring")


# ---------------------------------------------------------
# 2. Examine continuous variable distributions
# ---------------------------------------------------------

hist(TICS$BMI,
     main = "Distribution of BMI",
     xlab = "BMI")

hist(TICS$Years.of.Education,
     main = "Distribution of Years of Education",
     xlab = "Years of Education")

boxplot(Age01 ~ ptype,
        data = TICS,
        main = "Baseline Age by Study Group",
        xlab = "Study Group",
        ylab = "Age")

qqnorm(TICS$Age01)
qqline(TICS$Age01)


# ---------------------------------------------------------
# 3. Descriptive statistics
# ---------------------------------------------------------

vars_to_include <- c(
  "BMI",
  "Age01",
  "sex",
  "SH.Ever.Smoked.",
  "MC.Aspirin",
  "MC.Stroke",
  "MC.Diabetes.Mellitus",
  "MC.HTN",
  "MC.Coronary.Artery.Disease",
  "MC.Cancer",
  "MC.Heart.Attack",
  "Years.of.Education",
  "TICS01"
)

table1 <- TICS %>%
  tbl_summary(
    by = ptype,
    include = vars_to_include,
    type = list(
      BMI ~ "continuous",
      Age01 ~ "continuous",
      Years.of.Education ~ "continuous",
      sex ~ "categorical",
      SH.Ever.Smoked. ~ "categorical",
      MC.Aspirin ~ "categorical",
      MC.Stroke ~ "categorical",
      MC.Diabetes.Mellitus ~ "categorical",
      MC.HTN ~ "categorical",
      MC.Coronary.Artery.Disease ~ "categorical",
      MC.Cancer ~ "categorical",
      MC.Heart.Attack ~ "categorical"
    ),
    statistic = list(
      BMI ~ "{median} ({p25}, {p75})",
      Age01 ~ "{mean} ({sd})",
      Years.of.Education ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  add_p(
    test = list(
      BMI ~ "wilcox.test",
      Age01 ~ "t.test",
      Years.of.Education ~ "t.test",
      all_categorical() ~ "chisq.test"
    )
  ) %>%
  modify_header(
    p.value ~ "**p-value**"
  )

table1


# ---------------------------------------------------------
# 4. Additional descriptive comparisons
# ---------------------------------------------------------

# BMI
tapply(TICS$BMI, TICS$ptype, median, na.rm = TRUE)

tapply(
  TICS$BMI,
  TICS$ptype,
  quantile,
  probs = c(0.25, 0.75),
  na.rm = TRUE
)

wilcox.test(BMI ~ ptype, data = TICS)


# Years of education
tapply(
  TICS$Years.of.Education,
  TICS$ptype,
  mean,
  na.rm = TRUE
)

t.test(
  Years.of.Education ~ ptype,
  data = TICS
)


# ---------------------------------------------------------
# 5. Categorical comparisons
# ---------------------------------------------------------

table(TICS$MC.Heart.Attack, TICS$ptype)

prop.table(
  table(TICS$MC.Heart.Attack, TICS$ptype),
  2
) * 100

chisq.test(
  table(TICS$MC.Heart.Attack, TICS$ptype)
)


# ---------------------------------------------------------
# 6. Preliminary multivariable assessment
# ---------------------------------------------------------

baseline_model <- lm(
  TICS01 ~ ptype +
    BMI +
    Age01 +
    Years.of.Education +
    sex +
    SH.Ever.Smoked. +
    MC.Stroke +
    MC.Diabetes.Mellitus +
    MC.HTN +
    MC.Coronary.Artery.Disease +
    MC.Heart.Attack,
  data = TICS
)

summary(baseline_model)


# Optional multicollinearity assessment
library(car)

vif(baseline_model)
