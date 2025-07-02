# Pneumonia_in_babies
**Project Title**: Investigating Pneumonia Occurrence Rate in Babies 

---

## 📌 Project Overview

This project explores the relationship between several maternal and neonatal factors and the **occurrence of pneumonia in infants**. Specifically, the study investigates whether a mother's age, alcohol use during pregnancy, cigarette use, the baby's birthweight, and race affect how long a baby stays in the hospital for pneumonia.

I used **survival analysis** techniques such as the **Kaplan-Meier estimator**, **Nelson-Aalen estimator**, and the **Cox Proportional Hazards (PH) model** to analyze time-to-event data (time in hospital for pneumonia).

---

## 📊 Dataset

- **Source**: `pneumon` dataset from the R package `KMsurv`
- **Total Observations**: 3,470
- **Relevant Variables**:
  - `hospital`: Hospital ID
  - `mthage`: Age of the mother
  - `alcohol`: Alcohol use during pregnancy
  - `smoke`: Cigarette use during pregnancy
  - `bweight`: Whether the baby had a normal birthweight
  - `race`: Race of the mother
  - `agepn`: Age child is hospitalized for pneumonia (time-to-event)

---

## ❓ Research Question

> Does the mother’s age, alcohol use during pregnancy, baby’s birthweight, and race have an effect on pneumonia occurrence rate (i.e., duration in hospital) in babies?

---

## 🔬 Methodology

1. **Exploratory Data Analysis**: Identify variable types and missingness.
2. **Kaplan-Meier Estimate**: Estimate survival probabilities by covariate groups.
3. **Nelson-Aalen Estimate**: Estimate cumulative hazard.
4. **Cox Proportional Hazards Model**: Evaluate the effect of covariates on survival time.
5. **Backward Selection**: Model simplification using Wald statistics and AIC.
6. **Model Diagnostics**: 
   - Cox-Snell residuals
   - Martingale residuals
   - Assumption checks using time-dependent covariates

---

## 📈 Key Findings

- **Significant Variables**:
  - **Cigarette use during pregnancy** (HR = 2.14, *p* = 0.001)
  - **Normal birthweight** (HR = 0.02, *p* = 0.042)
  - **Mother's age at birth** (HR = 0.87, *p* = 0.048)
  - **Interaction (Mother’s age × Birthweight)** (HR = 1.24, *p* = 0.021)

- The final Cox PH model includes:
  ```r
  Surv(agepn, hospital) ~ mthage + factor(smoke) + factor(bweight) + mthage:factor(bweight)
