# Fundamental Frequency of Vowels – Statistical Modeling

Project for Applied Statistics (Feb 2025)  
Author: Georg Khella

## 🎯 Objective

To model the **fundamental frequency (f0)** of American English vowels using a **bi-lognormal distribution**.  
Approaches used:
- Expectation-Maximization (EM) algorithm with jittering
- Bayesian inference with MCMC (Gibbs sampling)
- Goodness-of-fit validation with bootstrap and statistical tests

## 📊 Dataset

- Based on Hillebrand et al. (1995)
- Binned acoustic measurements of vowel frequencies (f0) in Hz
- Data includes counts and percentages for children, male, and female speakers

## 🔍 Methods

- **EM Algorithm**: Handles binned data via jittering, estimates mixture model parameters.
- **Bayesian Approach**: Full posterior estimation using MCMC (Gibbs), accounts for parameter uncertainty.
- **Model Validation**:
  - Chi-square and Kolmogorov–Smirnov tests
  - Parametric bootstrap
  - Q–Q plots and histogram overlays

## ✅ Results

- Both methods fit the bi-lognormal model well.
- Bayesian inference shows **slightly better fit**, especially in tail behavior.
- Validates bi-lognormal as a solid model for vowel pitch distribution.

## 📁 Files

- `Georg_Khella_Project1.pdf`: Full report with theory, methods, and analysis

## 📚 References

Includes work from Gelman, Rubin, Silverman, McLachlan, Wasserman, and more (see report for full list).
