# Carica librerie necessarie
library(coda)
library(ggplot2)

# Carica il dataset
data <- read.csv("~/Downloads/lowest_freq.csv")

# Parametri iniziali
set.seed(123)
n_iter <- 5000  # Numero di iterazioni MCMC
burn_in <- 1000  # Iterazioni da scartare
lambda_chain <- numeric(n_iter)
mu1_chain <- numeric(n_iter)
sigma1_chain <- numeric(n_iter)
mu2_chain <- numeric(n_iter)
sigma2_chain <- numeric(n_iter)

# Inizializzazione
lambda_chain[1] <- 0.5
bin_centers <- (data$start_point + data$end_point) / 2
log_bin_centers <- log(bin_centers)



mu1_chain[1] <- mean(log_bin_centers) - 0.5
sigma1_chain[1] <- sd(log_bin_centers)
mu2_chain[1] <- mean(log_bin_centers) + 0.5
sigma2_chain[1] <- sd(log_bin_centers)

# Gibbs Sampling
for (i in 2:n_iter) {
  dens1 <- plnorm(data$end_point, meanlog = mu1_chain[i-1], sdlog = sigma1_chain[i-1]) - 
    plnorm(data$start_point, meanlog = mu1_chain[i-1], sdlog = sigma1_chain[i-1])
  dens2 <- plnorm(data$end_point, meanlog = mu2_chain[i-1], sdlog = sigma2_chain[i-1]) - 
    plnorm(data$start_point, meanlog = mu2_chain[i-1], sdlog = sigma2_chain[i-1])
  tau <- lambda_chain[i-1] * dens1 / (lambda_chain[i-1] * dens1 + (1 - lambda_chain[i-1]) * dens2)
  
  lambda_chain[i] <- rbeta(1, 1 + sum(tau * data$count), 1 + sum((1 - tau) * data$count))
  mu1_chain[i] <- rnorm(1, mean = sum(tau * log_bin_centers * data$count) / sum(tau * data$count), sd = sqrt(1 / sum(tau * data$count)))
  mu2_chain[i] <- rnorm(1, mean = sum((1 - tau) * log_bin_centers * data$count) / sum((1 - tau) * data$count), sd = sqrt(1 / sum((1 - tau) * data$count)))
  sigma1_chain[i] <- sqrt(1 / rgamma(1, shape = sum(tau * data$count) / 2, rate = sum(tau * data$count * (log_bin_centers - mu1_chain[i])^2) / 2))
  sigma2_chain[i] <- sqrt(1 / rgamma(1, shape = sum((1 - tau) * data$count) / 2, rate = sum((1 - tau) * data$count * (log_bin_centers - mu2_chain[i])^2) / 2))
}

# Rimuoviamo il burn-in
lambda_chain <- lambda_chain[(burn_in+1):n_iter]
mu1_chain <- mu1_chain[(burn_in+1):n_iter]
sigma1_chain <- sigma1_chain[(burn_in+1):n_iter]
mu2_chain <- mu2_chain[(burn_in+1):n_iter]
sigma2_chain <- sigma2_chain[(burn_in+1):n_iter]

# Stime finali
lambda_est <- mean(lambda_chain)
mu1_est <- mean(mu1_chain)
sigma1_est <- mean(sigma1_chain)
mu2_est <- mean(mu2_chain)
sigma2_est <- mean(sigma2_chain)

# Genera dati simulati
simulated <- ifelse(runif(sum(data$count)) < lambda_est,
                    rlnorm(sum(data$count), meanlog = mu1_est, sdlog = sigma1_est),
                    rlnorm(sum(data$count), meanlog = mu2_est, sdlog = sigma2_est))

# ---- Test Chi-quadrato sui dati reali ----
range_all <- range(simulated)
breaks <- seq(range_all[1], range_all[2], length.out = 10)
observed <- hist(simulated, breaks = breaks, plot = FALSE)$counts

expected_probs <- diff(plnorm(breaks, meanlog = mu1_est, sdlog = sigma1_est) * lambda_est +
                         plnorm(breaks, meanlog = mu2_est, sdlog = sigma2_est) * (1 - lambda_est))
expected <- expected_probs * length(simulated)
chi2_stat_original <- sum((observed - expected)^2 / expected)

# ---- Bootstrap con B = 1000 ----
set.seed(79)
B <- 1000
chi2_bootstrap <- numeric(B)

for (b in 1:B) {
  boot_sample <- ifelse(runif(length(simulated)) < lambda_est,
                        rlnorm(length(simulated), meanlog = mu1_est, sdlog = sigma1_est),
                        rlnorm(length(simulated), meanlog = mu2_est, sdlog = sigma2_est))
  breaks_boot <- seq(min(c(simulated, boot_sample)), max(c(simulated, boot_sample)), length.out = 10)
  boot_observed <- hist(boot_sample, breaks = breaks_boot, plot = FALSE)$counts
  chi2_bootstrap[b] <- sum((boot_observed - expected)^2 / expected)
}

p_value <- mean(chi2_bootstrap >= chi2_stat_original)

# Stampa risultati
cat("Statistica chi2 osservata:", chi2_stat_original, "\n")
cat("P-value bootstrap:", p_value, "\n")

# ---- Istogramma distribuzione bootstrap della statistica Chi2 ----
df_chi2 <- data.frame(chi2_bootstrap = chi2_bootstrap)

ggplot(df_chi2, aes(x = chi2_bootstrap)) +
  geom_histogram(binwidth = 1, fill = "black", alpha = 0.5, color = "black") +
  geom_vline(xintercept = chi2_stat_original, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Bootstrap Distribution of Chi-Square Statistic",
       x = "Chi-Square Statistic",
       y = "Frequency") +
  theme_minimal()


# ---- Test KS sui dati reali ----
ks_stat_original <- ks.test(simulated, 
                            function(x) lambda_est * plnorm(x, meanlog = mu1_est, sdlog = sigma1_est) + 
                              (1 - lambda_est) * plnorm(x, meanlog = mu2_est, sdlog = sigma2_est))$statistic

# ---- Bootstrap KS con B = 1000 ----
set.seed(42)
ks_bootstrap <- numeric(B)

for (b in 1:B) {
  boot_sample <- ifelse(runif(length(simulated)) < lambda_est,
                        rlnorm(length(simulated), meanlog = mu1_est, sdlog = sigma1_est),
                        rlnorm(length(simulated), meanlog = mu2_est, sdlog = sigma2_est))
  ks_bootstrap[b] <- ks.test(boot_sample, 
                             function(x) lambda_est * plnorm(x, meanlog = mu1_est, sdlog = sigma1_est) + 
                               (1 - lambda_est) * plnorm(x, meanlog = mu2_est, sdlog = sigma2_est))$statistic
}

# ---- Calcola p-value ----
p_value_ks <- mean(ks_bootstrap >= ks_stat_original)

# Stampa risultati
cat("Statistica KS osservata:", ks_stat_original, "\n")
cat("P-value bootstrap KS:", p_value_ks, "\n")


# ---- Istogramma distribuzione bootstrap della statistica KS ----
df_ks <- data.frame(ks_bootstrap = ks_bootstrap)

ggplot(df_ks, aes(x = ks_bootstrap)) +
  geom_histogram(binwidth = 0.0001, fill = "black", color = "black") +
  geom_vline(xintercept = ks_stat_original, color = "red", linetype = "dashed", size = 1.5) +
  labs(title = "Bootstrap Distribution of KS Statistic",
       x = "KS Statistic",
       y = "Frequency") +
  theme_minimal()

library(ggplot2)

# Genera quantili empirici
empirical_quantiles <- quantile(simulated, probs = ppoints(length(simulated)))

# Genera quantili teorici dalla mixture
set.seed(123)
theoretical_samples <- ifelse(runif(length(simulated)) < lambda_est,
                              rlnorm(length(simulated), meanlog = mu1_est, sdlog = sigma1_est),
                              rlnorm(length(simulated), meanlog = mu2_est, sdlog = sigma2_est))

theoretical_quantiles <- quantile(theoretical_samples, probs = ppoints(length(simulated)))

# Q-Q Plot con ggplot2
df_qq <- data.frame(Theoretical = theoretical_quantiles, Empirical = empirical_quantiles)

ggplot(df_qq, aes(x = Theoretical, y = Empirical)) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Q-Q Plot: Data vs Simulated Fitted Mixture",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_minimal()


# ---- Istogramma dati vs Densità stimata ----
x_seq <- seq(min(simulated), max(simulated), length.out = 1000)
y_density <- lambda_est * dlnorm(x_seq, meanlog = mu1_est, sdlog = sigma1_est) +
  (1 - lambda_est) * dlnorm(x_seq, meanlog = mu2_est, sdlog = sigma2_est)
df_curve <- data.frame(x = x_seq, y = y_density)

ggplot() +
  geom_histogram(data = data,
                 aes(x = (start_point + end_point) / 2, weight = count, y = ..density..),
                 binwidth = 10, fill = "blue", alpha = 0.5, color = "black") +
  geom_line(data = df_curve, aes(x = x, y = y), color = "red", size = 1) +
  labs(title = "Histogram of Data vs Estimated Bi-lognormal Density",
       x = "Fundamental Frequency (Hz)",
       y = "Density") +
  theme_minimal()


