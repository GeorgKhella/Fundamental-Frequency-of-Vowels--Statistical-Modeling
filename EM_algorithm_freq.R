library(ggplot2)

# Carica dataset
data <- read.csv("~/Downloads/lowest_freq.csv")

# Jittering
set.seed(123)
jittered <- unlist(lapply(1:nrow(data), function(i) {
  runif(data$count[i], min = data$start_point[i], max = data$end_point[i])
}))

# Funzione EM
EM_algorithm <- function(data, tol=1e-6, max_iter=1000) {
  n <- length(data)
  lambda <- 0.5
  mu1 <- mean(log(data)) - 0.5
  sigma1 <- sd(log(data))
  mu2 <- mean(log(data)) + 0.5
  sigma2 <- sd(log(data))
  loglik_old <- -Inf
  iter <- 0
  
  while(iter < max_iter) {
    iter <- iter + 1
    dens1 <- dlnorm(data, meanlog = mu1, sdlog = sigma1)
    dens2 <- dlnorm(data, meanlog = mu2, sdlog = sigma2)
    tau <- lambda * dens1 / (lambda * dens1 + (1 - lambda) * dens2)
    
    lambda_new <- mean(tau)
    mu1_new <- sum(tau * log(data)) / sum(tau)
    mu2_new <- sum((1 - tau) * log(data)) / sum(1 - tau)
    sigma1_new <- sqrt(sum(tau * (log(data) - mu1_new)^2) / sum(tau))
    sigma2_new <- sqrt(sum((1 - tau) * (log(data) - mu2_new)^2) / sum(1 - tau))
    
    dens_mix <- lambda_new * dens1 + (1 - lambda_new) * dens2
    loglik_new <- sum(log(dens_mix + 1e-10))
    
    if(abs(loglik_new - loglik_old) < tol) break
    
    lambda <- lambda_new
    mu1 <- mu1_new
    mu2 <- mu2_new
    sigma1 <- sigma1_new
    sigma2 <- sigma2_new
    loglik_old <- loglik_new
  }
  
  list(lambda = lambda, mu1 = mu1, sigma1 = sigma1, mu2 = mu2, sigma2 = sigma2, 
       iterations = iter, logLik = loglik_new)
}

# Stima EM
em_results <- EM_algorithm(jittered)

# ---- Test Chi-quadrato sui dati reali ----
# Espandiamo i break per evitare errori di binning
range_all <- range(jittered)
breaks <- seq(range_all[1], range_all[2], length.out = 10)
observed <- hist(jittered, breaks = breaks, plot = FALSE)$counts

# Frequenze attese
expected_probs <- diff(plnorm(breaks, meanlog = em_results$mu1, sdlog = em_results$sigma1) * em_results$lambda +
                         plnorm(breaks, meanlog = em_results$mu2, sdlog = em_results$sigma2) * (1 - em_results$lambda))
expected <- expected_probs * length(jittered)

# Statistica chi2 osservata
chi2_stat_original <- sum((observed - expected)^2 / expected)

# ---- Bootstrap con B = 1000 ----
set.seed(42)
B <- 1000
chi2_bootstrap <- numeric(B)

for (b in 1:B) {
  # Genera campione bootstrap
  boot_sample <- ifelse(runif(length(jittered)) < em_results$lambda,
                        rlnorm(length(jittered), meanlog = em_results$mu1, sdlog = em_results$sigma1),
                        rlnorm(length(jittered), meanlog = em_results$mu2, sdlog = em_results$sigma2))
  
  # Espandi breaks per includere range bootstrap
  breaks_boot <- seq(min(c(jittered, boot_sample)), max(c(jittered, boot_sample)), length.out = 10)
  
  # Frequenze bootstrap
  boot_observed <- hist(boot_sample, breaks = breaks_boot, plot = FALSE)$counts
  
  # Calcola la statistica chi2 nel bootstrap
  chi2_bootstrap[b] <- sum((boot_observed - expected)^2 / expected)
}

# ---- Calcola p-value ----
p_value <- mean(chi2_bootstrap >= chi2_stat_original)

# Stampa risultati
cat("Statistica chi2 osservata:", chi2_stat_original, "\n")
cat("P-value bootstrap:", p_value, "\n")

# ---- Istogramma distribuzione bootstrap della statistica chi2 ----
df_chi2 <- data.frame(chi2_bootstrap = chi2_bootstrap)

ggplot(df_chi2, aes(x = chi2_bootstrap)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.5, color = "black") +
  geom_vline(xintercept = chi2_stat_original, color = "red", linetype = "dashed", size = 1.2) +
  labs(title = "Bootstrap Distribution of Chi-Square Statistic",
       x = "Chi-Square Statistic",
       y = "Frequency") +
  theme_minimal()


x_seq <- seq(min(jittered), max(jittered), length.out = 1000)
y_density <- em_results$lambda * dlnorm(x_seq, meanlog = em_results$mu1, sdlog = em_results$sigma1) +
  (1 - em_results$lambda) * dlnorm(x_seq, meanlog = em_results$mu2, sdlog = em_results$sigma2)
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

# ---- Plot 2: Q-Q Plot ----
set.seed(123)
n <- length(jittered)
simulated <- ifelse(runif(n) < em_results$lambda,
                    rlnorm(n, meanlog = em_results$mu1, sdlog = em_results$sigma1),
                    rlnorm(n, meanlog = em_results$mu2, sdlog = em_results$sigma2))

df_qq <- data.frame(theoretical = sort(simulated),
                    sample = sort(jittered))

ggplot(df_qq, aes(x = theoretical, y = sample)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Q-Q Plot: Data vs Simulated Fitted Mixture",
       x = "Theoretical Quantiles",
       y = "Empirical Quantiles") +
  theme_minimal()

em_results



# ---- Test KS sui dati reali ----
ks_stat_original <- ks.test(jittered, 
                            function(x) em_results$lambda * plnorm(x, meanlog = em_results$mu1, sdlog = em_results$sigma1) + 
                              (1 - em_results$lambda) * plnorm(x, meanlog = em_results$mu2, sdlog = em_results$sigma2))$statistic

# ---- Bootstrap KS con B = 1000 ----
set.seed(42)
ks_bootstrap <- numeric(B)

for (b in 1:B) {
  # Genera campione bootstrap
  boot_sample <- ifelse(runif(length(jittered)) < em_results$lambda,
                        rlnorm(length(jittered), meanlog = em_results$mu1, sdlog = em_results$sigma1),
                        rlnorm(length(jittered), meanlog = em_results$mu2, sdlog = em_results$sigma2))
  
  # Calcola KS statistic nel bootstrap
  ks_bootstrap[b] <- ks.test(boot_sample, 
                             function(x) em_results$lambda * plnorm(x, meanlog = em_results$mu1, sdlog = em_results$sigma1) + 
                               (1 - em_results$lambda) * plnorm(x, meanlog = em_results$mu2, sdlog = em_results$sigma2))$statistic
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

