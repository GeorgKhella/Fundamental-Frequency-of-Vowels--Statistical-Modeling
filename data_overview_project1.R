# Caricamento delle librerie necessarie
library(ggplot2)
library(dplyr)
library(readr)

# Caricamento del dataset
data <- read.csv("~/Downloads/lowest_freq.csv")
summary(data)
# Definizione parametri per la miscela di due lognormali (esempio arbitrario)
p_mix   <- 0.3
mu1     <- log(130)
sigma1  <- 0.1
mu2     <- log(220)
sigma2  <- 0.1

# Funzione per la miscela di due lognormali
d2lnorm <- function(x, p, mu1, sigma1, mu2, sigma2) {
  p * dlnorm(x, meanlog = mu1, sdlog = sigma1) +
    (1 - p) * dlnorm(x, meanlog = mu2, sdlog = sigma2)
}

# Definizione del range e calcolo della densità
x_vals      <- seq(1, 300, length.out = 1000)
total_count <- sum(data$count)
binw        <- 10
y_vals      <- d2lnorm(x_vals, p_mix, mu1, sigma1, mu2, sigma2) * total_count * binw

# Data frame per la curva bi-lognormale
df_curve <- data.frame(x = x_vals, y = y_vals)

# Plot dell'istogramma e sovrapposizione della curva bi-lognormale,
# rimuovendo la parte prima del primo bin tramite i limiti sull'asse x
ggplot() +
  geom_histogram(data = data, 
                 aes(x = (start_point + end_point) / 2, weight = count),
                 binwidth = 10, fill = "blue", alpha = 0.5, color = "black") +
  geom_line(data = df_curve, aes(x = x, y = y), color = "red", size = 1) +
  scale_x_continuous(limits = c(60, 300)) +
  labs(title = "Histogram of Fundamental Frequency (f0) vs Bi-lognormal density",
       x = "Fundamental Frequency (Hz)",
       y = "Count") +
  theme_minimal()


# Istogramma delle frequenze fondamentali (f0) basato sui bin
ggplot(data, aes(x = (start_point + end_point) / 2, weight = count)) +
  geom_histogram(binwidth = 10, fill = "blue", alpha = 0.5, color = "black") +
  labs(title = "Histogram of Fundamental Frequency (f0)",
       x = "Fundamental Frequency (Hz)",
       y = "Count") +
  theme_minimal()