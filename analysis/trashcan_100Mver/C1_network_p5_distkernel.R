# ============================================================
# Movement Distribution Fitting
# Distance weight kernel for community contact network
# Data: Meta Data for Good - Movement Distribution Maps
# Proxy countries for DRC: AGO, BDI, COG, RWA, TZA, ZMB
# Model: Mixture of two exponentials (Stan)
# ============================================================

library(tidyverse)
library(rstan)
library(patchwork)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# ── 1. Load and prepare data ──────────────────────────────────

df <- read_csv("data/movement-distribution-maps_2026-04-01_2026-04-16.csv")

# DRC proxy countries (neighboring, similar mobility context)
drc_proxy <- c("COG", "AGO", "RWA", "BDI", "TZA", "ZMB")

# Average bin fractions across proxy countries
df_avg <- df %>%
  filter(country %in% drc_proxy) %>%
  group_by(home_to_ping_distance_category) %>%
  summarise(mean_fraction = mean(distance_category_ping_fraction, na.rm = TRUE),
            .groups = "drop")

# Exclude 0km bin (handled separately as household contacts)
# Normalize remaining 3 bins to sum to 1
df_avg_norm <- df_avg %>%
  filter(home_to_ping_distance_category != "0") %>%
  mutate(mean_fraction = mean_fraction / sum(mean_fraction))

# Observed proportions: [0-10km, 10-100km, 100+km]
obs <- df_avg_norm %>%
  arrange(match(home_to_ping_distance_category,
                c("(0, 10)", "[10, 100)", "100+"))) %>%
  pull(mean_fraction)

cat(sprintf("Observed proportions (community contacts only):\n"))
cat(sprintf("  0-10km   : %.3f\n", obs[1]))
cat(sprintf("  10-100km : %.3f\n", obs[2]))
cat(sprintf("  100+km   : %.3f\n", obs[3]))
cat(sprintf("  Sum      : %.3f\n", sum(obs)))


# ── 2. Stan model: mixture of two exponentials ────────────────
# f(d) = p * a1 * exp(-a1*d) + (1-p) * a2 * exp(-a2*d)
# Both components integrate to 1 over [0, inf)
# p     : mixing weight (local component)
# a1    : fast decay rate (local mobility)
# a2    : slow decay rate (long-distance mobility)

stan_code <- "
data {
  vector[3] obs;       // observed bin proportions (sum = 1)
  real N_pseudo;       // pseudo sample size for likelihood weighting
}

parameters {
  real<lower=0, upper=1> p;   // mixing weight (local component)
  real<lower=0, upper=1> a1;           // fast decay rate (local)
  real<lower=0, upper=1> a2;           // slow decay rate (long-distance)
}

transformed parameters {
  // Bin probabilities via integration of mixture PDF
  // int_{d1}^{d2} a*exp(-a*d) dd = exp(-a*d1) - exp(-a*d2)
  real p1 = p * (1 - exp(-10.0 * a1)) +
            (1-p) * (1 - exp(-10.0 * a2));                 // 0 to 10 km

  real p2 = p * (exp(-10.0 * a1) - exp(-100.0 * a1)) +
            (1-p) * (exp(-10.0 * a2) - exp(-100.0 * a2));  // 10 to 100 km

  real p3 = p * exp(-100.0 * a1) +
            (1-p) * exp(-100.0 * a2);                      // 100 km to inf

  simplex[3] bins = [p1, p2, p3]';
}

model {
  // Priors
  p  ~ beta(9, 1);                           // most weight on local component
  //a1 ~ lognormal(log(0.3), 0.5);             // fast decay: half-dist ~2 km
  //a2 ~ lognormal(log(0.01), 0.5);            // slow decay: half-dist ~70 km
  //a1 ~ normal(0.3, 0.5);             // fast decay: half-dist ~2 km
  //a2 ~ normal(0.01, 0.5);            // slow decay: half-dist ~70 km

  // Identifiability: a1 > a2 (local decay faster than long-distance)
  target += log(a1 > a2 ? 1 : 0);

  // Multinomial likelihood (pseudo-count weighted)
  target += N_pseudo * dot_product(obs, log(bins));
}

generated quantities {
  real half_dist1 = log(2) / a1;  // half-distance: local component (km)
  real half_dist2 = log(2) / a2;  // half-distance: long-distance component (km)
}
"

stan_data <- list(
  obs      = obs,
  N_pseudo = 1000
)

fit_mix <- stan(
  model_code = stan_code,
  data       = stan_data,
  iter       = 4000,
  warmup     = 2000,
  chains     = 4,
  seed       = 42,
  control    = list(adapt_delta = 0.95)
)


# ── 3. Results ────────────────────────────────────────────────

print(fit_mix,
      pars  = c("p", "a1", "a2", "p1", "p2", "p3", "half_dist1", "half_dist2"),
      probs = c(0.025, 0.5, 0.975))

post <- as.data.frame(fit_mix)

cat(sprintf("\nPosterior summary:\n"))
cat(sprintf("  p  (local weight)  : %.3f [%.3f, %.3f]\n",
            median(post$p), quantile(post$p, 0.025), quantile(post$p, 0.975)))
cat(sprintf("  a1 (local decay)   : %.4f [%.4f, %.4f]  half=%.1f km\n",
            median(post$a1), quantile(post$a1, 0.025), quantile(post$a1, 0.975),
            log(2)/median(post$a1)))
cat(sprintf("  a2 (long-dist)     : %.4f [%.4f, %.4f]  half=%.1f km\n",
            median(post$a2), quantile(post$a2, 0.025), quantile(post$a2, 0.975),
            log(2)/median(post$a2)))

cat(sprintf("\nPredicted: 0-10: %.3f | 10-100: %.3f | 100+: %.3f\n",
            mean(post$p1), mean(post$p2), mean(post$p3)))
cat(sprintf("Observed:  0-10: %.3f | 10-100: %.3f | 100+: %.3f\n",
            obs[1], obs[2], obs[3]))


# ── 4. Plots ──────────────────────────────────────────────────

dist_seq <- seq(0, 200, by = 0.5)
set.seed(42)
idx <- sample(nrow(post), 2000)

# Posterior predictive density curves
curves_df <- map_dfr(idx, function(i) {
  pp <- post$p[i]
  a1 <- post$a1[i]
  a2 <- post$a2[i]
  tibble(
    distance = dist_seq,
    density  = pp * a1 * exp(-a1 * dist_seq) +
      (1 - pp) * a2 * exp(-a2 * dist_seq)
  )
})

summary_df <- curves_df %>%
  group_by(distance) %>%
  summarise(
    median = median(density),
    lo     = quantile(density, 0.025),
    hi     = quantile(density, 0.975),
    .groups = "drop"
  )

# Observed density per bin (fraction / bin_width)
obs_df <- df_avg_norm %>%
  arrange(match(home_to_ping_distance_category,
                c("(0, 10)", "[10, 100)", "100+"))) %>%
  mutate(
    distance    = c(5, 55, 150),
    bin_width   = c(10, 90, 400),
    obs_density = mean_fraction / bin_width
  )



# Plot 1: posterior predictive check per bin
bin_post_df <- post[idx, ] %>%
  select(p1, p2, p3) %>%
  pivot_longer(everything(),
               names_to  = "bin",
               values_to = "prob") %>%
  mutate(
    bin = recode(bin, p1 = "0-10km", p2 = "10-100km", p3 = "100+km"),
    bin = factor(bin, levels = c("0-10km", "10-100km", "100+km"))
  )

obs_lines <- tibble(
  bin      = factor(c("0-10km", "10-100km", "100+km"),
                    levels = c("0-10km", "10-100km", "100+km")),
  observed = obs
)

p1_plot <- ggplot(bin_post_df, aes(x = prob)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = "#185FA5", alpha = 0.6) +
  geom_vline(data = obs_lines,
             aes(xintercept = observed),
             color = "#D85A30", linewidth = 1, linetype = "dashed") +
  facet_wrap(~ bin, scales = "free") +
  labs(
    x = "Bin probability",
    y = "Density"
  ) +
  theme_bw(base_size = 11) +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1))

# Plot 2: density curve
p2_plot <- ggplot() +
  geom_ribbon(data = summary_df,
              aes(x = distance, ymin = lo, ymax = hi),
              fill = "#185FA5", alpha = 0.3) +
  geom_line(data = summary_df,
            aes(x = distance, y = median),
            color = "#185FA5", linewidth = 1) +
  labs(
    x = "Distance from home (km)",
    y = ""
  ) +
  theme_bw(base_size = 11) +
  theme_minimal()+
  scale_x_continuous(
    trans  = "log10",
    limits = c(0.1, 125),
    breaks = c(0.1, 1, 10, 100),
    labels = c(0.1, 1, 10, 100)
  ) +
  scale_y_continuous(
    trans  = "log10",
    breaks = c(0.0001, 0.001, 0.01, 0.1, 0.5),
    labels = c("0.0001", "0.001", "0.01", "0.1", "0.5")
  )

p_merge <- p1_plot + p2_plot +
  plot_layout(widths = c(3, 2))

ggsave(file.path("figure/C1_p5_kernel/sar_by_w.png"),
       p_merge, width = 8, height = 3, dpi = 500)

# ── 5. Distance weight function for model use ─────────────────
# w(d) = p * a1 * exp(-a1*d) + (1-p) * a2 * exp(-a2*d)
# Use posterior median parameters

p_hat  <- median(post$p)
a1_hat <- median(post$a1)
a2_hat <- median(post$a2)

distance_weight <- function(d) {
  p_hat  * a1_hat * exp(-a1_hat * d) +
    (1 - p_hat) * a2_hat * exp(-a2_hat * d)
}

cat(sprintf("\nDistance weight function (posterior median):\n"))
cat(sprintf("  w(d) = %.3f * %.4f * exp(-%.4f*d) + %.3f * %.4f * exp(-%.4f*d)\n",
            p_hat, a1_hat, a1_hat,
            1 - p_hat, a2_hat, a2_hat))
cat(sprintf("\nExample weights:\n"))
for (d in c(1, 5, 10, 20, 50, 100, 200)) {
  cat(sprintf("  w(%3d km) = %.6f\n", d, distance_weight(d)))
}

# ── 6. Save distance kernel parameters ────────────────────────
# Posterior median parameters for use in community network construction
kernel_params <- list(
  p_hat  = p_hat,
  a1_hat = a1_hat,
  a2_hat = a2_hat,
  # Full posterior samples for uncertainty propagation if needed
  post_samples = data.frame(
    p  = post$p[idx],
    a1 = post$a1[idx],
    a2 = post$a2[idx]
  ),
  # Distance weight function (closure with baked-in params)
  distance_weight_fn = distance_weight
)

dir.create("output/kernel", showWarnings = FALSE, recursive = TRUE)
saveRDS(kernel_params, "output/kernel/community_distance_kernel.rds")

cat(sprintf("\nKernel saved to output/kernel/community_distance_kernel.rds\n"))
cat(sprintf("  p_hat  : %.4f\n", p_hat))
cat(sprintf("  a1_hat : %.4f\n", a1_hat))
cat(sprintf("  a2_hat : %.4f\n", a2_hat))
