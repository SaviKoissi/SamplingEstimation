# ==============================================================================
# Dynamic Spatial EVSI Engine: DDE-Coupled 3-Action Decision Framework
# Publication Pipeline (Low-Temperature Particle Filter Optimization)
# ==============================================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(scales)
  library(tidyr)
  library(future)
  library(furrr)
  library(parallel)
})

# ------------------------------------------------------------------------------
# 1. Parameter Configuration
# ------------------------------------------------------------------------------
national_params <- list(
  N_j         = 1000,     # Plant host population per farm
  R_regions   = 8,        # Total administrative regions
  rho_j       = 0.15,     # Within-farm ICC
  rho_R       = 0.08,     # Between-region ICC
  
  beta_v      = 0.35,     # Transmission rate constant
  gamma       = 0.005,    # Mortality/roguing rate (1/days)
  tau         = 14,       # Latent period delay (days)
  T_horizon   = 60,       # Projected decision horizon (days)
  
  eff_a0      = 0.00,     # Passive monitoring
  eff_a1      = 0.75,     # Local containment
  eff_a2      = 0.98,     # Regional eradication
  
  Se_vis      = 0.75, Sp_vis = 0.85,
  Se_pcr      = 0.80, Sp_pcr = 0.95,
  
  C0_base     = 100.00,   # Fixed regional deployment cost ($)
  C_travel    = 8.00,     # Travel cost per farm ($)
  C_vis       = 0.50,     # Visual inspection cost per plant ($)
  C_pcr       = 6.00,     # PCR assay cost per plant ($)
  B_max       = 250000.00,# Budget cap ($)
  
  L0_per_reg  = 2500000.0,# Total yield loss coefficient ($)
  C_act_a1    = 2500.00,  # Action a1 implementation cost ($)
  C_env_a1    = 500.00,   # Action a1 environmental penalty ($)
  C_act_a2    = 35000.00, # Action a2 implementation cost ($)
  C_env_a2    = 12000.00  # Action a2 environmental penalty ($)
)

# ------------------------------------------------------------------------------
# 2. DDE Integration & Probability Engines
# ------------------------------------------------------------------------------
sli_dde_system <- function(t, y, parms) {
  S <- y[1]; L <- y[2]; I <- y[3]
  Lambda <- parms$beta_v * I / parms$N_j
  Phi_t  <- Lambda * S
  
  if (t < parms$tau) {
    Phi_delayed <- 0.0
  } else {
    lag_y <- lagvalue(t - parms$tau)
    Phi_delayed <- (parms$beta_v * lag_y[3] / parms$N_j) * lag_y[1]
  }
  
  f_psi <- parms$current_eff
  dS <- -Phi_t + parms$gamma * parms$N_j * (1 - f_psi) - parms$gamma * S
  dL <- Phi_t - exp(-parms$gamma * parms$tau) * Phi_delayed - parms$gamma * L
  dI <- exp(-parms$gamma * parms$tau) * Phi_delayed - parms$gamma * I
  
  return(list(c(dS, dL, dI)))
}

solve_dde_terminal_prevalence <- function(pi_init, action_eff, p) {
  I_0 <- max(0.1, pi_init * p$N_j)
  y_init <- c(S = p$N_j - I_0, L = 0.0, I = I_0)
  times  <- seq(0, p$T_horizon, by = 1)
  
  p_run <- p; p_run$current_eff <- action_eff
  out <- tryCatch({
    dede(y = y_init, times = times, func = sli_dde_system, 
         parms = p_run, atol = 1e-6)
  }, error = function(e) NULL)
  
  if (is.null(out)) return(pi_init)
  return(pmax(0, pmin(1, tail(out[, "I"], 1) / p$N_j)))
}

compute_decision_thresholds <- function(p) {
  R <- p$R_regions
  mu_seq <- seq(0.0001, 0.15, length.out = 1500)
  
  I_a0 <- sapply(mu_seq, function(pi) solve_dde_terminal_prevalence(pi,
                                                                    p$eff_a0, p))
  I_a1 <- sapply(mu_seq, function(pi) solve_dde_terminal_prevalence(pi, 
                                                                    p$eff_a1, p))
  I_a2 <- sapply(mu_seq, function(pi) solve_dde_terminal_prevalence(pi, 
                                                                    p$eff_a2, p))
  
  loss_a0 <- (R * p$L0_per_reg) * I_a0
  loss_a1 <- R * (p$C_act_a1 + p$C_env_a1) + (R * p$L0_per_reg) * I_a1
  loss_a2 <- R * (p$C_act_a2 + p$C_env_a2) + (R * p$L0_per_reg) * I_a2
  
  idx_01 <- which(loss_a1 < loss_a0)[1]
  idx_12 <- which(loss_a2 < loss_a1)[1]
  
  mu_01 <- if (!is.na(idx_01)) mu_seq[idx_01] else 0.0025
  mu_12 <- if (!is.na(idx_12)) mu_seq[idx_12] else 0.0450
  
  return(c(mu_01 = mu_01, mu_12 = mu_12))
}

dbetabinom_log <- function(x, size, prob, rho) {
  if (rho <= 1e-6) return(dbinom(x, size = size, prob = prob, log = TRUE))
  a <- pmax(prob * (1 - rho) / rho, 1e-4)
  b <- pmax((1 - prob) * (1 - rho) / rho, 1e-4)
  lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b)
}

rbetabinom <- function(n, size, prob, rho) {
  if (rho <= 0) return(rbinom(n, size, prob))
  a <- pmax(prob * (1 - rho) / rho, 0.01)
  b <- pmax((1 - prob) * (1 - rho) / rho, 0.01)
  rbinom(n, size, rbeta(n, a, b))
}

compute_bayes_risk_precomputed <- function(weights, loss_a0, loss_a1, loss_a2) {
  min(sum(weights * loss_a0), sum(weights * loss_a1), sum(weights * loss_a2))
}

# ------------------------------------------------------------------------------
# 3. EVSI Particle Filter
# ------------------------------------------------------------------------------
evaluate_national_evsi_dde <- function(N, n, m, mu_j, p, N_particles = 2000, S_samples = 200, temp = 0.5) {
  R <- p$R_regions
  Total_farms <- N * R
  
  M_eff_region <- (1 - p$rho_R) / p$rho_R
  alpha_nat    <- pmax(mu_j * M_eff_region, 0.01)
  beta_nat     <- pmax((1 - mu_j) * M_eff_region, 0.01)
  
  pi_particles <- rbeta(N_particles, shape1 = alpha_nat, shape2 = beta_nat)
  
  I_term_a0 <- sapply(pi_particles, function(pi) solve_dde_terminal_prevalence(pi, p$eff_a0, p))
  I_term_a1 <- sapply(pi_particles, function(pi) solve_dde_terminal_prevalence(pi, p$eff_a1, p))
  I_term_a2 <- sapply(pi_particles, function(pi) solve_dde_terminal_prevalence(pi, p$eff_a2, p))
  
  loss_a0 <- (R * p$L0_per_reg) * I_term_a0
  loss_a1 <- R * (p$C_act_a1 + p$C_env_a1) + (R * p$L0_per_reg) * I_term_a1
  loss_a2 <- R * (p$C_act_a2 + p$C_env_a2) + (R * p$L0_per_reg) * I_term_a2
  
  w_prior <- rep(1 / N_particles, N_particles)
  v_prior <- compute_bayes_risk_precomputed(w_prior, loss_a0, loss_a1, loss_a2)
  
  v_post <- numeric(S_samples)
  N_vis_tot <- Total_farms * n
  N_pcr_tot <- Total_farms * m
  
  p_vis_part <- pmax(pmin(p$Se_vis * pi_particles + (1 - p$Sp_vis) * (1 - pi_particles), 0.999), 0.001)
  p_pcr_part <- pmax(pmin(p$Se_pcr * pi_particles + (1 - p$Sp_pcr) * (1 - pi_particles), 0.999), 0.001)
  
  for (s in 1:S_samples) {
    pi_true_s <- rbeta(1, shape1 = alpha_nat, shape2 = beta_nat)
    
    p_vis_s <- pmax(pmin(p$Se_vis * pi_true_s + (1 - p$Sp_vis) * (1 - pi_true_s), 0.999), 0.001)
    p_pcr_s <- pmax(pmin(p$Se_pcr * pi_true_s + (1 - p$Sp_pcr) * (1 - pi_true_s), 0.999), 0.001)
    
    Y_vis <- rbetabinom(1, size = N_vis_tot, prob = p_vis_s, rho = p$rho_j)
    Y_pcr <- rbetabinom(1, size = N_pcr_tot, prob = p_pcr_s, rho = p$rho_j)
    
    log_lik <- dbetabinom_log(Y_vis, size = N_vis_tot, prob = p_vis_part, rho = p$rho_j) +
      dbetabinom_log(Y_pcr, size = N_pcr_tot, prob = p_pcr_part, rho = p$rho_j)
    
    # Low-temperature scaling for weight stabilization
    max_log <- max(log_lik)
    weights <- exp((log_lik - max_log) * temp)
    
    if (sum(weights) <= 0 || is.na(sum(weights))) {
      weights <- rep(1 / N_particles, N_particles)
    } else {
      weights <- weights / sum(weights)
    }
    
    v_post[s] <- compute_bayes_risk_precomputed(weights, loss_a0, loss_a1, loss_a2)
  }
  
  gross_evsi <- v_prior - mean(v_post)
  return(max(0, gross_evsi))
}

compute_national_cost <- function(N, n, m, p) {
  p$R_regions * (p$C0_base + N * p$C_travel + N * n * p$C_vis + N * m * p$C_pcr)
}

run_parallel_sweep <- function(mu_j, params, n_workers = NULL) {
  if (is.null(n_workers)) n_workers <- max(1, parallel::detectCores() - 1)
  
  grid <- expand.grid(
    N_farms_per_reg = seq(10, 100, by = 15),
    n_plants        = c(10, 25, 50),
    m_pcr           = c(1, 3, 5)
  )
  grid$Cost <- compute_national_cost(grid$N_farms_per_reg, 
                                     grid$n_plants, grid$m_pcr, params)
  grid      <- grid %>% filter(Cost <= params$B_max & m_pcr <= n_plants)
  
  plan(multisession, workers = n_workers)
  evsi_results <- future_map_dfr(
    seq_len(nrow(grid)),
    function(i) {
      evsi_val <- evaluate_national_evsi_dde(
        grid$N_farms_per_reg[i], grid$n_plants[i], grid$m_pcr[i], mu_j, params
      )
      data.frame(grid_idx = i, EVSI = round(evsi_val, 2))
    },
    .options = furrr_options(seed = TRUE)
  )
  plan(sequential)
  
  grid$EVSI <- evsi_results$EVSI[match(1:nrow(grid), evsi_results$grid_idx)]
  grid$ENBS <- round(pmax(0, grid$EVSI - grid$Cost), 2)
  return(grid)
}

# ------------------------------------------------------------------------------
# 4. Pipeline Execution & Figure Generation
# ------------------------------------------------------------------------------
set.seed(42)

cat("\nSolving Decision Boundaries...\n")
thresholds <- compute_decision_thresholds(national_params)
cat(sprintf("Action Boundaries: μ*(a0->a1) = %.4f | μ*(a1->a2) = %.4f\n", thresholds["mu_01"], thresholds["mu_12"]))

cat("\nExecuting Parallel EVSI Sweeps near Action Boundaries...\n")
res_boundary_01 <- run_parallel_sweep(mu_j = thresholds["mu_01"], 
                                      params = national_params)
res_boundary_12 <- run_parallel_sweep(mu_j = thresholds["mu_12"], 
                                      params = national_params)

opt_01 <- res_boundary_01[which.max(res_boundary_01$ENBS), ]
opt_12 <- res_boundary_12[which.max(res_boundary_12$ENBS), ]

theme_q1 <- function() {
  theme_classic(base_size = 11) %+replace%
    theme(
      panel.grid.major = element_line(color = "#E5E5E5", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(face = "bold", size = 10),
      axis.text        = element_text(color = "black", size = 9),
      legend.position  = "bottom",
      legend.key.size  = unit(1, "lines"),
      legend.key.width = unit(2, "cm"),
      legend.text      = element_text(size = 6, margin = margin(t = 5, b = 5, 
                                                                unit = "pt"))
    )
}

# ------------------------------------------------------------------------------
# FIGURE 1: Dynamic Bayes Risk & Action Surfaces
# ------------------------------------------------------------------------------
mu_seq <- seq(0.0001, 0.08, length.out = 300)
R <- national_params$R_regions

I_a0 <- sapply(mu_seq, function(pi) solve_dde_terminal_prevalence(pi, national_params$eff_a0, national_params))
I_a1 <- sapply(mu_seq, function(pi) solve_dde_terminal_prevalence(pi, national_params$eff_a1, national_params))
I_a2 <- sapply(mu_seq, function(pi) solve_dde_terminal_prevalence(pi, national_params$eff_a2, national_params))

risk_df <- data.frame(
  mu = mu_seq,
  a0 = (R * national_params$L0_per_reg) * I_a0,
  a1 = R * (national_params$C_act_a1 + national_params$C_env_a1) + (R * national_params$L0_per_reg) * I_a1,
  a2 = R * (national_params$C_act_a2 + national_params$C_env_a2) + (R * national_params$L0_per_reg) * I_a2
) %>% pivot_longer(-mu, names_to = "Action", values_to = "Loss")

fig1_A <- ggplot(risk_df, aes(x = mu, y = Loss, color = Action)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = c(thresholds["mu_01"], thresholds["mu_12"]), 
             linetype = "dashed", color = "gray30") +
  scale_color_manual(values = c("a0" = "#D90429", "a1" = "#2A9D8F", 
                                "a2" = "#0D3B66"),
                     labels = c("a0: Passive", "a1: Containment", 
                                "a2: Eradication")) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = dollar_format()) +
  labs(title = "A. Dynamic Bayes Risk Lines & Boundaries", 
       x = "Prevalence (μ)", y = "Expected Posterior Loss ($)", color = NULL) +
  theme_q1()

fig1_B <- res_boundary_01 %>%
  filter(m_pcr == opt_01$m_pcr) %>%
  ggplot(aes(x = N_farms_per_reg, y = factor(n_plants), fill = ENBS)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "viridis", label = dollar_format()) +
  geom_point(data = opt_01, aes(x = N_farms_per_reg, y = factor(n_plants)), 
             color = "red", shape = 18, size = 4) +
  labs(title = sprintf("B. ENBS Surface at μ*(a0→a1) (m = %d PCR)", opt_01$m_pcr),
       x = "Farms Sampled / Region (N)", y = "Plants / Farm (n)", fill = "ENBS ($)") +
  theme_q1()

fig1_C <- res_boundary_12 %>%
  filter(m_pcr == opt_12$m_pcr) %>%
  ggplot(aes(x = N_farms_per_reg, y = factor(n_plants), fill = ENBS)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", label = dollar_format()) +
  geom_point(data = opt_12, aes(x = N_farms_per_reg, y = factor(n_plants)), 
             color = "cyan", shape = 18, size = 4) +
  labs(title = sprintf("C. ENBS Surface at μ*(a1→a2) (m = %d PCR)", opt_12$m_pcr),
       x = "Farms Sampled / Region (N)", y = "Plants / Farm (n)", fill = "ENBS ($)") +
  theme_q1()

fig1_D <- res_boundary_01 %>%
  filter(n_plants == opt_01$n_plants & m_pcr == opt_01$m_pcr) %>%
  ggplot(aes(x = N_farms_per_reg)) +
  geom_line(aes(y = EVSI, color = "Gross EVSI ($)"), linewidth = 1.1) +
  geom_line(aes(y = Cost, color = "Surveillance Cost ($)"), linewidth = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Gross EVSI ($)" = "#2A9D8F", 
                                "Surveillance Cost ($)" = "#E76F51")) +
  scale_y_continuous(labels = dollar_format()) +
  labs(title = "D. EVSI vs. Cost Trade-off Curve", x = "Farms Sampled / Region (N)", y = "Value / Cost ($)", color = NULL) +
  theme_q1()

fig1 <- grid.arrange(fig1_A, fig1_B, fig1_C, fig1_D, ncol = 2)
ggsave("Figure_1_EVSI_Framework_Q1.pdf", fig1, width = 11, height = 8.5, dpi = 300)
ggsave("Figure_1_EVSI_Framework_Q1.png", fig1, width = 11, height = 8.5, dpi = 300)

# ------------------------------------------------------------------------------
# FIGURE 2: Robust Sensitivity Analysis
# ------------------------------------------------------------------------------
c_act_seq <- seq(15000, 45000, length.out = 12)
c_env_seq <- seq(5000, 20000, length.out = 12)
cost_grid <- expand.grid(C_act_a2 = c_act_seq, C_env_a2 = c_env_seq)

cat("\nRunning Contour Sensitivity Sweep...\n")
cost_grid$mu_12 <- apply(cost_grid, 1, function(row) {
  p_temp <- national_params
  p_temp$C_act_a2 <- row[1]
  p_temp$C_env_a2 <- row[2]
  compute_decision_thresholds(p_temp)["mu_12"]
})

fig2_A <- ggplot(cost_grid, aes(x = C_act_a2 / 1000, y = C_env_a2 / 1000, z = mu_12 * 100)) +
  geom_contour_filled(bins = 8) +
  scale_fill_viridis_d(option = "cividis", name = "Threshold μ*12 (%)") +
  labs(title = "A. Eradication Switch Threshold (μ*12 %)",
       x = "Implementation Cost C_act(a2) ($k)", y = "Environmental Penalty C_env(a2) ($k)") +
  theme_q1()

icc_grid <- expand.grid(
  rho_R = c(0.04, 0.08, 0.12), 
  rho_j = c(0.05, 0.15, 0.25)
)
icc_results <- list()

cat("Running Robust Dynamic Spatial ICC Sensitivity Sweep...\n")
for (k in 1:nrow(icc_grid)) {
  p_temp <- national_params
  p_temp$rho_R <- icc_grid$rho_R[k]
  p_temp$rho_j <- icc_grid$rho_j[k]
  
  local_thresh <- compute_decision_thresholds(p_temp)["mu_01"]
  
  val <- evaluate_national_evsi_dde(
    N = opt_01$N_farms_per_reg, 
    n = opt_01$n_plants, 
    m = opt_01$m_pcr, 
    mu_j = local_thresh, 
    p = p_temp,
    N_particles = 2000,
    S_samples = 200,
    temp = 0.5
  )
  
  icc_results[[k]] <- data.frame(
    rho_R   = p_temp$rho_R, 
    rho_j   = p_temp$rho_j, 
    mu_star = local_thresh,
    EVSI    = val
  )
}

icc_df <- do.call(rbind, icc_results)

fig2_B <- ggplot(icc_df, aes(x = factor(rho_R), y = factor(rho_j), fill = EVSI / 1000)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(label = sprintf("$%.1fk\n(μ*=%.2f%%)", EVSI / 1000, mu_star * 100)), 
    color = "white", fontface = "bold", size = 3.0
  ) +
  scale_fill_viridis_c(option = "magma", labels = dollar_format(suffix = "k")) +
  labs(
    title = "B. Spatial ICC Impact on Gross EVSI ($k)", 
    x = "Between-Region ICC (ρ_R)", 
    y = "Within-Farm ICC (ρ_j)", 
    fill = "EVSI ($k)"
  ) +
  theme_q1()

fig2 <- grid.arrange(fig2_A, fig2_B, ncol = 2)
ggsave("Figure_2_Sensitivity_Suite_Q1.pdf", fig2, width = 11, height = 5, dpi = 300)
ggsave("Figure_2_Sensitivity_Suite_Q1.png", fig2, width = 11, height = 5, dpi = 300)

# ==============================================================================
# Dynamic Budget Saturation Sweep & Spatial Allocation Table
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(xtable)
})

# Define a spectrum of budget caps to test saturation
budget_spectrum <- c(2500, 5000, 10000, 25000, 50000, 100000, 250000)

# Function to extract optimal sampling strategy for a given budget cap
get_optimal_for_budget <- function(res_grid, B_cap) {
  valid_sub <- res_grid %>% filter(Cost <= B_cap)
  if (nrow(valid_sub) == 0) return(NULL)
  
  opt <- valid_sub[which.max(valid_sub$ENBS), ]
  return(opt)
}

# 1. Sweep across budget levels for both decision boundaries
saturation_list <- list()

for (b_cap in budget_spectrum) {
  opt_01_b <- get_optimal_for_budget(res_boundary_01, b_cap)
  opt_12_b <- get_optimal_for_budget(res_boundary_12, b_cap)
  
  if (!is.null(opt_01_b) && !is.null(opt_12_b)) {
    saturation_list[[as.character(b_cap)]] <- data.frame(
      B_max      = b_cap,
      # Boundary 01 (Containment)
      N_01       = opt_01_b$N_farms_per_reg,       n_01       = opt_01_b$n_plants,
      m_01       = opt_01_b$m_pcr,       Cost_01    = opt_01_b$Cost,
      EVSI_01    = opt_01_b$EVSI,       ENBS_01    = opt_01_b$ENBS,
      # Boundary 12 (Eradication)
      N_12       = opt_12_b$N_farms_per_reg,       n_12       = opt_12_b$n_plants,
      m_12       = opt_12_b$m_pcr,       Cost_12    = opt_12_b$Cost,
      EVSI_12    = opt_12_b$EVSI,       ENBS_12    = opt_12_b$ENBS
    )
  }
}

raw_saturation_df <- do.call(rbind, saturation_list)
