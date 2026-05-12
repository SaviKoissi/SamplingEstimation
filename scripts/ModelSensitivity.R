################################################################################
## JTB PIPELINE: SPATIAL DDE & CONSTRAINED HIERARCHICAL OPTIMIZATION
## Savi et al. 2026 (Under Review in JTB)
################################################################################

# Load required libraries 

library(deSolve)
library(tidyverse)
library(viridis)
library(patchwork)

# ==============================================================================
# 1. BIOLOGICAL ENGINE (DDE)
# ==============================================================================

M <- 100 # This is a place holder representing 100 farms)

# Attribute coordinate to each farm
coords <- cbind(runif(M, 0, 100), runif(M, 0, 100))

# Compute the distance matrix
D_mat  <- as.matrix(dist(coords)) # This is a square matrix of 100 \times 100 

# Parameters of the DDE model
# Here we use parameters from existing literature 
# User can load his own data and run a parameter search using ABC or POMP for e.g.


params_bio <- list(
  beta   = 0.50,          # local transmission
  alpha  = 0.33,          # between-farm transmission
  delta  = 0.50,          # distance decay
  gamma  = 0.05,          # removal / replanting rate
  tau    = 100,           # latency period
  kappa  = 8,             # behavioral response
  N      = rep(1000, M),  # plants per farm
  D      = D_mat
)


# The DDE model

cutting_model <- function(t, y, parms) {
  with(as.list(parms), {
    S <- y[1:M]; L <- y[(M + 1):(2 * M)]; I <- y[(2 * M + 1):(3 * M)]
    theta <- (L + I) / N
    W <- exp(-delta * D); diag(W) <- 0
    Lambda <- beta * theta + alpha * (W %*% theta)
    Phi <- as.vector(Lambda * S + gamma * N * (I/(N + kappa * I)))
    Phi_tau <- if (t > tau) lagvalue(t - tau, 3 * M + 1:M) * exp(-gamma * tau) 
    else rep(0, M)
    
    dS <- -Lambda * S + gamma * N * (1 - (I/(N + kappa * I))) - gamma * S
    dL <- Phi - Phi_tau - gamma * L
    dI <- Phi_tau - gamma * I
    return(list(c(dS, dL, dI, Phi), theta = theta))
  })
}

# Initial state 
state0 <- c(S=rep(980, M), L=rep(20, M), I=rep(0, M), Phi=rep(0, M))

## History function (used automatically by dede)
history <- function(t, parms) {
  c(
    S   = parms$N * (1 - theta0),
    L   = parms$N * theta0,
    I   = rep(0, M),
    Phi = rep(0, M)
  )
}

# solving the DDE
out_dde <- as.data.frame(dede(y=state0, times=seq(0, 150, by=2), 
                              func=cutting_model, parms=params_bio, 
                              initfunc = history))

# ==============================================================================
# 2. DECISION ENGINE (Automatic Thresholds & Optimizer)
# ==============================================================================

# The determination of the theta1 and theta2 are parameters dependent thus, with
# your own data, these metrics will change

get_thresholds <- function(p) {
  list(
    theta1 = p$gamma / p$beta, 
    # theta2: Spatial export threshold (Risk of secondary spread)
    theta2 = p$gamma / (p$alpha * exp(-p$delta * 10)) 
    
  )
}

# ==============================================================================
# 2. DECISION ENGINE: RESPONSIVE OPTIMIZER
# ==============================================================================

find_optimal_policy <- function(curr_theta, curr_phi, budget = 4000) {
  
  # --- Costs ---
  C_travel <- 150 # Cost per farm visit
  C_test   <- 5   # Cost per plant sampled
  
  # --- Value Scaling ---
  # Max_Value represents the societal/economic weight of the info
  Max_Value <- 8e5   
  
  # --- Threshold weighting (Focusing on the Decision Frontier) ---
  theta1 <- 0.08
  theta2 <- 0.20
  frontier_weight <- exp(-5 * abs(curr_theta - theta1)) + 
    exp(-5 * abs(curr_theta - theta2))
  
  # --- Effective Detectability ---
  sigma_I <- 0.85 # Symptomatic detection prob
  sigma_L <- 0.03 # Latent detection prob (low)
  p_eff   <- max(curr_theta * (sigma_I * curr_phi + sigma_L * (1 - curr_phi)),
                 1e-6)
  
  # --- Search Space ---
  expand.grid(
    n = seq(2, 600, by = 2), # Number of farms
    m = seq(1, 30, by = 1)   # Plants per farm (Constraint)
  ) %>%
    mutate(
      # Cost Logic: Fixed travel cost per farm + Variable test cost per farm
      Cost = n * C_travel + n * m * C_test
    ) %>%
    # filter(Cost <= budget) %>%
    # mutate(
    #   # Probability of detection at a single farm
    #   P_det_farm = 1 - (1 - p_eff)^m,
    #   
    #   # Expected Number of Detections (Utility of Information)
    #   # This ensures that more detections = more value, preventing flatlines.
    #   Expected_Detections = n * P_det_farm,
    #   
    #   # Responsive Value Function (Logarithmic diminishing returns)
    #   Value = Max_Value * frontier_weight * log(1 + Expected_Detections),
    #   
    #   ENGS = Value - Cost
    # ) %>%
    mutate(
      # 1. Cost Logic
      Cost = n * C_travel + n * m * C_test
    ) %>%
    filter(Cost <= budget) %>%
    mutate(
      # 2. Local success (at 1 farm)
      P_det_farm = 1 - (1 - p_eff)^m,
      
      # 3. Regional success (across n farms) - THE MISSING LINK
      P_det_total = 1 - (1 - P_det_farm)^n,
      
      # 4. Utility of Information
      # Use P_det_total so that increasing 'n' actually increases 'Value'
      Value = Max_Value * frontier_weight * log(1 + P_det_total * 100),
      
      ENGS = Value - Cost
    )
    arrange(desc(ENGS)) %>%
    slice(1)
}

#==============================================================================

library(dplyr)
library(purrr)

# =========================================================
# OPTIMAL SURVEILLANCE POLICY
# =========================================================

find_optimal_policy <- function(curr_theta,
                                curr_phi,
                                budget = 40000) {
  
  # =======================================================
  # COSTS
  # =======================================================
  
  C_travel <- 150   # cost per farm visited
  C_test   <- 5     # cost per plant sampled
  
  
  # =======================================================
  # VALUE SCALING
  # =======================================================
  
  Max_Value <- 5e4
  
  
  # =======================================================
  # THRESHOLD WEIGHTING
  # =======================================================
  
  theta1 <- 0.08
  theta2 <- 22.48
  
  frontier_weight <-
    exp(-5 * abs(curr_theta - theta1)) +
    exp(-5 * abs(curr_theta - theta2))
  
  
  # =======================================================
  # EFFECTIVE DETECTABILITY
  # =======================================================
  
  sigma_I <- 0.85
  sigma_L <- 0.03
  
  p_eff <- max(
    curr_theta *
      (sigma_I * curr_phi +
         sigma_L * (1 - curr_phi)),
    1e-6
  )
  
  
  # =======================================================
  # WITHIN-FARM CORRELATION
  # =======================================================
  # Plants within a farm are not independent.
  # alpha < 1 introduces diminishing returns
  # from additional plant samples.
  # =======================================================
  
  alpha <- 0.5
  
  
  # =======================================================
  # SEARCH SPACE
  # =======================================================
  
  results <- expand.grid(
    n = seq(2, 600, by = 2),   # farms visited
    m = seq(1, 30, by = 1)     # plants sampled per farm
  ) %>%
    
    mutate(
      
      # ---------------------------------------------------
      # TOTAL COST
      # ---------------------------------------------------
      
      Cost = n * C_travel + n * m * C_test
      
    ) %>%
    
    filter(Cost <= budget) %>%
    
    mutate(
      
      # ---------------------------------------------------
      # EFFECTIVE SAMPLE SIZE
      # ---------------------------------------------------
      # Reduces independence assumption
      # ---------------------------------------------------
      
      m_eff = m^alpha,
      
      
      # ---------------------------------------------------
      # DETECTION PROBABILITY AT ONE FARM
      # ---------------------------------------------------
      
      P_det_farm =
        1 - (1 - p_eff)^m_eff,
      
      
      # ---------------------------------------------------
      # EXPECTED DETECTIONS ACROSS FARMS
      # ---------------------------------------------------
      
      Expected_detections =
        n * P_det_farm,
      
      
      # ---------------------------------------------------
      # VALUE OF INFORMATION
      # ---------------------------------------------------
      
      Value =
        Max_Value *
        frontier_weight *
        log1p(Expected_detections),
      
      
      # ---------------------------------------------------
      # EXPECTED NET GAIN FROM SURVEILLANCE
      # ---------------------------------------------------
      
      ENGS = Value - Cost
    )
  
  
  # =======================================================
  # RETURN BEST POLICY
  # =======================================================
  
  best_policy <- results %>%
    arrange(desc(ENGS)) %>%
    slice(1)
  
  return(best_policy)
}


# =========================================================
# EXAMPLE EPIDEMIC PARAMETERS
# =========================================================

curr_theta <- 0.12
curr_phi   <- 0.25


# =========================================================
# BUDGET SCENARIOS
# STOP AT 80K
# =========================================================

test_budgets <- c(
  5000,
  15000,
  40000,
  80000
)


# =========================================================
# RUN SENSITIVITY ANALYSIS
# =========================================================

sensitivity_results <- map_df(test_budgets, function(b) {
  
  res <- find_optimal_policy(
    curr_theta = curr_theta,
    curr_phi   = curr_phi,
    budget     = b
  )
  
  res$TotalBudget <- b
  
  
  # Recompute p_eff for reporting
  
  res$p_eff <- max(
    curr_theta *
      (0.85 * curr_phi +
         0.03 * (1 - curr_phi)),
    1e-6
  )
  
  
  # Regional certainty diagnostic
  
  res$Regional_Certainty <-
    1 - (1 - res$P_det_farm)^res$n
  
  
  return(res)
})


# =========================================================
# FINAL OUTPUT TABLE
# =========================================================

final_table <- sensitivity_results %>%
  select(
    TotalBudget,
    n,
    m,
    m_eff,
    Cost,
    p_eff,
    P_det_farm,
    Expected_detections,
    Regional_Certainty,
    Value,
    ENGS
  )

print(final_table)






# ==============================================================================
# 3. BUDGET SENSITIVITY ANALYSIS
# ==============================================================================

# Extract a snapshot for testing (e.g., at t=80)
snap <- out_dde %>% filter(time == 80)
L_cols <- (M + 2):(2 * M + 1); I_cols <- (2 * M + 2):(3 * M + 1)
curr_theta <- mean((as.matrix(snap[, L_cols]) + as.matrix(snap[, I_cols])) / 1000)
curr_phi   <- mean(as.numeric(snap[, I_cols])) / (curr_theta * 1000 + 1e-6)

test_budgets <- c(5000, 15000, 40000, 80000, 150000)

sensitivity_results <- map_df(test_budgets, function(b) {
  res <- find_optimal_policy(curr_theta, curr_phi, budget = b)
  res$TotalBudget <- b
  return(res)
})
opt_policy <- find_optimal_policy(curr_theta, curr_phi, budget = 4000)
# ==============================================================================
# 4. VISUALIZATION & OUTPUT
# ==============================================================================

print(sensitivity_results %>% select(TotalBudget, n, m, Cost, ENGS))

# Figure: Policy Scaling with Budget
p_budget <- ggplot(sensitivity_results, aes(x = TotalBudget)) +
  geom_line(aes(y = n, color = "Number of Farms (n)"), linewidth = 1) +
  geom_point(aes(y = n, color = "Number of Farms (n)")) +
  geom_line(aes(y = m * 5, color = "Plants/Farm (m scaled)"), linewidth = 1, 
            linetype = "dashed") +
  scale_y_continuous(sec.axis = sec_axis(~./5, name = "Plants per Farm (m)")) +
  labs(title = "Budget Responsiveness of Surveillance Design",
       x = "Total Budget ($)", y = "Farms (n)", color = "Metric") +
  theme_minimal()

p_budget

# ==============================================================================
# 3. VISUALIZATION
# ==============================================================================

# FIGURE A: THE LATENCY TRAP
p1_trap <- out_dde %>%
  mutate(mean_L = rowMeans(.[, L_cols]), mean_I = rowMeans(.[, I_cols])) %>%
  ggplot(aes(x = time)) +
  geom_area(aes(y = mean_L + mean_I, fill = "Latent (Invisible)"), alpha = 0.4) +
  geom_area(aes(y = mean_I, fill = "Symptomatic (Visible)"), alpha = 0.8) +
  scale_fill_manual(values = c("Latent (Invisible)" = "grey70", "Symptomatic (Visible)" = "firebrick")) +
  labs(title = "A: Biological Latency Trap", y = "Population", x = "Days") + theme_minimal()

# FIGURE B: THE OPTIMIZATION FRONTIER (Corrected Axes)
p2_frontier <- expand.grid(n = seq(2, 250, 5), m = seq(1, 30, 1)) %>%
  mutate(Cost = n * 150 + n * m * 5) %>% filter(Cost <= 40000) %>%
  mutate(ENGS = (n * (1 - (1 - current_phi * current_theta)^m) * 2500) - Cost) %>%
  ggplot(aes(x = n, y = m, fill = ENGS)) +
  geom_tile() + scale_fill_viridis_c(option = "plasma") +
  geom_point(aes(x = opt_policy$n, y = opt_policy$m), color = "white", size = 3, shape = 18) +
  labs(title = "B: Optimal Allocation", x = "Number of Farms (n)", y = "Plants/Farm (m)") +
  theme_minimal()

# FIGURE C: SPATIAL FRONTIER
p3_map <- expand.grid(x = seq(0, 100, 2), y = seq(0, 100, 2)) %>%
  mutate(dist = sqrt((x - 50)^2 + (y - 50)^2),
         theta_local = 0.5 * exp(-0.05 * dist),
         W_j = theta_local * (1 - theta_local)) %>% 
  ggplot(aes(x, y, fill = W_j)) +
  geom_tile() + scale_fill_viridis_c(option = "mako") +
  labs(title = "C: Spatial Info Value", subtitle = "Wj Decision Curvature") + 
  theme_void()

# COMPOSITION
(p1_trap | p2_frontier) / (p3_map + plot_spacer()) + 
  plot_annotation(title = "Constrained Surveillance Optimization (m <= 30)")

# CONSOLE OUTPUT
cat("--- JTB THRESHOLD & POLICY RESULTS ---\n")
cat("Spatial Export Threshold (theta2):", round(thresholds$theta2, 4), "\n")
cat("Optimal Farms (n*):", opt_policy$n, "\n")
cat("Optimal Plants per Farm (m*):", opt_policy$m, "\n")

# ==============================================================================
# 5. GLOBAL SENSITIVITY ANALYSIS (LHS-Style Stress Test)
# ==============================================================================
library(reshape2)

# Define parameter ranges for sensitivity (+/- 20% of baseline)
# Define wider ranges and include Budget
param_ranges <- list(
  beta     = c(0.2, 0.8),    # Wider range
  tau      = c(10, 150),     # Wider range
  budget   = c(500, 5000), # Variable budget
  C_travel = c(50, 300)      # Variable costs
)

set.seed(123)
# Generate scenarios including Budget
sensitivity_data <- data.frame(
  scenario = 1:50,
  beta     = runif(50, param_ranges$beta[1], param_ranges$beta[2]),
  tau      = runif(50, param_ranges$tau[1], param_ranges$tau[2]),
  budget   = runif(50, param_ranges$budget[1], param_ranges$budget[2]),
  C_travel = runif(50, param_ranges$C_travel[1], param_ranges$C_travel[2])
)

# Run the loop using the scenario-specific budget
sensitivity_results_global <- sensitivity_data %>%
  rowwise() %>%
  mutate(
    n_star = find_optimal_policy(curr_theta, curr_phi, budget = budget)$n
  )


correlations <- sensitivity_results_global %>%
  gather(parameter, value, beta:C_travel) %>%
  group_by(parameter) %>%
  summarize(
    # Use handle_zero_sd logic: if sd is 0, correlation is 0
    correlation = if(sd(n_star) == 0) 0 else cor(value, n_star, method = "spearman")
  )

# FIGURE D: GLOBAL SENSITIVITY (TORNADO CHART)
p4_sensitivity <- ggplot(correlations, aes(x = reorder(parameter, correlation), y = correlation, fill = correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "firebrick", mid = "grey90", high = "steelblue") +
  labs(title = "D: Global Sensitivity Analysis", 
       subtitle = "Influence on Optimal Sample Size (n*)",
       x = "Parameter", y = "Spearman Correlation") +
  theme_minimal() + theme(legend.position = "none")

# ==============================================================================
# 6. FINAL COMPOSITION (JTB FULL FIGURE)
# ==============================================================================
(p1_trap | p2_frontier) / (p3_map | p4_sensitivity) + 
  plot_annotation(title = "Bio-Economic Surveillance Optimization Framework",
                  subtitle = "Integration of Latency Dynamics, Policy Frontiers, and Sensitivity Indices")



