############################################################
## Spatial DDE model for vegetatively propagated crops
## Seasonal dynamics with fixed latency, survival,
## behavioral recruitment, spatial coupling,
## and exact population conservation
############################################################

library(deSolve)
library(reshape2)
library(tidyverse)
library(lhs)
library(sensitivity)
library(geodata)
library(sf)
library(gstat)
library(viridis)
library(patchwork)
library(ggspatial) # For North Arrow and Scale Bar
library(scales)    # For better color labeling

############################################################
## 1. GLOBAL SETUP
############################################################

set.seed(123)

## Number of farms
M <- 40

## Season structure
season_length <- 270     # 9 months
n_seasons     <- 3
times_season  <- seq(0, season_length, by = 1)

############################################################
## 2. FARM LOCATIONS AND DISTANCES (BENIN-LIKE SETUP)
############################################################

coords <- cbind(
  x = runif(M, 0, 100),
  y = runif(M, 0, 100)
)

D <- as.matrix(dist(coords))   # distance matrix

############################################################
## 3. PARAMETERS
############################################################

params <- list(
  beta   = 0.50,          # local transmission
  alpha  = 0.33,          # between-farm transmission
  delta  = 0.50,          # distance decay
  gamma  = 0.05,          # removal / replanting rate
  tau    = 100,           # latency period
  kappa  = 8,             # behavioral response
  N      = rep(1000, M),  # plants per farm
  D      = D
)

############################################################
## 4. INITIAL CONDITIONS AND HISTORY
############################################################

## Initial latent prevalence (heterogeneous, low signal)
theta0 <- runif(M, 0.01, 0.20)

state0 <- c(
  S   = params$N * (1 - theta0),
  L   = params$N * theta0,
  I   = rep(0, M),
  Phi = rep(0, M)     # auxiliary variable for delay
)

## History function (used automatically by dede)
history <- function(t, parms) {
  c(
    S   = parms$N * (1 - theta0),
    L   = parms$N * theta0,
    I   = rep(0, M),
    Phi = rep(0, M)
  )
}

############################################################
## 5. MODEL FUNCTION
############################################################

cutting_model <- function(t, y, parms) {
  
  with(as.list(parms), {
    
    ## Unpack state
    S <- y[1:M]
    L <- y[(M + 1):(2 * M)]
    I <- y[(2 * M + 1):(3 * M)]
    
    ## Prevalence
    psi <- I / N
    
    ## Behavioral response
    fpsi <- psi / (1 + kappa * psi)
    
    ## Spatial kernel
    W <- exp(-delta * D)
    diag(W) <- 0
    
    ## Force of infection (local + spatial)
    Lambda <- beta * (L + I) / N +
      alpha * (W %*% ((L + I) / N))
    
    ## Latent inflow (infection + infected cuttings)
    Phi <- as.vector(Lambda * S + gamma * N * fpsi)
    
    ## Delayed progression with survival
    if (t > tau) {
      Phi_tau <- lagvalue(t - tau, 3 * M + 1:(M)) *
        exp(-gamma * tau)
    } else {
      Phi_tau <- rep(0, M)
    }
    
    ## DDE system
    dS <- -Lambda * S + gamma * N * (1 - fpsi) - gamma * S
    dL <- Phi - Phi_tau - gamma * L
    dI <- Phi_tau - gamma * I
    
    ## Return derivatives and diagnostics
    list(
      c(dS, dL, dI, Phi),
      Lambda = Lambda,
      theta  = (L + I) / N
    )
  })
}

############################################################
## 6. SEASONAL SIMULATION LOOP
############################################################

out_all <- list()
current_state <- state0
current_time_offset <- 0

for (s in 1:n_seasons) {
  
  cat("Simulating season", s, "\n")
  
  out_season <- dede(
    y        = current_state,
    times    = times_season,
    func     = cutting_model,
    parms    = params,
    initfunc = history
  )
  
  out_season <- as.data.frame(out_season)
  
  ## Make time continuous across seasons
  out_season$time <- out_season$time + current_time_offset
  
  ## Store output
  out_all[[s]] <- out_season
  
  ## Extract only the state variables (S, L, I, Phi) for next initial condition
  last_row <- tail(out_season, 1)
  current_state <- as.numeric(last_row[2:(4 * M + 1)])
  names(current_state) <- names(state0)
  
  current_time_offset <- current_time_offset + season_length
}


## Combine all seasons
out_df <- do.call(rbind, out_all)

############################################################
## 7. POST-PROCESSING: COMPARTMENT DYNAMICS
############################################################

## Identify compartment columns
S_cols <- 1:M
L_cols <- (M + 1):(2 * M)
I_cols <- (2 * M + 1):(3 * M)

## Aggregate across farms
compartments_df <- data.frame(
  time = out_df$time,
  S = rowSums(out_df[, S_cols]),
  L = rowSums(out_df[, L_cols]),
  I = rowSums(out_df[, I_cols])
)

## Long format for plotting
compartments_long <- melt(
  compartments_df,
  id.vars = "time",
  variable.name = "Compartment",
  value.name = "Plants"
)

############################################################
## 8. OPTIONAL PLOTTING 
############################################################

ggplot(compartments_long,
       aes(x = time, y = Plants, color = Compartment)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c(
      S = "blue",
      L = "green",
      I = "red"
    )) +
  labs(
    title = "Seasonal Evolution of S–L–I Compartments",
    x = "Time (days)",
    y = "Number of plants"
  ) +
  theme_minimal()

#==========================================================

# Combine all seasons into one data frame
combined_out <- do.call(rbind, out_all)

# Calculate prevalence per farm
M <- params$N %>% length()
S_cols <- 2:(M+1)
L_cols <- (M+2):(2*M+1)
I_cols <- (2*M+2):(3*M+1)
Phi_cols <- (3*M+2):(4*M+1)

# Compute prevalence for each row and farm: (L + I)/N
prevalence_mat <- (as.matrix(combined_out[, L_cols]) + as.matrix(combined_out[, I_cols])) / params$N

# Create a data.frame for prevalence with time and season info
prevalence_df <- data.frame(
  time = combined_out$time,
  season = rep(1:length(out_all), each = nrow(out_all[[1]]))
)

# Add farm columns
colnames(prevalence_mat) <- paste0("Farm_", 1:M)
prevalence_df <- cbind(prevalence_df, prevalence_mat)

# Melt for ggplot
prevalence_long <- melt(prevalence_df, id.vars = c("time", "season"),
                        variable.name = "Farm", value.name = "Prevalence")

# Plot prevalence by season
ggplot(prevalence_long, aes(x = time, y = Prevalence, group = Farm, color = Farm)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ season, scales = "free_x") +
  labs(title = "Farm-level Disease Prevalence Over Seasons",
       x = "Time (days)",
       y = "Prevalence (L + I) / N") +
  theme_minimal() +
  theme(legend.position = "none")

#=======================================================
# Effect of delay on the R0
#=======================================================


# 1. PARAMETERS (Matching your Theoretical Section)
params <- list(
  beta  = 0.5,    # Transmission rate
  gamma = 0.05,   # Removal rate
  alpha = 0.33,    # Spatial transmission scaling
  tau_base = 10,  # Base delay (days/weeks)
  delta = 0.05    # Spatial decay
)

# ----------------------------------------------------------
# 2. THRESHOLD IDENTIFICATION
# ----------------------------------------------------------

# Local Persistence Threshold (theta1): beta * theta1 = gamma
theta1 <- params$gamma / params$beta

# Spatial Export Threshold (theta2): Based on a neighboring farm at distance d
# We solve for theta_j such that alpha * theta_j * exp(-delta * d) = gamma - beta * theta_k
# Assuming a critical neighbor at 10km and subcritical theta_k = 0
d_crit <- 10
theta2 <- (params$gamma) / (params$alpha * exp(-params$delta * d_crit))

cat("Calculated Thresholds:\n")
cat("Local Persistence (theta1):", round(theta1, 4), "\n")
cat("Spatial Export (theta2):", round(theta2, 4), "\n")

# ----------------------------------------------------------
# 3. REPRODUCTION NUMBER (R0) VS DELAY (tau)
# ----------------------------------------------------------

# R0 formula from your section 5.5: 
# R0 = (beta/gamma)*(1 + exp(-gamma * tau)) + exp(-gamma * tau)

calculate_r0 <- function(tau, p = params) {
  term1 <- (p$beta / p$gamma) * (1 + exp(-p$gamma * tau))
  term2 <- exp(-p$gamma * tau)
  return(term1 + term2)
}

tau_range <- seq(0, 100, by = 1)
r0_data <- data.frame(
  tau = tau_range,
  R0 = sapply(tau_range, calculate_r0)
)

# 
# 4. VISUALIZATION
plot_r0 <- ggplot(r0_data, aes(x = tau, y = R0)) +
  geom_line(color = "#2E8B57", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 80, y = 1.5, label = "R0 = 1 Threshold", color = "red") +
  labs(
    title = expression("Effect of Latency Delay (" * tau * ") on Reproduction Potential"),
    subtitle = expression("Analysis of " * mathcal(R)[0] * " based on Next-Generation Operator"),
    x = expression("Delay Period " * tau * " (Time units)"),
    y = expression("Basic Reproduction Number " * mathcal(R)[0])
  ) +
  theme_minimal()

print(plot_r0)

#======================================================
# Information value frontier and Latency Penality
#======================================================

library(tidyverse)
library(extraDistr)
beta <- 0.5; gamma <- 0.05; kappa <- 8
# --- 1. SETTINGS & COSTS ---
costs <- list(C_epi = 12000, 
              eta = 0.4, 
              C_local = 200, 
              C_regional = 600, 
              c0 = 30, 
              c1 = 2)

# --- 2. THE OPTIMIZATION ENGINE ---
compute_n_star <- function(mu, tau) {
  # Observability and Risk Scaling
  rho <- exp(-gamma * tau)
  Rj_factor <- (beta/gamma) * (1 + exp(-gamma * tau))
  
  # Decision Grid
  n_grid <- seq(0, 200, by = 5)
  prec <- 6 # High uncertainty to drive sampling value
  a_prior <- mu * prec; b_prior <- (1 - mu) * prec
  
  # Calculate Expected Loss for Action A (A0, A1, A2)
  get_expected_loss <- function(act, a, b) {
    thetas <- seq(0.001, 0.999, length.out = 50)
    probs <- dbeta(thetas, a, b)
    losses <- sapply(thetas, function(th) {
      Rj <- Rj_factor * th
      epi <- if(act %in% c("A0", "A1") && th >= theta1) costs$C_epi * exp(costs$eta*(Rj-1)) else 0
      econ <- switch(act, "A0"=0, "A1"=200, "A2"=600)
      return(epi + econ)
    })
    return(sum(losses * probs) / sum(probs))
  }
  
  # Prior Decision
  L_prior <- min(sapply(c("A0", "A1", "A2"), function(act) get_expected_loss(act, a_prior, b_prior)))
  
  # ENGS over n
  engs <- sapply(n_grid, function(n) {
    if (n == 0) return(0)
    p_obs <- mu * rho
    x_vals <- round(seq(0, n, length.out = min(n+1, 20))) # Sampled for speed
    
    # Probability of observing x symptoms
    p_x <- dbbinom(x_vals, n, p_obs * prec, (1 - p_obs) * prec)
    
    # Posterior Loss
    L_post <- sum(p_x * sapply(x_vals, function(x) {
      # Adjust posterior for hidden burden (division by rho)
      mu_post <- ((p_obs * prec) + x) / (prec + n) / rho
      mu_post <- min(mu_post, 0.99)
      a_post <- mu_post * (prec + n); b_post <- (1 - mu_post) * (prec + n)
      min(sapply(c("A0", "A1", "A2"), function(act) get_expected_loss(act, a_post, b_post)))
    }))
    
    return((L_prior - L_post) - (costs$c0 + costs$c1 * n))
  })
  
  res <- n_grid[which.max(engs)]
  return(if(max(engs) <= 0) 0 else res)
}

# --- 3. RUN SCENARIOS ---
scenario_data <- expand.grid(
  Phase = factor(c("Onset", "Endemic"), levels = c("Onset", "Endemic")),
  Latency = factor(c("Early Delay", "Late Latency"), levels = c("Early Delay", "Late Latency"))
) %>%
  mutate(
    mu = ifelse(Phase == "Onset", 0.09, 0.45),
    tau = ifelse(Latency == "Early Delay", 5, 25)
  )

scenario_data$n_star <- mapply(compute_n_star, scenario_data$mu, scenario_data$tau)

# --- 4. VISUALIZATION ---
# 

ggplot(scenario_data, aes(x = Phase, y = n_star, fill = Latency)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  geom_text(aes(label = n_star), vjust = -0.5, position = position_dodge(width = 0.8), fontface = "bold") +
  scale_fill_manual(values = c("Early Delay" = "#56B4E9", "Late Latency" = "#D55E00")) +
  labs(title = "Optimization Results: Sampling n* by Epidemic Context",
       subtitle = "The 'Latency Penalty' significantly inflates effort only during the Onset phase.",
       y = expression("Optimal Sample Size (" * n^{"*"} * ")"),
       x = "Epidemic Phase") +
  theme_minimal()


# Generating the continuous Frontier for Early vs Late Latency
mu_grid <- seq(0.01, 0.45, length.out = 50)

frontier_results <- expand.grid(mu = mu_grid, tau = c(5, 25)) %>%
  mutate(Latency = ifelse(tau == 5, "Early Delay", "Late Latency"))

# Apply the compute_n_star function to every point on the grid
frontier_results$n_star <- mapply(compute_n_star, frontier_results$mu, frontier_results$tau)

# Visualization of the Frontier
ggplot(frontier_results, aes(x = mu, y = n_star, color = Latency)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = c(0.1, 0.25), linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.11, y = 140, label = "Onset Threshold", angle = 90, size = 3) +
  annotate("text", x = 0.26, y = 140, label = "Export Threshold", angle = 90, size = 3) +
  labs(title = "The Optimization Frontier",
       subtitle = "Comparing sampling effort requirements across biological delays",
       x = "Estimated Infectious Burden (Prevalence)",
       y = "Optimal Sample Size (n*)") +
  theme_classic()

#===============================================
# Sensitivity analysis on the Optimal effort
#===============================================

############################################################
## COMPLETE GLOBAL SENSITIVITY ANALYSIS: OPTIMAL EFFORT (n*)
############################################################

library(lhs)
library(tidyverse)
library(extraDistr)

# --- 1. SETTINGS & STRESS TEST RANGES ---
set.seed(456)
n_lhs <- 60  # Number of parameter sets to test

# Parameters to vary: 
# beta (trans.), gamma (removal), tau (latency), C_epi (risk), c1 (unit cost)
factors <- c("beta", "gamma", "tau", "C_epi", "c1")
lower <- c(0.20, 0.02, 5,   4000,  0.1)  # Low-risk / Cheap sampling
upper <- c(0.80, 0.10, 80, 25000, 10.0) # High-risk / Expensive sampling

# --- 2. GENERATE AND SCALE LHS DESIGN ---
lhs_design <- lhs::maximinLHS(n_lhs, length(factors))
lhs_params <- as.data.frame(lhs_design)
colnames(lhs_params) <- factors

for(i in 1:length(factors)) {
  lhs_params[,i] <- lower[i] + lhs_params[,i] * (upper[i] - lower[i])
}

# --- 3. DYNAMIC EVALUATION ENGINE ---
evaluate_n_star_sensitivity <- function(row) {
  
  # Extract parameters for this specific run
  b_curr     <- as.numeric(row["beta"])
  g_curr     <- as.numeric(row["gamma"])
  tau_curr   <- as.numeric(row["tau"])
  C_epi_curr <- as.numeric(row["C_epi"])
  c1_curr    <- as.numeric(row["c1"])
  
  # Decision Constants
  mu   <- 0.11     # Fixed at Onset Threshold to observe sensitivity
  prec <- 5        # Prior precision
  rho  <- exp(-g_curr * tau_curr)
  Rj_f <- (b_curr/g_curr) * (1 + exp(-g_curr * tau_curr))
  th1  <- g_curr / b_curr
  
  # Decision Grid for n*
  n_grid <- seq(0, 150, by = 5)
  
  # Prior Loss Calculation
  L_prior <- min(sapply(c("A0", "A1", "A2"), function(act) {
    Rj <- Rj_f * mu
    epi <- if(act %in% c("A0", "A1") && mu >= th1) C_epi_curr * exp(0.4*(Rj-1)) else 0
    econ <- switch(act, "A0"=0, "A1"=200, "A2"=600)
    return(epi + econ)
  }))
  
  # Expected Net Gain of Sampling (ENGS)
  engs <- sapply(n_grid, function(n) {
    if (n == 0) return(0)
    p_obs  <- mu * rho
    x_vals <- round(seq(0, n, length.out = min(n+1, 15))) 
    p_x    <- dbbinom(x_vals, n, p_obs * prec, (1 - p_obs) * prec)
    
    L_post <- sum(p_x * sapply(x_vals, function(x) {
      mu_p <- min(((p_obs * prec) + x) / (prec + n) / rho, 0.99)
      Rj_p <- Rj_f * mu_p
      min(
        (if(mu_p >= th1) C_epi_curr * exp(0.4*(Rj_p-1)) else 0),             
        (if(mu_p >= th1) C_epi_curr * exp(0.4*(Rj_p-1)) else 0) + 200,      
        600                                                                 
      )
    }))
    return((L_prior - L_post) - (30 + c1_curr * n))
  })
  
  return(n_grid[which.max(engs)])
}

# --- 4. EXECUTE SENSITIVITY RUNS ---
cat("Running Global Sensitivity Analysis on", n_lhs, "points...\n")
lhs_params$n_star <- apply(lhs_params, 1, evaluate_n_star_sensitivity)

# --- 5. ROBUST CORRELATION & VISUALIZATION ---
if(var(lhs_params$n_star) > 0) {
  
  # Calculate Spearman Rank Correlation (Standard for sensitivity)
  cor_results <- sapply(factors, function(f) {
    cor(lhs_params[[f]], lhs_params$n_star, method = "spearman")
  })
  
  prcc_df <- data.frame(
    Parameter = factors,
    Correlation = as.numeric(cor_results)
  )
  
  # TORNADO PLOT
  
  ggplot(prcc_df, aes(x = reorder(Parameter, Correlation), y = Correlation, fill = Correlation > 0)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#2E8B57", "FALSE" = "#D55E00"), guide = "none") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    labs(
      title = "Global Sensitivity of Optimal Sample Size (n*)",
      subtitle = "Spearman Rank Correlation: Biological & Economic Drivers",
      x = "Model Parameter",
      y = "Correlation with Optimal Effort (n*)"
    ) +
    theme_minimal()
  
} else {
  cat("No variance in n_star. Try broadening ranges further.")
}

# Print summary table
print(prcc_df %>% arrange(desc(abs(Correlation))))

#=================================================
# Spatial distribution 
#=================================================

# 1. Fetch Benin Map (Level 1)
benin_sf <- gadm(country = "BEN", level = 1, path = tempdir()) %>% 
  st_as_sf() %>% st_make_valid()

# 2. Border Anchors (Transboundary Entry Points)
border_coords <- matrix(c(2.65, 6.45, 2.70, 9.20, 1.65, 6.25, 1.70, 9.70, 2.50, 7.50), ncol = 2, byrow = TRUE)
confirmed_cases <- st_as_sf(as.data.frame(border_coords), coords = c("V1", "V2"), crs = 4326)

# 3. Scientific Map Plotting Function
plot_science_benin <- function(scenario_name, M_farms, decay_val, max_prev, palette_opt) {
  
  # Create high-res grid for IDW
  bbox <- st_bbox(benin_sf)
  grid <- expand.grid(
    x = seq(bbox[1], bbox[3], length.out = 200), # Increased resolution
    y = seq(bbox[2], bbox[4], length.out = 200)
  ) %>% st_as_sf(coords = c("x", "y"), crs = 4326) %>% st_intersection(st_union(benin_sf))
  
  # Calculate Biological Risk Surface
  dists_grid <- st_distance(grid, confirmed_cases)
  grid$prevalence <- exp(-as.numeric(apply(dists_grid, 1, min)) / decay_val) * max_prev
  
  # Optimization Weight (Entropy-based: P*(1-P))
  grid$opt_weight <- grid$prevalence * (1 - (grid$prevalence/max_prev))
  
  # Sample M_farms using Optimized Weights
  set.seed(123)
  optimized_farms <- grid[sample(1:nrow(grid), size = M_farms, prob = grid$opt_weight), ]
  
  # Convert to DF for geom_tile
  idw_df <- cbind(as.data.frame(st_coordinates(grid)), prevalence = grid$prevalence)
  
  ggplot() +
    # Layer 1: Prevalence Surface
    geom_tile(data = idw_df, aes(x = X, y = Y, fill = prevalence)) +
    scale_fill_viridis_c(option = palette_opt, name = "Prevalence", limits = c(0, max_prev), labels = label_percent()) +
    
    # Layer 2: Benin Department Borders
    geom_sf(data = benin_sf, fill = NA, color = "white", linewidth = 0.2, alpha = 0.5) +
    
    # Layer 3: Optimized Farm Locations (Strategic dots)
    geom_sf(data = optimized_farms, color = "black", size = 0.6, alpha = 0.9) +
    
    # Layer 4: Confirmed Entry Points (Materialized)
    geom_sf(data = confirmed_cases, shape = 21, fill = "white", color = "black", size = 2, stroke = 0.8) +
    
    # --- Scientific Additions ---
    annotation_scale(location = "bl", width_hint = 0.3, style = "ticks", text_family = "serif") +
    annotation_north_arrow(location = "tr", which_north = "true", 
                           pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
                           style = north_arrow_fancy_orienteering(text_family = "serif")) +
    
    labs(title = scenario_name, 
         subtitle = paste("Optimized Grid: M =", M_farms, "Farms")) +
    theme_minimal() +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", size = 12),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      legend.position = "right"
    )
}

# 4. Generate Final Figures
# Onset: Low prevalence (15% max), sparse entry-point clustering
p_onset <- plot_science_benin("Onset Scenario (Detection)", 50, 25000, 0.15, "magma")

# Endemic: High prevalence (80% max), national inland expansion
p_endemic <- plot_science_benin("Endemic Scenario (Mitigation)", 80, 160000, 0.80, "viridis")

# Side-by-Side Layout
final_map <- p_onset + p_endemic + plot_annotation(
  title = "Optimized Spatial Surveillance Allocation for Pathogen Incursion in Benin",
  theme = theme(plot.title = element_text(family = "serif", face = "bold", size = 16, hjust = 0.5))
)

final_map







