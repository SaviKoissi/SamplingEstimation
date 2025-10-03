library(deSolve)
library(ggplot2)
library(reshape2)
library(lhs)
library(sensitivity)


# Model definition
cutting_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Delayed latent
    L_tau <- if (t <= tau) L else lagvalue(t - tau, 2)
    
    # Behavioral function
    psi <- I / N
    f_psi <- psi / (1 + kappa * psi)
    
    # Force of infection
    lambda_local <- beta * I / N
    lambda_external <- alpha * I_ext / N_ext * exp(-delta * d_jk)
    lambda_total <- lambda_local + lambda_external
    
    # Differential equations
    dS <- -lambda_total * S + gamma * N * (1 - f_psi) - gamma * S
    dL <-  lambda_total * S - sigma * L + gamma * N * f_psi - gamma * L
    dI <-  sigma * L_tau - gamma * L
    
    return(list(c(dS, dL, dI)))
  })
}

# Initial state
state <- c(S = 1000, L = 1, I = 1)

# Parameters
params <- c(
  beta = 0.05,
  sigma = 0.5,
  gamma = 0.05,
  kappa = 3,
  tau = 100,
  alpha = 0.33,      # coupling from other regions
  I_ext = 1,        # external infection level
  N_ext = 1000,
  delta = 0.05,      # spatial decay
  d_jk = 5.25,          # distance from patch j to k
  N = sum(state)
)

# Time
times <- seq(0, 180, by = 1)

# Solve DDE
out <- dede(y = state, times = times, func = cutting_model, parms = params)
out_df <- as.data.frame(out)

# Plot
out_long <- melt(out_df, id = "time")
ggplot(out_long, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Cassava Disease Dynamics (with External Pressure)", y = "Plant Count") +
  theme_minimal()



# --- Sensitivity Analysis for Extended Model ---

# 1. Define parameter ranges
par_ranges <- data.frame(
  beta = c(0.01, 1),
  sigma = c(0.01, 0.5),
  gamma = c(0.01, 0.3),
  kappa = c(0.1, 10),
  alpha = c(0.001, 0.1),
  tau = c(0.001, 200),
  delta = c(0.01, 0.2)
)

param_names <- colnames(par_ranges)

# 2. Latin Hypercube Sampling
set.seed(123)
n_samples <- 1000
lhs_sample <- randomLHS(n_samples, length(param_names))

# 3. Scale LHS to actual ranges
param_matrix <- lhs_sample
for (i in 1:ncol(lhs_sample)) {
  param_matrix[, i] <- lhs_sample[, i] * 
    (par_ranges[2, i] - par_ranges[1, i]) + par_ranges[1, i]
}
colnames(param_matrix) <- param_names

# 4. Run model for each parameter set
output_I <- numeric(n_samples)

for (i in 1:n_samples) {
  current_params <- c(
    param_matrix[i, ],
    N = 1010,
    tau = 100,
    I_ext = 10,
    N_ext = 1000,
    d_jk = 1
  )
  
  tryCatch({
    out <- dede(y = state, times = times, func = cutting_model, parms = current_params)
    output_I[i] <- as.data.frame(out)[nrow(out), "I"]  # Final I
  }, error = function(e) {
    output_I[i] <- NA
  })
}

# 5. Clean up failed runs
valid <- !is.na(output_I)
param_matrix <- param_matrix[valid, ]
output_I <- output_I[valid]

# 6. PRCC Analysis
prcc_result <- pcc(X = as.data.frame(param_matrix), y = output_I, rank = TRUE)

# 7. Plot PRCC

prcc_df <- data.frame(
  Parameter = rownames(prcc_result$PRCC),
  PRCC = prcc_result$PRCC[, 1]
)

ggplot(prcc_df, aes(x = reorder(Parameter, PRCC), y = PRCC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "PRCC Sensitivity Analysis (Final I)", x = "Parameter", y = "PRCC") +
  theme_minimal()


#---==========stab


# Generate trajectories from different initial conditions
init_conditions <- list(
  c(S = 1000, L = 1, I = 1),
  c(S = 900, L = 50, I = 50),
  c(S = 800, L = 100, I = 100),
  c(S = 950, L = 20, I = 30)
)

trajectories <- lapply(seq_along(init_conditions), function(i) {
  out <- dede(y = init_conditions[[i]], times = times, func = cutting_model, 
              parms = params)
  as.data.frame(out)[, c("time", "S", "L", "I")]
})

# Add trajectory ID
for (i in seq_along(trajectories)) {
  trajectories[[i]]$trajectory <- paste0("traj_", i)
}
all_traj <- do.call(rbind, trajectories)

# Plot phase portrait: S vs I
ggplot(all_traj, aes(x = S, y = I, color = trajectory)) + 
  geom_path(size = 1.1) +
  facet_wrap(~ trajectory) +
  labs(title = "Phase Portrait: S vs I",
       x = "Susceptible (S)",
       y = "Infected (I)") +
  theme_minimal()

# Plot phase portrait: L vs I
ggplot(all_traj, aes(x = L, y = I, color = trajectory)) +
  geom_path(size = 1.1) +
  facet_wrap(~ trajectory) +
  labs(title = "Phase Portrait: L vs I",
       x = "Latent (L)",
       y = "Infected (I)") +
  theme_minimal()

#--------------========


# Fixed parameters
sigma <- 0.5
gamma <- 0.05
delta <- 0.05
d_jk <- 5.25
N_j <- 1000
N_k <- 1000

# Grid of beta and alpha
beta_vals <- seq(0.01, 1, length.out = 180)
alpha_vals <- seq(0.001, 0.5, length.out = 180)

# Create grid
grid <- expand.grid(beta = beta_vals, alpha = alpha_vals)

# Compute R0 for each (beta, alpha)
grid$R0 <- with(grid, {
  local_term <- beta * (N_j / N_k)
  spatial_term <- alpha * exp(-delta * d_jk) * (N_j / N_k)
  (local_term + spatial_term) / (sigma + gamma)
})

# Add stability label
grid$Stability <- ifelse(grid$R0 < 1, "Stable (R0 < 1)", "Unstable (R0 > 1)")

# Plot
ggplot(grid, aes(x = alpha, y = beta, fill = Stability)) +
  geom_tile() +
  scale_fill_manual(values = c("Stable (R0 < 1)" = "#91cf60", "Unstable (R0 > 1)" = "#d73027")) +
  labs(title = expression("Stability Region based on " ~ R[0]),
       x = expression(alpha),
       y = expression(beta)) +
  theme_minimal()



# --- Sample Size Calculation for Detection Probability ---

# Function to compute required sample size for detection
sample_size_calc <- function(p_detect, gamma_target = 0.95) {
  # n >= log(1 - gamma) / log(1 - p_detect)
  n <- ceiling(log(1 - gamma_target) / log(1 - p_detect))
  return(n)
}

# Example: compute detection probability from model output
# Here we approximate prevalence as I/N
out_df$prevalence <- out_df$I / params["N"]

# Assume diagnostic test sensitivity/specificity
Se <- 0.9
Sp <- 0.95

# Effective detection probability per plant
out_df$p_detect <- Se * out_df$prevalence + (1 - Sp) * (1 - out_df$prevalence)

# Required sample size per time step (95% confidence of detection)
out_df$sample_size <- sapply(out_df$p_detect, sample_size_calc, gamma_target = 0.95)

# Plot prevalence vs required sample size
ggplot(out_df, aes(x = prevalence, y = sample_size)) +
  geom_line(color = "darkred", linewidth = 1.2) +
  labs(title = "Optimal Sample Size vs Prevalence",
       x = "Infection Prevalence",
       y = "Required Sample Size (95% detection)") +
  theme_minimal()


# --- Scenario Comparison of Sample Size Requirements ---

scenarios <- data.frame(
  Scenario = c("Low prevalence", "Moderate prevalence", "High prevalence"),
  prevalence = c(0.01, 0.05, 0.2)
)

scenarios$p_detect <- Se * scenarios$prevalence + (1 - Sp) * (1 - scenarios$prevalence)
scenarios$sample_size <- sapply(scenarios$p_detect, sample_size_calc, gamma_target = 0.95)

# Bar plot mapping scenarios to sample sizes
ggplot(scenarios, aes(x = Scenario, y = sample_size, fill = Scenario)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sample_size), vjust = -0.5, size = 4.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Optimal Sample Size by Prevalence Scenario",
       x = NULL,
       y = "Sample Size (95% detection)") +
  theme_minimal() +
  theme(legend.position = "none")

#============Spatial Application ==========

# Required packages
library(sf)
library(sp)
library(raster)
library(gstat)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggnewscale)   # for separate color scales

set.seed(123)

# ------------ Parameters (easy to tweak) -----------------
n_farms <- 40
grid_res <- 0.005         # raster resolution (~5 km at equator)
idw_power <- 0.1            # IDW power (p); 0.001 ~ near-uniform, 1 = sensible gradients
kernel_eta <- 1         # exponent for power-law dispersal kernel
prior_strength <- 8       # sum a0+b0 effective prior "concentration"
pi_min <- 0.05            # detection threshold for "early detection" utility
batch_per_round <- 5      # number of samples allocated to each selected farm at each round
max_rounds <- 6           # maximum rounds of sequential sampling
n_select_per_round <- 5   # how many farms to sample each round
# per-farm sample target for scenarios:
target_samples_early <- 30
target_samples_established <- 50
# ---------------------------------------------------------

# 1. Load Benin boundary (naturalearth)
benin <- ne_countries(country = "Benin", scale = "medium", returnclass = "sf")

# 2. Generate random farm locations inside Benin
farms_pts <- st_sample(benin, size = n_farms, type = "random")
farms_sf <- st_sf(farm_id = 1:n_farms, geometry = farms_pts)

# 3. Assign synthetic 'true' prevalence (the unknown truth)
farms_sf$true_prev <- runif(n_farms, min = 0.01, max = 0.20)

# 4. Build prediction grid (raster) and perform IDW to get prior surface
# Convert to Spatial for gstat
farms_sp <- as(farms_sf, "Spatial")
r <- raster(benin, res = grid_res)                 # grid
grid_sp <- as(r, "SpatialPixelsDataFrame")
# If proj4string missing, use a default: safer to use st_crs, but keep for consistency with gstat
proj4string(grid_sp) <- proj4string(farms_sp)

# IDW interpolation (prevalence ~ 1)
idw_model <- gstat::idw(formula = true_prev ~ 1,
                        locations = farms_sp,
                        newdata = grid_sp,
                        idp = idw_power)          # set idw_power
idw_raster <- raster(idw_model)                    # raster of predicted prevalences
names(idw_raster) <- "idw_prev"

# Extract IDW-prevalence at farm locations (used to set prior mean per farm)
farms_sf$idw_prior_mean <- raster::extract(idw_raster, as(farms_sf, "Spatial"))

# sanitize any NA extraction (in case of edge issues)
if (any(is.na(farms_sf$idw_prior_mean))) {
  farms_sf$idw_prior_mean[is.na(farms_sf$idw_prior_mean)] <- mean(farms_sf$idw_prior_mean, na.rm = TRUE)
}

# 5. Build dispersal kernel K(d) between farms (rows sum to 1)
coords <- st_coordinates(farms_sf)
dist_mat <- as.matrix(dist(coords))                 # Euclidean distances between farms
# avoid exact zeros on diagonal for d^-eta: set diag to Inf to avoid self-influence
diag(dist_mat) <- Inf
K_raw <- dist_mat^(-kernel_eta)
# if any rows sum to zero (very unlikely), set them to small positive for numerical stability
row_sums_raw <- rowSums(K_raw, na.rm = TRUE)
row_sums_raw[row_sums_raw == 0] <- 1e-12
K <- K_raw / row_sums_raw                              # normalize rows

# 6. Helper: initialize priors (Beta) from IDW prior mean
init_priors <- function(farms, prior_strength) {
  a0 <- prior_strength * farms$idw_prior_mean
  b0 <- prior_strength * (1 - farms$idw_prior_mean)
  farms$a_post <- a0
  farms$b_post <- b0
  farms$posterior_mean <- farms$a_post / (farms$a_post + farms$b_post)
  farms$post_var <- (farms$a_post * farms$b_post) / ((farms$a_post + farms$b_post)^2 * (farms$a_post + farms$b_post + 1))
  farms$samples_taken <- integer(nrow(farms))
  farms$positives_total <- integer(nrow(farms))
  return(farms)
}

# 7. Utilities
# early detection utility: prob(pi > thr) under Beta(a,b) = 1 - CDF
detection_prob <- function(a, b, thr) {
  # ensure numeric vectors and positive shapes
  a <- pmax(a, .Machine$double.eps)
  b <- pmax(b, .Machine$double.eps)
  1 - pbeta(thr, a, b)
}

# monitoring utility: posterior variance (we want to reduce it)
posterior_variance <- function(a, b) {
  (a * b) / ((a + b)^2 * (a + b + 1))
}

# 8. Sequential sampling routine (one function) ------------------------------
sequential_sampling <- function(farms_init, scenario_name = "A",
                                target_samples = 30,
                                n_rounds = max_rounds,
                                batch = batch_per_round,
                                n_select = n_select_per_round,
                                regime = c("detection", "monitoring"),
                                K_mat = K,
                                true_prev_col = "true_prev",
                                pi_threshold = NULL) {
  regime <- match.arg(regime)
  farms <- farms_init
  farms$scenario <- scenario_name
  
  # if threshold not provided, use global pi_min
  if (is.null(pi_threshold)) pi_threshold <- pi_min
  
  # record history
  history <- vector("list", n_rounds)
  
  for (r in seq_len(n_rounds)) {
    # compute utilities
    if (regime == "detection") {
      farms$util <- detection_prob(farms$a_post, farms$b_post, pi_threshold)
    } else {
      farms$util <- posterior_variance(farms$a_post, farms$b_post)
    }
    
    # choose candidate farms that still need samples (haven't reached target)
    need_idx <- which(farms$samples_taken < target_samples)
    if (length(need_idx) == 0) break
    
    # rank by utility within the 'need' set
    candidate_order <- order(farms$util[need_idx], decreasing = TRUE)
    selected_idx <- need_idx[candidate_order][1:min(n_select, length(need_idx))]
    
    # allocate batch samples to each selected farm
    for (i in selected_idx) {
      to_take <- min(batch, target_samples - farms$samples_taken[i])
      if (to_take <= 0) next
      # simulate observations: binomial with true prevalence
      y <- rbinom(1, size = to_take, prob = farms[[true_prev_col]][i])
      # update posterior Beta counts
      farms$a_post[i] <- farms$a_post[i] + y
      farms$b_post[i] <- farms$b_post[i] + (to_take - y)
      farms$positives_total[i] <- farms$positives_total[i] + y
      farms$samples_taken[i] <- farms$samples_taken[i] + to_take
    }
    
    # recompute posterior mean and var
    farms$posterior_mean <- farms$a_post / (farms$a_post + farms$b_post)
    farms$post_var <- posterior_variance(farms$a_post, farms$b_post)
    
    # --- propagate information via kernel K: smooth posterior means across neighbors ---
    neighbor_mean <- as.numeric(K_mat %*% farms$posterior_mean)
    
    # weight for kernel influence (tuneable)
    w_kernel <- 0.25
    smoothed_mean <- (1 - w_kernel) * farms$posterior_mean + w_kernel * neighbor_mean
    
    # convert smoothed mean back to a_post/b_post while preserving concentration (a+b)
    total_counts <- farms$a_post + farms$b_post
    # avoid zero total_counts - set min to 1e-8 for safety
    total_counts <- pmax(total_counts, 1e-8)
    farms$a_post <- smoothed_mean * total_counts
    farms$b_post <- (1 - smoothed_mean) * total_counts
    
    # recompute posterior stats after smoothing
    farms$posterior_mean <- farms$a_post / (farms$a_post + farms$b_post)
    farms$post_var <- posterior_variance(farms$a_post, farms$b_post)
    
    # store round snapshot
    history[[r]] <- list(round = r,
                         selected = farms$farm_id[selected_idx],
                         samples_taken = farms$samples_taken,
                         posterior_mean = farms$posterior_mean,
                         post_var = farms$post_var)
    
    # termination if all farms reached target
    if (all(farms$samples_taken >= target_samples)) break
  }
  
  return(list(farms = farms, history = history))
}

# 9. Initialize priors (based on IDW prior mean)
farms_init <- init_priors(farms_sf, prior_strength = prior_strength)

# 10. Run Scenario A: Early detection regime (target_samples_early)
resA <- sequential_sampling(farms_init,
                            scenario_name = "A_Early",
                            target_samples = target_samples_early,
                            n_rounds = max_rounds,
                            batch = batch_per_round,
                            n_select = n_select_per_round,
                            regime = "detection",
                            K_mat = K,
                            pi_threshold = pi_min)

# 11. Run Scenario B: Established epidemic (target_samples_established)
# re-init priors fresh (independent run)
farms_init2 <- init_priors(farms_sf, prior_strength = prior_strength)
resB <- sequential_sampling(farms_init2,
                            scenario_name = "B_Established",
                            target_samples = target_samples_established,
                            n_rounds = max_rounds,
                            batch = batch_per_round,
                            n_select = n_select_per_round,
                            regime = "monitoring",
                            K_mat = K,
                            pi_threshold = pi_min)

# 12. Prepare plotting objects ----------------------------------------------
# Convert idw raster to data.frame for ggplot
idw_df <- as.data.frame(rasterToPoints(idw_raster))
colnames(idw_df) <- c("x", "y", "idw_prev")

# farms dataframes for plotting
farmsA <- resA$farms
farmsB <- resB$farms

# choose next sampling candidates to mark on map: those with highest current utility
get_next_candidates <- function(res_obj, n = 5, regime = "detection", pi_threshold = pi_min) {
  farms <- res_obj$farms
  if (regime == "detection") {
    farms$util <- detection_prob(farms$a_post, farms$b_post, pi_threshold)
  } else {
    farms$util <- posterior_variance(farms$a_post, farms$b_post)
  }
  top_idx <- order(farms$util, decreasing = TRUE)[1:min(n, nrow(farms))]
  return(farms[top_idx, ])
}
nextA <- get_next_candidates(resA, n = 5, regime = "detection", pi_threshold = pi_min)
nextB <- get_next_candidates(resB, n = 5, regime = "monitoring", pi_threshold = pi_min)

# 13. Plotting function (IDW transparent + farms posterior colored)
plot_scenario_map <- function(idw_df, benin_sf, farms_plot, next_points, title, target_samples) {
  ggplot() +
    # IDW prevalence raster (transparent)
    geom_raster(data = idw_df, aes(x = x, y = y, fill = idw_prev), alpha = 0.45) +
    scale_fill_viridis_c(name = "IDW prior\nprevalence", option = "A") +
    ggnewscale::new_scale_fill() +
    
    # country outline on top
    geom_sf(data = benin_sf, fill = NA, color = "black", size = 0.5) +
    
    # farms: posterior mean colored with a separate scale
    geom_sf(data = farms_plot, aes(color = posterior_mean), size = 2) +
    scale_color_viridis_c(name = "Posterior\nmean", option = "C") +
    
    # next sampling locations
    geom_sf(data = next_points, color = "red", size = 3, shape = 8, stroke = 1.2) +
    
    # annotate sample target in subtitle and caption
    labs(title = title,
         subtitle = paste0("Per-farm target samples = ", target_samples,
                           " | Red stars = next sampling locations"),
         caption = "IDW layer= spatial prior; colored points = posterior mean after updates") +
    theme_minimal()
}

# 14. Create maps
mapA <- plot_scenario_map(idw_df, benin, farmsA, nextA,
                          title = "Scenario A: Early detection ",
                          target_samples = target_samples_early)
mapB <- plot_scenario_map(idw_df, benin, farmsB, nextB,
                          title = "Scenario B: Established epidemic",
                          target_samples = target_samples_established)

# 15. Display maps (separately)
print(mapA)
print(mapB)

library(patchwork)

# Combine side by side
combined_map <- mapA + mapB + plot_layout(ncol = 2)

# Or if you prefer stacked vertically:
# combined_map <- mapA / mapB

# Display
print(combined_map)

# Save to file
ggsave("combined_map.png", combined_map, width = 12, height = 6)

# 16. Optional: Summarize final sample counts and uncertainties
summaryA <- farmsA %>% st_drop_geometry() %>%
  select(farm_id, true_prev, samples_taken, positives_total, posterior_mean, post_var)
summaryB <- farmsB %>% st_drop_geometry() %>%
  select(farm_id, true_prev, samples_taken, positives_total, posterior_mean, post_var)

# show a quick table in console (first 8 farms)
print(head(summaryA, 8))
print(head(summaryB, 8))


